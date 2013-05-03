// KDEMeshTally.cpp

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"

#include "KDEKernel.hpp"
#include "KDEMeshTally.hpp"
#include "KDENeighborhood.hpp"
#include "TallyEvent.hpp"

// initialize static variables
bool KDEMeshTally::seed_is_set = false;
const char* const KDEMeshTally::kde_tally_names[] = {"collision",
                                                     "integral-track",
                                                     "sub-track"};

// quadrature points and weights for the integral_track_score function
const double quad_points[4] = {0.339981043585, -0.339981043585,
                               0.861136311594, -0.861136311594};

const double quad_weights[4] = {0.652145154863, 0.652145154863,
                                0.347854845137, 0.347854845137};

//-----------------------------------------------------------------------------
static double parse_bandwidth_value(const std::string& key,
                                    const std::string& value,
                                    double default_value = 0.01)
{
    char* end;
    double bandwidth_value = strtod(value.c_str(), &end);

    if( value.c_str() == end || bandwidth_value <= 0 )
    {
        std::cerr << "Warning: invalid bandwidth value for " << key
                  << " = " << value << std::endl;
        std::cerr << "    using default value " << key << " = "
                  << default_value << std::endl;
        bandwidth_value = default_value;
    }

    return bandwidth_value;
}
//-----------------------------------------------------------------------------
KDEMeshTally::KDEMeshTally(const MeshTallyInput& input,
                           KDEMeshTally::TallyType type)
    : MeshTally(input),
      mbi(new moab::Core()),
      bandwidth(moab::CartVect(0.01, 0.01, 0.01)),
      kde_tally(type),
      kernel(NULL),
      num_subtracks(3)
{
    std::cout << "Creating KDE " << kde_tally_names[kde_tally]
              << " mesh tally " << input.tally_id << std::endl;

    std::cout << "    for input mesh: " << input.input_filename
              << ", output file: " << output_filename << std::endl;

    // set up KDEMeshTally member variables from MeshTallyInput
    parse_tally_options();

    // create second-order epanechnikov kernel if user did not specify type
    if (kernel == NULL)
    {
        kernel = KDEKernel::createKernel("epanechnikov");
    }
 
    std::cout << "    using " << kernel->get_kernel_name()
              << " kernel and bandwidth " << bandwidth << std::endl;

    if (kde_tally == SUB_TRACK)
    {
        std::cout << "    splitting full tracks into "
                  << num_subtracks << " sub-tracks" << std::endl;

        // set seed value if not already set by another instance
        if (!seed_is_set)
        {
            srand(time(NULL));
            seed_is_set = true;
        }
    }

    // initialize MeshTally member variables representing the mesh data
    moab::ErrorCode rval = initialize_mesh_data();

    if (rval != moab::MB_SUCCESS)
    {
        std::cout << "Error: Could not load mesh data for KDE mesh tally "
                  << input_data.tally_id << " from input file "
                  << input_data.input_filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // initialize running variance variables
    max_collisions = false;
    numCollisions = 0;
    Mn = moab::CartVect(0, 0, 0);
    Sn = moab::CartVect(0, 0, 0);
}
//-----------------------------------------------------------------------------
KDEMeshTally::~KDEMeshTally()
{
    delete kd_tree;
    delete kernel;
    delete mbi;  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::compute_score(const TallyEvent& event, int ebin)
{
    // initialize common weighting factor for this tally event
    double weight = event.get_weighting_factor();

    // set up tally event based on KDE mesh tally type
    TrackData track;
    CollisionData collision;
    std::vector<moab::CartVect> subtrack_points;
    bool event_is_set = false;

    if (kde_tally == INTEGRAL_TRACK || kde_tally == SUB_TRACK)
    {
        event_is_set = event.get_track_data(track);

        if (kde_tally == SUB_TRACK && event_is_set == true)
        {
            // multiply weight by track length and set up sub-track points
            weight *= track.track_length;
            subtrack_points = choose_points(track, num_subtracks);
        }
    }
    else // kde_tally == COLLISION
    {
        event_is_set = event.get_collision_data(collision);

        if (event_is_set)
        {
            // divide weight by cross section and update optimal bandwidth
            weight /= collision.total_cross_section;
            update_variance(collision.collision_point);
        }
    }

    // check that a valid tally event has been set for this KDE mesh tally
    if (!event_is_set)
    {
        std::cerr << "Error: Tally event is not valid for KDE mesh tally ";
        std::cerr << input_data.tally_id << std::endl;
        exit(EXIT_FAILURE);
    }

    // create the neighborhood region and find all of the calculations points
    KDENeighborhood region(event, bandwidth, *kd_tree, kd_tree_root);

    std::vector<moab::EntityHandle> calculation_points;
    moab::ErrorCode rval = region.get_points(calculation_points);

    assert(moab::MB_SUCCESS == rval);

    // iterate through the calculation points
    std::vector<moab::EntityHandle>::iterator i;
    double coords[3];

    for (i = calculation_points.begin(); i != calculation_points.end(); ++i)
    {
        // get coordinates of this point
        moab::EntityHandle point = *i;
        rval = mbi->get_coords(&point, 1, coords);

        assert(moab::MB_SUCCESS == rval);

        // compute the final contribution to the tally for this point
        double score = weight;

        if (kde_tally == INTEGRAL_TRACK)
        {
            score *= integral_track_score(track, moab::CartVect(coords));
        }
        else if (kde_tally == SUB_TRACK)
        {
            score *= subtrack_score(subtrack_points, moab::CartVect(coords));
        }
        else // kde_tally == COLLISION
        {
            score *= collision_score(collision, moab::CartVect(coords));
        }

        // add score to KDE mesh tally for the current history
        add_score_to_tally(point, score, ebin);
    }
}
//-----------------------------------------------------------------------------
void KDEMeshTally::end_history()
{
  
  std::set<moab::EntityHandle>::iterator i;
 
  // add results from current history to the tally for each calculation point
  for ( i = visited_this_history.begin() ; i != visited_this_history.end() ; ++i ) {
    
    for( unsigned int j = 0; j < num_energy_bins; ++j ){
      double& history = get_data( temp_tally_data, *i, j );
      double& tally =   get_data( tally_data, *i, j );
      double& error =   get_data( error_data, *i, j );
      
      tally += history;
      error += ( history * history );
      
      // reset temp_tally for the next particle history
      history = 0;
      
    }
  }
  visited_this_history.clear();

}
//-----------------------------------------------------------------------------
void KDEMeshTally::print( double sp_norm, double fmesh_fact )
{

  // tags tally/error results to the nodes and writes mesh to output file
  write_results( sp_norm, fmesh_fact );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::write_results( double sp_norm, double fmesh_fact )
{

  double tally = 0;
  double error = 0, rel_err = 0;
  
  moab::ErrorCode rval = moab::MB_SUCCESS;

  // print the optimal bandwidth if it was computed
  if ( kde_tally == COLLISION ) {
  
    std::cout << std::endl << "optimal bandwidth for " << numCollisions;
    std::cout  << " collisions is: " << get_optimal_bandwidth() << std::endl;

  }

  // tag tally and relative error results to the mesh for each entity
  moab::Range::iterator i;
  
  for ( i = tally_points.begin() ; i != tally_points.end() ; ++i ) {

    moab::EntityHandle point = *i;

    for ( unsigned int j = 0; j < num_energy_bins; ++ j){

      tally = get_data( tally_data, point, j);
      error = get_data( error_data, point, j );
      
      // compute relative error for the tally
      // Use 0 as the rel_err value if nothing has been computed for this tally point;
      // this reflects MCNP's approach to avoiding a divide-by-zero situation.
      rel_err = 0; 
      if( error != 0 ){
        rel_err = sqrt( error / ( tally * tally ) - 1.0 / sp_norm );
      }
      
      // normalizing mesh tally results by the number of source particles
      tally /= sp_norm;
      
      // applying the fmesh multiplication FACTOR to the mesh tally results
      tally *= fmesh_fact;
      
      // set tally and error tag values for this entity
      rval = mbi->tag_set_data( tally_tags[j], &point, 1, &tally );
      assert( moab::MB_SUCCESS == rval );
      
      rval = mbi->tag_set_data( error_tags[j], &point, 1, &rel_err );
      assert( moab::MB_SUCCESS == rval ); 
      
    } 
    
  }

  // create a global tag to store the bandwidth value
  moab::Tag bandwidth_tag;
  rval = mbi->tag_get_handle( "BANDWIDTH_TAG", 3, moab::MB_TYPE_DOUBLE, bandwidth_tag,
                             moab::MB_TAG_MESH|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == rval );

  // add bandwidth tag to the root set
  moab::EntityHandle bandwidth_set = mbi->get_root_set();
  rval = mbi->tag_set_data( bandwidth_tag, &bandwidth_set, 1, &bandwidth );
  assert( moab::MB_SUCCESS == rval );

  // define list of tags to include and write mesh to output file
  std::vector<moab::Tag> output_tags = tally_tags;
  output_tags.insert( output_tags.end(), error_tags.begin(), error_tags.end() );
  output_tags.push_back( bandwidth_tag );  

  rval = mbi->write_file( output_filename.c_str(),
                         NULL, NULL, &tally_mesh_set, 1, &(output_tags[0]), output_tags.size() );
  assert( moab::MB_SUCCESS == rval );
  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::parse_tally_options()
{
    const MeshTallyInput::TallyOptions& options = input_data.options;  
    MeshTallyInput::TallyOptions::const_iterator it;

    for (it = options.begin(); it != options.end(); ++it)
    {
        std::string key = it->first;
        std::string value = it->second;

        // process tally option according to key
        if      (key == "hx") bandwidth[0] = parse_bandwidth_value(key, value);
        else if (key == "hy") bandwidth[1] = parse_bandwidth_value(key, value);
        else if (key == "hz") bandwidth[2] = parse_bandwidth_value(key, value);
        else if (key == "kernel")
        {
            kernel = KDEKernel::createKernel(value);
        }
        else if (key == "seed" && kde_tally == SUB_TRACK)
        {
            // override random number seed if requested by user
            unsigned long int seed = strtol(value.c_str(), NULL, 10);
            srand(seed);
            seed_is_set = true;
            std::cout << "    setting random seed to " << seed
                      << " for choosing sub-track points" << std::endl;
        }
        else if (key == "subtracks" && kde_tally == SUB_TRACK)
        {
            char* end;
            int subtracks = strtol(value.c_str(), &end, 10);

            if (value.c_str() == end || subtracks <= 0)
            {
                std::cerr << "Warning: '" << value << "' is an invalid value"
                          << " for the number of subtracks" << std::endl;
                std::cerr << "    using default value " << key << " = 3\n";
                subtracks = 3;
            }
            
            num_subtracks = subtracks;
        }
        else // invalid tally option
        {
            std::cerr << "Warning: input data for KDE mesh tally "
                      << input_data.tally_id
                      << " has unknown key '" << key << "'" << std::endl;
        }
    }
}
//-----------------------------------------------------------------------------
moab::ErrorCode KDEMeshTally::initialize_mesh_data()
{
    // load the MOAB mesh data from the input file for this KDE mesh tally
    moab::EntityHandle moab_mesh_set;
    moab::ErrorCode rval = load_moab_mesh(mbi, moab_mesh_set);

    if (rval != moab::MB_SUCCESS) return rval;

    // get all of the mesh nodes from the MOAB mesh set
    moab::Range mesh_nodes;
    rval = mbi->get_entities_by_type(moab_mesh_set, moab::MBVERTEX, mesh_nodes);

    if (rval != moab::MB_SUCCESS) return rval;

    // initialize MeshTally::tally_points to include all mesh nodes
    set_tally_points(mesh_nodes);

    // build a kd-tree from all of the mesh nodes
    kd_tree = new moab::AdaptiveKDTree(mbi);
    rval = kd_tree->build_tree(mesh_nodes, kd_tree_root);

    if (rval != moab::MB_SUCCESS) return rval;

    // reduce the loaded MOAB mesh set to include only 3D elements
    moab::Range mesh_cells;
    rval = reduce_meshset_to_3D(mbi, moab_mesh_set, mesh_cells);  

    if (rval != moab::MB_SUCCESS) return rval;

    // initialize MeshTally::tally_mesh_set and set up tags for energy bins
    tally_mesh_set = moab_mesh_set;
    rval = setup_tags(mbi, "KDE_");

    if (rval != moab::MB_SUCCESS) return rval;

    return moab::MB_SUCCESS; 
}
//-----------------------------------------------------------------------------
void KDEMeshTally::update_variance(const moab::CartVect& collision_point)
{
 
  if ( numCollisions != LLONG_MAX ) {

    ++numCollisions;
    
    // obtain previous value for the mean
    moab::CartVect Mn_prev = Mn;
  
    // compute new values for the mean and variance
    if ( numCollisions == 1 )
      Mn = collision_point;
    else {

      Mn += (collision_point - Mn_prev) / numCollisions;
    
      for ( int i = 0 ; i < 3 ; ++i )
        Sn[i] += (collision_point[i] - Mn_prev[i]) * (collision_point[i] - Mn[i]);
    
    }

  }
  else if ( !max_collisions ) {
  
    std::cerr << "Warning: number of collisions exceeds maximum\n"
              << "    optimal bandwidth will be based on " << numCollisions
              << " collisions.\n";

    max_collisions = true;

  }

}
//-----------------------------------------------------------------------------
moab::CartVect KDEMeshTally::get_optimal_bandwidth()
{
  
  double stdev = 0;
  moab::CartVect optimal_bandwidth;
  
  for ( int i = 0 ; i < 3 ; ++i ) {

    stdev = sqrt( Sn[i] / ( numCollisions - 1 ) );
    optimal_bandwidth[i] = 0.968625 * stdev * pow( numCollisions, -1.0/7.0 );

  }

  return optimal_bandwidth;

}
//-----------------------------------------------------------------------------
void KDEMeshTally::add_score_to_tally( moab::EntityHandle mesh_point,
                                       double score,
                                       int ebin )
{

  get_data( temp_tally_data, mesh_point, ebin ) += score;

  // tally the total energy bin if requested
  if ( input_data.total_energy_bin )
    get_data( temp_tally_data, mesh_point, (num_energy_bins-1) ) += score;

  visited_this_history.insert( mesh_point );

}
//-----------------------------------------------------------------------------
// NOTE: integral_track_estimator uses the 4-point gaussian quadrature method
double KDEMeshTally::integral_track_score(const TrackData& data,
                                          const moab::CartVect& X)
{
    // determine the limits of integration
    std::pair<double, double> limits;
    
    bool valid_limits = set_integral_limits(data, X, limits);

    // compute value of the integral only if valid limits exist
    if (valid_limits)
    {
        // define scaling constants
        double c1 = 0.5 * (limits.second - limits.first);
        double c2 = 0.5 * (limits.second + limits.first);

        // sum contributions for all quadrature points
        double sum = 0;

        for (int i = 0; i < 4; ++i)
        {
            // define scaled quadrature point
            double s = c1 * quad_points[i] + c2;

            // compute the value of the kernel function K(X, s)
            double kernel_value = 1;

            for (int j = 0; j < 3; ++j)
            {
                double u = X[j] - data.start_point[j] - s * data.direction[j];
                u /= bandwidth[j];
                kernel_value *= kernel->evaluate(u) / bandwidth[j];
            }
        
            // multiply by quadrature weight and add to sum
            sum += quad_weights[i] * kernel_value;
        }

        // return value of the integral
        return c1 * sum;
    }
    else
    {
        // integration limits are not valid so no score is computed
        return 0;
    }
}
//-----------------------------------------------------------------------------
bool KDEMeshTally::set_integral_limits(const TrackData& data,
                                       const moab::CartVect& X,
                                       std::pair<double, double>& limits)
{
    bool valid_limits = false;

    // set initial integral limits to the full track length (default values)
    limits = std::make_pair(0, data.track_length);

    // check limits against the valid path length interval for each dimension
    for (int i = 0; i < 3; ++i)
    {
        double path_min = limits.first;
        double path_max = limits.second;

        // compute valid path length interval Si = [path_min, path_max]
        if (data.direction[i] > 0)
        {
            path_min = X[i] - data.start_point[i] - bandwidth[i];
            path_min /= data.direction[i];

            path_max = X[i] - data.start_point[i] + bandwidth[i];
            path_max /= data.direction[i];
        }
        else if (data.direction[i] < 0)
        {
            path_min = X[i] - data.start_point[i] + bandwidth[i];
            path_min /= data.direction[i];

            path_max = X[i] - data.start_point[i] - bandwidth[i];
            path_max /= data.direction[i];
        }

        // set lower limit to highest minimum
        if (path_min > limits.first)
        {
            limits.first = path_min;
        }

        // set upper limit to lowest maximum
        if (path_max < limits.second)
        {
            limits.second = path_max;
        }
    }
  
    // limits are only valid if upper limit is greater than lower limit
    if (limits.first < limits.second)
    {
        valid_limits = true;
    }

    return valid_limits;
}
//-----------------------------------------------------------------------------
double KDEMeshTally::subtrack_score(const std::vector<moab::CartVect>& points,
                                    const moab::CartVect& X)
{
    // iterate through the sub-track points
    std::vector<moab::CartVect>::const_iterator i;
    double score = 0;

    for (i = points.begin(); i != points.end(); ++i)
    {
        // Compute the value of the kernel function K(X)
        double kernel_value = 1;

        for (int j = 0; j < 3; ++j)
        {
            double u = (X[j] - (*i)[j]) / bandwidth[j];
            kernel_value *= kernel->evaluate(u) / bandwidth[j];
        }

        // add kernel contribution for sub-track point to sum
        score += kernel_value;
    }

    // normalize by the total number of sub-track points
    score /= points.size();

    return score;
}
//-----------------------------------------------------------------------------
std::vector<moab::CartVect> KDEMeshTally::choose_points(const TrackData& data,
                                                        int p)
{
    // make sure the number of sub-tracks is valid
    assert(p > 0);

    // compute sub-track length, assumed to be equal for all sub-tracks
    double sub_track_length = data.track_length / p;

    // set the starting point to the beginning of the track segment
    moab::CartVect start_point = data.start_point;

    // choose a random position along each sub-track
    std::vector<moab::CartVect> random_points;

    for (int i = 0; i < p; ++i)
    {
        double path_length = rand() * sub_track_length / RAND_MAX;
        
        // add the coordinates of the corresponding point
        random_points.push_back(start_point + path_length * data.direction);

        // shift starting point to the next sub-track
        start_point += sub_track_length * data.direction;
    }
 
    return random_points;
}
//-----------------------------------------------------------------------------
double KDEMeshTally::collision_score(const CollisionData& data,
                                     const moab::CartVect& X)
{
    // compute the value of the kernel function K(X)
    double score = 1;

    for (int i = 0; i < 3; ++i)
    {
        double u = (X[i] - data.collision_point[i]) / bandwidth[i];
        score *= kernel->evaluate(u) / bandwidth[i];
    }

    return score;
}
//-----------------------------------------------------------------------------
