#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/FileOptions.hpp"
#include "Internals.hpp"

#ifdef MOAB_HAVE_CGM
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "CubitVector.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <set>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#ifndef M_PI  /* windows */
# define M_PI 3.14159265358979323846
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define FACETING_TOL_TAG_NAME "FACETING_TOL"

namespace moab {

/* Tolerance Summary

   Facet Tolerance:
   Maximum distance between continuous solid model surface and faceted surface.
     Performance:  increasing tolerance increased performance (fewer triangles)
     Robustness:   should not be affected
     Knowledge:    user must understand how coarser faceting influences accuracy
                   of results

   Overlap Thickness:
   This tolerance is the maximum distance across an overlap. It should be zero
   unless the geometry has overlaps. The overlap thickness is set using the dagmc
   card. Overlaps must be small enough not to significantly affect physics.
     Performance: increasing tolerance decreases performance
     Robustness:  increasing tolerance increases robustness
     Knowledge:   user must have intuition of overlap thickness

   Numerical Precision:
   This tolerance is used for obb.intersect_ray, finding neighborhood of
   adjacent triangle for edge/node intersections, and error in advancing
   geometric position of particle (x' ~= x + d*u). When determining the
   neighborhood of adjacent triangles for edge/node intersections, the facet
   based model is expected to be watertight.
     Performance: increasing tolerance decreases performance (but not very much)
     Robustness:  increasing tolerance increases robustness
     Knowledge:   user should not change this tolerance

*/

  const bool debug    = false; /* controls print statements */
  const bool counting = false; /* controls counts of ray casts and pt_in_vols */

// Empty synonym map for DagMC::parse_metadata()
const std::map<std::string, std::string> DagMC::no_synonyms;

DagMC::DagMC(Interface *mb_impl) {
  moab_instance_created = false;
  // if we arent handed a moab instance create one
  if (NULL == mb_impl) {
    mb_impl = new moab::Core();
    moab_instance_created = true;
  }
  if(moab_instance_created){
    std::cout << "Making new moab instance" << std::endl;
  } else {
    std::cout << "Using old moab instance" << std::endl;
  }

  // other wise take the existing one
  MBI = mb_impl;

  // make new obbtree
  obbTree = new moab::OrientedBoxTreeTool(MBI,MB_OBB_TREE_TAG_NAME,true);

  // This is the correct place to uniquely define default values for the dagmc settings
  overlapThickness = 0; // must be nonnegative
  defaultFacetingTolerance = .001;
  numericalPrecision = .001;
  useCAD = false;

  memset( implComplName, 0, NAME_TAG_SIZE );
  strcpy( implComplName , "impl_complement" );
}

// Destructor
DagMC::~DagMC(){
  // delete the obb tree
  obbTree->~OrientedBoxTreeTool();
  // if we created the moab instance
  // clear it
  if(moab_instance_created) {
    MBI->delete_mesh();
    MBI->~Interface();
 }
}

float DagMC::version(std::string *version_string) {
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

unsigned int DagMC::interface_revision() {
  unsigned int result = 0;
  const char* interface_string = DAGMC_INTERFACE_REVISION;
  if( strlen(interface_string) >= 5 ){
    // start looking for the revision number after "$Rev: "
    result = strtol( interface_string+5, NULL, 10 );
  }
  return result;
}

DagMC::DagMC(Interface *mb_impl)
  : mbImpl(mb_impl), obbTree(mb_impl), impl_compl_handle(0),
    obbTag(0), geomTag(0), idTag(0), nameTag(0), senseTag(0), facetingTolTag(0),
    setOffset(0), facetingTolerance(0.0), have_cgm_geom(false),
    n_pt_in_vol_calls(0), n_ray_fire_calls(0)
{
    // This is the correct place to uniquely define default values for the dagmc settings
  overlapThickness = 0; // must be nonnegative
  defaultFacetingTolerance = .001;
  numericalPrecision = .001;
  useCAD = false;

  memset( implComplName, 0, NAME_TAG_SIZE );
  strcpy( implComplName , "impl_complement" );

}

/* SECTION I: Geometry Initialization and problem setup */

// the standard DAGMC load file method
ErrorCode DagMC::load_file(const char* cfile,
                           const double facet_tolerance) {
  ErrorCode rval;
  
  std::cout << "Requested faceting tolerance: " << facet_tolerance << std::endl;
  
#ifdef MOAB_HAVE_CGM
  // cgm must be initialized so we can check it for CAD data after the load
  InitCGMA::initialize_cgma();
#endif

  facetingTolerance = defaultFacetingTolerance;
    // override default value of facetingTolerance with passed value
  if (facet_tolerance > 0 )
    facetingTolerance = facet_tolerance;

  char facetTolStr[16];

  sprintf(facetTolStr,"%g",facetingTolerance);

  // load options 
  char options[120] = {0};
  char file_ext[4] = "" ; // file extension

  // get the last 4 chars of file .i.e .h5m .sat etc
  memcpy(file_ext, &cfile[strlen(cfile) - 4] ,4);
  // these options only needed if faceting a sat file at DAG runtime
  // not recommended as load_file overloads to ReadCGM load_file
  if ( strstr(file_ext,".sat") != NULL )
    {
      strcat(options,"CGM_ATTRIBS=yes;FACET_DISTANCE_TOLERANCE=");
      strcat(options,facetTolStr);
    }

  EntityHandle file_set;
  rval = MBI->create_meshset( MESHSET_SET, file_set );
  if (MB_SUCCESS != rval)
    return rval;

  rval = MBI->load_file(cfile, &file_set, options, NULL, 0, 0);

  if( MB_UNHANDLED_OPTION == rval ){
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end in '.h5m'
    std::string filename(cfile);
    if( filename.length() < 4 || filename.substr(filename.length()-4) != ".h5m"){
      std::cerr << "DagMC warning: unhandled file loading options." << std::endl;
    }
  }
  else if (MB_SUCCESS != rval) {
    std::cerr << "DagMC Couldn't read file " << cfile << std::endl;
    std::string message;
    if (MB_SUCCESS == MBI->get_last_error(message) && !message.empty())
        std::cerr << "Error message: " << message << std::endl;

    return rval;
  }

#ifdef MOAB_HAVE_CGM
  // check to see if CGM has data; if so, assume it corresponds to the data we loaded in.
  if( GeometryQueryTool::instance()->num_ref_volumes() > 0 ){
    have_cgm_geom = true;
  }
#endif

  return finish_loading();

}

// helper function to load the existing contents of a MOAB instance into DAGMC
ErrorCode DagMC::load_existing_contents( ){

  return finish_loading();
}

// helper function to finish setting up required tags.
ErrorCode DagMC::finish_loading()
{

  ErrorCode rval;

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, NULL, false);

  idTag = get_tag( GLOBAL_ID_TAG_NAME, 1, MB_TAG_DENSE, MB_TYPE_INTEGER );

  geomTag = get_tag( GEOM_DIMENSION_TAG_NAME, 1, MB_TAG_DENSE, MB_TYPE_INTEGER );

  obbTag = get_tag( MB_OBB_TREE_TAG_NAME, 1, MB_TAG_DENSE, MB_TYPE_HANDLE );

  facetingTolTag = get_tag(FACETING_TOL_TAG_NAME, 1, MB_TAG_SPARSE, MB_TYPE_DOUBLE );

    // get sense of surfaces wrt volumes
  senseTag = get_tag( "GEOM_SENSE_2", 2, MB_TAG_SPARSE, MB_TYPE_HANDLE );

  // search for a tag that has the faceting tolerance
  Range tagged_sets;
  double facet_tol_tagvalue = 0;
  bool other_set_tagged = false, root_tagged = false;

  // get list of entity sets that are tagged with faceting tolerance
  // (possibly empty set)
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &facetingTolTag,
                                            NULL, 1, tagged_sets );
  // if NOT empty set
  if (MB_SUCCESS == rval && !tagged_sets.empty()) {
    rval = MBI->tag_get_data( facetingTolTag, &(*tagged_sets.begin()), 1, &facet_tol_tagvalue );
    if (MB_SUCCESS != rval) return rval;
    other_set_tagged = true;
  }
  else if (MB_SUCCESS == rval) {
    // check to see if interface is tagged
    EntityHandle root = 0;
    rval = MBI->tag_get_data( facetingTolTag, &root, 1, &facet_tol_tagvalue );
    if (MB_SUCCESS == rval) root_tagged = true;
    else rval = MB_SUCCESS;
  }

  if ( (root_tagged || other_set_tagged) && facet_tol_tagvalue > 0) {
    facetingTolerance = facet_tol_tagvalue;
  }

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;

  return MB_SUCCESS;

}

// setup the implicit compliment
ErrorCode DagMC::setup_impl_compl()
{
  // If it doesn't already exist, create implicit complement
  // Create data structures for implicit complement
  ErrorCode rval = get_impl_compl();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle." << std::endl;
    return rval;
  }
  return MB_SUCCESS;
}

// gets the entity sets tagged with geomtag 2 and 3
// surfaces and volumes respectively
ErrorCode DagMC::setup_geometry(Range &surfs, Range &vols)
{
  ErrorCode rval;

  const int three = 3;
  const void* const three_val[] = {&three};
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag,
                                            three_val, 1, vols );
  if (MB_SUCCESS != rval)
    return rval;

  const int two = 2;
  const void* const two_val[] = {&two};
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag,
                                            two_val, 1, surfs );
  if (MB_SUCCESS != rval)
    return rval;

  return MB_SUCCESS;
}

// sets up the obb tree for the problem
ErrorCode DagMC::setup_obbs()
{
  ErrorCode rval;
  Range surfs,vols;
  rval = setup_geometry(surfs,vols);
  if(MB_SUCCESS != rval)
    {
      std::cerr << "Failed to setup the geometry" << std::endl;
      return rval;
    }

  // Build OBB trees for everything, but only if we only read geometry
  // Changed to build obb tree if tree does not already exist. -- JK
  if (!have_obb_tree()) {
    rval = build_obbs(surfs, vols);MB_CHK_SET_ERR(rval, "Failed to build obb.");
  }
  return MB_SUCCESS;
}

// setups of the indices for the problem, builds a list of
ErrorCode DagMC::setup_indices()
{
  Range surfs, vols;
  ErrorCode rval = setup_geometry(surfs,vols);

  // If we haven't got the implicit compliment it would be silly to add it
  if(have_impl_compl())
    {
      // build_indices expects the implicit complement to be in vols.
      if( vols.find(impl_compl_handle) == vols.end() )
	     {
	        vols.insert( vols.end(), impl_compl_handle );
	     }
    }

  // build the various index vectors used for efficiency
  rval = build_indices(surfs, vols);MB_CHK_SET_ERR(rval, "Failed to build surface/volume indices.");
  return MB_SUCCESS;
}

// initialise the obb tree
ErrorCode DagMC::init_OBBTree()
{
  ErrorCode rval;
  // implicit compliment
  rval = setup_impl_compl();MB_CHK_SET_ERR(rval, "Failed to setup the implicit compliment");

  // build obbs
  rval = setup_obbs();MB_CHK_SET_ERR(rval, "Failed to setup the OBBs");

  // setup indices
  rval = setup_indices();MB_CHK_SET_ERR(rval, "Failed to setup problem indices");

  return MB_SUCCESS;
}


/* SECTION I (private) */

bool DagMC::have_obb_tree()
{
  Range entities;
  ErrorCode rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &obbTag, 0, 1,
                                                           entities );
  return MB_SUCCESS == rval && !entities.empty();
}

bool DagMC::have_impl_compl()
{
  Range entities;
  const void* const tagdata[] = {implComplName};
  ErrorCode rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &nameTag, tagdata, 1,
                                                           entities );MB_CHK_ERR(rval);
  if (!entities.empty())
    return true;
  else
    return false;
}

ErrorCode DagMC::get_impl_compl()
{
  Range entities;
  const void* const tagdata[] = {implComplName};
  ErrorCode rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &nameTag, tagdata, 1,
                                                           entities );
  // query error
  if (MB_SUCCESS != rval) {
    std::cerr << "Unable to query for implicit complement." << std::endl;
    return rval;
  }

  // found too many
  if (entities.size() > 1) {
    std::cerr << "Too many implicit complement sets." << std::endl;
    return MB_MULTIPLE_ENTITIES_FOUND;
  }

  // found none
  if (entities.empty()) {
    rval= MBI->create_meshset(MESHSET_SET,impl_compl_handle);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to create mesh set for implicit complement." << std::endl;
      return rval;
    }
      // tag this entity with name for implicit complement
    rval = MBI->tag_set_data(nameTag,&impl_compl_handle,1,&implComplName);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to tag new entity as implicit complement." << std::endl;
    }

    return rval;

  } else {
    // found a single implicit complement
    impl_compl_handle = entities.front();
    return MB_SUCCESS;
  }


}

ErrorCode DagMC::build_obbs(Range &surfs, Range &vols)
{
  ErrorCode rval = MB_SUCCESS;

  for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    EntityHandle root;
    Range tris;
    rval = MBI->get_entities_by_dimension( *i, 2, tris );
    if (MB_SUCCESS != rval)
      return rval;
    if (tris.empty())
      std::cerr << "WARNING: Surface " << get_entity_id(*i) << " has no facets." << std::endl;
    rval = obbTree->build( tris, root );
    if (MB_SUCCESS != rval)
      return rval;
    rval = MBI->add_entities( root, &*i, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    rval = MBI->tag_set_data( obbTag, &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;
  }

  for (Range::iterator i = vols.begin(); i != vols.end(); ++i) {
      // get all surfaces in volume
    moab::Range tmp_surfs;
    rval = MBI->get_child_meshsets( *i, tmp_surfs );
    if (MB_SUCCESS != rval)
      return rval;

      // get OBB trees for each surface
    moab::EntityHandle root;
    moab::Range trees;
    trees.clear();
    std::cout << trees[0] << " a " << trees.size() << std::endl;
    //std::vector<moab::EntityHandle> trees_v;
    for (Range::iterator j = tmp_surfs.begin();  j != tmp_surfs.end(); ++j) {
      // skip any surfaces that are non-manifold in the volume
      // because point containment code will get confused by them
      std::cout << TYPE_FROM_HANDLE(*j) << std::endl;
      std::cout << trees[0] << " b " << trees.size() << std::endl;
      int sense = 0;
      std::cout << trees[0] << " c " << trees.size() << std::endl;
      rval = surface_sense( *i, *j, sense );
      std::cout << trees[0] << " d " << trees.size() << std::endl;
      if (MB_SUCCESS != rval) {
        std::cerr << "Surface/Volume sense data missing." << std::endl;
        return rval;
      }
      std::cout << trees[0] << " e " << trees.size() << std::endl;
      if (!sense)
        continue;
      std::cout << trees[0] << " f " << trees.size() << std::endl;
      rval = MBI->tag_get_data( obbTag, &*j, 1, &root );
      std::cout << trees[0] << " g " << trees.size() << std::endl;
      if (MB_SUCCESS != rval || !root) return MB_FAILURE;
      //if(!root) return moab::MB_FAILURE;
      std::cout << trees[0] << " h  " << trees.size() << std::endl;
      //      trees_v.push_back(root);
      trees.insert( root );
      std::cout << trees[0] << " i " << trees.size() << std::endl;
    }
    for ( moab::Range::iterator it = trees.begin() ; it != trees.end() ; ++it ) {
      std::cout << *it << " " << TYPE_FROM_HANDLE(*it) << std::endl;
    }
    //for ( int j = 0 ; j < trees_v.size() ; j++ ) {
      // trees.insert(trees_v[j]);
      //    }
    
    // build OBB tree for volume
    rval = obbTree->join_trees( trees, root );
    std::cout << "rval join " << rval << std::endl;
    if (MB_SUCCESS != rval) return rval;

    rval = MBI->tag_set_data( obbTag, &*i, 1, &root );
    std::cout << "rval set " << rval << std::endl;
    if (MB_SUCCESS != rval) return rval;

  }

  if ( !(have_impl_compl()) )
    {
      std::cerr << "Warning, there is no implicit compliment" << std::endl;
    }
  else
    {
      rval = build_obb_impl_compl(surfs);
      if (MB_SUCCESS != rval) {
	std::cerr << "Unable to build OBB tree for implicit complement." << std::endl;
	return rval;
      }
    }

  return MB_SUCCESS;
}

ErrorCode DagMC::build_obb_impl_compl(Range &surfs)
{
  EntityHandle comp_root, surf_obb_root;
  Range comp_tree;
  ErrorCode rval;
  std::vector<EntityHandle> parent_vols;

  int impl_compl_surf_count = 0;
  double impl_compl_surf_area = 0.0;

    // search through all surfaces
  for (Range::iterator surf_i = surfs.begin(); surf_i != surfs.end(); ++surf_i) {

    parent_vols.clear();
      // get parents of each surface
    rval = MBI->get_parent_meshsets( *surf_i, parent_vols );
    if (MB_SUCCESS != rval)
      return rval;

      // if only one parent, get the OBB root for this surface
    if (parent_vols.size() == 1 ) {

      double a;
      measure_area( *surf_i, a );
      impl_compl_surf_count += 1;
      impl_compl_surf_area  += a;

      rval = MBI->tag_get_data( obbTag, &*surf_i, 1, &surf_obb_root );
      if (MB_SUCCESS != rval)
        return rval;
      if (!surf_obb_root)
        return MB_FAILURE;

        // add obb root to list of obb roots
      comp_tree.insert( surf_obb_root );

      // add this surf to the topology of the implicit complement volume
      rval = MBI->add_parent_child(impl_compl_handle,*surf_i);
      if (MB_SUCCESS != rval)
        return rval;

      // get the surface sense wrt original volume
      EntityHandle sense_data[2] = {0,0};
      rval = MBI->tag_get_data( sense_tag(), &(*surf_i), 1, sense_data );
      if (MB_SUCCESS != rval) return rval;

      // set the surface sense wrt implicit complement volume
      if(0==sense_data[0] && 0==sense_data[1]) return MB_FAILURE;
      if(0==sense_data[0])
        sense_data[0] = impl_compl_handle;
      else if(0==sense_data[1])
        sense_data[1] = impl_compl_handle;
      else
        return MB_FAILURE;
      rval = MBI->tag_set_data( sense_tag(), &(*surf_i), 1, sense_data );
      if (MB_SUCCESS != rval)  return rval;

    }
  }
  // print info about the implicit complement if one was created
  if( impl_compl_surf_count ){
    bool one = (impl_compl_surf_count == 1);
    std::cout << "The implicit complement bounds " << impl_compl_surf_count
              << (one ? " surface" : " surfaces") << std::endl;
    std::cout << "The implicit complement's total surface area = "
              << impl_compl_surf_area << std::endl;
  }

    // join surface trees to make OBB tree for implicit complement
  rval = obbTree->join_trees( comp_tree, comp_root );
  if (MB_SUCCESS != rval)
    return rval;

    // tag the implicit complement handle with the handle for its own OBB tree
  rval = MBI->tag_set_data( obbTag, &impl_compl_handle, 1, &comp_root );
  if (MB_SUCCESS != rval)
    return rval;

  // following ReadCGM, assign dimension and category tags
  int three = 3;
  rval = MBI->tag_set_data(geomTag, &impl_compl_handle, 1, &three );
  if (MB_SUCCESS != rval)
    return rval;

  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                             MB_TAG_SPARSE, MB_TYPE_OPAQUE);
  static const char volume_category[CATEGORY_TAG_SIZE] = "Volume\0";
  rval = MBI->tag_set_data(category_tag, &impl_compl_handle, 1, volume_category );
  if (MB_SUCCESS != rval)
    return rval;

  return MB_SUCCESS;

}

  /* SECTION II: Fundamental Geometry Operations/Queries */
void DagMC::RayHistory::reset() {
  prev_facets.clear();
}

void DagMC::RayHistory::reset_to_last_intersection() {

  if( prev_facets.size() > 1 ){
    prev_facets[0] = prev_facets.back();
    prev_facets.resize( 1 );
  }

}

void DagMC::RayHistory::rollback_last_intersection() {
  if( prev_facets.size() )
    prev_facets.pop_back();
}

ErrorCode DagMC::ray_fire(const EntityHandle vol,
                          const double point[3], const double dir[3],
                          EntityHandle& next_surf, double& next_surf_dist,
                          RayHistory* history, double user_dist_limit,
			  int ray_orientation,
                          OrientedBoxTreeTool::TrvStats* stats ) {

  // take some stats that are independent of nps
  if(counting) {
    ++n_ray_fire_calls;
    if(0==n_ray_fire_calls%10000000) {
      std::cout << "n_ray_fires="   << n_ray_fire_calls
                << " n_pt_in_vols=" << n_pt_in_vol_calls << std::endl;
    }
  }

  if (debug) {
    std::cout << "ray_fire:"
              << " xyz=" << point[0] << " " << point[1] << " " << point[2]
              << " uvw=" << dir[0] << " " << dir[1] << " " << dir[2]
              << " vol_id=" << id_by_index(3, index_by_handle(vol)) << std::endl;
    }

  const double huge_val = std::numeric_limits<double>::max();
  double dist_limit = huge_val;
  if( user_dist_limit > 0 )
    dist_limit = user_dist_limit;

  // don't recreate these every call
  std::vector<double>       &dists       = distList;
  std::vector<EntityHandle> &surfs       = surfList;
  std::vector<EntityHandle> &facets      = facetList;
  dists.clear();
  surfs.clear();
  facets.clear();

  assert(vol - setOffset < rootSets.size());
  const EntityHandle root = rootSets[vol - setOffset];
  ErrorCode rval;

  // check behind the ray origin for intersections
  double neg_ray_len;
  if(0 == overlapThickness) {
    neg_ray_len = -numericalPrecision;
  } else {
    neg_ray_len = -overlapThickness;
  }

  // optionally, limit the nonneg_ray_len with the distance to next collision.
  double nonneg_ray_len = dist_limit;

  // the nonneg_ray_len should not be less than -neg_ray_len, or an overlap
  // may be missed due to optimization within ray_intersect_sets
  if(nonneg_ray_len < -neg_ray_len) nonneg_ray_len = -neg_ray_len;
  assert(0 <= nonneg_ray_len);
  assert(0 >     neg_ray_len);

  // min_tolerance_intersections is passed but not used in this call
  const int min_tolerance_intersections = 0;

  // numericalPrecision is used for box.intersect_ray and find triangles in the
  // neighborhood of edge/node intersections.
  rval = obbTree->ray_intersect_sets( dists, surfs, facets,
                                     root, numericalPrecision,
                                     min_tolerance_intersections,
                                     point, dir, &nonneg_ray_len,
                                     stats, &neg_ray_len, &vol, &senseTag,
                                     &ray_orientation,
                                     history ? &(history->prev_facets) : NULL );
  assert( MB_SUCCESS == rval );
  if(MB_SUCCESS != rval) return rval;

  // if useCAD is true at this point, then we know we can call CGM's ray casting function.
  if (useCAD) {
    rval = CAD_ray_intersect( point, dir, huge_val, dists, surfs, nonneg_ray_len );
    if (MB_SUCCESS != rval) return rval;
  }

  // If no distances are returned, the particle is lost unless the physics limit
  // is being used. If the physics limit is being used, there is no way to tell
  // if the particle is lost. To avoid ambiguity, DO NOT use the distance limit
  // unless you know lost particles do not occur.
  if( dists.empty() ) {
    next_surf = 0;
    if(debug) {
      std::cout << "          next_surf=0 dist=(undef)" << std::endl;
    }
    return MB_SUCCESS;
  }

  // Assume that a (neg, nonneg) pair of RTIs could be returned,
  // however, only one or the other may exist. dists[] may be populated, but
  // intersections are ONLY indicated by nonzero surfs[] and facets[].
  assert(2 == dists.size());
  assert(2 == facets.size());
  assert(0.0 >= dists[0]);
  assert(0.0 <= dists[1]);

  // If both negative and nonnegative RTIs are returned, the negative RTI must
  // closer to the origin.
  if(0!=facets[0] && 0!=facets[1]) {
    assert(-dists[0] <= dists[1]);
  }

  // If an RTI is found at negative distance, perform a PMT to see if the
  // particle is inside an overlap.
  int exit_idx = -1;
  if(0!=facets[0]) {
    // get the next volume
    std::vector<EntityHandle> vols;
    EntityHandle nx_vol;
    rval = MBI->get_parent_meshsets( surfs[0], vols );
    if(MB_SUCCESS != rval) return rval;
    assert(2 == vols.size());
    if(vols.front() == vol) {
      nx_vol = vols.back();
    } else {
      nx_vol = vols.front();
    }
    // Check to see if the point is actually in the next volume.
    // The list of previous facets is used to topologically identify the
    // "on_boundary" result of the PMT. This avoids a test that uses proximity
    // (a tolerance).
    int result;
    rval = point_in_volume( nx_vol, point, result, dir, history );
    if(MB_SUCCESS != rval) return rval;
    if(1==result) exit_idx = 0;

  }

  // if the negative distance is not the exit, try the nonnegative distance
  if(-1==exit_idx && 0!=facets[1]) exit_idx = 1;

  // if the exit index is still unknown, the particle is lost
  if(-1 == exit_idx) {
    next_surf = 0;
    if (debug) {
      std::cout << "next surf hit = 0, dist = (undef)" << std::endl;
    }
    return MB_SUCCESS;
  }

  // return the intersection
  next_surf = surfs[exit_idx];
  next_surf_dist = ( 0>dists[exit_idx] ? 0 : dists[exit_idx]);

  if( history ){
    history->prev_facets.push_back( facets[exit_idx] );
  }

  if (debug) {
    if( 0 > dists[exit_idx] ){
      std::cout << "          OVERLAP track length=" << dists[exit_idx] << std::endl;
    }
    std::cout << "          next_surf = " <<  id_by_index(2, index_by_handle(next_surf))
              << ", dist = " << next_surf_dist << " new_pt=";
    for( int i = 0; i < 3; ++i ){
      std::cout << point[i]+dir[i]*next_surf_dist << " ";
    }
    std::cout << std::endl;
  }

  return MB_SUCCESS;
}

ErrorCode DagMC::point_in_volume(const EntityHandle volume,
                                 const double xyz[3],
                                 int& result,
                                 const double *uvw,
                                 const RayHistory *history) {
  // take some stats that are independent of nps
  if(counting) ++n_pt_in_vol_calls;

  // get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  EntityHandle root = rootSets[volume - setOffset];

  // Don't recreate these every call. These cannot be the same as the ray_fire
  // vectors because both are used simultaneously.
  std::vector<double>       &dists = disList;
  std::vector<EntityHandle> &surfs = surList;
  std::vector<EntityHandle> &facets= facList;
  std::vector<int>          &dirs  = dirList;
  dists.clear();
  surfs.clear();
  facets.clear();
  dirs.clear();

  // if uvw is not given or is full of zeros, use a random direction
  double u = 0, v = 0, w = 0;

  if( uvw ){
    u = uvw[0]; v=uvw[1], w=uvw[2];
  }

  if( u == 0 && v == 0 && w == 0 )
  {
    u = rand();
    v = rand();
    w = rand();
    const double magnitude = sqrt( u*u + v*v + w*w );
    u /= magnitude;
    v /= magnitude;
    w /= magnitude;
  }

  const double ray_direction[] = { u, v, w };

  // if overlaps, ray must be cast to infinity and all RTIs must be returned
  const double   large       = 1e15;
  const double   ray_length  = large;

  // If overlaps occur, the pt is inside if traveling along the ray from the
  // origin, there are ever more exits than entrances. In lieu of implementing
  // that, all intersections to infinity are required if overlaps occur (expensive)
  int min_tolerance_intersections;
  if(0 != overlapThickness) {
    min_tolerance_intersections = -1;
  // only the first intersection is needed if overlaps do not occur (cheap)
  } else {
    min_tolerance_intersections = 1;
  }

  // Get intersection(s) of forward and reverse orientation. Do not return
  // glancing intersections or previous facets.
  ErrorCode rval = obbTree->ray_intersect_sets( dists, surfs, facets, root,
                                               numericalPrecision,
                                               min_tolerance_intersections,
                                               xyz, ray_direction,
                                               &ray_length, NULL, NULL, &volume,
                                               &senseTag, NULL,
                                               history ? &(history->prev_facets) : NULL );
  if(MB_SUCCESS != rval) return rval;

  // determine orientation of all intersections
  // 1 for entering, 0 for leaving, -1 for tangent
  // Tangent intersections are not returned from ray_tri_intersect.
  dirs.resize(dists.size());
  for(unsigned i=0; i<dists.size(); ++i) {
    rval = boundary_case( volume, dirs[i], u, v, w, facets[i], surfs[i] );
    if(MB_SUCCESS != rval) return rval;
  }

  // count all crossings
  if(0 != overlapThickness) {
    int sum = 0;
    for(unsigned i=0; i<dirs.size(); ++i) {
      if     ( 1==dirs[i]) sum+=1; // +1 for entering
      else if( 0==dirs[i]) sum-=1; // -1 for leaving
      else if(-1==dirs[i]) {       //  0 for tangent
        std::cout << "direction==tangent" << std::endl;
        sum+=0;
      } else {
        std::cout << "error: unknown direction" << std::endl;
        return MB_FAILURE;
      }
    }

    // inside/outside depends on the sum
    if(0<sum)                          result = 0; // pt is outside (for all vols)
    else if(0>sum)                     result = 1; // pt is inside  (for all vols)
    else if(impl_compl_handle==volume) result = 1; // pt is inside  (for impl_compl_vol)
    else                               result = 0; // pt is outside (for all other vols)

  // Only use the first crossing
  } else {
      if( dirs.empty() ) {
      result = 0; // pt is outside
    } else {
      int smallest = std::min_element( dists.begin(), dists.end() ) - dists.begin();
      if     ( 1==dirs[smallest] ) result = 0; // pt is outside
      else if( 0==dirs[smallest] ) result = 1; // pt is inside
      else if(-1==dirs[smallest] ) {
        // Should not be here because Plucker ray-triangle test does not
        // return coplanar rays as intersections.
        std::cout << "direction==tangent" << std::endl;
        result = -1;
      } else {
        std::cout << "error: unknown direction" << std::endl;
        return MB_FAILURE;
      }
    }
  }

  if(debug)
    std::cout << "pt_in_vol: result=" << result
              << " xyz=" << xyz[0] << " " << xyz[1] << " " << xyz[2] << " uvw=" << u << " " << v << " " << w
              << " vol_id=" << id_by_index(3, index_by_handle(volume)) << std::endl;

  return MB_SUCCESS;
}

ErrorCode DagMC::test_volume_boundary( const EntityHandle volume, const EntityHandle surface,
                                       const double xyz[3], const double uvw[3], int& result,
                                       const RayHistory* history )
{
  ErrorCode rval;
  int dir;

  if( history && history->prev_facets.size() ){
    // the current facet is already available
    rval = boundary_case( volume, dir, uvw[0], uvw[1], uvw[2], history->prev_facets.back(), surface );
    if (MB_SUCCESS != rval) return rval;
  }
  else{
    // look up nearest facet

    // Get OBB Tree for surface
    assert(volume - setOffset < rootSets.size());
    EntityHandle root = rootSets[volume - setOffset];

    // Get closest triangle on surface
    const CartVect point(xyz);
    CartVect nearest;
    EntityHandle facet_out;
    rval = obbTree->closest_to_location( point.array(), root, nearest.array(), facet_out );
    if (MB_SUCCESS != rval) return rval;

    rval = boundary_case( volume, dir, uvw[0], uvw[1], uvw[2], facet_out, surface );
    if (MB_SUCCESS != rval) return rval;

  }

  result = dir;

  return MB_SUCCESS;

}

// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMC::point_in_volume_slow( EntityHandle volume, const double xyz[3], int& result )
{
  ErrorCode rval;
  Range faces;
  std::vector<EntityHandle> surfs;
  std::vector<int> senses;
  double sum = 0.0;
  const CartVect point(xyz);

  rval = MBI->get_child_meshsets( volume, surfs );
  if (MB_SUCCESS != rval)
    return rval;

  senses.resize( surfs.size() );
  rval = surface_sense( volume, surfs.size(), &surfs[0], &senses[0] );
  if (MB_SUCCESS != rval)
    return rval;

  for (unsigned i = 0; i < surfs.size(); ++i) {
    if (!senses[i])  // skip non-manifold surfaces
      continue;

    double surf_area = 0.0, face_area;
    faces.clear();
    rval = MBI->get_entities_by_dimension( surfs[i], 2, faces );
    if (MB_SUCCESS != rval)
      return rval;

    for (Range::iterator j = faces.begin(); j != faces.end(); ++j) {
      rval = poly_solid_angle( *j, point, face_area );
      if (MB_SUCCESS != rval)
        return rval;

      surf_area += face_area;
    }

    sum += senses[i] * surf_area;
  }

  result = fabs(sum) > 2.0*M_PI;
  return MB_SUCCESS;
}



// detemine distance to nearest surface
ErrorCode DagMC::closest_to_location( EntityHandle volume, const double coords[3], double& result)
{
    // Get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  EntityHandle root = rootSets[volume - setOffset];

    // Get closest triangles in volume
  const CartVect point(coords);
  CartVect nearest;
  EntityHandle facet_out;
  ErrorCode rval = obbTree->closest_to_location( point.array(), root, nearest.array(), facet_out );
  if (MB_SUCCESS != rval) return rval;

  // calculate distance between point and nearest facet
  result = (point-nearest).length();

  return MB_SUCCESS;

}

// calculate volume of polyhedron
ErrorCode DagMC::measure_volume( EntityHandle volume, double& result )
{
  ErrorCode rval;
  std::vector<EntityHandle> surfaces;
  result = 0.0;

   // don't try to calculate volume of implicit complement
  if (volume == impl_compl_handle) {
    result = 1.0;
    return MB_SUCCESS;
  }

    // get surfaces from volume
  rval = MBI->get_child_meshsets( volume, surfaces );
  if (MB_SUCCESS != rval) return rval;

    // get surface senses
  std::vector<int> senses( surfaces.size() );
  rval = surface_sense( volume, surfaces.size(), &surfaces[0], &senses[0] );
  if (MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }

  for (unsigned i = 0; i < surfaces.size(); ++i) {
      // skip non-manifold surfaces
    if (!senses[i])
      continue;

      // get triangles in surface
    Range triangles;
    rval = MBI->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (MB_SUCCESS != rval)
      return rval;
    if (!triangles.all_of_type(MBTRI)) {
      std::cout << "WARNING: Surface " << get_entity_id(surfaces[i])
                << " contains non-triangle elements. Volume calculation may be incorrect."
                << std::endl;
      triangles.clear();
      rval = MBI->get_entities_by_type( surfaces[i], MBTRI, triangles );
      if (MB_SUCCESS != rval) return rval;
    }

      // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const EntityHandle *conn;
    int len;
    CartVect coords[3];
    for (Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = MBI->get_connectivity( *j, conn, len, true );
      if (MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = MBI->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;

      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }

  result /= 6.0;
  return MB_SUCCESS;
}

// sum area of elements in surface
ErrorCode DagMC::measure_area( EntityHandle surface, double& result )
{
    // get triangles in surface
  Range triangles;
  ErrorCode rval = MBI->get_entities_by_dimension( surface, 2, triangles );
  if (MB_SUCCESS != rval)
    return rval;
  if (!triangles.all_of_type(MBTRI)) {
    std::cout << "WARNING: Surface " << get_entity_id(surface)
              << " contains non-triangle elements. Area calculation may be incorrect."
              << std::endl;
    triangles.clear();
    rval = MBI->get_entities_by_type( surface, MBTRI, triangles );
    if (MB_SUCCESS != rval) return rval;
  }

    // calculate sum of area of triangles
  result = 0.0;
  const EntityHandle *conn;
  int len;
  CartVect coords[3];
  for (Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
    rval = MBI->get_connectivity( *j, conn, len, true );
    if (MB_SUCCESS != rval) return rval;
    assert(3 == len);
    rval = MBI->get_coords( conn, 3, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    coords[0] = coords[1] * coords[2];
    result += coords[0].length();
  }
  result *= 0.5;
  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense( EntityHandle volume,
                           int num_surfaces,
                           const EntityHandle* surfaces,
                           int* senses_out )
{

  /* The sense tags do not reference the implicit complement handle.
     All surfaces that interact with the implicit complement should have
     a null handle in the direction of the implicit complement. */
  //if (volume == impl_compl_handle)
  //  volume = (EntityHandle) 0;

  std::vector<EntityHandle> surf_volumes( 2*num_surfaces );
  ErrorCode rval = MBI->tag_get_data( sense_tag(), surfaces, num_surfaces, &surf_volumes[0] );
  if (MB_SUCCESS != rval)  return rval;

  const EntityHandle* end = surfaces + num_surfaces;
  std::vector<EntityHandle>::const_iterator surf_vols = surf_volumes.begin();
  while (surfaces != end) {
    EntityHandle forward = *surf_vols; ++surf_vols;
    EntityHandle reverse = *surf_vols; ++surf_vols;
    if (volume == forward)
      *senses_out = (volume != reverse); // zero if both, otherwise 1
    else if (volume == reverse)
      *senses_out = -1;
    else
      return MB_ENTITY_NOT_FOUND;

    ++surfaces;
    ++senses_out;
  }

  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense( EntityHandle volume,
                                  EntityHandle surface,
                                  int& sense_out )
{
  /* The sense tags do not reference the implicit complement handle.
     All surfaces that interact with the implicit complement should have
     a null handle in the direction of the implicit complement. */
  //if (volume == impl_compl_handle)
  //  volume = (EntityHandle) 0;

    // get sense of surfaces wrt volumes
  EntityHandle surf_volumes[2];
  ErrorCode rval = MBI->tag_get_data( sense_tag(), &surface, 1, surf_volumes );
  if (MB_SUCCESS != rval)  return rval;

  if (surf_volumes[0] == volume)
    sense_out = (surf_volumes[1] != volume); // zero if both, otherwise 1
  else if (surf_volumes[1] == volume)
    sense_out = -1;
  else
    return MB_ENTITY_NOT_FOUND;

  return MB_SUCCESS;
}

ErrorCode DagMC::get_angle(EntityHandle surf, const double in_pt[3], double angle[3], const RayHistory* history )
{
  EntityHandle root = rootSets[surf - setOffset];
  ErrorCode rval;

  std::vector<EntityHandle> facets;

  // if no history or history empty, use nearby facets
  if( !history || (history->prev_facets.size() == 0) ){
    rval = obbTree->closest_to_location( in_pt, root, numericalPrecision, facets );
    assert(MB_SUCCESS == rval);
    if (MB_SUCCESS != rval) return rval;
  }
  // otherwise use most recent facet in history
  else{
    facets.push_back( history->prev_facets.back() );
  }

  CartVect coords[3], normal(0.0);
  const EntityHandle *conn;
  int len;
  for (unsigned i = 0; i < facets.size(); ++i) {
    rval = MBI->get_connectivity( facets[i], conn, len );
    assert( MB_SUCCESS == rval );
    assert( 3 == len );

    rval = MBI->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == rval);

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal += coords[1] * coords[2];
  }

  normal.normalize();
  normal.get( angle );

  return MB_SUCCESS;
}

ErrorCode DagMC::next_vol( EntityHandle surface, EntityHandle old_volume,
                           EntityHandle& new_volume )
{
  std::vector<EntityHandle> parents;
  ErrorCode rval = MBI->get_parent_meshsets( surface, parents );

  if (MB_SUCCESS == rval) {
    if (parents.size() != 2)
      rval = MB_FAILURE;
    else if (parents.front() == old_volume)
      new_volume = parents.back();
    else if( parents.back() == old_volume )
      new_volume = parents.front();
    else
      rval = MB_FAILURE;
  }

  if( rval != MB_SUCCESS ){
    std::cerr << "DAGMC: mesh error in next_vol for surf " << get_entity_id(surface) << std::endl;
  }

  return rval;

}

/* SECTION II (private) */

ErrorCode DagMC::CAD_ray_intersect(
#if defined(CGM) && defined(HAVE_CGM_FIRE_RAY)
    const double *point,
    const double *dir,
    const double huge_val,
    std::vector<double> &distances,
    std::vector<EntityHandle> &surfaces,
    double &len
#else
    const double *,
    const double *,
    const double ,
    std::vector<double> &,
    std::vector<EntityHandle> &,
    double &
#endif
)
{
#ifdef MOAB_HAVE_CGM
#ifndef HAVE_CGM_FIRE_RAY
  return MB_NOT_IMPLEMENTED;
#else
  std::vector<double>::iterator dit = distances.begin();
  std::vector<EntityHandle>::iterator sit = surfaces.begin();
  static DLIList<double> ray_params;

  for (; dit != distances.end(); ++dit, ++sit) {
      // get the RefFace
    RefEntity *this_face = geomEntities[*sit - setOffset];
      // get the ray distance to this face
    ray_params.clean_out();
    int result = GeometryQueryTool::instance()->
      fire_ray(dynamic_cast<RefFace*>(this_face), CubitVector(point),
               CubitVector(dir), ray_params);
    assert(CUBIT_SUCCESS == result);
    if(CUBIT_SUCCESS != result) return MB_FAILURE;
    if (ray_params.size() != 0) {
      ray_params.reset();
      *dit = ray_params.get();
    }
    else *dit = huge_val;
  }

    // now bubble sort list
  bool done = false;
  while (!done) {
    dit = distances.begin();
    sit = surfaces.begin();
    done = true;
    for (; dit != distances.end(); ++dit, ++sit) {
      if (dit+1 != distances.end() && *dit > *(dit+1)) {
        double tmp_dist = *dit;
        *dit = *(dit+1);
        *(dit+1) = tmp_dist;
        EntityHandle tmp_hand = *sit;
        *sit = *(sit+1);
        *(sit+1) = tmp_hand;
        done = false;
      }
    }
  }

  if (!distances.empty()) len = distances[0];

  return MB_SUCCESS;
#endif
#else
  return MB_FAILURE;
#endif
}

// If point is on boundary, then this function is called to
// discriminate cases in which the ray is entering or leaving.
// result= 1 -> inside volume or entering volume
// result= 0 -> outside volume or leaving volume
// result=-1 -> on boundary with null or tangent uvw
ErrorCode DagMC::boundary_case( EntityHandle volume, int& result,
                                double u, double v, double w,
                                EntityHandle facet,
                                EntityHandle surface)
{
  ErrorCode rval;

  // test to see if uvw is provided
  if ( u <= 1.0 && v <= 1.0 && w <= 1.0 ) {

    const CartVect ray_vector(u, v, w);
    CartVect coords[3], normal(0.0);
    const EntityHandle *conn;
    int len, sense_out;

    rval = MBI->get_connectivity( facet, conn, len );
    assert( MB_SUCCESS == rval );
    if(MB_SUCCESS != rval) return rval;
    assert( 3 == len );

    rval = MBI->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == rval);
    if(MB_SUCCESS != rval) return rval;

    rval = surface_sense( volume, surface, sense_out );
    assert( MB_SUCCESS == rval);
    if(MB_SUCCESS != rval) return rval;

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal = sense_out * (coords[1] * coords[2]);

    double sense = ray_vector % normal;

    if ( sense < 0.0 ) {
      result = 1;     // inside or entering
    } else  if ( sense > 0.0 ) {
      result = 0;     // outside or leaving
    } else  if ( sense == 0.0 ) {
      result = -1;    // tangent, therefore on boundary
    } else {
      result = -1;    // failure
      return MB_FAILURE;
    }

  // if uvw not provided, return on_boundary.
  } else {
    result = -1;      // on boundary
    return MB_SUCCESS;

  }

  return MB_SUCCESS;
}

// point_in_volume_slow, including poly_solid_angle helper subroutine
// are adapted from "Point in Polyhedron Testing Using Spherical Polygons", Paulo Cezar
// Pinto Carvalho and Paulo Roma Cavalcanti, _Graphics Gems V_, pg. 42.  Original algorithm
// was described in "An Efficient Point In Polyhedron Algorithm", Jeff Lane, Bob Magedson,
// and Mike Rarick, _Computer Vision, Graphics, and Image Processing 26_, pg. 118-225, 1984.

// helper function for point_in_volume_slow.  calculate area of a polygon
// projected into a unit-sphere space
   ErrorCode DagMC::poly_solid_angle( EntityHandle face, const CartVect& point, double& area )
{
  ErrorCode rval;

    // Get connectivity
  const EntityHandle* conn;
  int len;
  rval = MBI->get_connectivity( face, conn, len, true );
  if (MB_SUCCESS != rval)
    return rval;

    // Allocate space to store vertices
  CartVect coords_static[4];
  std::vector<CartVect> coords_dynamic;
  CartVect* coords = coords_static;
  if ((unsigned)len > (sizeof(coords_static)/sizeof(coords_static[0]))) {
    coords_dynamic.resize(len);
    coords = &coords_dynamic[0];
  }

    // get coordinates
  rval = MBI->get_coords( conn, len, coords->array() );
  if (MB_SUCCESS != rval)
    return rval;

    // calculate normal
  CartVect norm(0.0), v1, v0 = coords[1] - coords[0];
  for (int i = 2; i < len; ++i) {
    v1 = coords[i] - coords[0];
    norm += v0 * v1;
    v0 = v1;
  }

    // calculate area
  double s, ang;
  area = 0.0;
  CartVect r, n1, n2, b, a = coords[len-1] - coords[0];
  for (int i = 0; i < len; ++i) {
    r = coords[i] - point;
    b = coords[(i+1)%len] - coords[i];
    n1 = a * r; // = norm1 (magnitude is important)
    n2 = r * b; // = norm2 (magnitude is important)
    s = (n1 % n2) / (n1.length() * n2.length()); // = cos(angle between norm1,norm2)
    ang = s <= -1.0 ? M_PI : s >= 1.0 ? 0.0 : acos(s); // = acos(s)
    s = (b * a) % norm; // =orientation of triangle wrt point
    area += s > 0.0 ? M_PI - ang : M_PI + ang;
    a = -b;
  }

  area -= M_PI * (len - 2);
  if ((norm % r) > 0)
    area = -area;
  return MB_SUCCESS;
}

/* SECTION III */

EntityHandle DagMC::entity_by_id( int dimension, int id )
{
  assert(0 <= dimension && 3 >= dimension);
  const Tag tags[] = { idTag, geomTag };
  const void* const vals[] = { &id, &dimension };
  ErrorCode rval;

  Range results;
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, results );

  if ( MB_SUCCESS != rval )
      return 0;

  if ( results.empty() ){
    // old versions of dagmc did not set tags correctly on the implicit complement 'volume',
    // causing it to not be found by the call above.  This check allows this function to work
    // correctly, even on reloaded files from older versions.
    if( dimension == 3 && get_entity_id(impl_compl_handle) == id )
      return impl_compl_handle;
    else
      return 0;
  }

  return results.front();
}

int DagMC::id_by_index( int dimension, int index )
{
  EntityHandle h = entity_by_index( dimension, index );
  if (!h)
    return 0;

  int result = 0;
  MBI->tag_get_data( idTag, &h, 1, &result );
  return result;
}

int DagMC::get_entity_id(EntityHandle this_ent)
{
  int id = 0;
  ErrorCode result = MBI->tag_get_data(idTag, &this_ent, 1, &id);
  if (MB_TAG_NOT_FOUND == result)
    id = MBI->id_from_handle(this_ent);

  return id;
}

ErrorCode DagMC::build_indices(Range &surfs, Range &vols)
{
  ErrorCode rval = MB_SUCCESS;

    // surf/vol offsets are just first handles
  setOffset = std::min( *surfs.begin(), *vols.begin() );

    // max
  EntityHandle tmp_offset = std::max( surfs.back(), vols.back() );

    // set size
  rootSets.resize(tmp_offset - setOffset + 1);
  entIndices.resize(rootSets.size());

    // store surf/vol handles lists (surf/vol by index) and
    // index by handle lists
  surf_handles().resize( surfs.size() + 1 );
  std::vector<EntityHandle>::iterator iter = surf_handles().begin();
  *(iter++) = 0;
  std::copy( surfs.begin(), surfs.end(), iter );
  int idx = 1;
  for (Range::iterator rit = surfs.begin(); rit != surfs.end(); ++rit)
    entIndices[*rit-setOffset] = idx++;

  vol_handles().resize( vols.size() + 1 );
  iter = vol_handles().begin();
  *(iter++) = 0;
  std::copy( vols.begin(), vols.end(), iter );

  idx = 1;
  int max_id = -1;
  for (Range::iterator rit = vols.begin(); rit != vols.end(); ++rit)    {
    entIndices[*rit-setOffset] = idx++;

    if( *rit != impl_compl_handle ){
      int result=0;
      MBI->tag_get_data( idTag, &*rit, 1, &result );
      max_id = std::max( max_id, result );
    }
  }
    // assign ID to implicit complement
    // for consistency with earlier versions of DagMC, make sure it always has the highest ID
  max_id++;
  MBI->tag_set_data(idTag, &impl_compl_handle, 1, &max_id);

#ifdef MOAB_HAVE_CGM
  if ( have_cgm_geom ) {
    // TODO: this block should only execute if the user has explicitly requested useCAD for ray firing.
    // however, this function currently executes before we know if useCAD will be specified, so do it every time.

    geomEntities.resize(rootSets.size());
      // get geometry entities by id and cache in this vector
    std::vector<int> ids;
    ids.resize(surfs.size());
    rval = MBI->tag_get_data(id_tag(), surfs, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    int i = 0;
    Range::iterator rit = surfs.begin();
    for (; rit != surfs.end(); ++rit, i++) {
      RefEntity *this_surf = GeometryQueryTool::instance()->
        get_ref_face(ids[i]);
      assert(NULL != this_surf);
      geomEntities[*rit - setOffset] = this_surf;
    }
    ids.resize(vols.size());
    rval = MBI->tag_get_data(id_tag(), vols, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    i = 0;
    rit = vols.begin();
    for (; rit != vols.end(); ++rit, i++) {
      if( is_implicit_complement( *rit ) ) continue;
      RefEntity *this_vol = GeometryQueryTool::instance()->
        get_ref_volume(ids[i]);
      assert(NULL != this_vol);
      geomEntities[*rit - setOffset] = this_vol;
    }
  }
#endif

    // get group handles
  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                               MB_TAG_SPARSE, MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category+CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};
  Range groups;
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &category_tag,
                                           group_val, 1, groups);
  if (MB_SUCCESS != rval)
    return rval;
  group_handles().resize(groups.size()+1);
  group_handles()[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles()[1]);

    // populate root sets vector
  std::vector<EntityHandle> rsets;
  rsets.resize(surfs.size());
  rval = MBI->tag_get_data(obb_tag(), surfs, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  Range::iterator rit;
  int i;
  for (i = 0, rit = surfs.begin(); rit != surfs.end(); ++rit, i++)
    rootSets[*rit-setOffset] = rsets[i];

  rsets.resize(vols.size());
  rval = MBI->tag_get_data(obb_tag(), vols, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  for (i = 0, rit = vols.begin(); rit != vols.end(); ++rit, i++)
    rootSets[*rit-setOffset] = rsets[i];

  return MB_SUCCESS;
}



/* SECTION IV */

void DagMC::set_overlap_thickness( double new_thickness ){

  if (new_thickness < 0 || new_thickness > 100) {
    std::cerr << "Invalid overlap_thickness = " << new_thickness << std::endl;
  }
  else{
    overlapThickness = new_thickness;
  }
  std::cout << "Set overlap thickness = " << overlapThickness << std::endl;

}

void DagMC::set_numerical_precision( double new_precision ){

  if ( new_precision <= 0 || new_precision > 1) {
    std::cerr << "Invalid numerical_precision = " << numericalPrecision << std::endl;
  }
  else{
    numericalPrecision = new_precision;
  }

  std::cout << "Set numerical precision = " << numericalPrecision << std::endl;

}

void DagMC::set_use_CAD( bool use_cad ){
  useCAD = use_cad;
  if( useCAD ){
    if( !have_cgm_geom ){
      std::cerr << "Warning: CAD-based ray tracing not available, because CGM has no data." << std::endl;
      std::cerr << "         your input file was probably not a CAD format." << std::endl;
      useCAD = false;
    }

#ifndef HAVE_CGM_FIRE_RAY
    {
      std::cerr << "Warning: use_cad = 1 not supported with this build of CGM/DagMC." << std:: endl;
      std::cerr << "         Required ray-fire query not available. (Cubit-based CGM?)" <<  std::endl;
      useCAD = false;
    }
#endif
  }

  std::cout << "Turned " << (useCAD?"ON":"OFF") << " ray firing on full CAD model." << std::endl;

}

ErrorCode DagMC::write_mesh(const char* ffile,
                            const int flen)
{
  ErrorCode rval;

    // write out a mesh file if requested
  if (ffile && 0 < flen) {
    rval = MBI->write_mesh(ffile);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return rval;
    }
  }

  return MB_SUCCESS;
}

/* SECTION V: Metadata handling */

ErrorCode DagMC::get_group_name( EntityHandle group_set, std::string& name )
{
  ErrorCode rval;
  const void* v = NULL;
  int ignored;
  rval = MBI->tag_get_by_ptr(name_tag(), &group_set, 1, &v, &ignored);
  if( MB_SUCCESS != rval ) return rval;
  name = static_cast<const char*>(v);
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_group_name( EntityHandle group_set, prop_map& result, const char *delimiters )
{
  ErrorCode rval;
  std::string group_name;
  rval = get_group_name( group_set, group_name );
  if( rval != MB_SUCCESS ) return rval;

  std::vector< std::string > group_tokens;
  tokenize( group_name, group_tokens, delimiters );

  // iterate over all the keyword positions
  // keywords are even indices, their values (optional) are odd indices
  for( unsigned int i = 0; i < group_tokens.size(); i += 2 ){
    std::string groupkey = group_tokens[i];
    std::string groupval;
    if( i < group_tokens.size() - 1 )
      groupval = group_tokens[i+1];
    result[groupkey] = groupval;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::detect_available_props( std::vector<std::string>& keywords_list, const char *delimiters )
{
  ErrorCode rval;
  std::set< std::string > keywords;
  for( std::vector<EntityHandle>::const_iterator grp=group_handles().begin();
       grp != group_handles().end(); ++grp )
  {
    std::map< std::string, std::string > properties;
    rval = parse_group_name( *grp, properties, delimiters );
    if( rval == MB_TAG_NOT_FOUND ) continue;
    else if( rval != MB_SUCCESS ) return rval;

    for( prop_map::iterator i = properties.begin();
         i != properties.end(); ++i )
    {
      keywords.insert( (*i).first );
    }
  }
  keywords_list.assign( keywords.begin(), keywords.end() );
  return MB_SUCCESS;
}

ErrorCode DagMC::append_packed_string( Tag tag, EntityHandle eh,
                                       std::string& new_string )
{
    // When properties have multiple values, the values are tagged in a single character array
    // with the different values separated by null characters
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = MBI->tag_get_by_ptr( tag, &eh, 1, &p, &len );
  if( rval == MB_TAG_NOT_FOUND ){
    // This is the first entry, and can be set directly
    p = new_string.c_str();
    return MBI->tag_clear_data( tag, &eh, 1, p, new_string.length()+1);
  }
  else if( rval != MB_SUCCESS ) return rval;
  else{
    str = static_cast<const char*>(p);
  }

  // append a new value for the property to the existing property string
  unsigned int tail_len = new_string.length() + 1;
  char* new_packed_string = new char[ len + tail_len ];
  memcpy( new_packed_string, str, len );
  memcpy( new_packed_string + len, new_string.c_str(), tail_len );

  int new_len = len + tail_len;
  p = new_packed_string;
  rval = MBI->tag_set_by_ptr( tag, &eh, 1, &p, &new_len );
  delete[] new_packed_string;
  return rval;
}

ErrorCode DagMC::unpack_packed_string( Tag tag, EntityHandle eh,
                                       std::vector< std::string >& values )
{
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = MBI->tag_get_by_ptr( tag, &eh, 1, &p, &len );
  if( rval != MB_SUCCESS ) return rval;
  str = static_cast<const char*>(p);
  int idx = 0;
  while( idx < len ){
    std::string item(str + idx);
    values.push_back( item );
    idx += item.length() + 1;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_properties( const std::vector<std::string>& keywords,
                                   const std::map<std::string, std::string>& keyword_synonyms,
				   const char *delimiters)
{
  ErrorCode rval;

  // master keyword map, mapping user-set words in cubit to canonical property names
  std::map< std::string, std::string > keyword_map( keyword_synonyms );

  for( std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i )
  {
    keyword_map[*i] = *i;
  }

  // the set of all canonical property names
  std::set< std::string > prop_names;
  for( prop_map::iterator i = keyword_map.begin();
       i != keyword_map.end(); ++i )
  {
    prop_names.insert((*i).second);
  }

  // set up DagMC's property tags based on what's been requested
  for( std::set<std::string>::iterator i = prop_names.begin();
       i != prop_names.end(); ++i )
  {
    std::string tagname("DAGMCPROP_");
    tagname += (*i);

    Tag new_tag;
    rval = MBI->tag_get_handle( tagname.c_str(), 0, MB_TYPE_OPAQUE, new_tag,
                                MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_CREAT );
    if( MB_SUCCESS != rval ) return rval;
    property_tagmap[(*i)] = new_tag;
  }

  // now that the keywords and tags are ready, iterate over all the actual geometry groups
  for( std::vector<EntityHandle>::iterator grp=group_handles().begin();
       grp != group_handles().end(); ++grp )
  {

    prop_map properties;
    rval = parse_group_name( *grp, properties, delimiters );
    if( rval == MB_TAG_NOT_FOUND ) continue;
    else if( rval != MB_SUCCESS ) return rval;

    Range grp_sets;
    rval = MBI->get_entities_by_type( *grp, MBENTITYSET, grp_sets);
    if( MB_SUCCESS != rval ) return rval;
    if( grp_sets.size() == 0 ) continue;

    for( prop_map::iterator i = properties.begin();
         i != properties.end(); ++i )
    {
      std::string groupkey = (*i).first;
      std::string groupval = (*i).second;

      if( property_tagmap.find( groupkey ) != property_tagmap.end() ){
        Tag proptag = property_tagmap[groupkey];
        const unsigned int groupsize = grp_sets.size();
        for( unsigned int j = 0; j < groupsize; ++j){
            rval = append_packed_string( proptag, grp_sets[j], groupval );
        }
      }
    }
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_value( EntityHandle eh, const std::string& prop, std::string& value )
{
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return MB_TAG_NOT_FOUND;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr( proptag, &eh, 1, &data, &ignored );
  if( rval != MB_SUCCESS ) return rval;
  value = static_cast<const char*>(data);
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_values( EntityHandle eh, const std::string& prop,
                              std::vector< std::string >& values )
{

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  return unpack_packed_string( proptag, eh, values );

}

bool DagMC::has_prop( EntityHandle eh, const std::string& prop )
{
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return false;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr( proptag, &eh, 1, &data, &ignored );
  return ( rval == MB_SUCCESS );

}


ErrorCode DagMC::get_all_prop_values( const std::string& prop, std::vector<std::string>& return_list )
{
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &proptag, NULL, 1, all_ents );
  if( MB_SUCCESS != rval ) return rval;

  std::set<std::string> unique_values;
  for( Range::iterator i = all_ents.begin(); i!= all_ents.end(); ++i){
    std::vector<std::string> values;
    rval = prop_values( *i, prop, values );
    if( MB_SUCCESS != rval ) return rval;
    unique_values.insert( values.begin(), values.end() );
  }

  return_list.assign( unique_values.begin(), unique_values.end() );
  return MB_SUCCESS;
}

ErrorCode DagMC::entities_by_property( const std::string& prop, std::vector<EntityHandle>& return_list,
                                       int dimension, const std::string* value )
{
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  // Note that we cannot specify values for proptag here-- the passed value,
  // if it exists, may be only a subset of the packed string representation
  // of this tag.
  Tag tags[2] = {proptag, geomTag };
  void* vals[2] = {NULL, (dimension!=0) ? &dimension : NULL };
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, all_ents );
  if( MB_SUCCESS != rval ) return rval;

  std::set<EntityHandle> handles;
  for( Range::iterator i = all_ents.begin(); i!=all_ents.end(); ++i){
    std::vector<std::string> values;
    rval = prop_values( *i, prop, values );
    if( MB_SUCCESS != rval ) return rval;
    if( value ){
      if( std::find(values.begin(), values.end(), *value) != values.end() ){
        handles.insert(*i);
      }
    }
    else{
      handles.insert(*i);
    }
  }

  return_list.assign( handles.begin(), handles.end() );
  return MB_SUCCESS;
}

bool DagMC::is_implicit_complement(EntityHandle volume)
{
  return volume == impl_compl_handle;
}

void DagMC::tokenize( const std::string& str,
                      std::vector<std::string>& tokens,
                      const char* delimiters ) const
{
  std::string::size_type last = str.find_first_not_of( delimiters, 0 );
  std::string::size_type pos  = str.find_first_of( delimiters, last );
  if ( std::string::npos == pos )
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back( str.substr( last, pos - last ) );
      last = str.find_first_not_of( delimiters, pos );
      pos  = str.find_first_of( delimiters, last );
      if(std::string::npos == pos)
        pos = str.size();
    }
}

Tag DagMC::get_tag( const char* name, int size, TagType store,
                    DataType type, const void* def_value,
                    bool create_if_missing)
{
  Tag retval = 0;
  unsigned flags = store|MB_TAG_CREAT;
  // NOTE: this function seems to be broken in that create_if_missing has
  // the opposite meaning from what its name implies.  However, changing the
  // behavior causes tests to fail, so I'm leaving the existing behavior
  // in place.  -- j.kraftcheck.
  if (!create_if_missing)
    flags |= MB_TAG_EXCL;
  ErrorCode result = MBI->tag_get_handle(name, size, type, retval, flags, def_value);
  if (create_if_missing && MB_SUCCESS != result)
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;

  return retval;
}


ErrorCode DagMC::getobb(EntityHandle volume, double minPt[3],
                          double maxPt[3])
{
  double center[3], axis1[3], axis2[3], axis3[3];

    // get center point and vectors to OBB faces
  ErrorCode rval = getobb(volume, center, axis1, axis2, axis3);
  if (MB_SUCCESS != rval)
    return rval;

    // compute min and max vertices
  for (int i=0; i<3; i++)
  {
    double sum = fabs(axis1[i]) + fabs(axis2[i]) + fabs(axis3[i]);
    minPt[i] = center[i] - sum;
    maxPt[i] = center[i] + sum;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::getobb(EntityHandle volume, double center[3], double axis1[3],
                          double axis2[3], double axis3[3])
{
    //find EntityHandle node_set for use in box
  EntityHandle root = rootSets[volume - setOffset];

    // call box to get center and vectors to faces
  return obbTree->box(root, center, axis1, axis2, axis3);

}


} // namespace moab
