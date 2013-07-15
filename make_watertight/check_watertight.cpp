// ********************************************************************
// Brandon Smith
// August 2009

// This is a function to test DagMC-style mesh for watertightness. For
// now this will be a stand-alone code that uses MOAB. For volumes to 
// be watertight, the facet edges of each surface must be matched 
// one-to-one. By default this checks for topological watertightness.
// To instead use a geometrical tolerance between two vertices that are
// considered the same, pass in a tolerance.
//
// input:  h5m file name, tolerance(optional)
// output: list of unmatched facet edges and their parent surfaces, by volume.

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

/* Assume no MBEDGEs exist for fastest skinning
   For each volume:
     For each child surface:
       skin surface tris
       enter data into coords_and_id
     }
     match edges
   }
   Each surface is skinned twice, but the logic is simple and the memory handling 
   is easy.
*/   

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

MBInterface *MBI();

// struct to hold coordinates of skin edge, it's surface id, and a matched flag
struct coords_and_id {
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  int  surf_id;
  bool matched;
  MBEntityHandle vert1;
  MBEntityHandle vert2;
};

/* qsort struct comparision function */
int compare_by_handle(const void *a, const void *b)
{
  struct coords_and_id *ia = (struct coords_and_id *)a;
  struct coords_and_id *ib = (struct coords_and_id *)b;
  if(ia->vert1 == ib->vert1) 
  {
    return (int)(ia->vert2 - ib->vert2);
  } 
  else 
  {
    return (int)(ia->vert1 - ib->vert1);
  }
  /* float comparison: returns negative if b > a 
     and positive if a > b. We multiplied result by 100.0
     to preserve decimal fraction */
} 

/* qsort struct comparision function */
// This is tricky because doubles always get rounded down to ints.
int compare_by_coords(const void *a, const void *b)
{
  struct coords_and_id *ia = (struct coords_and_id *)a;
  struct coords_and_id *ib = (struct coords_and_id *)b;
  if(ia->x1 == ib->x1) {
    if(ia->y1 == ib->y1) {
      if(ia->z1 == ib->z1) {
        if(ia->x2 == ib->x2) {
          if(ia->y2 == ib->y2) {
            if(ia->z2 == ib->z2) {
              return ia->surf_id - ib->surf_id;
            } else {
              return (ia->z2 > ib->z2) - (ia->z2 < ib->z2);
            }
          } else {
            return (ia->y2 > ib->y2) - (ia->y2 < ib->y2);
          }
        } else {
          return (ia->x2 > ib->x2) - (ia->x2 < ib->x2);
        }
      } else {
        return (ia->z1 > ib->z1) - (ia->z1 < ib->z1);
      }
    } else {
      return (ia->y1 > ib->y1) - (ia->y1 < ib->y1);;
    }
  } else {
    return (ia->x1 > ib->x1) - (ia->x1 < ib->x1);
  }
  /* float comparison: returns negative if b > a 
     and positive if a > b. We multiplied result by 100.0
     to preserve decimal fraction */
} 
 
int main(int argc, char **argv) 
{

  // ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time;
  start_time = clock();

  // check input args
  if( argc < 2 || argc > 5) 
    {
    std::cout << "To check using topology of facet points:              " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)>" << std::endl;
    std::cout << "To check using geometry tolerance of facet points:    " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)> <tolerance>" << std::endl;
    return 1;
    }

  // load file and get tolerance from input argument
  MBErrorCode result;
  std::string filename = argv[1]; //set filename
  MBEntityHandle input_set;
  result = MBI()->create_meshset( MESHSET_SET, input_set ); //create handle to meshset
  if(MB_SUCCESS != result) 
    {
      return result;
    }

  result = MBI()->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(MB_SUCCESS != result) 
    {
      // failed to load the file
      std::cout << "could not load file" << std::endl;
      return result;
    }

  double tol; // tolerance for two verts to be considered the same
  bool check_topology, verbose;

  if(2 == argc) // set topological check
    {
      std::cout << "topology check" << std::endl;
      check_topology = true;
      verbose = false;
    } 
  else if (3 == argc)  // set topological check with different tolerance
    {
      std::cout << "topology check" << std::endl;
      check_topology = true;
      const std::string verbose_string = argv[2];
      verbose = ( 0==verbose_string.compare("true") );
    } 
  else // otherwise do geometry check
    {
      std::cout << "geometry check";
      check_topology = false;
      tol = atof( argv[3] );
      std::cout<< " tolerance=" << tol << std::endl;
      const std::string verbose_string = argv[2];
      verbose = ( 0==verbose_string.compare("true") );
    }

  // create tags on geometry
  MBTag geom_tag, id_tag;
  //result = MBI()->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                            //MB_TYPE_INTEGER, geom_tag, 0, true );
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, sizeof(int), 
                            MB_TYPE_INTEGER, geom_tag, MB_TAG_DENSE, 0, 0 );
  if(MB_SUCCESS != result) 
    {
      return result;  
    }

  //result = MBI()->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                            //MB_TYPE_INTEGER, id_tag, 0, true );
  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, sizeof(int), 
                            MB_TYPE_INTEGER, id_tag, MB_TAG_DENSE, 0, 0 );
  if(MB_SUCCESS != result) 
    {
      return result;
    }
  
  // get surface and volume sets
  MBRange surf_sets, vol_sets; // MBRange of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, surf_sets);
  if(MB_SUCCESS != result) 
    {
      return result;
    }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, vol_sets);
  if(MB_SUCCESS != result)
    {
      return result;
    }

  std::cout<< "number of surfaces=" << surf_sets.size() << std::endl;
  std::cout<< "number of volumes="  << vol_sets.size() << std::endl;

  // ******************************************************************
  // Watertightness is a property of volumes. Check each surface for
  // watertightness.
  // ******************************************************************
  // counted leaky surfaces
  int total_counter = 0, unmatched_counter=0;
  std::set<int> leaky_surfs, leaky_vols;
  MBSkinner tool(MBI());

  // remove all edges for fast skinning
  MBRange edges;

  result = MBI()->get_entities_by_type( 0, MBEDGE, edges ); // get all edges
  if(MB_SUCCESS != result) // failed to get edge data
    {
      return result; // failed
    }

  result = MBI()->delete_entities( edges ); //otherwise delete all edge
  
  if(MB_SUCCESS != result) // failed to delete edge data
    {
      return result; // failed
    }
  
  // loop over each volume meshset
  int vol_counter = 0;
  for(MBRange::iterator i=vol_sets.begin(); i!=vol_sets.end(); ++i) 
    {
      ++vol_counter;
      int surf_counter=0;
      MBRange child_sets;

      result = MBI()->get_child_meshsets( *i, child_sets ); // get child set
      if(MB_SUCCESS != result)  
	{
	  return result; // failed
	}

      // get the volume id of the volume meshset to print a status message
      int vol_id=0;
      // i is the iterator, so &(*i) is a pointer to the first element of MBRange
      result = MBI()->tag_get_data( id_tag, &(*i), 1, &vol_id );

      if(MB_SUCCESS != result)
	{
	  return result;
	}

      if(verbose) 
	{
	  std::cout << "checking volume " << vol_counter << "/" << vol_sets.size()
			  << " id=" << vol_id << std::endl;
	}

      // determine how many skin edges are in each volume
      int n_tris = 0;

      for(MBRange::iterator j=child_sets.begin(); j!=child_sets.end(); ++j) 
	{
	  result = MBI()->get_number_entities_by_type( *j, MBTRI, n_tris ); // for each child set get number of triangles
	  if(MB_SUCCESS != result) 
	    {
	      return result;
	    }
	}

      // save the edges in a vector that is large enough to avoid resizing
      // presumably some kind of efficiency thing?? ad ??
      std::vector<coords_and_id> the_coords_and_id;
      the_coords_and_id.reserve(n_tris);

      // loop over the surface meshsets of each volume meshset
      for(MBRange::iterator j=child_sets.begin(); j!=child_sets.end(); ++j) 
	{

	  // get the surface id of the surface meshset
	  surf_counter++;
	  int surf_id=0;
	  result = MBI()->tag_get_data( id_tag, &(*j), 1, &surf_id );
	  if(MB_SUCCESS != result) 
	    {
	      return result;
	    }

	  // get the range of facets of the surface meshset
	  MBRange facets;
	  result = MBI()->get_entities_by_type( *j, MBTRI, facets );
	  if(MB_SUCCESS != result) 
	    {
	      return result;
	    }

	  // get the range of skin edges from the range of facets
	  // Fiasco: Jason wrote an optimized function (find_skin_vertices) that performed
	  // almost as well as my specialized version (gen::find_skin). When he made then
	  // generalized find_skin_vertices for MOAB it killed performance. As it stands,
	  // gen::find_skin is ~7x faster (January 29, 2010).
	  MBRange skin_edges;
	  if(!facets.empty()) 
	    {
	      result = gen::find_skin( facets, 1, skin_edges, false );
	      if(MB_SUCCESS != result) 
		{
		  return result;
		}
	    }

      // count the number of skin edges in the range
      if(verbose) 
	{
	  std::cout << "surface " << surf_counter << "/" << child_sets.size()
		    << " id=" << surf_id << " contains " << facets.size() 
		    << " facets and " << skin_edges.size() << " skin edges" << std::endl;
	}
     
      for(MBRange::const_iterator k=skin_edges.begin(); k!=skin_edges.end(); ++k) {
	// get the endpoint vertices of the facet edge
	MBRange verts;
	result = MBI()->get_adjacencies( &(*k), 1, 0, false, verts );
        if(MB_SUCCESS != result) return result;
        if(2 != verts.size()) {
          std::cout << "  WARNING: verts.size()=" << verts.size() << std::endl;
	  continue;
        }  

        // Save the range of verts to an array of verts that can store duplicates.
        coords_and_id temp;
        if(check_topology) {
          temp.vert1 = verts[0];
          temp.vert2 = verts[1];
        } else {
          // get the coordinates of endpoint vertices
          double coords0[3], coords1[3];
	  result = MBI()->get_coords( &(verts.front()), 1, coords0 );
          if(MB_SUCCESS != result) return result;
	  result = MBI()->get_coords( &(verts.back()), 1, coords1 );
          if(MB_SUCCESS != result) return result;
          // orient the edge by endpoint coords
          if(!check_topology && 
             coords1[0]< coords0[0] ||
            (coords1[0]==coords0[0] && coords1[1]< coords0[1]) ||
            (coords1[0]==coords0[0] && coords1[1]==coords0[1] && coords1[2]< coords0[2])) {
            temp.x1 = coords1[0];
            temp.y1 = coords1[1];
            temp.z1 = coords1[2];
            temp.x2 = coords0[0];
            temp.y2 = coords0[1];
            temp.z2 = coords0[2];
            temp.vert1 = verts[1];
  	    temp.vert2 = verts[0];
          } else {
            temp.x1 = coords0[0];
            temp.y1 = coords0[1];
            temp.z1 = coords0[2];
            temp.x2 = coords1[0];
            temp.y2 = coords1[1];
            temp.z2 = coords1[2];
            temp.vert1 = verts[0];
  	    temp.vert2 = verts[1];
          }
        }
        temp.surf_id = surf_id;
	temp.matched = false;
        the_coords_and_id.push_back(temp);
      }

      // clean up the edges for the next find_skin call
      result = MBI()->delete_entities( skin_edges );
      if(MB_SUCCESS != result) return result;
      int n_edges;
      result = MBI()->get_number_entities_by_type(0, MBEDGE, n_edges );
      if(MB_SUCCESS != result) return result;
      if(0 != n_edges) return MB_MULTIPLE_ENTITIES_FOUND;
    }

    // sort the edges by the first vert. The first vert has a lower handle than the second.
    int n_edges = the_coords_and_id.size();
    total_counter += n_edges;
    if(check_topology) {
      qsort( &the_coords_and_id[0], n_edges, sizeof(struct coords_and_id), compare_by_handle);
    } else {
      qsort( &the_coords_and_id[0], n_edges, sizeof(struct coords_and_id), compare_by_coords);
    }

    // ******************************************************************
    // Iterate through each facet edge, looking for its match. If a match
    // is found set the edge's flag to 'matched' so that we do not check
    // it again.
    // WARNING: The logic is different for checking by topology vs. proximity.
    // ******************************************************************
    // loop over each facet edge in the volume
    for(int j=0; j!=n_edges; ++j) {

      // if the edge has already been matched, skip it
      if (the_coords_and_id[j].matched) continue;
    
      // try to match the edge with another facet edge:
      for(int k=j+1; k!=n_edges+1; ++k) {

        // look for a match
        if(check_topology) {
	  if( the_coords_and_id[j].vert1==the_coords_and_id[k].vert1 &&
	      the_coords_and_id[j].vert2==the_coords_and_id[k].vert2 ) {
            the_coords_and_id[j].matched = true;
            the_coords_and_id[k].matched = true;
  	    //std::cout<< "matched by handle" << std::endl;
            break;
	  }
        } else {
          // When matching by proximity, it is possible that the k edge has already
          // been matched. If so, skip it.
	  if (the_coords_and_id[k].matched) continue;

          // see if the edge matches
          MBCartVect diff0(the_coords_and_id[j].x1-the_coords_and_id[k].x1,
                           the_coords_and_id[j].y1-the_coords_and_id[k].y1,
                           the_coords_and_id[j].z1-the_coords_and_id[k].z1);
          MBCartVect diff1(the_coords_and_id[j].x2-the_coords_and_id[k].x2,
			   the_coords_and_id[j].y2-the_coords_and_id[k].y2,
			   the_coords_and_id[j].z2-the_coords_and_id[k].z2);
          double d0 = diff0.length_squared();
          double d1 = diff1.length_squared();
          if( d0<tol*tol && d1<tol*tol ) {
            the_coords_and_id[j].matched = true;
            the_coords_and_id[k].matched = true;
  	    //std::cout<< "matched by proximity" << std::endl;
            break;
	  }
          // Due to the sort, once the x-coords are out of tolerance, a match 
          // cannot exist.
          if( the_coords_and_id[k].x1 - the_coords_and_id[j].x1 <= tol) continue;
        }

        // If no break or continue has been hit, the edge is unmatched.    
        // if we have a new leaky surface, save it
	std::set<int>::iterator found;
        found = leaky_surfs.find( the_coords_and_id[j].surf_id );
	if(found == leaky_surfs.end()) {
	  leaky_surfs.insert( the_coords_and_id[j].surf_id );
	}
	found = leaky_vols.find( vol_id );
        if(found == leaky_vols.end()) {
          leaky_vols.insert( vol_id );
        }
        // print info for unmatched edge
        if(verbose) {
          // get the coordinates if we don't already have them
          if(check_topology) {
  	    double endpt_coords[3];
	    result = MBI()->get_coords( &the_coords_and_id[j].vert1, 1, endpt_coords );
            if(MB_SUCCESS != result) return result;
	    the_coords_and_id[j].x1 = endpt_coords[0]; 
	    the_coords_and_id[j].y1 = endpt_coords[1];
	    the_coords_and_id[j].z1 = endpt_coords[2];
            result = MBI()->get_coords( &the_coords_and_id[j].vert2, 1, endpt_coords );
            if(MB_SUCCESS != result) return result;
	    the_coords_and_id[j].x2 = endpt_coords[0]; 
	    the_coords_and_id[j].y2 = endpt_coords[1];
	    the_coords_and_id[j].z2 = endpt_coords[2];
          }
          std::cout << "  edge of surf " << the_coords_and_id[j].surf_id 
                 << " unmatched: " << " (" 
                 << the_coords_and_id[j].x1 << "," 
	         << the_coords_and_id[j].y1 << ","
	         << the_coords_and_id[j].z1 << ") (" 
                 << the_coords_and_id[j].x2 << ","
	         << the_coords_and_id[j].y2 << ","
	         << the_coords_and_id[j].z2 << ")" 
		 << " v0=" << the_coords_and_id[j].vert1 
                 << " v1=" << the_coords_and_id[j].vert2
                 << " j=" << j << " k=" << k <<std::endl;
        }
        unmatched_counter++;
        break;        
      }  // k loop
    }    // j loop
  }      // volume loop

  // print time and summary
  clock_t end_time = clock();
  std::cout << std::endl << unmatched_counter << "/" << total_counter << " ("
            << double(100.0*unmatched_counter)/total_counter 
            << "%) unmatched edges" << std::endl;
  std::cout << leaky_surfs.size() << "/" << surf_sets.size() << " ("
            << double(100.0*leaky_surfs.size())/surf_sets.size() 
            << "%) unsealed surfaces" << std::endl;
  std::cout << leaky_vols.size() << "/" << vol_sets.size() << " ("
            << double(100.0*leaky_vols.size())/vol_sets.size() 
            << "%) unsealed volumes" << std::endl;
  std::cout << (double) (end_time-start_time)/CLOCKS_PER_SEC << " seconds" << std::endl;

  // list the leaky surface and volume ids
  std::cout << "leaky surface ids=";
  for( std::set<int>::iterator i=leaky_surfs.begin(); i!=leaky_surfs.end(); i++) {
    std::cout << *i << " ";
  }
  std::cout << std::endl;
  std::cout << "leaky volume ids=";
  for( std::set<int>::iterator i=leaky_vols.begin(); i!=leaky_vols.end(); i++) {
    std::cout << *i << " ";
  }
  std::cout << std::endl;

}

MBInterface* MBI() {
  static MBCore instance;
  return &instance;
}
