// *********************************

// Patrick Shriwise
// October, 2013

// functions needed to check models for watertightness

#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <set>
#include <algorithm>

// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "CheckWatertight.hpp"

moab::ErrorCode CheckWatertight::check_mesh_for_watertightness( moab::EntityHandle input_set, double tol, bool &sealed, bool test, bool verbose, bool check_topology )
{

  moab::ErrorCode result;

// create tags on geometry
  moab::Tag geom_tag, id_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
                                  moab::MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  MB_CHK_SET_ERR(result,"could not get GEOM_DIMENSION_TAG_NAME handle");

  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1,
                                  moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );

  MB_CHK_SET_ERR(result, "could not get GLOBAL_ID_TAG_NAME handle");


  // get surface and volume sets
  moab::Range surf_sets, vol_sets; // moab::Range of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, moab::MBENTITYSET, &geom_tag,
           input_dim, 1, surf_sets);
  if(moab::MB_SUCCESS != result) {
    return result;
  }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, moab::MBENTITYSET, &geom_tag,
           input_dim, 1, vol_sets);
  if(moab::MB_SUCCESS != result) {
    return result;
  }
  if(!test) {
    std::cout<< "number of surfaces=" << surf_sets.size() << std::endl;
    std::cout<< "number of volumes="  << vol_sets.size() << std::endl;
  }
  // ******************************************************************
  // Watertightness is a property of volumes. Check each surface for
  // watertightness.
  // ******************************************************************
  // counted leaky surfaces
  int total_counter = 0, unmatched_counter=0;
  std::set<int> leaky_surfs, leaky_vols;
  moab::Skinner tool(MBI());

  // remove all edges for fast skinning
  moab::Range edges;

  result = MBI()->get_entities_by_type( 0, moab::MBEDGE, edges ); // get all edges
  if(moab::MB_SUCCESS != result) { // failed to get edge data
    return result; // failed
  }
//10/11/13
//removed as a result of the change from the gen::find_skin function to the moab::Skinner:find_skin function
  /*
    result = MBI()->delete_entities( edges ); //otherwise delete all edge

    if(moab::MB_SUCCESS != result) // failed to delete edge data
      {
        return result; // failed
      }
  */
  // loop over each volume meshset
  int vol_counter = 0;
  for(moab::Range::iterator i=vol_sets.begin(); i!=vol_sets.end(); ++i) {
    ++vol_counter;
    int surf_counter=0;
    moab::Range child_sets;

    result = MBI()->get_child_meshsets( *i, child_sets ); // get child set
    if(moab::MB_SUCCESS != result) {
      return result; // failed
    }

    // get the volume id of the volume meshset to print a status message
    int vol_id=0;
    // i is the iterator, so &(*i) is a pointer to the first element of moab::Range
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &vol_id );

    if(moab::MB_SUCCESS != result) {
      return result;
    }

    if(verbose) {
      std::cout << "checking volume " << vol_counter << "/" << vol_sets.size()
                << " id=" << vol_id << std::endl;
    }

    // determine how many skin edges are in each volume
    int n_tris = 0;

    for(moab::Range::iterator j=child_sets.begin(); j!=child_sets.end(); ++j) {
      result = MBI()->get_number_entities_by_type( *j, moab::MBTRI, n_tris ); // for each child set get number of triangles
      if(moab::MB_SUCCESS != result) {
        return result;
      }
    }

    // save the edges in a vector that is large enough to avoid resizing
    // presumably some kind of efficiency thing?? ad ??
    std::vector<coords_and_id> the_coords_and_id;
    the_coords_and_id.reserve(n_tris);

    // loop over the surface meshsets of each volume meshset
    for(moab::Range::iterator j=child_sets.begin(); j!=child_sets.end(); ++j) {

      // get the surface id of the surface meshset
      surf_counter++;
      int surf_id=0;
      result = MBI()->tag_get_data( id_tag, &(*j), 1, &surf_id );
      if(moab::MB_SUCCESS != result) {
        return result;
      }

      // get the range of facets of the surface meshset
      moab::Range facets;
      result = MBI()->get_entities_by_type( *j, moab::MBTRI, facets );
      if(moab::MB_SUCCESS != result) {
        return result;
      }

      // get the range of skin edges from the range of facets
      // Fiasco: Jason wrote an optimized function (find_skin_vertices) that performed
      // almost as well as my specialized version (gen::find_skin). When he made then
      // generalized find_skin_vertices for MOAB it killed performance. As it stands,
      // gen::find_skin is ~7x faster (January 29, 2010).
      moab::Range skin_edges;
      moab::Skinner sk(MBI());
      if(!facets.empty()) {
        result = sk.find_skin( 0 , facets, 1, skin_edges, false );
        if(moab::MB_SUCCESS != result) {
          return result;
        }
      }

      // count the number of skin edges in the range
      if(verbose) {
        std::cout << "surface " << surf_counter << "/" << child_sets.size()
                  << " id=" << surf_id << " contains " << facets.size()
                  << " facets and " << skin_edges.size() << " skin edges" << std::endl;
      }

      for(moab::Range::const_iterator k=skin_edges.begin(); k!=skin_edges.end(); ++k) {
        // get the endpoint vertices of the facet edge
        moab::Range verts;
        result = MBI()->get_adjacencies( &(*k), 1, 0, false, verts );
        if(moab::MB_SUCCESS != result) return result;
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
          if(moab::MB_SUCCESS != result) return result;
          result = MBI()->get_coords( &(verts.back()), 1, coords1 );
          if(moab::MB_SUCCESS != result) return result;
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
      //10/10/13
      // Removed the following to avoid altering the data set at all
      // -No need to delete skin_edges with the moab:Skinner find_skin function
      // -skin_edges size will be reset to zero upon new moab::Range skin_edges; call
      // clean up the edges for the next find_skin call
      //result = MBI()->delete_entities( skin_edges );
      //if(moab::MB_SUCCESS != result) return result;

      //10/10/13
      // - No ned to ensure edges aren't in the meshset with moab::Skinner find_skin function
      //int n_edges;
      //result = MBI()->get_number_entities_by_type(0, moab::MBEDGE, n_edges );
      //if(moab::MB_SUCCESS != result) return result;
      //if(0 != n_edges) MB_CHK_SET_ERR(moab::MB_MULTIPLE_ENTITIES_FOUND,"n_edges not equal to zero");
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
          moab::CartVect diff0(the_coords_and_id[j].x1-the_coords_and_id[k].x1,
                               the_coords_and_id[j].y1-the_coords_and_id[k].y1,
                               the_coords_and_id[j].z1-the_coords_and_id[k].z1);
          moab::CartVect diff1(the_coords_and_id[j].x2-the_coords_and_id[k].x2,
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
            if(moab::MB_SUCCESS != result) return result;
            the_coords_and_id[j].x1 = endpt_coords[0];
            the_coords_and_id[j].y1 = endpt_coords[1];
            the_coords_and_id[j].z1 = endpt_coords[2];
            result = MBI()->get_coords( &the_coords_and_id[j].vert2, 1, endpt_coords );
            if(moab::MB_SUCCESS != result) return result;
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






  if(!test) {
    // print time and summary

    std::cout << std::endl << unmatched_counter << "/" << total_counter << " ("
              << double(100.0*unmatched_counter)/total_counter
              << "%) unmatched edges" << std::endl;
    std::cout << leaky_surfs.size() << "/" << surf_sets.size() << " ("
              << double(100.0*leaky_surfs.size())/surf_sets.size()
              << "%) unsealed surfaces" << std::endl;
    std::cout << leaky_vols.size() << "/" << vol_sets.size() << " ("
              << double(100.0*leaky_vols.size())/vol_sets.size()
              << "%) unsealed volumes" << std::endl;


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


  if(unmatched_counter == 0 && leaky_vols.size() == 0 && leaky_surfs.size() == 0 ) {
    sealed=true;
  } else {
    sealed=false;
  }
  return moab::MB_SUCCESS;
}

