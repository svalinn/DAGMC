// ********************************************************************
// Brandon Smith
// August, 2009

/* _curve_to_be_tested_for_watertightness_
      vert1 X X vert1
            | |
      vert2 X |
  surf1     | |    surf2
            | |
      vert3 X X vert2
            | |
      vert4 X X vert3                   */

// input:  h5m filename, tolerance
// output: watertight h5m

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <cassert>
#include <cmath>
#include <ctime>
#include <vector>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"

MBInterface *MBI();

MBErrorCode delete_all_edges() {
  // delete all of the edges. Never keep edges. They are too hard to track and use
  // due to orientation and multiple entities errors when merging.
  MBErrorCode result;
  MBRange edges;
  result = MBI()->get_entities_by_type( 0, MBEDGE, edges );
  assert(MB_SUCCESS == result);
  result = MBI()->delete_entities( edges );
  assert(MB_SUCCESS == result);
  return MB_SUCCESS;
}

// Input: unordered sets of curves that do not track ownership
// Output: ordered sets of verts that do track ownership. All edges are deleted.
MBErrorCode prepare_curves(MBRange &curve_sets, 
                           MBTag geom_tag, MBTag id_tag, MBTag merge_tag,
                           double const MERGE_TOL, double const FACET_TOL) {
  MBErrorCode result;

  // process each curve
  for(MBRange::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++ ) {
    // get the curve id of the curve meshset
    int id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &id );
    assert(MB_SUCCESS == result);
    std::cout << "curve " << id << std::endl;

    // get the range of edges of the curve meshset
    std::vector<MBEntityHandle> curve_edges;
    result = MBI()->get_entities_by_type( *i, MBEDGE, curve_edges );
    assert( MB_SUCCESS == result );

    /* Merge the endpoints of the curve and remove its edges if it is too small.
    Use the MERGE_TOL because these edges will be merged with the MERGE_TOL 
    during zipping anyhow. Doing this now removes small curves from zipping and
    reduces ambiguity. */
    if(MERGE_TOL > gen::length(curve_edges)) {
      std::cout << "  removed curve with length=" << gen::length(curve_edges) 
                << " n_verts=" << curve_edges.size()+1 << std::endl;

      // get the endpoints of the curve
      MBRange endpt_sets;
      result = MBI()->get_child_meshsets( *i, endpt_sets );
      assert(MB_SUCCESS==result);
      std::cout << "  endpt_sets.size()=" << endpt_sets.size() << std::endl;
      if(endpt_sets.empty()) {
        assert(false);
      } else if(1 == endpt_sets.size()) {
        // nothing
      } else if(2 == endpt_sets.size()) {
        MBRange front_endpt, back_endpt;
        result = MBI()->get_entities_by_type( endpt_sets.front(), MBVERTEX, front_endpt);
        assert(MB_SUCCESS == result);
        assert(1 == front_endpt.size());
        result = MBI()->get_entities_by_type( endpt_sets.back(), MBVERTEX, back_endpt);
        assert(MB_SUCCESS == result);
        assert(1 == back_endpt.size());
        // merge the endpoints-ALWAYS CHECK TO AVOID MERGING THE SAME ENTITY!!!
        if(front_endpt[0] != back_endpt[0]) {
          result = MBI()->merge_entities( front_endpt[0], back_endpt[0], false, true);
          assert(MB_SUCCESS == result);
          // check for and remove degenerate edges caused by the merge
          MBRange edges;
          MBEntityHandle temp = front_endpt[0];
          result = MBI()->get_adjacencies( &temp, 1, 1, false, edges);
          assert(MB_SUCCESS == result);
          for(MBRange::iterator j=edges.begin(); j!=edges.end(); j++) {
            const MBEntityHandle *conn;
            int n_verts;
            result = MBI()->get_connectivity( *j, conn, n_verts); 
            assert(MB_SUCCESS == result);
            if(conn[0] == conn[1]) {
              result = MBI()->delete_entities( &(*j), 1 );
              assert(MB_SUCCESS == result);
            }
          }
	}
      } else {
        assert(false);
      }
      // It is possible that the endpoints themselves are orphaned. Should these
      // be deleted?
      
      // Remove the curve set. This also removes parent-child relationships.
      result = MBI()->delete_entities( &(*i), 1);
      assert(MB_SUCCESS == result);
      i = curve_sets.erase(i) - 1;
      continue;
    }

    // convert the curve of edges into a curve of verts 
    std::vector<MBEntityHandle> ordered_verts;
    result = gen::ordered_verts_from_ordered_edges( curve_edges, ordered_verts);
    assert(MB_SUCCESS == result);

    // the edges are no longer needed
    result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
    assert(MB_SUCCESS == result);

    // replace the unordered edges with the ordered verts
    result = arc::set_meshset( *i, ordered_verts );
    assert(MB_SUCCESS == result);      
  }

  return MB_SUCCESS;
}

/* Isolate the failure by removing the curve and loop that failed. The zip_loop
   function will be called again on the remaining loops and curves. */
MBErrorCode remove_failed_loop_and_curve( std::vector<std::vector<MBEntityHandle> > &skin,
                                          std::vector<std::vector<MBEntityHandle> > &curves,
                                          std::vector<int> &curve_ids,
                                          std::vector<MBEntityHandle> &curve_sets,
                                          //MBRange &curve_sets,
					  const unsigned int loop,
                                          const unsigned int curve ) {
  skin.erase( skin.begin()+loop );
  curves.erase( curves.begin()+curve );
  curve_ids.erase( curve_ids.begin()+curve ); 
  curve_sets.erase( curve_sets.begin()+curve );
  std::cout << "remove_failed_loop: removed loop " << loop << std::endl;
  return MB_SUCCESS;
}

  // input: surface sets, ordered curve sets,
  // output: skin arcs corresponding to curves are added to parent surface sets
MBErrorCode prepare_surfaces(MBRange &surface_sets,
                             MBTag geom_tag, MBTag id_tag, MBTag normal_tag, MBTag merge_tag,
                             const double SME_RESABS_TOL, const double FACET_TOL, 
                               const double MERGE_TOL) {
    
    MBErrorCode result;

    // loop over each surface meshset
    for(MBRange::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) {

      // get the surf id of the surface meshset
      int surf_id;
      result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );
      assert(MB_SUCCESS == result);
      std::cout << "  surf id=" << surf_id << std::endl;

      // get facets of the surface meshset
      MBRange tris;
      result = MBI()->get_entities_by_type( *i, MBTRI, tris );
      assert(MB_SUCCESS == result);

      // Get the curves sets
      std::vector<MBEntityHandle> curve_sets, unmerged_curve_sets;
      result = MBI()->get_child_meshsets( *i, curve_sets );
      assert(MB_SUCCESS==result);

      // Update the curve_sets with that contain entity_to_delete curves with their
      // entity_to_keep curves. Unmerged_curve_sets will end up holding the curves
      // of this surface that are not merged with another curve in this surface.
      for(std::vector<MBEntityHandle>::iterator j=curve_sets.begin();
	  j!=curve_sets.end(); j++) {
        MBEntityHandle merged_curve, curve;
        result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );
        assert(MB_TAG_NOT_FOUND==result || MB_SUCCESS==result);     
        if(MB_TAG_NOT_FOUND==result) {
          curve = *j;
        } else if(MB_SUCCESS == result) {
	  std::cout << "  curve " << gen::geom_id_by_handle(*j) 
                    << " is entity_to_delete" << std::endl;
          curve = merged_curve;
          // should parent-childs be updated for the entity_to_keep?
        } else {
	  std::cout << "prepare_surfaces: result=" << result << std::endl;
          return result;        
        }
      
        // If the curve is in unmerged_curve_sets, then remove it. Otherwise add it.
	std::vector<MBEntityHandle>::iterator k=find(unmerged_curve_sets.begin(), 
	  unmerged_curve_sets.end(), curve);
        if(unmerged_curve_sets.end() == k) {
  	  //std::cout << "  curve " << gen::geom_id_by_handle(*k) 
          //          << " is entity_to_keep" << std::endl;
          unmerged_curve_sets.push_back(curve);
        } else {
          unmerged_curve_sets.erase(k);
        }
      }

      // If all of the curves are merged, remove the surfaces facets.
      if(unmerged_curve_sets.empty()) {
        result = MBI()->remove_entities( *i, tris);                                           
	assert(MB_SUCCESS == result);                                                         
	std::cout << "  removed " << tris.size() << " facets and deleted surface" << std::endl;
	result = MBI()->delete_entities( tris );                                              
	assert(MB_SUCCESS == result);
        // remove the surface set itself
        result = MBI()->delete_entities( &(*i), 1);
        assert(MB_SUCCESS == result);
        i = surface_sets.erase(i) - 1;
        continue;
      }

      // Try zipping without curves that are merged with each other
      curve_sets.swap(unmerged_curve_sets);

      // Save the normals of the facets. These will later be used to determine if
      // the tri became inverted.
      result = gen::save_normals( tris, normal_tag );
      assert(MB_SUCCESS == result);
  
      // get the range of skin edges from the range of facets
      MBSkinner tool(MBI());
      MBRange skin_edges;

      // merge the vertices of the skin
      // BRANDON: For some reason cgm2moab does not do this? This was the 
      // problem with mod13 surf 881. Two skin verts were coincident. A tol=1e-10
      // found the verts, but tol=0 did not.
      MBRange skin_verts;
      result = MBI()->get_adjacencies( skin_edges, 0, false, skin_verts, 
                                       MBInterface::UNION );
      assert(MB_SUCCESS == result);
      result = gen::merge_vertices( skin_verts, SME_RESABS_TOL );         
      if (MB_SUCCESS != result) {                                             
	std::cout << "result= " << result << std::endl;                  
	std::cout << "SURFACE_ZIPPING_FAILURE: could not merge vertices, surf_id="   
		  << surf_id << std::endl;  
	continue;
      }  
						      
      // Create loops with the skin edges.  
      std::vector< std::vector<MBEntityHandle> > skin_loops_of_edges;
      if(MB_SUCCESS != result) {
	std::cout << "SURFACE_ZIPPING_FAILURE: could not create loops for surf_id=" 
		  << surf_id << std::endl;
	continue;
      }
      std::cout << "    surf has " << skin_loops_of_edges.size() 
                << " skin loop(s)." << std::endl;
    
      // Convert the loops of skin edges to loops of skin verts.
      std::vector< std::vector<MBEntityHandle> > skin(skin_loops_of_edges.size());
      for(unsigned int j=0; j<skin_loops_of_edges.size(); j++) {
        result = gen::ordered_verts_from_ordered_edges( skin_loops_of_edges[j], skin[j] );
        assert(MB_SUCCESS == result);
	// check to make sure that the loop is closed
	assert(skin[j].front() == skin[j].back());
      }

      // edges are no longer needed       
      result = delete_all_edges();
      assert(MB_SUCCESS == result);

      /* Get the curves that are part of the surface. Use vectors to store all curve
	 stuff so that we can remove curves from the set as they are zipped. */
      //curve_sets.clear();

      //result = MBI()->get_child_meshsets( *i, curve_sets );
      //assert(MB_SUCCESS==result);
      std::vector<int> curve_ids;
      int curve_id;
      std::vector<std::vector<MBEntityHandle> > curves;
      //for(MBRange::iterator j=curve_sets.begin(); j!=curve_sets.end(); j++) {
      //for(unsigned int j=0; j<curve_sets.size(); j++) {
      for(std::vector<MBEntityHandle>::iterator j=curve_sets.begin(); 
        j!=curve_sets.end(); j++) {

        // If a delete_curve, replace it with the keep_curve. This approach allows
        // for duplicates because we are using vectors instead of ranges. Note that
        // parent-child links also cannot store duplicate handles.
        MBEntityHandle merged_curve;
	//MBEntityHandle temp = curve_sets[j];
        result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );     
        assert(MB_TAG_NOT_FOUND==result || MB_SUCCESS==result);
        if(MB_SUCCESS == result) *j = merged_curve;

	// do not add a curve if it contains nothing
	//temp = curve_sets[j];
	result = MBI()->tag_get_data( id_tag, &(*j), 1, &curve_id );
	assert(MB_SUCCESS == result);
	std::cout << "  curve_id=" << curve_id << " handle=" << *j << std::endl;
	curve_ids.push_back(curve_id);
	std::vector<MBEntityHandle> curve;
	result = arc::get_meshset( *j, curve );
	assert(MB_SUCCESS == result);
	curves.push_back( curve );
      }

      // Keep zipping loops until each is either zipped or failed. This function
      // returns only after all loops are zipped or a failure occurs.
      while(!skin.empty()) {
        //result = zip_loop( normal_tag, FACET_TOL, MERGE_TOL, 
	//                 curves, skin, curve_ids, *i, curve_sets );
        if(MB_SUCCESS != result) {
	  std::cout << "SURFACE_ZIPPING_FAILURE: could not zip surf_id=" << surf_id << std::endl;
        }
      }

      // mod13surf2996, 3028 and 2997 are adjacent to the same bad geometry (figure 8 loop)
      //assert(MB_SUCCESS==result || 2996==surf_id || 2997==surf_id || 3028==surf_id);
    }
    return MB_SUCCESS;
  }

  MBErrorCode test_edges() {
    MBErrorCode result;
    MBRange edges;
    result = MBI()->get_entities_by_dimension( 0, 1, edges );
    assert(MB_SUCCESS == result);
    MBI()->list_entities( edges );
    return MB_SUCCESS;
  }


  int main(int argc, char **argv) {

    // ******************************************************************
    // Load the h5m file and create tags.
    // ******************************************************************

    clock_t start_time = clock(), prep_time, zip_time;
    if(2 > argc) {
      std::cout << "usage: do not use" << std::endl;
      std::cout << "./post_process <input_file>" << std::endl;
      return 1;
    }
    MBErrorCode result;
    std::string input_name = argv[1];

    // The root name does not have an extension
    std::string root_name = argv[1];
    int len = root_name.length();
    root_name.erase(len - 4);
    const double MERGE_TOL = 1e-3; // should this depend on FACET_TOL? 

    // load the input file
    MBEntityHandle input_meshset;
    result = MBI()->create_meshset( MESHSET_SET, input_meshset );
    assert(MB_SUCCESS == result);
    if(std::string::npos != input_name.find("h5m")) {
      result = MBI()->load_file( input_name.c_str(), &input_meshset );
      assert( MB_SUCCESS == result );
    } else {
      std::cout << "invalid input file: must be h5m" << std::endl;
      return 1;
    }

    // create tags
    MBTag geom_tag, id_tag, sense_tag, normal_tag, merge_tag;
    result = MBI()->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE,
				MB_TYPE_INTEGER, geom_tag, 0, true );
    assert( MB_SUCCESS == result );
    result = MBI()->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
				MB_TYPE_INTEGER, id_tag, 0, true );
    assert( MB_SUCCESS == result );
    result = MBI()->tag_create( "GEOM_SENSE_2", 2*sizeof(MBEntityHandle), MB_TAG_DENSE,
                                MB_TYPE_HANDLE, sense_tag, 0, true );
    assert( MB_SUCCESS == result );
    result = MBI()->tag_create( "NORMAL", sizeof(MBCartVect), MB_TAG_DENSE,
                                MB_TYPE_OPAQUE, normal_tag, 0, true );
    assert( MB_SUCCESS == result );
    result = MBI()->tag_create( "MERGE", sizeof(MBEntityHandle), MB_TAG_SPARSE,
                                MB_TYPE_HANDLE, merge_tag, 0, true );
    assert( MB_SUCCESS == result );  

    // get all geometry sets
    MBRange geom_sets[4];
    for(unsigned dim=0; dim<4; dim++) {
      void *val[] = {&dim};
      result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
      assert(MB_SUCCESS == result);
      // make sure that sets TRACK membership and curves are ordered
      // MESHSET_TRACK_OWNER=0x1, MESHSET_SET=0x2, MESHSET_ORDERED=0x4
      for(MBRange::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
        unsigned int options;
        result = MBI()->get_meshset_options(*i, options );
        assert(MB_SUCCESS == result);
    
        // if options are wrong change them
        if(dim==1) {
          if( !(MESHSET_TRACK_OWNER&options) || !(MESHSET_ORDERED&options) ) {
	    result = MBI()->set_meshset_options(*i, MESHSET_TRACK_OWNER|MESHSET_ORDERED);
            assert(MB_SUCCESS == result);
          }
        } else {
          if( !(MESHSET_TRACK_OWNER&options) ) {        
	    result = MBI()->set_meshset_options(*i, MESHSET_TRACK_OWNER);
            assert(MB_SUCCESS == result);
          }
        }
      }
    }
    std::cout << geom_sets[3].size() << " volumes, " 
              << geom_sets[2].size() << " surfaces, and "
              << geom_sets[1].size() << " curves" << std::endl;  

    result = cleanup::delete_small_edges(geom_sets[2], MERGE_TOL);
    assert(MB_SUCCESS == result);
  
    std::string output_filename = root_name + "_tri.h5m";
    // PROBLEM: If I write the input meshset the writer returns MB_FAILURE.
    // This happens only if I delete vertices when merging.
    // result = MBI()->write_mesh( filename_new.c_str(), &input_meshset, 1);
    result = MBI()->write_mesh( output_filename.c_str() );
    if (MB_SUCCESS != result) std::cout << "result= " << result << std::endl;
    assert(MB_SUCCESS == result);

    zip_time = clock();
    std::cout << "zipping took " << (double) (zip_time-prep_time)/CLOCKS_PER_SEC 
	      << " sec." << std::endl;
  
    return 0;  
  }

  MBInterface *MBI() {
    static MBCore instance;
    return &instance;
  }

