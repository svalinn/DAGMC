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
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"

moab::Interface *MBI();

moab::ErrorCode delete_all_edges()
{
  // delete all of the edges. Never keep edges. They are too hard to track and use
  // due to orientation and multiple entities errors when merging.
  moab::ErrorCode result;
  moab::Range edges;
  result = MBI()->get_entities_by_type( 0, moab::MBEDGE, edges );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->delete_entities( edges );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

// Input: unordered sets of curves that do not track ownership
// Output: ordered sets of verts that do track ownership. All edges are deleted.
moab::ErrorCode prepare_curves(moab::Range &curve_sets,
                               moab::Tag geom_tag, moab::Tag id_tag, moab::Tag merge_tag,
                               double const MERGE_TOL, double const FACET_TOL)
{
  moab::ErrorCode result;

  // process each curve
  for(moab::Range::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++ ) {
    // get the curve id of the curve meshset
    int id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &id );
    assert(moab::MB_SUCCESS == result);
    std::cout << "curve " << id << std::endl;

    // get the range of edges of the curve meshset
    std::vector<moab::EntityHandle> curve_edges;
    result = MBI()->get_entities_by_type( *i, moab::MBEDGE, curve_edges );
    assert( moab::MB_SUCCESS == result );

    /* Merge the endpoints of the curve and remove its edges if it is too small.
    Use the MERGE_TOL because these edges will be merged with the MERGE_TOL
    during zipping anyhow. Doing this now removes small curves from zipping and
    reduces ambiguity. */
    if(MERGE_TOL > gen::length(curve_edges)) {
      std::cout << "  removed curve with length=" << gen::length(curve_edges)
                << " n_verts=" << curve_edges.size()+1 << std::endl;

      // get the endpoints of the curve
      moab::Range endpt_sets;
      result = MBI()->get_child_meshsets( *i, endpt_sets );
      assert(moab::MB_SUCCESS==result);
      std::cout << "  endpt_sets.size()=" << endpt_sets.size() << std::endl;
      if(endpt_sets.empty()) {
        assert(false);
      } else if(1 == endpt_sets.size()) {
        // nothing
      } else if(2 == endpt_sets.size()) {
        moab::Range front_endpt, back_endpt;
        result = MBI()->get_entities_by_type( endpt_sets.front(), moab::MBVERTEX, front_endpt);
        assert(moab::MB_SUCCESS == result);
        assert(1 == front_endpt.size());
        result = MBI()->get_entities_by_type( endpt_sets.back(), moab::MBVERTEX, back_endpt);
        assert(moab::MB_SUCCESS == result);
        assert(1 == back_endpt.size());
        // merge the endpoints-ALWAYS CHECK TO AVOID MERGING THE SAME ENTITY!!!
        if(front_endpt[0] != back_endpt[0]) {
          result = MBI()->merge_entities( front_endpt[0], back_endpt[0], false, true);
          assert(moab::MB_SUCCESS == result);
          // check for and remove degenerate edges caused by the merge
          moab::Range edges;
          moab::EntityHandle temp = front_endpt[0];
          result = MBI()->get_adjacencies( &temp, 1, 1, false, edges);
          assert(moab::MB_SUCCESS == result);
          for(moab::Range::iterator j=edges.begin(); j!=edges.end(); j++) {
            const moab::EntityHandle *conn;
            int n_verts;
            result = MBI()->get_connectivity( *j, conn, n_verts);
            assert(moab::MB_SUCCESS == result);
            if(conn[0] == conn[1]) {
              result = MBI()->delete_entities( &(*j), 1 );
              assert(moab::MB_SUCCESS == result);
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
      assert(moab::MB_SUCCESS == result);
      i = curve_sets.erase(i) - 1;
      continue;
    }

    // convert the curve of edges into a curve of verts
    std::vector<moab::EntityHandle> ordered_verts;
    result = gen::ordered_verts_from_ordered_edges( curve_edges, ordered_verts);
    assert(moab::MB_SUCCESS == result);

    // the edges are no longer needed
    result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
    assert(moab::MB_SUCCESS == result);

    // replace the unordered edges with the ordered verts
    result = arc::set_meshset( *i, ordered_verts );
    assert(moab::MB_SUCCESS == result);
  }

  return moab::MB_SUCCESS;
}

/* Isolate the failure by removing the curve and loop that failed. The zip_loop
   function will be called again on the remaining loops and curves. */
moab::ErrorCode remove_failed_loop_and_curve( std::vector<std::vector<moab::EntityHandle> > &skin,
    std::vector<std::vector<moab::EntityHandle> > &curves,
    std::vector<int> &curve_ids,
    std::vector<moab::EntityHandle> &curve_sets,
    //moab::Range &curve_sets,
    const unsigned int loop,
    const unsigned int curve )
{
  skin.erase( skin.begin()+loop );
  curves.erase( curves.begin()+curve );
  curve_ids.erase( curve_ids.begin()+curve );
  curve_sets.erase( curve_sets.begin()+curve );
  std::cout << "remove_failed_loop: removed loop " << loop << std::endl;
  return moab::MB_SUCCESS;
}

// input: surface sets, ordered curve sets,
// output: skin arcs corresponding to curves are added to parent surface sets
moab::ErrorCode prepare_surfaces(moab::Range &surface_sets,
                                 moab::Tag geom_tag, moab::Tag id_tag, moab::Tag normal_tag, moab::Tag merge_tag,
                                 const double SME_RESABS_TOL, const double FACET_TOL,
                                 const double MERGE_TOL)
{

  moab::ErrorCode result;

  // loop over each surface meshset
  for(moab::Range::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) {

    // get the surf id of the surface meshset
    int surf_id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );
    assert(moab::MB_SUCCESS == result);
    std::cout << "  surf id=" << surf_id << std::endl;

    // get facets of the surface meshset
    moab::Range tris;
    result = MBI()->get_entities_by_type( *i, moab::MBTRI, tris );
    assert(moab::MB_SUCCESS == result);

    // Get the curves sets
    std::vector<moab::EntityHandle> curve_sets, unmerged_curve_sets;
    result = MBI()->get_child_meshsets( *i, curve_sets );
    assert(moab::MB_SUCCESS==result);

    // Update the curve_sets with that contain entity_to_delete curves with their
    // entity_to_keep curves. Unmerged_curve_sets will end up holding the curves
    // of this surface that are not merged with another curve in this surface.
    for(std::vector<moab::EntityHandle>::iterator j=curve_sets.begin();
        j!=curve_sets.end(); j++) {
      moab::EntityHandle merged_curve, curve;
      result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );
      assert(moab::MB_TAG_NOT_FOUND==result || moab::MB_SUCCESS==result);
      if(moab::MB_TAG_NOT_FOUND==result) {
        curve = *j;
      } else if(moab::MB_SUCCESS == result) {
        std::cout << "  curve " << gen::geom_id_by_handle(*j)
                  << " is entity_to_delete" << std::endl;
        curve = merged_curve;
        // should parent-childs be updated for the entity_to_keep?
      } else {
        std::cout << "prepare_surfaces: result=" << result << std::endl;
        return result;
      }

      // If the curve is in unmerged_curve_sets, then remove it. Otherwise add it.
      std::vector<moab::EntityHandle>::iterator k=find(unmerged_curve_sets.begin(),
          unmerged_curve_sets.end(), curve);
      if(unmerged_curve_sets.end() == k) {
        unmerged_curve_sets.push_back(curve);
      } else {
        unmerged_curve_sets.erase(k);
      }
    }

    // If all of the curves are merged, remove the surfaces facets.
    if(unmerged_curve_sets.empty()) {
      result = MBI()->remove_entities( *i, tris);
      assert(moab::MB_SUCCESS == result);
      std::cout << "  removed " << tris.size() << " facets and deleted surface" << std::endl;
      result = MBI()->delete_entities( tris );
      assert(moab::MB_SUCCESS == result);
      // remove the surface set itself
      result = MBI()->delete_entities( &(*i), 1);
      assert(moab::MB_SUCCESS == result);
      i = surface_sets.erase(i) - 1;
      continue;
    }

    // Try zipping without curves that are merged with each other
    curve_sets.swap(unmerged_curve_sets);

    // Save the normals of the facets. These will later be used to determine if
    // the tri became inverted.
    result = gen::save_normals( tris, normal_tag );
    assert(moab::MB_SUCCESS == result);

    // get the range of skin edges from the range of facets
    moab::Skinner tool(MBI());
    moab::Range skin_edges;

    // merge the vertices of the skin
    // BRANDON: For some reason cgm2moab does not do this? This was the
    // problem with mod13 surf 881. Two skin verts were coincident. A tol=1e-10
    // found the verts, but tol=0 did not.
    moab::Range skin_verts;
    result = MBI()->get_adjacencies( skin_edges, 0, false, skin_verts,
                                     moab::Interface::UNION );
    assert(moab::MB_SUCCESS == result);
    result = gen::merge_vertices( skin_verts, SME_RESABS_TOL );
    if (moab::MB_SUCCESS != result) {
      std::cout << "result= " << result << std::endl;
      std::cout << "SURFACE_ZIPPING_FAILURE: could not merge vertices, surf_id="
                << surf_id << std::endl;
      continue;
    }

    // Create loops with the skin edges.
    std::vector< std::vector<moab::EntityHandle> > skin_loops_of_edges;
    if(moab::MB_SUCCESS != result) {
      std::cout << "SURFACE_ZIPPING_FAILURE: could not create loops for surf_id="
                << surf_id << std::endl;
      continue;
    }
    std::cout << "    surf has " << skin_loops_of_edges.size()
              << " skin loop(s)." << std::endl;

    // Convert the loops of skin edges to loops of skin verts.
    std::vector< std::vector<moab::EntityHandle> > skin(skin_loops_of_edges.size());
    for(unsigned int j=0; j<skin_loops_of_edges.size(); j++) {
      result = gen::ordered_verts_from_ordered_edges( skin_loops_of_edges[j], skin[j] );
      assert(moab::MB_SUCCESS == result);
      // check to make sure that the loop is closed
      assert(skin[j].front() == skin[j].back());
    }

    // edges are no longer needed
    result = delete_all_edges();
    assert(moab::MB_SUCCESS == result);

    /* Get the curves that are part of the surface. Use vectors to store all curve
    stuff so that we can remove curves from the set as they are zipped. */
    std::vector<int> curve_ids;
    int curve_id;
    std::vector<std::vector<moab::EntityHandle> > curves;
    for(std::vector<moab::EntityHandle>::iterator j=curve_sets.begin();
        j!=curve_sets.end(); j++) {

      // If a delete_curve, replace it with the keep_curve. This approach allows
      // for duplicates because we are using vectors instead of ranges. Note that
      // parent-child links also cannot store duplicate handles.
      moab::EntityHandle merged_curve;
      result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );
      assert(moab::MB_TAG_NOT_FOUND==result || moab::MB_SUCCESS==result);
      if(moab::MB_SUCCESS == result) *j = merged_curve;

      // do not add a curve if it contains nothing
      result = MBI()->tag_get_data( id_tag, &(*j), 1, &curve_id );
      assert(moab::MB_SUCCESS == result);
      std::cout << "  curve_id=" << curve_id << " handle=" << *j << std::endl;
      curve_ids.push_back(curve_id);
      std::vector<moab::EntityHandle> curve;
      result = arc::get_meshset( *j, curve );
      assert(moab::MB_SUCCESS == result);
      curves.push_back( curve );
    }

    // Keep zipping loops until each is either zipped or failed. This function
    // returns only after all loops are zipped or a failure occurs.
    while(!skin.empty()) {
      if(moab::MB_SUCCESS != result) {
        std::cout << "SURFACE_ZIPPING_FAILURE: could not zip surf_id=" << surf_id << std::endl;
      }
    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode test_edges()
{
  moab::ErrorCode result;
  moab::Range edges;
  result = MBI()->get_entities_by_dimension( 0, 1, edges );
  assert(moab::MB_SUCCESS == result);
  MBI()->list_entities( edges );
  return moab::MB_SUCCESS;
}


int main(int argc, char **argv)
{

  // ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time = clock(), prep_time, zip_time;
  if(2 > argc) {
    std::cout << "usage: do not use" << std::endl;
    std::cout << "./post_process <input_file>" << std::endl;
    return 1;
  }
  moab::ErrorCode result;
  std::string input_name = argv[1];

  // The root name does not have an extension
  std::string root_name = argv[1];
  int len = root_name.length();
  root_name.erase(len - 4);
  const double MERGE_TOL = 1e-3; // should this depend on FACET_TOL?

  // load the input file
  moab::EntityHandle input_meshset;
  result = MBI()->create_meshset( moab::MESHSET_SET, input_meshset );
  assert(moab::MB_SUCCESS == result);
  if(std::string::npos != input_name.find("h5m")) {
    result = MBI()->load_file( input_name.c_str(), &input_meshset );
    assert( moab::MB_SUCCESS == result );
  } else {
    std::cout << "invalid input file: must be h5m" << std::endl;
    return 1;
  }

  // create tags
  moab::Tag geom_tag, id_tag, sense_tag, normal_tag, merge_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, sizeof(int),
                                  moab::MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE, 0, 0 );
  assert( moab::MB_SUCCESS == result );
  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, sizeof(int),
                                  moab::MB_TYPE_INTEGER, id_tag,moab::MB_TAG_DENSE, 0, 0 );
  assert( moab::MB_SUCCESS == result );
  result = MBI()->tag_get_handle( "GEOM_SENSE_2", 2*sizeof(moab::EntityHandle),
                                  moab::MB_TYPE_HANDLE, sense_tag, moab::MB_TAG_DENSE, 0, 0 );
  assert( moab::MB_SUCCESS == result );
  result = MBI()->tag_get_handle( "NORMAL", sizeof(moab::CartVect),
                                  moab::MB_TYPE_OPAQUE, normal_tag, moab::MB_TAG_DENSE, 0, 0 );
  assert( moab::MB_SUCCESS == result );
  result = MBI()->tag_get_handle( "MERGE", sizeof(moab::EntityHandle),
                                  moab::MB_TYPE_HANDLE, merge_tag, moab::MB_TAG_SPARSE, 0, 0 );
  assert( moab::MB_SUCCESS == result );

  // get all geometry sets
  moab::Range geom_sets[4];
  for(unsigned dim=0; dim<4; dim++) {
    void *val[] = {&dim};
    result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &geom_tag,
             val, 1, geom_sets[dim] );
    assert(moab::MB_SUCCESS == result);
    for(moab::Range::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
      unsigned int options;
      result = MBI()->get_meshset_options(*i, options );
      assert(moab::MB_SUCCESS == result);

      // if options are wrong change them
      if(dim==1) {
        if( !(moab::MESHSET_TRACK_OWNER&options) || !(moab::MESHSET_ORDERED&options) ) {
          result = MBI()->set_meshset_options(*i, moab::MESHSET_TRACK_OWNER|moab::MESHSET_ORDERED);
          assert(moab::MB_SUCCESS == result);
        }
      } else {
        if( !(moab::MESHSET_TRACK_OWNER&options) ) {
          result = MBI()->set_meshset_options(*i, moab::MESHSET_TRACK_OWNER);
          assert(moab::MB_SUCCESS == result);
        }
      }
    }
  }
  std::cout << geom_sets[3].size() << " volumes, "
            << geom_sets[2].size() << " surfaces, and "
            << geom_sets[1].size() << " curves" << std::endl;

  result = cleanup::delete_small_edges(geom_sets[2], MERGE_TOL);
  assert(moab::MB_SUCCESS == result);

  std::string output_filename = root_name + "_tri.h5m";
  result = MBI()->write_mesh( output_filename.c_str() );
  if (moab::MB_SUCCESS != result) std::cout << "result= " << result << std::endl;
  assert(moab::MB_SUCCESS == result);

  zip_time = clock();
  std::cout << "zipping took " << (double) (zip_time-prep_time)/CLOCKS_PER_SEC
            << " sec." << std::endl;

  return 0;
}

moab::Interface *MBI()
{
  static moab::Core instance;
  return &instance;
}

