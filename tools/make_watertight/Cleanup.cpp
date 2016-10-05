#include <iostream>

// moab includes
#include "moab/OrientedBoxTreeTool.hpp"

#include "Cleanup.hpp"

// The obbtrees are no longer valid because the triangles have been altered.
//  -Surface and volume sets are tagged with tags holding the obb tree
//   root handles.
//  -Surface/volume set handles are added to the root meshset.
// Somehow, delete the old tree without deleting the
// surface and volume sets, then build a new tree.
moab::ErrorCode Cleanup::remove_obb_tree(bool verbose)
{
  moab::ErrorCode result;
  moab::Range obb_entities;
  moab::Tag obbTag;
  result = MBI()->tag_get_handle( "OBB_TREE", sizeof(moab::EntityHandle),
                                  moab::MB_TYPE_HANDLE, obbTag, moab::MB_TAG_DENSE, NULL, 0 );
  if(verbose) {
    MB_CHK_SET_ERR(result,"could not get OBB tree handle");
  } else {
    MB_CHK_SET_ERR(result,"");
  }
  // This gets the surface/volume sets. I don't want to delete the sets.
  // I want to remove the obbTag that contains the tree root handle and
  // delete the tree.
  result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET,
           &obbTag, 0, 1, obb_entities );
  assert(moab::MB_SUCCESS == result);
  std::cout << "  found " << obb_entities.size() << " OBB entities" << std::endl;

  // find tree roots
  moab::Range trees;
  moab::OrientedBoxTreeTool tool( MBI() );
  moab::Tag rootTag;
  for(moab::Range::iterator i=obb_entities.begin(); i!=obb_entities.end(); i++) {
    moab::EntityHandle root;
    result = MBI()->tag_get_data( obbTag, &(*i), 1, &root );
    MB_CHK_SET_ERR(result, "coule not get OBB tree data");
    tool.delete_tree( root );
  }
  result = MBI()->tag_delete( obbTag ); // use this for DENSE tags
  assert(moab::MB_SUCCESS == result);

  bool created;
  result = MBI()->tag_get_handle ( "OBB", sizeof(double),
                                   moab::MB_TYPE_DOUBLE, rootTag, moab::MB_TAG_SPARSE, 0, &created);
  assert(moab::MB_SUCCESS==result || moab::MB_ALREADY_ALLOCATED==result);
  // Were all of the trees deleted? Perhaps some of the roots we found were
  // child roots that got deleted with their parents.
  trees.clear();
  result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &rootTag,
           NULL, 1, trees );
  assert(moab::MB_SUCCESS == result);
  std::cout << "  " << trees.size() << " OBB tree(s) contained in file" << std::endl;
  return moab::MB_SUCCESS;
}

moab::ErrorCode Cleanup::delete_small_edge_and_tris( const moab::EntityHandle vert0,
    moab::EntityHandle &vert1,
    const double tol )
{
  // If the verts are the same, this is not meaningful.
  if(vert0 == vert1) return moab::MB_SUCCESS;
  moab::ErrorCode result;

  // If the edge is small, delete it and the adjacent tris.
  if(tol > gen->dist_between_verts(vert0, vert1)) {
    // get tris to delete
    moab::Range tris;
    moab::EntityHandle verts[2] = {vert0, vert1};
    result = MBI()->get_adjacencies( verts, 2, 2, false, tris );
    assert(moab::MB_SUCCESS == result);
    result = MBI()->delete_entities( tris );
    assert(moab::MB_SUCCESS == result);
    std::cout << "delete_small_edge_and_tris: deleted " << tris.size()
              << " tris." << std::endl;
    // now merge the verts, keeping the first one
    // IN FUTURE, AVERAGE THE LOCATIONS???????????
    result = MBI()->merge_entities( vert0, vert1, false, true);
    assert(moab::MB_SUCCESS == result);
    vert1 = vert0;
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Cleanup::delete_small_edges(const moab::Range &surfaces, const double FACET_TOL)
{
  // PROBLEM: THIS IS INVALID BECAUSE TRIS CAN HAVE LONG EDGES BUT
  // SMALL AREA. All three pts are in a line. This is the nature of
  // faceting vs. meshing.
  /* Remove small triangles by removing edges that are too small.
  Remove small edges by merging their endpoints together, creating
  degenerate triangles. Delete the degenerate triangles. */
  moab::ErrorCode result;
  for(moab::Range::const_iterator i=surfaces.begin(); i!=surfaces.end(); i++) {
    std::cout << "surf_id=" << gen->geom_id_by_handle(*i) << std::endl;

    // get all tris
    moab::Range tris;
    result = MBI()->get_entities_by_type( *i, moab::MBTRI, tris );
    assert(moab::MB_SUCCESS == result);

    // Check to ensure there area no degenerate tris
    for(moab::Range::iterator j=tris.begin(); j!=tris.end(); j++) {
      // get endpts
      const moab::EntityHandle *endpts;
      int n_verts;
      result = MBI()->get_connectivity( *j, endpts, n_verts);
      assert(moab::MB_SUCCESS == result);
      assert(3 == n_verts);
      assert( endpts[0]!=endpts[1] && endpts[1]!=endpts[2] );
    }


    // get the skin first, because my find_skin does not check before creating edges.
    moab::Range skin_edges;
    //result = gen->find_skin( tris, 1, skin_edges, false );
    moab::Skinner tool(MBI());
    result = tool.find_skin( 0 , tris, 1, skin_edges, false );
    assert(moab::MB_SUCCESS == result);

    // create the edges
    moab::Range edges;
    result = MBI()->get_adjacencies( tris, 1, true, edges, moab::Interface::UNION );
    if(moab::MB_SUCCESS != result) {
      std::cout << "result=" << result << std::endl;
    }
    assert(moab::MB_SUCCESS == result);

    // get the internal edges
    moab::Range internal_edges = subtract(edges, skin_edges);

    for(moab::Range::iterator j=internal_edges.begin(); j!=internal_edges.end(); j++) {
      int n_internal_edges = internal_edges.size();
      std::cout << "edge=" << *j << std::endl;
      MBI()->list_entity( *j );
      assert(moab::MB_SUCCESS == result);

      // get endpts
      const moab::EntityHandle *endpts;
      int n_verts;
      result = MBI()->get_connectivity( *j, endpts, n_verts);
      assert(moab::MB_SUCCESS == result);
      assert(2 == n_verts);

      // does another edge exist w the same endpts? Why would it?
      moab::Range duplicate_edges;
      result = MBI()->get_adjacencies( endpts, 2, 1, true, duplicate_edges );
      assert(moab::MB_SUCCESS == result);
      if(1 < duplicate_edges.size()) MBI()->list_entities( duplicate_edges );
      assert(1 == duplicate_edges.size());

      // if the edge length is less than MERGE_TOL do nothing
      if(FACET_TOL < gen->dist_between_verts( endpts[0], endpts[1] )) continue;

      // quick check
      for(moab::Range::iterator k=internal_edges.begin(); k!=internal_edges.end(); k++) {
        const moab::EntityHandle *epts;
        int n_vts;
        result = MBI()->get_connectivity( *k, epts, n_vts);
        assert(moab::MB_SUCCESS == result);
        assert(2 == n_vts);
        // The skin edges/verts cannot be moved, therefore both endpoints cannot
        // be on the skin. If they are, continue.
        moab::Range adj_edges0;
        result = MBI()->get_adjacencies( &epts[0], 1, 1, true, adj_edges0 );
        assert(moab::MB_SUCCESS == result);
        if(3 > adj_edges0.size()) {
          std::cout << "adj_edges0.size()=" << adj_edges0.size()
                    << " epts[0]=" << epts[0] << std::endl;
          MBI()->list_entity( epts[0] );
          assert(moab::MB_SUCCESS == result);
        }
        assert(3 <= adj_edges0.size());
        moab::Range adj_skin_edges0 = intersect( adj_edges0, skin_edges );
        bool endpt0_is_skin;
        if(adj_skin_edges0.empty()) endpt0_is_skin = false;
        else endpt0_is_skin = true;

        moab::Range adj_edges1;
        result = MBI()->get_adjacencies( &epts[1], 1, 1, true, adj_edges1 );
        assert(moab::MB_SUCCESS == result);
        if(3 > adj_edges1.size()) {
          std::cout << "adj_edges1.size()=" << adj_edges1.size()
                    << " epts[1]=" << epts[1] << std::endl;
          MBI()->list_entity( epts[1] );
          assert(moab::MB_SUCCESS == result);
        }
        assert(3 <= adj_edges1.size());
      }


      // The skin edges/verts cannot be moved, therefore both endpoints cannot
      // be on the skin. If they are, continue.
      moab::Range adj_edges0;
      result = MBI()->get_adjacencies( &endpts[0], 1, 1, true, adj_edges0 );
      assert(moab::MB_SUCCESS == result);
      if(3 > adj_edges0.size()) {
        std::cout << "adj_edges0.size()=" << adj_edges0.size()
                  << " endpts[0]=" << endpts[0] << std::endl;
        MBI()->list_entity( endpts[0] );
        assert(moab::MB_SUCCESS == result);
      }
      assert(3 <= adj_edges0.size());
      moab::Range adj_skin_edges0 = intersect( adj_edges0, skin_edges );
      bool endpt0_is_skin;
      if(adj_skin_edges0.empty()) endpt0_is_skin = false;
      else endpt0_is_skin = true;

      moab::Range adj_edges1;
      result = MBI()->get_adjacencies( &endpts[1], 1, 1, true, adj_edges1 );
      assert(moab::MB_SUCCESS == result);
      if(3 > adj_edges1.size()) {
        std::cout << "adj_edges1.size()=" << adj_edges1.size()
                  << " endpts[1]=" << endpts[1] << std::endl;
        MBI()->list_entity( endpts[1] );
        assert(moab::MB_SUCCESS == result);
      }
      assert(3 <= adj_edges1.size());
      moab::Range adj_skin_edges1 = intersect( adj_edges1, skin_edges );
      bool endpt1_is_skin;
      if(adj_skin_edges1.empty()) endpt1_is_skin = false;
      else endpt1_is_skin = true;
      if(endpt0_is_skin && endpt1_is_skin) continue;

      // Keep the skin endpt, and delete the other endpt
      moab::EntityHandle keep_endpt, delete_endpt;
      if(endpt0_is_skin) {
        keep_endpt   = endpts[0];
        delete_endpt = endpts[1];
      } else {
        keep_endpt   = endpts[1];
        delete_endpt = endpts[0];
      }

      // get the adjacent tris
      std::vector<moab::EntityHandle> adj_tris;
      result = MBI()->get_adjacencies( &(*j), 1, 2, false, adj_tris );
      assert(moab::MB_SUCCESS == result);
      if(2 != adj_tris.size()) {
        std::cout << "adj_tris.size()=" << adj_tris.size() << std::endl;
        for(unsigned int i=0; i<adj_tris.size(); i++) gen->print_triangle( adj_tris[i], true );
      }
      assert(2 == adj_tris.size());

      // When merging away an edge, a tri and 2 edges will be deleted.
      // Get each triangle's edge other edge the will be deleted.
      moab::Range tri0_delete_edge_verts;
      result = MBI()->get_adjacencies( &adj_tris[0], 1, 0, true, tri0_delete_edge_verts );
      assert(moab::MB_SUCCESS == result);
      assert(3 == tri0_delete_edge_verts.size());
      tri0_delete_edge_verts.erase( keep_endpt );
      moab::Range tri0_delete_edge;
      result = MBI()->get_adjacencies( tri0_delete_edge_verts, 1, true, tri0_delete_edge );
      assert(moab::MB_SUCCESS == result);
      assert(1 == tri0_delete_edge.size());

      moab::Range tri1_delete_edge_verts;
      result = MBI()->get_adjacencies( &adj_tris[1], 1, 0, true, tri1_delete_edge_verts );
      assert(moab::MB_SUCCESS == result);
      assert(3 == tri1_delete_edge_verts.size());
      tri1_delete_edge_verts.erase( keep_endpt );
      moab::Range tri1_delete_edge;
      result = MBI()->get_adjacencies( tri1_delete_edge_verts, 1, true, tri1_delete_edge );
      assert(moab::MB_SUCCESS == result);
      assert(1 == tri1_delete_edge.size());

      // *********************************************************************
      // Test to see if the merge would create inverted tris. Several special
      // cases to avoid all result in inverted tris.
      // *********************************************************************
      // get all the tris adjacent to the point the will be moved.
      moab::Range altered_tris;
      result = MBI()->get_adjacencies( &delete_endpt, 1, 2, false, altered_tris );
      assert(moab::MB_SUCCESS == result);
      bool inverted_tri = false;
      for(moab::Range::const_iterator k=altered_tris.begin(); k!=altered_tris.end(); ++k) {
        const moab::EntityHandle *conn;
        int n_verts;
        result = MBI()->get_connectivity( *k, conn, n_verts );
        assert(moab::MB_SUCCESS == result);
        assert(3 == tris.size());
        moab::EntityHandle new_conn[3];
        for(unsigned int i=0; i<3; ++i) {
          new_conn[i] = (conn[i]==delete_endpt) ? keep_endpt : conn[i];
        }
        double area;
        result = gen->triangle_area( new_conn, area );
        assert(moab::MB_SUCCESS == result);
        if(0 > area) {
          std::cout << "inverted tri detected, area=" << area << std::endl;
          inverted_tri = true;
          break;
        }
      }
      if(inverted_tri) continue;

      // If we got here, then the merge will occur. Delete the edge the will be
      // made degenerate and the edge that will be made duplicate.
      std::cout << "A merge will occur" << std::endl;
      gen->print_triangle( adj_tris[0], true );
      gen->print_triangle( adj_tris[1], true );
      internal_edges.erase( tri0_delete_edge.front() );
      internal_edges.erase( tri1_delete_edge.front() );
      std::cout << "merged verts=" << keep_endpt << " " << delete_endpt << std::endl;
      MBI()->list_entity( keep_endpt );
      MBI()->list_entity( delete_endpt );
      result = MBI()->merge_entities( keep_endpt, delete_endpt, false, true );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->delete_entities( tri0_delete_edge );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->delete_entities( tri1_delete_edge );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->delete_entities( &(*j), 1 );
      assert(moab::MB_SUCCESS == result);
      std::cout << "deleted edges=" << *j << " " << tri0_delete_edge.front()
                << " " << tri1_delete_edge.front() << std::endl;

      // delete degenerate tris
      result = MBI()->delete_entities( &adj_tris[0], 2 );
      assert(moab::MB_SUCCESS == result);
      MBI()->list_entity( keep_endpt );

      // remove the edge from the range
      j = internal_edges.erase(*j) - 1;
      std::cout << "next iter=" << *j << std::endl;

      moab::Range new_tris;
      result = MBI()->get_entities_by_type( *i, moab::MBTRI, new_tris );
      assert(moab::MB_SUCCESS == result);
      moab::Range new_skin_edges;
      result = tool.find_skin( *i, new_tris, 1, new_skin_edges, false );
      assert(moab::MB_SUCCESS == result);
      assert(skin_edges.size() == new_skin_edges.size());
      for(unsigned int k=0; k<skin_edges.size(); k++) {
        if(skin_edges[k] != new_skin_edges[k]) {
          MBI()->list_entity( skin_edges[k] );
          MBI()->list_entity( new_skin_edges[k] );
        }
        assert(skin_edges[k] == new_skin_edges[k]);
      }
      assert(n_internal_edges = internal_edges.size()+3);
    }

    // cleanup edges
    result = MBI()->get_entities_by_type( 0, moab::MBEDGE, edges );
    assert(moab::MB_SUCCESS == result);
    result = MBI()->delete_entities( edges );
    assert(moab::MB_SUCCESS == result);

  }
  return moab::MB_SUCCESS;
}

// Lots of edges have been created but are no longer needed.
// Delete edges that are not in curves. These should be the only edges
// that remain. This incredibly speeds up the watertight_check tool (100x?).
moab::ErrorCode Cleanup::cleanup_edges( moab::Range curve_meshsets )
{
  moab::ErrorCode result;
  moab::Range edges, edges_to_keep;
  for(moab::Range::iterator i=curve_meshsets.begin(); i!=curve_meshsets.end(); i++) {
    result = MBI()->get_entities_by_dimension( *i, 1, edges );
    assert(moab::MB_SUCCESS == result);
    edges_to_keep.merge( edges );
  }

  moab::Range all_edges;
  result = MBI()->get_entities_by_dimension( 0, 1, all_edges );
  assert(moab::MB_SUCCESS == result);

  // delete the edges that are not in curves.
  moab::Range edges_to_delete = subtract( all_edges, edges_to_keep );
  std::cout << "deleting " << edges_to_delete.size() << " unused edges" << std::endl;
  result = MBI()->delete_entities( edges_to_delete );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

