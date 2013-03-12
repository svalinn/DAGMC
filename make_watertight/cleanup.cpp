#include <iostream>
#include "cleanup.hpp"
//#include "MBAdaptiveKDTree.hpp"
#include "MBOrientedBoxTreeTool.hpp"

namespace cleanup {
  // The obbtrees are no longer valid because the triangles have been altered.
  //  -Surface and volume sets are tagged with tags holding the obb tree
  //   root handles.
  //  -Surface/volume set handles are added to the root meshset.
  // Somehow, delete the old tree without deleting the
  // surface and volume sets, then build a new tree.
  MBErrorCode remove_obb_tree() {
    MBErrorCode result;
    MBRange obb_entities;
    MBTag obbTag;
    result = MBI()->tag_create( "OBB_TREE", sizeof(MBEntityHandle),
	       		  MB_TAG_DENSE, MB_TYPE_HANDLE, obbTag, NULL, true );
    assert(MB_SUCCESS == result);
    // This gets the surface/volume sets. I don't want to delete the sets.
    // I want to remove the obbTag that contains the tree root handle and
    // delete the tree.
    result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET,
						    &obbTag, 0, 1, obb_entities );
    assert(MB_SUCCESS == result);
    std::cout << "  found " << obb_entities.size() << " OBB entities" << std::endl;
    //gen::print_range( obb_entities );
    //result = MBI()->delete_entities( obb_entities );
 
    // find tree roots
    MBRange trees;
    MBOrientedBoxTreeTool tool( MBI() );
    MBTag rootTag;
    for(MBRange::iterator i=obb_entities.begin(); i!=obb_entities.end(); i++) {
      MBEntityHandle root;
      result = MBI()->tag_get_data( obbTag, &(*i), 1, &root );
      //if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
      //assert(MB_SUCCESS == result);
      tool.delete_tree( root );
    }
    result = MBI()->tag_delete( obbTag ); // use this for DENSE tags
    assert(MB_SUCCESS == result);


    result = MBI()->tag_create( "OBB", sizeof(double), MB_TAG_SPARSE,
     				  MB_TYPE_DOUBLE, rootTag, 0, false);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
    /*    result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &rootTag, 
                                                   NULL, 1, trees );
    if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
    assert(MB_SUCCESS == result);
    //tool.find_all_trees( trees );
    std::cout << trees.size() << " tree(s) contained in file" << std::endl;
    //gen::print_range( trees );
  
    // delete the trees
    for (MBRange::iterator i = trees.begin(); i != trees.end(); ++i) {
      result = tool.delete_tree( *i );
      //if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
      //assert(MB_SUCCESS == result);
    }
    */ 
    // Were all of the trees deleted? Perhaps some of the roots we found were
    // child roots that got deleted with their parents.
    trees.clear();
    result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &rootTag,
						    NULL, 1, trees );
    assert(MB_SUCCESS == result);
    std::cout << "  " << trees.size() << " OBB tree(s) contained in file" << std::endl;
    return MB_SUCCESS;  
  }

  MBErrorCode delete_small_edge_and_tris( const MBEntityHandle vert0, 
                                          MBEntityHandle &vert1,
                                          const double tol ) {
    // If the verts are the same, this is not meaningful.
    if(vert0 == vert1) return MB_SUCCESS;
    MBErrorCode result;

    // If the edge is small, delete it and the adjacent tris.
    if(tol > gen::dist_between_verts(vert0, vert1)) {
      // get tris to delete
      MBRange tris;              
      MBEntityHandle verts[2] = {vert0, vert1};                      
      result = MBI()->get_adjacencies( verts, 2, 2, false, tris );                   
      assert(MB_SUCCESS == result);
      result = MBI()->delete_entities( tris );
      assert(MB_SUCCESS == result);
      std::cout << "delete_small_edge_and_tris: deleted " << tris.size() 
                << " tris." << std::endl;
      // now merge the verts, keeping the first one
      // IN FUTURE, AVERAGE THE LOCATIONS???????????
      result = MBI()->merge_entities( vert0, vert1, false, true);
      assert(MB_SUCCESS == result);
      vert1 = vert0;
    }
    return MB_SUCCESS;
  }

  MBErrorCode delete_small_edges(const MBRange &surfaces, const double FACET_TOL) {
    // PROBLEM: THIS IS INVALID BECAUSE TRIS CAN HAVE LONG EDGES BUT
    // SMALL AREA. All three pts are in a line. This is the nature of
    // faceting vs. meshing.
    /* Remove small triangles by removing edges that are too small. 
    Remove small edges by merging their endpoints together, creating
    degenerate triangles. Delete the degenerate triangles. */
    MBErrorCode result;
    for(MBRange::const_iterator i=surfaces.begin(); i!=surfaces.end(); i++) {
      std::cout << "surf_id=" << gen::geom_id_by_handle(*i) << std::endl;

      // get all tris
      MBRange tris;
      result = MBI()->get_entities_by_type( *i, MBTRI, tris );
      assert(MB_SUCCESS == result);

      // Check to ensure there area no degenerate tris
      for(MBRange::iterator j=tris.begin(); j!=tris.end(); j++) {
        // get endpts
        const MBEntityHandle *endpts;                                
        int n_verts;                                      
        result = MBI()->get_connectivity( *j, endpts, n_verts);
        assert(MB_SUCCESS == result);
        assert(3 == n_verts);
        assert( endpts[0]!=endpts[1] && endpts[1]!=endpts[2] );
      }


      // get the skin first, because my find_skin does not check before creating edges.
      MBRange skin_edges;
      //result = gen::find_skin( tris, 1, skin_edges, false );
      MBSkinner tool(MBI());
      result = tool.find_skin( tris, 1, skin_edges, false );
      assert(MB_SUCCESS == result);

      // create the edges
      MBRange edges;
      result = MBI()->get_adjacencies( tris, 1, true, edges, MBInterface::UNION );
      if(MB_SUCCESS != result) {
	std::cout << "result=" << result << std::endl;
      }
      assert(MB_SUCCESS == result);

      // get the internal edges
      MBRange internal_edges = subtract(edges, skin_edges);

      for(MBRange::iterator j=internal_edges.begin(); j!=internal_edges.end(); j++) {
        int n_internal_edges = internal_edges.size();
	std::cout << "edge=" << *j << std::endl;
        MBI()->list_entity( *j );
        assert(MB_SUCCESS == result);

        // get endpts
        const MBEntityHandle *endpts;                                
        int n_verts;                                      
        result = MBI()->get_connectivity( *j, endpts, n_verts);
        assert(MB_SUCCESS == result);
        assert(2 == n_verts);

        // does another edge exist w the same endpts? Why would it?
        MBRange duplicate_edges;
        result = MBI()->get_adjacencies( endpts, 2, 1, true, duplicate_edges );
        assert(MB_SUCCESS == result);
        if(1 < duplicate_edges.size()) MBI()->list_entities( duplicate_edges );
        assert(1 == duplicate_edges.size());

        // if the edge length is less than MERGE_TOL do nothing
        if(FACET_TOL < gen::dist_between_verts( endpts[0], endpts[1] )) continue; 
 
        // quick check
        for(MBRange::iterator k=internal_edges.begin(); k!=internal_edges.end(); k++) {
          const MBEntityHandle *epts;                                
          int n_vts;                                      
          result = MBI()->get_connectivity( *k, epts, n_vts);
          assert(MB_SUCCESS == result);
          assert(2 == n_vts);
	  // The skin edges/verts cannot be moved, therefore both endpoints cannot 
	  // be on the skin. If they are, continue.
	  MBRange adj_edges0;
	  result = MBI()->get_adjacencies( &epts[0], 1, 1, true, adj_edges0 );
	  assert(MB_SUCCESS == result);
	  if(3 > adj_edges0.size()) {
	    std::cout << "adj_edges0.size()=" << adj_edges0.size() 
		      << " epts[0]=" << epts[0] << std::endl;
	    MBI()->list_entity( epts[0] );
	    //MBI()->write_mesh( "test_output.h5m" );
	    assert(MB_SUCCESS == result);
	  }
	  assert(3 <= adj_edges0.size());
	  MBRange adj_skin_edges0 = intersect( adj_edges0, skin_edges );
	  bool endpt0_is_skin;
	  if(adj_skin_edges0.empty()) endpt0_is_skin = false;
	  else endpt0_is_skin = true;

	  MBRange adj_edges1;
	  result = MBI()->get_adjacencies( &epts[1], 1, 1, true, adj_edges1 );
	  assert(MB_SUCCESS == result);
	  if(3 > adj_edges1.size()) {
	    std::cout << "adj_edges1.size()=" << adj_edges1.size() 
		      << " epts[1]=" << epts[1] << std::endl;
	    MBI()->list_entity( epts[1] );
	    //MBI()->write_mesh( "test_output.h5m" );
	    assert(MB_SUCCESS == result);
	  }
	  assert(3 <= adj_edges1.size());
        }


        // The skin edges/verts cannot be moved, therefore both endpoints cannot 
        // be on the skin. If they are, continue.
        MBRange adj_edges0;
        result = MBI()->get_adjacencies( &endpts[0], 1, 1, true, adj_edges0 );
        assert(MB_SUCCESS == result);
        if(3 > adj_edges0.size()) {
          std::cout << "adj_edges0.size()=" << adj_edges0.size() 
                    << " endpts[0]=" << endpts[0] << std::endl;
          MBI()->list_entity( endpts[0] );
          //MBI()->write_mesh( "test_output.h5m" );
          assert(MB_SUCCESS == result);
        }
        assert(3 <= adj_edges0.size());
        MBRange adj_skin_edges0 = intersect( adj_edges0, skin_edges );
        bool endpt0_is_skin;
        if(adj_skin_edges0.empty()) endpt0_is_skin = false;
        else endpt0_is_skin = true;

        MBRange adj_edges1;
        result = MBI()->get_adjacencies( &endpts[1], 1, 1, true, adj_edges1 );
        assert(MB_SUCCESS == result);
        if(3 > adj_edges1.size()) {
          std::cout << "adj_edges1.size()=" << adj_edges1.size() 
                    << " endpts[1]=" << endpts[1] << std::endl;
          MBI()->list_entity( endpts[1] );
          //MBI()->write_mesh( "test_output.h5m" );
          assert(MB_SUCCESS == result);
        }
        assert(3 <= adj_edges1.size());
        MBRange adj_skin_edges1 = intersect( adj_edges1, skin_edges );
        bool endpt1_is_skin;
        if(adj_skin_edges1.empty()) endpt1_is_skin = false;
        else endpt1_is_skin = true;
        if(endpt0_is_skin && endpt1_is_skin) continue;
        
        // Keep the skin endpt, and delete the other endpt
        MBEntityHandle keep_endpt, delete_endpt;
        if(endpt0_is_skin) {
          keep_endpt   = endpts[0];
          delete_endpt = endpts[1];
        } else {
          keep_endpt   = endpts[1];
          delete_endpt = endpts[0];
        }

        // get the adjacent tris
	std::vector<MBEntityHandle> adj_tris;
        result = MBI()->get_adjacencies( &(*j), 1, 2, false, adj_tris );
        assert(MB_SUCCESS == result);
	        if(2 != adj_tris.size()) {
	std::cout << "adj_tris.size()=" << adj_tris.size() << std::endl;
	for(unsigned int i=0; i<adj_tris.size(); i++) gen::print_triangle( adj_tris[i], true );
        } 
        assert(2 == adj_tris.size());
 
        // When merging away an edge, a tri and 2 edges will be deleted.
        // Get each triangle's edge other edge the will be deleted.
        MBRange tri0_delete_edge_verts;
        result = MBI()->get_adjacencies( &adj_tris[0], 1, 0, true, tri0_delete_edge_verts );
        assert(MB_SUCCESS == result);
        assert(3 == tri0_delete_edge_verts.size());
        tri0_delete_edge_verts.erase( keep_endpt );
        MBRange tri0_delete_edge;
        result = MBI()->get_adjacencies( tri0_delete_edge_verts, 1, true, tri0_delete_edge );
        assert(MB_SUCCESS == result);
        assert(1 == tri0_delete_edge.size());      
 
        MBRange tri1_delete_edge_verts;
        result = MBI()->get_adjacencies( &adj_tris[1], 1, 0, true, tri1_delete_edge_verts );
        assert(MB_SUCCESS == result);
        assert(3 == tri1_delete_edge_verts.size());
        tri1_delete_edge_verts.erase( keep_endpt );
        MBRange tri1_delete_edge;
        result = MBI()->get_adjacencies( tri1_delete_edge_verts, 1, true, tri1_delete_edge );
        assert(MB_SUCCESS == result);
        assert(1 == tri1_delete_edge.size());      

        // When an edge is merged, it will be deleted and its to adjacent tris
        // will be deleted because they are degenerate. We cannot alter the skin.
        // How many skin edges does tri0 have?
	/*        MBRange tri_edges;
        result = MBI()->get_adjacencies( &adj_tris[0], 1, 1, false, tri_edges );
        assert(MB_SUCCESS == result);
        assert(3 == tri_edges.size());
        MBRange tri0_internal_edges = intersect(tri_edges, internal_edges);

        // Cannot merge the edge away if the tri has more than one skin edge.
        // Otherwise we would delete a skin edge. We already know the edges in 
        // edge_set are not on the skin.
        if(2 > tri0_internal_edges.size()) continue;

        // check the other tri
        tri_edges.clear();
        result = MBI()->get_adjacencies( &adj_tris[1], 1, 1, false, tri_edges );
        assert(MB_SUCCESS == result);
        assert(3 == tri_edges.size());
        MBRange tri1_internal_edges = intersect(tri_edges, internal_edges);
        if(2 > tri1_internal_edges.size()) continue;

        // Check to make sure that the internal edges are not on the skin
        MBRange temp;
        temp = intersect( tri0_internal_edges, skin_edges );
        assert(temp.empty());
        temp = intersect( tri1_internal_edges, skin_edges );
        assert(temp.empty());

        // We know that the edge will be merged away. Find the keep_vert and
        // delete_vert. The delete_vert should never be a skin vertex because the
        // skin must not move.
        MBRange delete_vert;
        result = MBI()->get_adjacencies( tri0_internal_edges, 0, false, delete_vert);
        assert(MB_SUCCESS == result);
        assert(1 == delete_vert.size());        
	*/

        // *********************************************************************
        // Test to see if the merge would create inverted tris. Several special
        // cases to avoid all result in inverted tris.
        // *********************************************************************
        // get all the tris adjacent to the point the will be moved.
        MBRange altered_tris;
        result = MBI()->get_adjacencies( &delete_endpt, 1, 2, false, altered_tris );
        assert(MB_SUCCESS == result);
        bool inverted_tri = false;
        for(MBRange::const_iterator k=altered_tris.begin(); k!=altered_tris.end(); ++k) {
          const MBEntityHandle *conn;
          int n_verts;
          result = MBI()->get_connectivity( *k, conn, n_verts );
          assert(MB_SUCCESS == result);
          assert(3 == tris.size());
          MBEntityHandle new_conn[3];
          for(unsigned int i=0; i<3; ++i) {
            new_conn[i] = (conn[i]==delete_endpt) ? keep_endpt : conn[i];
          }
          double area;
          result = gen::triangle_area( new_conn, area );
          assert(MB_SUCCESS == result);
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
      	gen::print_triangle( adj_tris[0], true );
	gen::print_triangle( adj_tris[1], true );
        //tri0_internal_edges.erase( *j );
        //tri1_internal_edges.erase( *j );
        internal_edges.erase( tri0_delete_edge.front() );
        internal_edges.erase( tri1_delete_edge.front() );
	std::cout << "merged verts=" << keep_endpt << " " << delete_endpt << std::endl;
	MBI()->list_entity( keep_endpt );
	MBI()->list_entity( delete_endpt );
        result = MBI()->merge_entities( keep_endpt, delete_endpt, false, true );          
	assert(MB_SUCCESS == result);
        result = MBI()->delete_entities( tri0_delete_edge );
        assert(MB_SUCCESS == result);
        result = MBI()->delete_entities( tri1_delete_edge );
        assert(MB_SUCCESS == result);        
        result = MBI()->delete_entities( &(*j), 1 );
        assert(MB_SUCCESS == result);
	std::cout << "deleted edges=" << *j << " " << tri0_delete_edge.front()
                  << " " << tri1_delete_edge.front() << std::endl;
	          
	// delete degenerate tris                          
	result = MBI()->delete_entities( &adj_tris[0], 2 );                 
	assert(MB_SUCCESS == result); 
	MBI()->list_entity( keep_endpt );

        // remove the edge from the range
        j = internal_edges.erase(*j) - 1; 
	std::cout << "next iter=" << *j << std::endl;        

        MBRange new_tris;
        result = MBI()->get_entities_by_type( *i, MBTRI, new_tris );
        assert(MB_SUCCESS == result);
        MBRange new_skin_edges;
        result = tool.find_skin( new_tris, 1, new_skin_edges, false );
        assert(MB_SUCCESS == result);
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
      result = MBI()->get_entities_by_type( 0, MBEDGE, edges );
      assert(MB_SUCCESS == result);
      result = MBI()->delete_entities( edges ); 
      assert(MB_SUCCESS == result); 
   
    }
    return MB_SUCCESS;
  } 
  
  // Lots of edges have been created but are no longer needed.
  // Delete edges that are not in curves. These should be the only edges
  // that remain. This incredibly speeds up the watertight_check tool (100x?).
  MBErrorCode cleanup_edges( MBRange curve_meshsets ) {
    MBErrorCode result;
    MBRange edges, edges_to_keep;
    for(MBRange::iterator i=curve_meshsets.begin(); i!=curve_meshsets.end(); i++) {
      result = MBI()->get_entities_by_dimension( *i, 1, edges );
      assert(MB_SUCCESS == result);
      edges_to_keep.merge( edges );
    }

    MBRange all_edges;
    result = MBI()->get_entities_by_dimension( 0, 1, all_edges );
    assert(MB_SUCCESS == result);

    // delete the edges that are not in curves.
    //MBRange edges_to_delete = all_edges.subtract( edges_to_keep );
    MBRange edges_to_delete = subtract( all_edges, edges_to_keep );
    std::cout << "deleting " << edges_to_delete.size() << " unused edges" << std::endl;
    result = MBI()->delete_entities( edges_to_delete );
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }
 
  
}
