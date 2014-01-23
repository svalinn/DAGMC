#include <iostream>
#include "zip.hpp"
#include "MBOrientedBoxTreeTool.hpp"

namespace zip {
  MBErrorCode t_joint( MBTag normal_tag,
                       const MBEntityHandle vert0,
                       const MBEntityHandle vert1,
                       const MBEntityHandle vert2 ) {
    struct triangles {
      MBEntityHandle before_tri;
      const MBEntityHandle *before;
      MBCartVect     before_norm;
      MBEntityHandle after0[3];
      MBEntityHandle after1[3];
      MBCartVect     after0_norm;
      MBCartVect     after1_norm;
      MBEntityHandle surf_set;
    };   

    // Get all of the old information before changing anything. 
    // This is important because once the
    // new connectivity is set stuff becomes stale.
    // get the edge

    // get endpoints of the edge
    MBErrorCode result;
    MBEntityHandle endpts[2] = { vert0, vert2 }; 
    MBRange tris;
    result = MBI()->get_adjacencies( endpts, 2, 2, true, tris );
    assert(MB_SUCCESS == result);
    //std::cout << "t_joint: tris.size()=" << tris.size() << std::endl;
    //MBI()->list_entities( tris );

    triangles joints[tris.size()];
    for(unsigned int i=0; i<tris.size(); i++) {
      joints[i].before_tri = tris[i];

      // Find the surface set that the tri is in.
      MBRange surf_sets;
      result = MBI()->get_adjacencies( &joints[i].before_tri, 1, 4, false, surf_sets);
      assert(MB_SUCCESS == result);
      //std::cout << "t_joint: " << surf_sets.size() << " surface sets found for triangle" 
      //        << std::endl;
  
      // Check to make sure we found a set     
      if(1 != surf_sets.size()) {
	std::cout << "    t_joint: " << surf_sets.size() << " surface sets found for triangle " 
                  << joints[i].before_tri << std::endl;
        assert(1 == surf_sets.size());
	//if(1!=surf_sets.size()) return MB_FAILURE;
      }
      joints[i].surf_set = surf_sets.front();
      //std::cout << "t_joint: surf id=" << gen::geom_id_by_handle( joints[i].surf_set )
      //        << std::endl;  
      //gen::print_triangle( joints[i].before_tri, false );

      // get old  connectivity
      int n_verts;
      result = MBI()->get_connectivity( joints[i].before_tri, joints[i].before, n_verts);
      if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
      assert(MB_SUCCESS == result);
      if(3 != n_verts) std::cout << "n_verts=" << n_verts << std::endl;
      assert(3 == n_verts);
   
      // test to make sure not degenerate
      //if(conn[0]==conn[1] || conn[1]==conn[2] || conn[2]==conn[0]) {
      if( gen::triangle_degenerate( joints[i].before_tri )) {
        std::cout << "    t_joint: degenerate input triangle" << std::endl;
        gen::print_triangle( joints[i].before_tri, false);
        return MB_FAILURE;
      }

      // make new connectivity
      for(int j=0; j<3; j++) {
	joints[i].after0[j] = (joints[i].before[j]==endpts[0]) ? vert1 : joints[i].before[j];
	joints[i].after1[j] = (joints[i].before[j]==endpts[1]) ? vert1 : joints[i].before[j]; 
      }

      // test to make sure not degenerate
      //if(conn0[0]==conn0[1] || conn0[1]==conn0[2] || conn0[2]==conn0[0]) {
      if(gen::triangle_degenerate( joints[i].after0[0], joints[i].after0[1], 
                                                        joints[i].after0[2])) {
	std::cout << "    t_joint: degenerate output triangle 1" << std::endl;	
        gen::print_triangle( joints[i].before_tri, false );
	//	return MB_FAILURE;
      }
      // test to make sure not degenerate
      //if(conn1[0]==conn1[1] || conn1[1]==conn1[2] || conn1[2]==conn1[0]) {
      if(gen::triangle_degenerate( joints[i].after1[0], joints[i].after1[1], 
                                                        joints[i].after1[2])) {
	std::cout << "    t_joint: degenerate output triangle 2" << std::endl;
        gen::print_triangle( joints[i].before_tri, false );
	//gen::print_triangle( *i, true );
	//return MB_FAILURE;
      }

      // set the new connectivity on the original triangle
      result = MBI()->set_connectivity( joints[i].before_tri, joints[i].after0, 3 );
      assert(MB_SUCCESS == result);
      // set the new connectivity on the new triangle
      MBEntityHandle new_tri;
      result = MBI()->create_element( MBTRI, joints[i].after1, 3, new_tri );
      assert(MB_SUCCESS == result);

      // copy the original normal to the new triangle
      MBCartVect normal;
      result = MBI()->tag_get_data( normal_tag, &joints[i].before_tri, 1, &normal);
      assert(MB_SUCCESS == result);
      result = MBI()->tag_set_data( normal_tag, &new_tri, 1, &normal);
      assert(MB_SUCCESS == result);
    
      // add the new triangle to the same surface set as the original
      result = MBI()->add_entities( joints[i].surf_set, &new_tri, 1);
      assert(MB_SUCCESS == result);

      // catch-all to remove degenerate tris
      result = zip::delete_degenerate_tris( joints[i].before_tri );
      if(gen::error(MB_SUCCESS!=result,"could not delete degenerate tri")) return result;
      result = zip::delete_degenerate_tris( new_tri );
      if(gen::error(MB_SUCCESS!=result,"could not delete degenerate tri")) return result;

      //gen::print_triangle( tri, false );
      //gen::print_triangle( new_tri, false );
    }
    return MB_SUCCESS;
  }

  // Delete degenerate triangles in the range.
  MBErrorCode delete_degenerate_tris( MBEntityHandle tri ) {
    MBErrorCode result;
    const MBEntityHandle *con;
    int n_verts;
    result = MBI()->get_connectivity( tri, con, n_verts);
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    if(con[0]==con[1] || con[1]==con[2] || con[2]==con[0]) {
      //std::cout << "delete_degenerate_tris: degenerate triangle=" << tri << " deleted" << std::endl;
      result = MBI()->delete_entities( &tri, 1 );
      assert(MB_SUCCESS == result);
    }
    return MB_SUCCESS;
  }
  MBErrorCode delete_degenerate_tris( MBRange tris ) {
    MBErrorCode result;
    for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
      result = delete_degenerate_tris( *i );
      assert(MB_SUCCESS == result);
    }
    return MB_SUCCESS;
  }

  MBErrorCode delete_adj_degenerate_tris( const MBEntityHandle adj_vert ) {
    // get the adjacent triangles
    MBErrorCode result;
    MBRange tris;
    result = MBI()->get_adjacencies( &adj_vert, 1, 2, false, tris );
    assert(MB_SUCCESS == result);
    result = delete_degenerate_tris( tris );
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }

  // Problem: MOAB can check for degenerate tris and delete them.
  // MOAB check the curves and arcs passed in to update the merged vertex, but
  // still ends up with degenerate edges. The curves that are in MOAB as sets
  // are updated but also contain degenerate edges due to merging.
  MBErrorCode merge_verts( const MBEntityHandle keep_vert, 
                           const MBEntityHandle delete_vert,
                           std::vector<MBEntityHandle> &arc0,
                           std::vector<MBEntityHandle> &arc1 ) {

    MBErrorCode rval;
    // first update the arcs with the keep_vert
    for(std::vector<MBEntityHandle>::iterator i=arc0.begin(); i!=arc0.end(); ++i) {
      if(delete_vert == *i) *i = keep_vert;
    }
    for(std::vector<MBEntityHandle>::iterator i=arc1.begin(); i!=arc1.end(); ++i) {
      if(delete_vert == *i) *i = keep_vert;
    }

    // UPDATE: This is slower than using the O(n) linear search above.
    // Let moab update adjacencies. Unless moab stores data is must be 
    // merge-updated manually to prevent stale handles.
    /*    MBEntityHandle arc0_set, arc1_set;
    rval = MBI()->create_meshset( MESHSET_TRACK_OWNER|MESHSET_ORDERED, arc0_set );
    if(gen::error(MB_SUCCESS!=rval,"creating arc0_set failed")) return rval;
    rval = MBI()->create_meshset( MESHSET_TRACK_OWNER|MESHSET_ORDERED, arc1_set );
    if(gen::error(MB_SUCCESS!=rval,"creating arc1_set failed")) return rval;
    rval = arc::set_meshset( arc0_set, arc0 );                                          
    if(gen::error(MB_SUCCESS!=rval,"setting arc0_set failed")) return rval; 
    rval = arc::set_meshset( arc1_set, arc1 );                                          
    if(gen::error(MB_SUCCESS!=rval,"setting arc1_set failed")) return rval; 
    */

    // get adjacent tris
    MBRange tris;
    MBEntityHandle verts[2]={keep_vert, delete_vert};
    rval = MBI()->get_adjacencies( verts, 2, 2, false, tris, MBInterface::UNION );
    if(gen::error(MB_SUCCESS!=rval,"getting adjacent tris failed")) return rval;
    //if(0 == tris.size()) {
    //  std::cout << "merge_verts: cannot find any triangles adjacent to vertices" << std::endl;
    //  return MB_ENTITY_NOT_FOUND;
    //}
    //gen::print_triangles( tris );

    // actually do the merge
    rval = MBI()->merge_entities( keep_vert, delete_vert, false, true );
    if(gen::error(MB_SUCCESS!=rval,"merge entities failed")) return rval;
    //std::cout << "  tris.size()=" << tris.size() << std::endl;

    // delete degenerate tris
    rval = delete_degenerate_tris( tris );
    if(gen::error(MB_SUCCESS!=rval,"deleting degenerate tris failed")) return rval;

    // Get the merge-updated arcs back.
    /*    rval = arc::get_meshset( arc0_set, arc0 ); 
    if(gen::error(MB_SUCCESS!=rval,"getting arc0 set failed")) return rval;
    rval = arc::get_meshset( arc1_set, arc1 ); 
    if(gen::error(MB_SUCCESS!=rval,"getting arc1 set failed")) return rval;
    rval = MBI()->delete_entities( &arc0_set, 1 );                   
    if(gen::error(MB_SUCCESS!=rval,"deleting arc0_set failed")) return rval;
    rval = MBI()->delete_entities( &arc1_set, 1 );                   
    if(gen::error(MB_SUCCESS!=rval,"deleting arc1_set failed")) return rval;
    */

    return MB_SUCCESS;
  }                           

  // Test to make sure the triangle normal vectors have not been inverted.
  MBErrorCode test_normals( const std::vector<MBCartVect> norms0,
                            const std::vector<MBCartVect> norms1,
                            std::vector<int> &inverted_tri_indices ) {
    assert(norms0.size() == norms1.size());
    for(unsigned int i=0; i<norms0.size(); i++) {
      MBErrorCode result = test_normals( norms0[i], norms1[i]);
      if(MB_SUCCESS != result) {
        //std::cout << "test_normals: failed on i=" << i << std::endl;
        inverted_tri_indices.push_back(i);
      }
    }
    return MB_SUCCESS;
  }
  MBErrorCode test_normals( const MBCartVect norm0, const MBCartVect norm1 ) {
    if(0 > norm0 % norm1) {
      //std::cout << "test_normals: tri is inverted, dot product=" 
      //          << norm0 % norm1 << std::endl;
      return MB_FAILURE;
    } else {
      return MB_SUCCESS;
    }
  }

  /* Accepts a range of inverted tris. Refacets affected surface so that no tris
     are inverted. */
  MBErrorCode remove_inverted_tris(MBTag normal_tag, MBRange tris, const bool debug ) {
      
    MBErrorCode result;
    bool failures_occur = false;
    while(!tris.empty()) {

      /* Get a group of triangles to re-facet. They must be adjacent to each other
	 and in the same surface. */
      MBRange tris_to_refacet;
      tris_to_refacet.insert( tris.front() );
      MBRange surf_set;
      result = MBI()->get_adjacencies( tris_to_refacet, 4, false, surf_set );
      assert(MB_SUCCESS == result);
      if(1 != surf_set.size()) {
	std::cout << "remove_inverted_tris: tri is in " << surf_set.size() 
                  << " surfaces" << std::endl;
        return MB_FAILURE;
      }

      // get all tris in the surface
      MBRange surf_tris;
      result = MBI()->get_entities_by_type( surf_set.front(), MBTRI, surf_tris );
      assert(MB_SUCCESS == result); 

      /* Find all of the adjacent inverted triangles of the same surface. Keep
	 searching until a search returns no new triangles. */
      bool search_again = true;
      while(search_again) {

        // Here edges are being created. Remember to delete them. Outside of this
        // function. Skinning gets bogged down if unused MBEdges (from other 
        // surfaces) accumulate.
        MBRange tri_edges;
        result = MBI()->get_adjacencies( tris_to_refacet, 1, true, tri_edges,
                                         MBInterface::UNION );
        assert(MB_SUCCESS == result);
        MBRange connected_tris;
        result = MBI()->get_adjacencies( tri_edges, 2, false, connected_tris, 
                                         MBInterface::UNION );
        assert(MB_SUCCESS == result);
        result = MBI()->delete_entities( tri_edges );
        assert(MB_SUCCESS == result);
        MBRange tris_to_refacet2 = intersect( tris_to_refacet, connected_tris );
        tris_to_refacet2 = intersect( tris_to_refacet, surf_tris );

        if(tris_to_refacet.size() == tris_to_refacet2.size()) search_again = false;
        tris_to_refacet.swap( tris_to_refacet2 );
      }

      // Remove the inverted tris that will be refaceted.
      tris = subtract( tris, tris_to_refacet );

        // do edges already exist?
	MBRange temp;
          result = MBI()->get_entities_by_type(0, MBEDGE, temp );
          assert(MB_SUCCESS == result);
          if(!temp.empty()) MBI()->list_entities( temp );
	  assert(temp.empty());


      // keep enlarging patch until we have tried to refacet the entire surface
      int counter=0;
      while(true) {
        // do edges already exist?
	temp.clear();
          result = MBI()->get_entities_by_type(0, MBEDGE, temp );
          assert(MB_SUCCESS == result);
          if(!temp.empty()) MBI()->list_entities( temp );
	  assert(temp.empty());


        counter++;
        // Only try enlarging each patch a few times
        if(48 == counter) {
          failures_occur = true;
	  if(debug) std::cout << "remove_inverted_tris: ear clipping failed, counter="
		              << counter << std::endl;
	  break;
        }
        // THIS PROVIDES A BAD EXIT. MUST FIX

        // get the edges of the patch of inverted tris
	MBRange tri_edges;
	result = MBI()->get_adjacencies( tris_to_refacet, 1, true, tri_edges,
                                         MBInterface::UNION );
	assert(MB_SUCCESS == result);

	// get all adjacent tris to the patch of inverted tris in the surface
	MBRange adj_tris;
	result = MBI()->get_adjacencies( tri_edges, 2, false, adj_tris, 
                                         MBInterface::UNION );
	assert(MB_SUCCESS == result);
        result = MBI()->delete_entities( tri_edges );
        assert(MB_SUCCESS == result);
	tris_to_refacet = intersect( surf_tris, adj_tris );
        if(tris_to_refacet.empty()) continue;
	//gen::print_triangles( tris_to_refacet );    
	
	// get an area-weighted normal of the adj_tris
	MBCartVect plane_normal(0,0,0);
	//for(unsigned int i=0; i<tris_to_refacet.size(); i++) {
	for(MBRange::iterator i=tris_to_refacet.begin(); i!=tris_to_refacet.end(); i++) {
	  MBCartVect norm;
	  result = MBI()->tag_get_data( normal_tag, &(*i), 1, &norm);
	  assert(MB_SUCCESS == result);
	  double area;
          result = gen::triangle_area( *i, area );
          assert(MB_SUCCESS == result);
	  if(debug) std::cout << "norm=" << norm << " area=" << area << std::endl;
	  //plane_normal += norm*area;
	  plane_normal += norm;
	}
	plane_normal.normalize();

        // do edges already exist?
	temp.clear();
          result = MBI()->get_entities_by_type(0, MBEDGE, temp );
          assert(MB_SUCCESS == result);
          if(!temp.empty()) MBI()->list_entities( temp );
	  assert(temp.empty());
 
	// skin the tris
	MBRange unordered_edges;
	//MBSkinner tool(MBI());
	//result = tool.find_skin( tris_to_refacet, 1, unordered_edges, false );
	result = gen::find_skin( tris_to_refacet, 1, unordered_edges, false );
	assert(MB_SUCCESS == result);
        if(unordered_edges.empty()) {
        // do edges already exist?
          MBRange temp;
          result = MBI()->get_entities_by_type(0, MBEDGE, temp );
          assert(MB_SUCCESS == result);
          if(!temp.empty()) MBI()->list_entities( temp );
	  assert(temp.empty());
          continue;
        }

	//std::cout << "remove_inverted_tris: surf_id=" 
	//  << gen::geom_id_by_handle(surf_set.front()) << std::endl;
	//result = MBI()->list_entities( tris_to_refacet );
	//assert(MB_SUCCESS == result);

	// assemble into a polygon
	std::vector<MBEntityHandle> polygon_of_verts;
	result = arc::order_verts_by_edge( unordered_edges, polygon_of_verts );
	if(debug) gen::print_loop( polygon_of_verts ); 
	//assert(MB_SUCCESS == result);
	if(MB_SUCCESS != result) {
	  if(debug) std::cout << "remove_inverted_tris: couldn't order polygon by edge" << std::endl;
	  return MB_FAILURE;
	}

        // remember to remove edges
        result = MBI()->delete_entities( unordered_edges );
        assert(MB_SUCCESS == result);

	// remove the duplicate endpt
	polygon_of_verts.pop_back();

	// the polygon should have at least 3 verts
	if(3 > polygon_of_verts.size()) {
	  if(debug) std::cout << "remove_inverted_tris: polygon has too few points" << std::endl;
	  return MB_FAILURE;
	}

	// orient the polygon with the triangles (could be backwards)
	// get the first adjacent tri
	MBEntityHandle edge[2] = { polygon_of_verts[0], polygon_of_verts[1] };
	MBRange one_tri;
	result = MBI()->get_adjacencies( edge, 2, 2, false, one_tri );
	assert(MB_SUCCESS == result);
	one_tri = intersect( tris_to_refacet, one_tri );
	assert(1 == one_tri.size());
	const MBEntityHandle *conn;
	int n_conn;
	result = MBI()->get_connectivity( one_tri.front(), conn, n_conn );
	assert(MB_SUCCESS == result);
	assert(3 == n_conn);
	if( (edge[0]==conn[1] && edge[1]==conn[0]) ||
	    (edge[0]==conn[2] && edge[1]==conn[1]) ||
	    (edge[0]==conn[0] && edge[1]==conn[2]) ) {
	  reverse( polygon_of_verts.begin(), polygon_of_verts.end() );
	  if(debug) std::cout << "remove_inverted_tris: polygon reversed" << std::endl;
	}

	/* facet the polygon. Returns MB_FAILURE if it fails to facet the polygon. */
	MBRange new_tris;
	result = gen::ear_clip_polygon( polygon_of_verts, plane_normal, new_tris );

        // break if the refaceting is successful
	if(MB_SUCCESS == result) {
          // summarize tri area
          for(MBRange::iterator i=new_tris.begin(); i!=new_tris.end(); i++) {
            double area;
            result = gen::triangle_area( *i, area );
            assert(MB_SUCCESS == result);
	    if(debug) std::cout << "  new tri area=" << area << std::endl;
          }

  	  // check the new normals
	  std::vector<MBCartVect> new_normals;
	  result = gen::triangle_normals( new_tris, new_normals );
	  if(MB_SUCCESS != result) return result;

	  // test the new triangles
	  std::vector<int> inverted_tri_indices;
	  std::vector<MBCartVect> normals ( new_normals.size(), plane_normal );
	  result = zip::test_normals( normals, new_normals, inverted_tri_indices );
	  assert(MB_SUCCESS == result);
	  if(inverted_tri_indices.empty()) {
  	    // remove the tris that were re-faceted
            tris = subtract( tris, tris_to_refacet );
  	    result = MBI()->remove_entities( surf_set.front(), tris_to_refacet );
	    assert(MB_SUCCESS == result);
	    result = MBI()->delete_entities( tris_to_refacet ); 
	    assert(MB_SUCCESS == result);

	    // add the new tris to the surf set
	    result = MBI()->add_entities( surf_set.front(), new_tris );
	    assert(MB_SUCCESS == result);

            // put the new normals on the new tris
            result = gen::save_normals( new_tris, normal_tag );
            assert(MB_SUCCESS == result);
	    if(debug) std::cout << "remove_inverted_tris: success fixing a patch" << std::endl;
            break;
          }
        }

        // remember to delete the tris that were created from the failed ear clipping
        else {
          result = MBI()->delete_entities( new_tris );
          assert(MB_SUCCESS == result);
        }

        // If the entire surface could not be ear clipped, give up
        if (tris_to_refacet.size() == surf_tris.size()) {
	  if(debug) std::cout << "remove_inverted_tris: ear clipping entire surface failed"
			    << std::endl;
	  return MB_FAILURE;
	}

      } // loop until the entire surface has attempted to be refaceted
    }   // loop over each patch of inverted tris
   
    if(failures_occur) {
      if(debug) std::cout << "remove_inverted_facets: at least one failure occured" << std::endl;
      return MB_FAILURE;
    } else {
      return MB_SUCCESS;
    }
  }


    // we do not merge edges, just vert. check the verts
  MBErrorCode test_zipping(const double FACET_TOL,
                           const std::vector< std::vector<MBEntityHandle> > arcs ) {
      MBErrorCode result;

      // make sure each arc has the same number of edges
      for(unsigned int i=1; i<arcs.size(); i++) {
	if(arcs[0].size() != arcs[i].size()) {
	  std::cout << "The curve has " << arcs[0].size() << " edges but arc "
		    << i << " has " << arcs[i].size() << " edges." << std::endl;
	  gen::print_arcs( arcs );
	  return MB_FAILURE;
	}
      }
  
      // loop over every edge of the curve (first arc)
      for(unsigned int i=0; i<arcs[0].size()-1; i++) {
	// check for degenerate edge
	if(arcs[0][i] == arcs[0][i+1]) {
	  std::cout << "degenerate edge at pos " << i << " and " << i+1 << " with verts "
		    << arcs[0][i] << " and " << arcs[0][i+1] << std::endl;
	  return MB_FAILURE;
	}
      
	// check for edge of zero dist
	double d = gen::dist_between_verts( arcs[0][i], arcs[0][i+1] );
	if(FACET_TOL >= d) {
	  std::cout << "edge length=" << d << " betwee pos " << i << " and " << i+1
		    << " with verts " << arcs[0][i] << " and " << arcs[0][i+1] << std::endl;
	  return MB_FAILURE;
	}
     
	// loop over every arc
	for( unsigned int j=0; j<arcs.size(); j++) {      
	  // make sure vertices match
	  if(arcs[0][i]!=arcs[j][i] || arcs[0][i+1]!=arcs[j][i+1]) {
	    std::cout << "arc " << j << " vertices do not match curve vertices, pos= " 
		      << i << "/" << arcs[j].size() << std::endl;
	    return MB_FAILURE;
	  }
	}  	
        
	// make sure triangles have area
	MBRange tris;
	result = MBI()->get_adjacencies( &(arcs[0][i]), 2, 2, false, tris );
	assert(MB_SUCCESS == result);
	for(MBRange::iterator k=tris.begin(); k!=tris.end(); k++) {
	  // We know that there are not degenerate edges along the curve.
	  // Sometimes degenerate tris are created due to merging curve endpts.
	  // here we do not remove tri from the surf meshset, but we should
	  if( gen::triangle_degenerate(*k) ) {
	    //result = MBI()->delete_entities( &(*k), 1);
	    //assert(MB_SUCCESS == result);
	    std::cout << "  arc=" << 0 << " pos=" << i << " vert=" << arcs[0][i] 
		      << " degenerate triangle" << std::endl;
	    gen::print_triangle(*k, false);
	    //print_edge( edge );
	    //continue;
	    return MB_FAILURE;
	  }

	  double area;
          result = gen::triangle_area( *k, area );
          assert(MB_SUCCESS == result);
	  // I found a valid tri on a curve with only one edge (1e-5 long)
	  // that had an area of 1e-11.
	  if(1e-8 > area) {
	    std::cout << "    arc=" << 0 << " pos=" << i << " vert=" << arcs[0][i]
		      << " small triangle " << std::endl;
	    gen::print_triangle(*k, false);            
	    //print_edge( edge );
  	    gen::print_arcs( arcs );
	    //if(0.0 >= area) return MB_FAILURE;
	  } 
	}
      }
      return MB_SUCCESS;
    }          
  
  
  }
