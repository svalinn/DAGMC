#include <iostream>
#include "Zip.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/Skinner.hpp"

moab::ErrorCode Zip::t_joint( moab::Tag normal_tag,
                              const moab::EntityHandle vert0,
                              const moab::EntityHandle vert1,
                              const moab::EntityHandle vert2,
                              bool debug )
{

  // Get all of the old information before changing anything.
  // This is important because once the
  // new connectivity is set stuff becomes stale.
  // get the edge

  // get endpoints of the edge
  moab::ErrorCode result;
  moab::EntityHandle endpts[2] = { vert0, vert2 };
  moab::Range tris;
  result = MBI()->get_adjacencies( endpts, 2, 2, true, tris );
  assert(moab::MB_SUCCESS == result);

  std::vector<triangles> joints(tris.size());
  for(unsigned int i=0; i<tris.size(); i++) {
    joints[i].before_tri = tris[i];

    // Find the surface set that the tri is in.
    moab::Range surf_sets;
    result = MBI()->get_adjacencies( &joints[i].before_tri, 1, 4, false, surf_sets);
    assert(moab::MB_SUCCESS == result);

    // Check to make sure we found a set
    if(1 != surf_sets.size()) {
      if(debug) std::cout << "    t_joint: " << surf_sets.size() << " surface sets found for triangle "
                            << joints[i].before_tri << std::endl;
      assert(1 == surf_sets.size());
    }
    joints[i].surf_set = surf_sets.front();

    // get old  connectivity
    int n_verts;
    result = MBI()->get_connectivity( joints[i].before_tri, joints[i].before, n_verts);
    if(moab::MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
    assert(moab::MB_SUCCESS == result);
    if(3 != n_verts) std::cout << "n_verts=" << n_verts << std::endl;
    assert(3 == n_verts);

    // test to make sure not degenerate
    if( gen->triangle_degenerate( joints[i].before_tri ) && debug ) {
      std::cout << "    t_joint: degenerate input triangle" << std::endl;
      gen->print_triangle( joints[i].before_tri, false);
      return moab::MB_FAILURE;
    }

    // make new connectivity
    for(int j=0; j<3; j++) {
      joints[i].after0[j] = (joints[i].before[j]==endpts[0]) ? vert1 : joints[i].before[j];
      joints[i].after1[j] = (joints[i].before[j]==endpts[1]) ? vert1 : joints[i].before[j];
    }

    // test to make sure not degenerate
    if(gen->triangle_degenerate( joints[i].after0[0], joints[i].after0[1],
                                 joints[i].after0[2]) && debug ) {
      std::cout << "    t_joint: degenerate output triangle 1" << std::endl;
      gen->print_triangle( joints[i].before_tri, false );
    }
    // test to make sure not degenerate
    if(gen->triangle_degenerate( joints[i].after1[0], joints[i].after1[1],
                                 joints[i].after1[2]) && debug) {
      std::cout << "    t_joint: degenerate output triangle 2" << std::endl;
      gen->print_triangle( joints[i].before_tri, false );
    }

    // set the new connectivity on the original triangle
    result = MBI()->set_connectivity( joints[i].before_tri, joints[i].after0, 3 );
    assert(moab::MB_SUCCESS == result);
    // set the new connectivity on the new triangle
    moab::EntityHandle new_tri;
    result = MBI()->create_element( moab::MBTRI, joints[i].after1, 3, new_tri );
    assert(moab::MB_SUCCESS == result);

    // copy the original normal to the new triangle
    moab::CartVect normal;
    result = MBI()->tag_get_data( normal_tag, &joints[i].before_tri, 1, &normal);
    assert(moab::MB_SUCCESS == result);
    result = MBI()->tag_set_data( normal_tag, &new_tri, 1, &normal);
    assert(moab::MB_SUCCESS == result);

    // add the new triangle to the same surface set as the original
    result = MBI()->add_entities( joints[i].surf_set, &new_tri, 1);
    assert(moab::MB_SUCCESS == result);

    // catch-all to remove degenerate tris
    result = delete_degenerate_tris( joints[i].before_tri );
    MB_CHK_SET_ERR(result,"could not delete degenerate tri");
    result = delete_degenerate_tris( new_tri );
    MB_CHK_SET_ERR(result,"could not delete degenerate tri");
  }
  return moab::MB_SUCCESS;
}

// Delete degenerate triangles in the range.
moab::ErrorCode Zip::delete_degenerate_tris( moab::EntityHandle tri )
{
  moab::ErrorCode result;
  const moab::EntityHandle *con;
  int n_verts;
  result = MBI()->get_connectivity( tri, con, n_verts);
  assert(moab::MB_SUCCESS == result);
  assert(3 == n_verts);
  if(con[0]==con[1] || con[1]==con[2] || con[2]==con[0]) {
    result = MBI()->delete_entities( &tri, 1 );
    assert(moab::MB_SUCCESS == result);
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Zip::delete_degenerate_tris( moab::Range tris )
{
  moab::ErrorCode result;
  for(moab::Range::iterator i=tris.begin(); i!=tris.end(); i++) {
    result = delete_degenerate_tris( *i );
    assert(moab::MB_SUCCESS == result);
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Zip::delete_adj_degenerate_tris( const moab::EntityHandle adj_vert )
{
  // get the adjacent triangles
  moab::ErrorCode result;
  moab::Range tris;
  result = MBI()->get_adjacencies( &adj_vert, 1, 2, false, tris );
  assert(moab::MB_SUCCESS == result);
  result = delete_degenerate_tris( tris );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

// Problem: MOAB can check for degenerate tris and delete them.
// MOAB check the curves and arcs passed in to update the merged vertex, but
// still ends up with degenerate edges. The curves that are in MOAB as sets
// are updated but also contain degenerate edges due to merging.

// Test to make sure the triangle normal vectors have not been inverted.
moab::ErrorCode Zip::test_normals( const std::vector<moab::CartVect> norms0,
                                   const std::vector<moab::CartVect> norms1,
                                   std::vector<int> &inverted_tri_indices )
{
  assert(norms0.size() == norms1.size());
  for(unsigned int i=0; i<norms0.size(); i++) {
    moab::ErrorCode result = test_normals( norms0[i], norms1[i]);
    if(moab::MB_SUCCESS != result) {
      inverted_tri_indices.push_back(i);
    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Zip::test_normals( const moab::CartVect norm0, const moab::CartVect norm1 )
{
  if(0 > norm0 % norm1) {
    return moab::MB_FAILURE;
  } else {
    return moab::MB_SUCCESS;
  }
}

/* Accepts a range of inverted tris. Refacets affected surface so that no tris
   are inverted. */
moab::ErrorCode Zip::remove_inverted_tris(moab::Tag normal_tag, moab::Range tris, const bool debug )
{

  moab::ErrorCode result;
  bool failures_occur = false;
  while(!tris.empty()) {

    /* Get a group of triangles to re-facet. They must be adjacent to each other
    and in the same surface. */
    moab::Range tris_to_refacet;
    tris_to_refacet.insert( tris.front() );
    moab::Range surf_set;
    result = MBI()->get_adjacencies( tris_to_refacet, 4, false, surf_set );
    assert(moab::MB_SUCCESS == result);
    if(1 != surf_set.size()) {
      std::cout << "remove_inverted_tris: tri is in " << surf_set.size()
                << " surfaces" << std::endl;
      return moab::MB_FAILURE;
    }

    // get all tris in the surface
    moab::Range surf_tris;
    result = MBI()->get_entities_by_type( surf_set.front(), moab::MBTRI, surf_tris );
    assert(moab::MB_SUCCESS == result);

    /* Find all of the adjacent inverted triangles of the same surface. Keep
    searching until a search returns no new triangles. */
    bool search_again = true;
    while(search_again) {

      // Here edges are being created. Remember to delete them. Outside of this
      // function. Skinning gets bogged down if unused MBEdges (from other
      // surfaces) accumulate.
      moab::Range tri_edges;
      result = MBI()->get_adjacencies( tris_to_refacet, 1, true, tri_edges,
                                       moab::Interface::UNION );
      assert(moab::MB_SUCCESS == result);
      moab::Range connected_tris;
      result = MBI()->get_adjacencies( tri_edges, 2, false, connected_tris,
                                       moab::Interface::UNION );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->delete_entities( tri_edges );
      assert(moab::MB_SUCCESS == result);
      moab::Range tris_to_refacet2 = intersect( tris_to_refacet, connected_tris );
      tris_to_refacet2 = intersect( tris_to_refacet, surf_tris );

      if(tris_to_refacet.size() == tris_to_refacet2.size()) search_again = false;
      tris_to_refacet.swap( tris_to_refacet2 );
    }

    // Remove the inverted tris that will be refaceted.
    tris = subtract( tris, tris_to_refacet );

    // do edges already exist?
    moab::Range temp;
    result = MBI()->get_entities_by_type(0, moab::MBEDGE, temp );
    assert(moab::MB_SUCCESS == result);
    if(!temp.empty()) MBI()->list_entities( temp );
    assert(temp.empty());


    // keep enlarging patch until we have tried to refacet the entire surface
    int counter=0;
    while(true) {
      // do edges already exist?
      temp.clear();
      result = MBI()->get_entities_by_type(0, moab::MBEDGE, temp );
      assert(moab::MB_SUCCESS == result);
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
      moab::Range tri_edges;
      result = MBI()->get_adjacencies( tris_to_refacet, 1, true, tri_edges,
                                       moab::Interface::UNION );
      assert(moab::MB_SUCCESS == result);

      // get all adjacent tris to the patch of inverted tris in the surface
      moab::Range adj_tris;
      result = MBI()->get_adjacencies( tri_edges, 2, false, adj_tris,
                                       moab::Interface::UNION );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->delete_entities( tri_edges );
      assert(moab::MB_SUCCESS == result);
      tris_to_refacet = intersect( surf_tris, adj_tris );
      if(tris_to_refacet.empty()) continue;
      // get an area-weighted normal of the adj_tris
      moab::CartVect plane_normal(0,0,0);
      for(moab::Range::iterator i=tris_to_refacet.begin(); i!=tris_to_refacet.end(); i++) {
        moab::CartVect norm;
        result = MBI()->tag_get_data( normal_tag, &(*i), 1, &norm);
        assert(moab::MB_SUCCESS == result);
        double area;
        result = gen->triangle_area( *i, area );
        assert(moab::MB_SUCCESS == result);
        if(debug) std::cout << "norm=" << norm << " area=" << area << std::endl;
        plane_normal += norm;
      }
      plane_normal.normalize();

      // do edges already exist?
      temp.clear();
      result = MBI()->get_entities_by_type(0, moab::MBEDGE, temp );
      assert(moab::MB_SUCCESS == result);
      if(!temp.empty()) MBI()->list_entities( temp );
      assert(temp.empty());

      // skin the tris
      moab::Range unordered_edges;
      result = gen->find_skin( tris_to_refacet, 1, unordered_edges, false );
      assert(moab::MB_SUCCESS == result);
      if(unordered_edges.empty()) {
        // do edges already exist?
        moab::Range temp;
        result = MBI()->get_entities_by_type(0, moab::MBEDGE, temp );
        assert(moab::MB_SUCCESS == result);
        if(!temp.empty()) MBI()->list_entities( temp );
        assert(temp.empty());
        continue;
      }

      // assemble into a polygon
      std::vector<moab::EntityHandle> polygon_of_verts;
      result = order_verts_by_edge( unordered_edges, polygon_of_verts );
      if(debug) gen->print_loop( polygon_of_verts );
      if(moab::MB_SUCCESS != result) {
        if(debug) std::cout << "remove_inverted_tris: couldn't order polygon by edge" << std::endl;
        return moab::MB_FAILURE;
      }

      // remember to remove edges
      result = MBI()->delete_entities( unordered_edges );
      assert(moab::MB_SUCCESS == result);

      // remove the duplicate endpt
      polygon_of_verts.pop_back();

      // the polygon should have at least 3 verts
      if(3 > polygon_of_verts.size()) {
        if(debug) std::cout << "remove_inverted_tris: polygon has too few points" << std::endl;
        return moab::MB_FAILURE;
      }

      // orient the polygon with the triangles (could be backwards)
      // get the first adjacent tri
      moab::EntityHandle edge[2] = { polygon_of_verts[0], polygon_of_verts[1] };
      moab::Range one_tri;
      result = MBI()->get_adjacencies( edge, 2, 2, false, one_tri );
      assert(moab::MB_SUCCESS == result);
      one_tri = intersect( tris_to_refacet, one_tri );
      assert(1 == one_tri.size());
      const moab::EntityHandle *conn;
      int n_conn;
      result = MBI()->get_connectivity( one_tri.front(), conn, n_conn );
      assert(moab::MB_SUCCESS == result);
      assert(3 == n_conn);
      if( (edge[0]==conn[1] && edge[1]==conn[0]) ||
          (edge[0]==conn[2] && edge[1]==conn[1]) ||
          (edge[0]==conn[0] && edge[1]==conn[2]) ) {
        reverse( polygon_of_verts.begin(), polygon_of_verts.end() );
        if(debug) std::cout << "remove_inverted_tris: polygon reversed" << std::endl;
      }

      /* facet the polygon. Returns moab::MB_FAILURE if it fails to facet the polygon. */
      moab::Range new_tris;
      result = gen->ear_clip_polygon( polygon_of_verts, plane_normal, new_tris );

      // break if the refaceting is successful
      if(moab::MB_SUCCESS == result) {
        // summarize tri area
        for(moab::Range::iterator i=new_tris.begin(); i!=new_tris.end(); i++) {
          double area;
          result = gen->triangle_area( *i, area );
          assert(moab::MB_SUCCESS == result);
          if(debug) std::cout << "  new tri area=" << area << std::endl;
        }

        // check the new normals
        std::vector<moab::CartVect> new_normals;
        result = gen->triangle_normals( new_tris, new_normals );
        if(moab::MB_SUCCESS != result) return result;

        // test the new triangles
        std::vector<int> inverted_tri_indices;
        std::vector<moab::CartVect> normals ( new_normals.size(), plane_normal );
        result = test_normals( normals, new_normals, inverted_tri_indices );
        assert(moab::MB_SUCCESS == result);
        if(inverted_tri_indices.empty()) {
          // remove the tris that were re-faceted
          tris = subtract( tris, tris_to_refacet );
          result = MBI()->remove_entities( surf_set.front(), tris_to_refacet );
          assert(moab::MB_SUCCESS == result);
          result = MBI()->delete_entities( tris_to_refacet );
          assert(moab::MB_SUCCESS == result);

          // add the new tris to the surf set
          result = MBI()->add_entities( surf_set.front(), new_tris );
          assert(moab::MB_SUCCESS == result);

          // put the new normals on the new tris
          result = gen->save_normals( new_tris, normal_tag );
          assert(moab::MB_SUCCESS == result);
          if(debug) std::cout << "remove_inverted_tris: success fixing a patch" << std::endl;
          break;
        }
      }

      // remember to delete the tris that were created from the failed ear clipping
      else {
        result = MBI()->delete_entities( new_tris );
        assert(moab::MB_SUCCESS == result);
      }

      // If the entire surface could not be ear clipped, give up
      if (tris_to_refacet.size() == surf_tris.size()) {
        if(debug) std::cout << "remove_inverted_tris: ear clipping entire surface failed"
                              << std::endl;
        return moab::MB_FAILURE;
      }

    } // loop until the entire surface has attempted to be refaceted
  }   // loop over each patch of inverted tris

  if(failures_occur) {
    if(debug) std::cout << "remove_inverted_facets: at least one failure occured" << std::endl;
    return moab::MB_FAILURE;
  } else {
    return moab::MB_SUCCESS;
  }
}

// we do not merge edges, just vert. check the verts
moab::ErrorCode Zip::test_zipping(const double FACET_TOL,
                                  const std::vector< std::vector<moab::EntityHandle> > arcs )
{
  moab::ErrorCode result;

  // make sure each arc has the same number of edges
  for(unsigned int i=1; i<arcs.size(); i++) {
    if(arcs[0].size() != arcs[i].size()) {
      std::cout << "The curve has " << arcs[0].size() << " edges but arc "
                << i << " has " << arcs[i].size() << " edges." << std::endl;
      gen->print_arcs( arcs );
      return moab::MB_FAILURE;
    }
  }

  // loop over every edge of the curve (first arc)
  for(unsigned int i=0; i<arcs[0].size()-1; i++) {
    // check for degenerate edge
    if(arcs[0][i] == arcs[0][i+1]) {
      std::cout << "degenerate edge at pos " << i << " and " << i+1 << " with verts "
                << arcs[0][i] << " and " << arcs[0][i+1] << std::endl;
      return moab::MB_FAILURE;
    }

    // check for edge of zero dist
    double d = gen->dist_between_verts( arcs[0][i], arcs[0][i+1] );
    if(FACET_TOL >= d) {
      std::cout << "edge length=" << d << " betwee pos " << i << " and " << i+1
                << " with verts " << arcs[0][i] << " and " << arcs[0][i+1] << std::endl;
      return moab::MB_FAILURE;
    }

    // loop over every arc
    for( unsigned int j=0; j<arcs.size(); j++) {
      // make sure vertices match
      if(arcs[0][i]!=arcs[j][i] || arcs[0][i+1]!=arcs[j][i+1]) {
        std::cout << "arc " << j << " vertices do not match curve vertices, pos= "
                  << i << "/" << arcs[j].size() << std::endl;
        return moab::MB_FAILURE;
      }
    }

    // make sure triangles have area
    moab::Range tris;
    result = MBI()->get_adjacencies( &(arcs[0][i]), 2, 2, false, tris );
    assert(moab::MB_SUCCESS == result);
    for(moab::Range::iterator k=tris.begin(); k!=tris.end(); k++) {
      // We know that there are not degenerate edges along the curve.
      // Sometimes degenerate tris are created due to merging curve endpts.
      // here we do not remove tri from the surf meshset, but we should
      if( gen->triangle_degenerate(*k) ) {
        std::cout << "  arc=" << 0 << " pos=" << i << " vert=" << arcs[0][i]
                  << " degenerate triangle" << std::endl;
        gen->print_triangle(*k, false);
        return moab::MB_FAILURE;
      }

      double area;
      result = gen->triangle_area( *k, area );
      assert(moab::MB_SUCCESS == result);
      // I found a valid tri on a curve with only one edge (1e-5 long)
      // that had an area of 1e-11.
      if(1e-8 > area) {
        std::cout << "    arc=" << 0 << " pos=" << i << " vert=" << arcs[0][i]
                  << " small triangle " << std::endl;
        gen->print_triangle(*k, false);
        gen->print_arcs( arcs );
      }
    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Zip::order_verts_by_edge( moab::Range unordered_edges,
    std::vector<moab::EntityHandle> &ordered_verts )
{
  if(unordered_edges.empty()) return moab::MB_SUCCESS;

  // get the endpoints of the curve. It should have 2 endpoints, unless is it a circle.
  moab::Range end_verts;
  moab::Skinner tool(MBI());
  moab::ErrorCode result;
  result = tool.find_skin( 0 , unordered_edges, 0, end_verts, false );
  if(moab::MB_SUCCESS != result) gen->print_range_of_edges( unordered_edges );
  assert(moab::MB_SUCCESS == result);

  // start with one endpoint
  moab::EntityHandle vert, edge;
  if(2 == end_verts.size()) {
    vert = end_verts.front();
  } else if (0 == end_verts.size()) {
    result = MBI()->get_adjacencies( &unordered_edges.front(), 1, 0, false, end_verts );
    assert(moab::MB_SUCCESS == result);
    assert(2 == end_verts.size());
    vert = end_verts.front();
  } else return moab::MB_FAILURE;

  // build the ordered set of verts. It will be as large as the number
  // of edges, plus one extra endpoint.
  ordered_verts.clear();
  ordered_verts.push_back( vert );

  // this cannot be used if multiple loops exist
  while(!unordered_edges.empty()) {
    // get an edge of the vert
    moab::Range adj_edges;
    result = MBI()->get_adjacencies( &vert, 1, 1, false, adj_edges );
    assert(moab::MB_SUCCESS == result);
    adj_edges = intersect( adj_edges, unordered_edges );
    if(adj_edges.empty()) {
      std::cout << "    order_verts_by_edgs: adj_edges is empty" << std::endl;
      return moab::MB_FAILURE;
    }
    edge = adj_edges.front();
    unordered_edges.erase( edge );

    // get the next vert
    end_verts.clear();
    result = MBI()->get_adjacencies( &edge, 1, 0, false, end_verts );
    assert(moab::MB_SUCCESS == result);
    if(2 != end_verts.size()) {
      std::cout << "end_verts.size()=" << end_verts.size() << std::endl;
      gen->print_edge( edge );
    }
    assert(2 == end_verts.size());
    vert = end_verts.front()==vert ? end_verts.back() : end_verts.front();
    ordered_verts.push_back( vert );
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Zip::merge_verts( const moab::EntityHandle keep_vert,
                                  const moab::EntityHandle delete_vert,
                                  std::vector<moab::EntityHandle> &arc0,
                                  std::vector<moab::EntityHandle> &arc1 )
{

  moab::ErrorCode rval;
  // first update the arcs with the keep_vert
  for(std::vector<moab::EntityHandle>::iterator i=arc0.begin(); i!=arc0.end(); ++i) {
    if(delete_vert == *i) *i = keep_vert;
  }
  for(std::vector<moab::EntityHandle>::iterator i=arc1.begin(); i!=arc1.end(); ++i) {
    if(delete_vert == *i) *i = keep_vert;
  }

  // get adjacent tris
  moab::Range tris;
  moab::EntityHandle verts[2]= {keep_vert, delete_vert};
  rval = MBI()->get_adjacencies( verts, 2, 2, false, tris, moab::Interface::UNION );
  MB_CHK_SET_ERR(rval,"getting adjacent tris failed");

  // actually do the merge
  rval = MBI()->merge_entities( keep_vert, delete_vert, false, true );
  MB_CHK_SET_ERR(rval,"merge entities failed");

  // delete degenerate tris
  rval = delete_degenerate_tris( tris );
  MB_CHK_SET_ERR(rval,"deleting degenerate tris failed");

  return moab::MB_SUCCESS;
}

