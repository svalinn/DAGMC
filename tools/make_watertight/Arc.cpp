#include <iostream>
#include <iomanip> // for setprecision
#include <limits>  // for double min/max
#include <assert.h>
#include <vector>

// moab includes
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/FileOptions.hpp"

#include "Arc.hpp"

moab::ErrorCode Arc::orient_edge_with_tri( const moab::EntityHandle edge, const moab::EntityHandle tri )
{
  moab::ErrorCode result;
  // get the connected vertices, properly ordered
  const moab::EntityHandle *tri_conn;
  int n_verts;
  result = MBI()->get_connectivity( tri, tri_conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert( 3 == n_verts );

  // get the endpoints of the edge
  const moab::EntityHandle *edge_conn;
  result = MBI()->get_connectivity( edge, edge_conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert( 2 == n_verts );

  // if the edge is backwards, reverse it
  if (( edge_conn[0]==tri_conn[0] && edge_conn[1]==tri_conn[2] ) ||
      ( edge_conn[0]==tri_conn[1] && edge_conn[1]==tri_conn[0] ) ||
      ( edge_conn[0]==tri_conn[2] && edge_conn[1]==tri_conn[1] ) ) {
    moab::EntityHandle new_conn[2];
    new_conn[0] = edge_conn[1];
    new_conn[1] = edge_conn[0];
    result = MBI()->set_connectivity( edge, new_conn, 2 );
    assert(moab::MB_SUCCESS == result);
  }
  return moab::MB_SUCCESS;
}

// Degenerate edges (same topological endpts) are caused by a prior step in which
// coincident verts are merged.
moab::ErrorCode Arc::remove_degenerate_edges( moab::Range &edges, const bool debug )
{
  moab::Range::iterator i = edges.begin();
  while (i!=edges.end()) {
    // get the endpoints of the edge
    moab::ErrorCode rval;
    const moab::EntityHandle *endpts;
    int n_verts;
    rval = MBI()->get_connectivity( *i, endpts, n_verts );
    MB_CHK_SET_ERR(rval,"could not get connectivity");


    // remove the edge if degenerate
    if(2==n_verts && endpts[0]!=endpts[1]) {
      ++i;
    }	else if( (2==n_verts && endpts[0]==endpts[1]) ||
               (1==n_verts                        ) ) {
      if(debug) {
        std::cout << "remove_degenerate_edges: deleting degenerate edge and tris "
                  << std::endl;
      }
      rval = zip->delete_adj_degenerate_tris( endpts[0] );
      MB_CHK_SET_ERR(rval,"could not delete degenerate tris");
      rval = MBI()->delete_entities( &(*i), 1 );
      MB_CHK_SET_ERR(rval,"could not delete degenerate edge");
      i = edges.erase(i);
    } else {
      std::cout << "remove_degenerate_edge: wrong edge connectivity size" << std::endl;
      return moab::MB_FAILURE;
    }

  }
  return moab::MB_SUCCESS;
}


// Given a range of edges, remove pairs that have vertices (a,b) (b,a)
moab::ErrorCode Arc::remove_opposite_pairs_of_edges( moab::Range &edges, const bool debug )
{

  // do this in O(n) by using adjacencies instead of O(n^2)
  moab::ErrorCode result;
  for(unsigned int i=0; i<edges.size(); i++) {
    moab::EntityHandle the_edge = edges[i];

    // get endpoint verts
    moab::Range two_verts;
    result = MBI()->get_adjacencies( &the_edge, 1, 0, false, two_verts);
    if(moab::MB_SUCCESS != result) {
      std::cout << "result=" << result << " could not get adjacencies of edge" << std::endl;
      return result;
    }

    // get adjacent edges, but only keep the edges adjacent to both verts
    moab::Range adj_edges;
    result = MBI()->get_adjacencies( two_verts, 1, false, adj_edges, moab::Interface::INTERSECT);
    assert(moab::MB_SUCCESS == result);

    // if any other edges exist, they are opposite the original edge and should be
    // removed from the skin
    if ( 1<adj_edges.size() ) {
      if(debug) {
        std::cout << adj_edges.size()
                  << " opposite edges will be removed from the surface skin "
                  << adj_edges[0] << " " << adj_edges[1] << std::endl;
      }
      edges = subtract( edges, adj_edges );
      result = MBI()->delete_entities( adj_edges );
      assert(moab::MB_SUCCESS == result);
      i--;
    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::remove_opposite_pairs_of_edges_fast( moab::Range &edges, const bool debug)
{
  // special case
  moab::ErrorCode rval;
  if(1==edges.size()) {
    std::cout << "cannot remove pairs: only one input edge" << std::endl;
    return moab::MB_FAILURE;
  }

  // populate edge array, used only for searching
  unsigned n_orig_edges = edges.size();
  edge *my_edges = new edge[n_orig_edges];
  unsigned j = 0;
  for(moab::Range::const_iterator i=edges.begin(); i!=edges.end(); ++i) {
    // get the endpoints of the edge
    const moab::EntityHandle *endpts;
    int n_verts;
    rval = MBI()->get_connectivity( *i, endpts, n_verts );
    if(moab::MB_SUCCESS!=rval || 2!=n_verts) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"could not get connectivity");
    }

    // store the edges
    my_edges[j].edge = *i;
    my_edges[j].v0   = endpts[0];
    my_edges[j].v1   = endpts[1];

    // sort edge by handle
    if(my_edges[j].v1 < my_edges[j].v0) {
      moab::EntityHandle temp = my_edges[j].v0;
      my_edges[j].v0 = my_edges[j].v1;
      my_edges[j].v1 = temp;
    }
    ++j;
  }

  // sort edge array
  qsort(my_edges, n_orig_edges, sizeof(struct edge), compare_edge);

  // find duplicate edges
  j=0;
  moab::Range duplicate_edges;
  for(unsigned i=1; i<n_orig_edges; ++i) {
    // delete edge if a match exists
    if(my_edges[j].v0==my_edges[i].v0 && my_edges[j].v1==my_edges[i].v1) {
      duplicate_edges.insert( my_edges[j].edge );
      // find any remaining matches
      while( my_edges[j].v0==my_edges[i].v0 &&
             my_edges[j].v1==my_edges[i].v1 &&
             i<n_orig_edges) {
        duplicate_edges.insert( my_edges[i].edge );
        ++i;
      }
      // delete the matches
      edges = subtract( edges, duplicate_edges );
      rval = MBI()->delete_entities( duplicate_edges );
      if(moab::MB_SUCCESS!=rval) {
        delete[] my_edges;
        MB_CHK_SET_ERR(moab::MB_FAILURE,"cannot delete edge");
      }
      if(debug) {
        std::cout << "remove_opposite_edges: deleting " << duplicate_edges.size()
                  << " edges" << std::endl;
      }
      duplicate_edges.clear();
    }
    j = i;
  }

  delete[] my_edges;
  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::get_next_oriented_edge( const moab::Range edges,
    const moab::EntityHandle edge,
    moab::EntityHandle &next_edge )
{

  // get the back vertex
  moab::ErrorCode result;
  const moab::EntityHandle *end_verts;
  int n_verts;
  result = MBI()->get_connectivity( edge, end_verts, n_verts );
  assert(moab::MB_SUCCESS==result);
  assert( 2 == n_verts );

  // get the edges adjacent to the back vertex
  moab::Range adj_edges;
  result = MBI()->get_adjacencies( &(end_verts[1]), 1, 1, false, adj_edges );
  assert(moab::MB_SUCCESS==result);

  // keep the edges that are part of the input range
  adj_edges = intersect( adj_edges, edges );
  // don't want the input edge
  adj_edges.erase( edge );

  // make sure the edge is oriented correctly
  for(moab::Range::iterator i=adj_edges.begin(); i!=adj_edges.end(); i++) {
    const moab::EntityHandle *adj_end_verts;
    result = MBI()->get_connectivity( *i, adj_end_verts, n_verts );
    if(moab::MB_SUCCESS != result) {
      MBI()->list_entity(*i);
      std::cout << "result=" << result
                << " could not get connectivity of edge" << std::endl;
      return result;
      //print_edge( *i );
    }
    assert(moab::MB_SUCCESS==result);
    assert( 2 == n_verts );
    if ( end_verts[1]!=adj_end_verts[0] ) i = adj_edges.erase(i) - 1;
  }

  /* The next edge could be ambiguous if more than one exists.This happens in
     surfaces that are ~1D, and in surfaces that have pinch points
     (mod13surf881). Although I didn't handle this case yet, if it occurs:
     -Remember that a pinch point could have not just 2, but multiple ears.
     -Select the next edge as the edge that shares the same triangle.
     -This approach will only work if the pinch point also exists in the
      geometric curves.
     -If the pinch point does not exist in the geometric curves (~1D surfs),
      there is no robust way to handle it. There should be an input assumption
      that this never happens. "The faceting skin cannot have pinch points
      unless they also occur in the surface's geometric curves."
   */
  if ( 0==adj_edges.size() ) {
    next_edge = 0;
  } else if ( 1==adj_edges.size() ) {
    next_edge = adj_edges.front();
  } else {
    std::cout << "get_next_oriented_edge: " << adj_edges.size() <<
              " possible edges indicates a pinch point." << std::endl;
    result = MBI()->list_entity( end_verts[1] );
    //assert(moab::MB_SUCCESS == result);
    //return moab::MB_MULTIPLE_ENTITIES_FOUND;
    next_edge = adj_edges.front();
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::create_loops_from_oriented_edges_fast( moab::Range edges,
    std::vector< std::vector<moab::EntityHandle> > &loops_of_edges,
    const bool debug )
{
  // place all edges in map
  std::multimap<moab::EntityHandle,edge> my_edges;
  moab::ErrorCode rval;
  for(moab::Range::const_iterator i=edges.begin(); i!=edges.end(); ++i) {
    // get the endpoints of the edge
    const moab::EntityHandle *endpts;
    int n_verts;
    rval = MBI()->get_connectivity( *i, endpts, n_verts );
    if(moab::MB_SUCCESS!=rval || 2!=n_verts) {
      MB_CHK_SET_ERR(moab::MB_FAILURE, "could not get connectivity");
    }
    // store the edges
    edge temp;
    temp.edge = *i;
    temp.v0   = endpts[0];
    temp.v1   = endpts[1];
    my_edges.insert( std::pair<moab::EntityHandle,edge>(temp.v0,temp) );
  }
  std::cout << "error: function not complete" << std::endl;
  return moab::MB_FAILURE;

  return moab::MB_SUCCESS;
}

// This function should be rewritten using multimaps or something to avoid
// upward adjacency searching. Vertices are searched for their adjacent edges.
moab::ErrorCode Arc::create_loops_from_oriented_edges( moab::Range edges,
    std::vector< std::vector<moab::EntityHandle> > &loops_of_edges,
    const bool debug )
{

  // conserve edges
  moab::ErrorCode result;
  unsigned int n_edges_in  = edges.size();
  unsigned int n_edges_out = 0;
  // there could be several arcs for each surface
  while ( 0!= edges.size() ) {
    std::vector<moab::EntityHandle> loop_of_edges;
    // pick initial edge and point
    moab::EntityHandle edge = edges.front();

    // 20091201 Update: Pinch points may not be important. If not, there is no
    // purpose detecting them. Instead assume that pinch points coincide with
    // the endpoints of geometric curves. Also assume that the loop creation at
    // pinch points does not matter. Pinch points can result in one or more
    // loops, depending upon the path of traversal through the point.

    // Check to make sure the beginning endpt of the first edge is not a pinch
    // point. If it is a pinch point the loop is ambiguous. Maybe--see watertightness notes for 20091201
    {
      const moab::EntityHandle *end_verts;
      int n_verts;
      result = MBI()->get_connectivity( edge, end_verts, n_verts );
      assert(moab::MB_SUCCESS==result);
      assert( 2 == n_verts );
      // get the edges adjacent to the back vertex
      moab::Range adj_edges;
      result = MBI()->get_adjacencies( &(end_verts[0]), 1, 1, false, adj_edges );
      assert(moab::MB_SUCCESS==result);
      // keep the edges that are part of the input range
      adj_edges = intersect( adj_edges, edges );
      if(2!=adj_edges.size() && debug) {
        std::cout << "  create_loops: adj_edges.size()=" << adj_edges.size() << std::endl;
        std::cout << "  create_loops: pinch point exists" << std::endl;
        result = MBI()->list_entity( end_verts[0] );
        assert(moab::MB_SUCCESS == result);
      }
    }

    // add it to the loop
    loop_of_edges.push_back( edge );
    if(debug) std::cout << "push_back: " << edge << std::endl;
    n_edges_out++;
    edges.erase( edge );

    // find connected edges and add to the loop
    moab::EntityHandle next_edge = 0;
    while (true) {

      // get the next vertex and next edge
      result = get_next_oriented_edge( edges, edge, next_edge );
      if(moab::MB_ENTITY_NOT_FOUND == result) {
        return result;
      } else if(moab::MB_SUCCESS != result) {
        gen->print_arc_of_edges( loop_of_edges );
        return result;
      }

      // if the next edge was found
      if ( 0!=next_edge ) {
        // add it to the loop
        loop_of_edges.push_back( next_edge );
        if(debug) std::cout << "push_back: " << next_edge << std::endl;
        n_edges_out++;

        // remove the edge from the possible edges
        edges.erase( next_edge );

        // set the new reference vertex
        edge = next_edge;

        // if another edge was not found
      } else {
        break;

      }
    }

    // check to ensure the arc is closed
    moab::Range first_edge;
    first_edge.insert( loop_of_edges.front() );
    result = get_next_oriented_edge( first_edge, loop_of_edges.back(), next_edge );
    assert(moab::MB_SUCCESS == result);
    if(next_edge != first_edge.front()) {
      std::cout << "create_loops: loop is not closed" << std::endl;
      gen->print_arc_of_edges(loop_of_edges);
      return moab::MB_FAILURE;
    }

    // add the current arc to the vector of arcs
    loops_of_edges.push_back(loop_of_edges);
  }

  // check to make sure that we have the same number of verts as we started with
  if(n_edges_in!=n_edges_out) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"edges not conserved");
  }
  assert( n_edges_in == n_edges_out );

  return moab::MB_SUCCESS;
}

// return a set of ordered_verts and remaining unordered_edges
moab::ErrorCode Arc::order_verts_by_edge( moab::Range unordered_edges,
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

moab::ErrorCode Arc::get_meshset( const moab::EntityHandle set, std::vector<moab::EntityHandle> &vec)
{
  moab::ErrorCode result;
  vec.clear();
  result = MBI()->get_entities_by_handle( set, vec );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::set_meshset( const moab::EntityHandle set, const std::vector<moab::EntityHandle> vec)
{
  moab::ErrorCode result;
  result = MBI()->clear_meshset( &set, 1 );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->add_entities( set, &vec[0], vec.size() );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::merge_verts( const moab::EntityHandle keep_vert,
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
  rval = zip->delete_degenerate_tris( tris );
  MB_CHK_SET_ERR(rval,"deleting degenerate tris failed");

  return moab::MB_SUCCESS;
}

moab::ErrorCode Arc::merge_curves( moab::Range curve_sets, const double facet_tol,
                                   moab::Tag id_tag, moab::Tag merge_tag, const bool debug )
{
  // find curve endpoints to add to kd tree
  moab::ErrorCode result;
  const double SQR_TOL = facet_tol*facet_tol;
  moab::Range endpoints;
  for(moab::Range::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++) {
    std::vector<moab::EntityHandle> curve;
    result = get_meshset( *i, curve );
    assert(moab::MB_SUCCESS == result);
    assert(1 < curve.size());
    moab::EntityHandle front_endpt = curve[0];
    moab::EntityHandle back_endpt  = curve[curve.size()-1];
    // ADD CODE TO HANDLE SPECIAL CASES!!
    if(front_endpt == back_endpt) { // special case
      if(0 == gen->length(curve)) { // point curve
      } else {                      // circle
      }
    } else {                        // normal curve
      endpoints.insert( front_endpt );
      endpoints.insert( back_endpt );
    }
  }

  // add endpoints to kd-tree. Tree must track ownership to know when verts are
  // merged away (deleted).

  moab::AdaptiveKDTree kdtree(MBI()); //, 0, MESHSET_TRACK_OWNER);
  moab::EntityHandle root;

  //set tree options
  const char settings[]="MAX_PER_LEAF=1;SPLITS_PER_DIR=1;PLANE_SET=0;MESHSET_FLAGS=0x1;TAG_NAME=0";
  moab::FileOptions fileopts(settings);

  // initialize the tree and pass the root entity handle back into root
  result = kdtree.build_tree( endpoints, &root, &fileopts);
  assert(moab::MB_SUCCESS == result);
  // create tree iterator
  moab::AdaptiveKDTreeIter tree_iter;
  kdtree.get_tree_iterator( root, tree_iter );

  // search for other endpoints that match each curve's endpoints
  for(moab::Range::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++) {
    std::vector<moab::EntityHandle> curve_i_verts;
    result = get_meshset( *i, curve_i_verts );
    assert(moab::MB_SUCCESS == result);
    double curve_length = gen->length( curve_i_verts );
    if(curve_i_verts.empty()) continue;
    moab::EntityHandle endpts[2] = { curve_i_verts.front(), curve_i_verts.back() };
    moab::CartVect endpt_coords;
    std::vector<moab::EntityHandle> leaves;

    // initialize an array which will contain matched of front points in [0] and
    // matches for back points in [1]
    moab::Range adj_curves[2];
    // match the front then back endpts
    for(unsigned int j=0; j<2; j++) {
      result = MBI()->get_coords( &endpts[j], 1, endpt_coords.array());
      assert(moab::MB_SUCCESS == result);
      // takes all leaves of the tree within a distance (facet_tol) of the coordinates
      // passed in by endpt_coords.array() and returns them in leaves
      result = kdtree.distance_search( endpt_coords.array(),
                                       facet_tol, leaves, root );
      assert(moab::MB_SUCCESS == result);
      for(unsigned int k=0; k<leaves.size(); k++) {
        // retrieve all information about vertices in leaves
        std::vector<moab::EntityHandle> leaf_verts;
        result = MBI()->get_entities_by_type( leaves[k], moab::MBVERTEX, leaf_verts);
        assert(moab::MB_SUCCESS == result);
        for(unsigned int l=0; l<leaf_verts.size(); l++) {
          double sqr_dist;
          result = gen->squared_dist_between_verts( endpts[j], leaf_verts[l], sqr_dist);
          assert(moab::MB_SUCCESS == result);
          /* Find parent curves. There will be no parent curves if the curve has
          already been merged and no longer exists. */
          if(SQR_TOL >= sqr_dist) {
            moab::Range temp_curves;
            // get the curves for all vertices that are within the squared distance of each other
            result = MBI()->get_adjacencies( &leaf_verts[l], 1, 4, false, temp_curves);
            assert(moab::MB_SUCCESS == result);
            // make sure the sets are curve sets before adding them to the list
            // of candidates.
            temp_curves = intersect(temp_curves, curve_sets);
            adj_curves[j].merge( temp_curves );
          }
        }
      }
    }

    // now find the curves that have matching endpts
    moab::Range candidate_curves;
    // select curves that do not have coincident front AND back points
    // place them into candidate curves
    candidate_curves =intersect( adj_curves[0], adj_curves[1] );
    if(candidate_curves.empty()) continue;

    // subtract the current curve
    int n_before = candidate_curves.size();
    candidate_curves.erase( *i );
    int n_after = candidate_curves.size();


    // now find curves who's interior vertices are also coincident and merge them
    for(moab::Range::iterator j=candidate_curves.begin(); j!=candidate_curves.end(); j++) {
      std::vector<moab::EntityHandle> curve_j_verts;
      result = get_meshset( *j, curve_j_verts );
      assert(moab::MB_SUCCESS == result);
      double j_curve_length = gen->length( curve_j_verts );

      int i_id, j_id;
      result = MBI()->tag_get_data( id_tag, &(*i), 1, &i_id );
      assert(moab::MB_SUCCESS == result);
      result = MBI()->tag_get_data( id_tag, &(*j), 1, &j_id );
      assert(moab::MB_SUCCESS == result);
      if(debug) {
        std::cout << "curve i_id=" << i_id << " j_id=" << j_id
                  << " leng0=" << curve_length << " leng1=" << j_curve_length << std::endl;
      }

      // reject curves with significantly different length (for efficiency)
      if( facet_tol < abs(curve_length - j_curve_length)) continue;

      // align j_curve to be the same as i_curve
      bool reversed;
      if(gen->dist_between_verts(curve_i_verts.front(), curve_j_verts.front()) >
         gen->dist_between_verts(curve_i_verts.front(), curve_j_verts.back())) {
        reverse( curve_j_verts.begin(), curve_j_verts.end() );
        reversed = true;
      } else {
        reversed = false;
      }

      // Reject curves if the average distance between them is greater than
      // facet_tol.
      double dist;
      result = gen->dist_between_arcs( debug, curve_i_verts, curve_j_verts, dist );
      assert(moab::MB_SUCCESS == result);
      if( facet_tol < dist ) continue;

      // THE CURVE WILL BE MERGED
      if (debug) {
        std::cout << "  merging curve " << j_id << " to curve " << i_id
                  << ", dist_between_curves=" << dist << " cm" << std::endl;
      }

      // Merge the endpts of the curve to preserve topology. Merging (and deleting)
      // the endpoints will also remove them from the KDtree so that the merged
      // curve cannot be selected again. Update the curves when merging to avoid
      // stale info.
      if(curve_i_verts.front() != curve_j_verts.front()) {
        result = zip->merge_verts( curve_i_verts.front(), curve_j_verts.front(),
                                   curve_i_verts, curve_j_verts );
        if(moab::MB_SUCCESS != result) std::cout << result << std::endl;
        assert(moab::MB_SUCCESS == result);
      }
      if(curve_i_verts.back() != curve_j_verts.back()) {
        result = zip->merge_verts( curve_i_verts.back(), curve_j_verts.back(),
                                   curve_i_verts, curve_j_verts );
        if(moab::MB_SUCCESS != result) std::cout << result << std::endl;
        assert(moab::MB_SUCCESS == result);
      }

      // Tag the curve that is merged away. We do not delete it so that its
      // parent-child links are preserved. Later, a surface's facets will be
      // deleted if all of its curves are deleted or merged with themselves.
      // Only after this can the merged away curves be deleted.
      result = MBI()->tag_set_data( merge_tag, &(*j), 1, &(*i) );
      assert(moab::MB_SUCCESS == result);

      // clear the sets contents
      curve_j_verts.clear();
      result = set_meshset( *j, curve_j_verts );
      assert(moab::MB_SUCCESS == result);

      // reverse the curve-surf senses if the merged curve was found to be opposite the
      // curve we keep
      if(reversed) {
        std::vector<moab::EntityHandle> surfs;
        std::vector<int> senses;
        moab::GeomTopoTool gt(MBI(), false);
        result = gt.get_senses( *j, surfs, senses );
        MB_CHK_SET_ERR(result,"failed to get senses");
        for(unsigned k=0; k<surfs.size(); ++k) {
          //forward to reverse
          if(moab::SENSE_FORWARD==senses[k])
            senses[k] = moab::SENSE_REVERSE;
          //reverse to forward
          else if(moab::SENSE_REVERSE==senses[k])
            senses[k] = moab::SENSE_FORWARD;
          //unknown to unknown
          else if(moab::SENSE_BOTH==senses[k])
            senses[k] = moab::SENSE_BOTH;
          else {
            MB_CHK_SET_ERR(moab::MB_FAILURE,"unrecognized sense");
          }
        }
        result = gt.set_senses( *j, surfs, senses );
        MB_CHK_SET_ERR(result,"failed to set senses");
      }

    }
  }
  return moab::MB_SUCCESS;
}

