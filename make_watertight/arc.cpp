#include <iostream>
#include <iomanip> // for setprecision
#include <limits>  // for double min/max
#include <assert.h>
#include <vector>
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "gen.hpp"

#include "moab/GeomTopoTool.hpp"

namespace arc {

  MBErrorCode orient_edge_with_tri( const MBEntityHandle edge, const MBEntityHandle tri ) {
    MBErrorCode result;
    // get the connected vertices, properly ordered
    const MBEntityHandle *tri_conn;
    int n_verts;
    result = MBI()->get_connectivity( tri, tri_conn, n_verts );
    assert(MB_SUCCESS == result);
    assert( 3 == n_verts );

    // get the endpoints of the edge
    const MBEntityHandle *edge_conn;
    result = MBI()->get_connectivity( edge, edge_conn, n_verts );
    assert(MB_SUCCESS == result);
    assert( 2 == n_verts );

    // if the edge is backwards, reverse it
    if (( edge_conn[0]==tri_conn[0] && edge_conn[1]==tri_conn[2] ) ||
	( edge_conn[0]==tri_conn[1] && edge_conn[1]==tri_conn[0] ) ||
	( edge_conn[0]==tri_conn[2] && edge_conn[1]==tri_conn[1] ) ) {
      MBEntityHandle new_conn[2];
      new_conn[0] = edge_conn[1];
      new_conn[1] = edge_conn[0];
      result = MBI()->set_connectivity( edge, new_conn, 2 );
      assert(MB_SUCCESS == result);
    }
    return MB_SUCCESS;
  } 

  // Degenerate edges (same topological endpts) are caused by a prior step in which
  // coincident verts are merged.
  MBErrorCode remove_degenerate_edges( MBRange &edges, const bool debug ) {
    MBRange::iterator i = edges.begin();
    while (i!=edges.end()) {
      // get the endpoints of the edge
      MBErrorCode rval;
      const MBEntityHandle *endpts;
      int n_verts;
      rval = MBI()->get_connectivity( *i, endpts, n_verts );
      if(gen::error(MB_SUCCESS!=rval,"could not get connectivity")) 
        return rval;

      // remove the edge if degenerate
      if(2==n_verts && endpts[0]!=endpts[1]) {
        ++i;
      }	else if( (2==n_verts && endpts[0]==endpts[1]) ||
                 (1==n_verts                        ) ) {
	if(debug) {
          std::cout << "remove_degenerate_edges: deleting degenerate edge and tris " 
                    << std::endl;
        }
        rval = zip::delete_adj_degenerate_tris( endpts[0] );
        if(gen::error(MB_SUCCESS!=rval,"could not delete degenerate tris")) return rval;
        rval = MBI()->delete_entities( &(*i), 1 );
        if(gen::error(MB_SUCCESS!=rval,"could not delete degenerate edge")) return rval;
        i = edges.erase(i);
      } else {
	std::cout << "remove_degenerate_edge: wrong edge connectivity size" << std::endl;
        return MB_FAILURE;
      }
    
    }
    return MB_SUCCESS;
  }


  // Given a range of edges, remove pairs that have vertices (a,b) (b,a)
  MBErrorCode remove_opposite_pairs_of_edges( MBRange &edges, const bool debug ) {

    // do this in O(n) by using adjacencies instead of O(n^2)
    MBErrorCode result;
    //for(MBRange::iterator i=edges.begin(); i!=edges.end(); i++ ) {
    for(unsigned int i=0; i<edges.size(); i++) {
      MBEntityHandle the_edge = edges[i];

      // get endpoint verts
      MBRange two_verts;
      result = MBI()->get_adjacencies( &the_edge, 1, 0, false, two_verts);
      if(MB_SUCCESS != result) {
        std::cout << "result=" << result << " could not get adjacencies of edge" << std::endl;
        return result;
      }

      // get adjacent edges, but only keep the edges adjacent to both verts
      MBRange adj_edges;
      result = MBI()->get_adjacencies( two_verts, 1, false, adj_edges, MBInterface::INTERSECT);
      assert(MB_SUCCESS == result);
      // remove the original edge
      //adj_edges.erase( *i );

      // if any other edges exist, they are opposite the original edge and should be
      // removed from the skin
      if ( 1<adj_edges.size() ) {
     	if(debug) {
          std::cout << adj_edges.size() 
                    << " opposite edges will be removed from the surface skin " 
                    << adj_edges[0] << " " << adj_edges[1] << std::endl;
	}
	//gen::print_range_of_edges( adj_edges );
	//gen::print_range_of_edges( edges );
	//edges = edges.subtract( adj_edges );
        //edges.erase( *i );
        edges = subtract( edges, adj_edges );
        result = MBI()->delete_entities( adj_edges );
        assert(MB_SUCCESS == result);
        i--;
      }
    }
    return MB_SUCCESS;
  }

  MBErrorCode remove_opposite_pairs_of_edges_fast( MBRange &edges, const bool debug) {
    // special case
    MBErrorCode rval;
    if(1==edges.size()) {
      std::cout << "cannot remove pairs: only one input edge" << std::endl;
      return MB_FAILURE;
    }

    // populate edge array, used only for searching
    unsigned n_orig_edges = edges.size();
    gen::edge *my_edges = new gen::edge[n_orig_edges];
    unsigned j = 0;
    for(MBRange::const_iterator i=edges.begin(); i!=edges.end(); ++i) {
      // get the endpoints of the edge
      const MBEntityHandle *endpts;
      int n_verts;
      rval = MBI()->get_connectivity( *i, endpts, n_verts );
      if(gen::error(MB_SUCCESS!=rval || 2!=n_verts,"could not get connectivity")) 
        return rval;
    
      // store the edges
      my_edges[j].edge = *i;
      my_edges[j].v0   = endpts[0];
      my_edges[j].v1   = endpts[1];

      // sort edge by handle
      if(my_edges[j].v1 < my_edges[j].v0) {
        MBEntityHandle temp = my_edges[j].v0;
        my_edges[j].v0 = my_edges[j].v1;
        my_edges[j].v1 = temp;
      }
      ++j;
    }

    // sort edge array
    qsort(my_edges, n_orig_edges, sizeof(struct gen::edge), gen::compare_edge);

    // find duplicate edges
    j=0;
    MBRange duplicate_edges;
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
        if(gen::error(MB_SUCCESS!=rval,"cannot delete edge")) {
          delete[] my_edges;
          return rval;
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
    return MB_SUCCESS;
  }

  MBErrorCode get_next_oriented_edge( const MBRange edges, const MBEntityHandle edge,
				      MBEntityHandle &next_edge ) {

    // get the back vertex
    MBErrorCode result;
    const MBEntityHandle *end_verts;
    int n_verts;
    result = MBI()->get_connectivity( edge, end_verts, n_verts );
    assert(MB_SUCCESS==result);
    assert( 2 == n_verts );

    // get the edges adjacent to the back vertex
    MBRange adj_edges;
    result = MBI()->get_adjacencies( &(end_verts[1]), 1, 1, false, adj_edges );
    assert(MB_SUCCESS==result);

    // keep the edges that are part of the input range
    adj_edges = intersect( adj_edges, edges );
    // don't want the input edge
    adj_edges.erase( edge );

    // make sure the edge is oriented correctly
    for(MBRange::iterator i=adj_edges.begin(); i!=adj_edges.end(); i++) {
      const MBEntityHandle *adj_end_verts;
      result = MBI()->get_connectivity( *i, adj_end_verts, n_verts );
      if(MB_SUCCESS != result) {
        MBI()->list_entity(*i);
        std::cout << "result=" << result 
                  << " could not get connectivity of edge" << std::endl;
        return result;
	//gen::print_edge( *i );
      }
      assert(MB_SUCCESS==result);
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
      //assert(MB_SUCCESS == result);
      //return MB_MULTIPLE_ENTITIES_FOUND;
      next_edge = adj_edges.front();
    }
    return MB_SUCCESS;
  }

  MBErrorCode create_loops_from_oriented_edges_fast( MBRange edges,
						     std::vector< std::vector<MBEntityHandle> > &loops_of_edges,
                                                     const bool debug ) {
    // place all edges in map
    std::multimap<MBEntityHandle,gen::edge> my_edges;
    MBErrorCode rval;
    for(MBRange::const_iterator i=edges.begin(); i!=edges.end(); ++i) {
      // get the endpoints of the edge
      const MBEntityHandle *endpts;
      int n_verts;
      rval = MBI()->get_connectivity( *i, endpts, n_verts );
      if(gen::error(MB_SUCCESS!=rval || 2!=n_verts,"could not get connectivity")) 
        return rval;
    
      // store the edges
      gen::edge temp;
      temp.edge = *i;
      temp.v0   = endpts[0];
      temp.v1   = endpts[1];
      my_edges.insert( std::pair<MBEntityHandle,gen::edge>(temp.v0,temp) );  
    }
    std::cout << "error: function not complete" << std::endl;
    return MB_FAILURE;

    return MB_SUCCESS;
  }

  // This function should be rewritten using multimaps or something to avoid
  // upward adjacency searching. Vertices are searched for their adjacent edges.
  MBErrorCode create_loops_from_oriented_edges( MBRange edges,
						std::vector< std::vector<MBEntityHandle> > &loops_of_edges,
                                                const bool debug ) {

    // conserve edges
    MBErrorCode result;
    unsigned int n_edges_in  = edges.size();
    unsigned int n_edges_out = 0;
    //gen::print_range_of_edges( edges );
    // there could be several arcs for each surface
    while ( 0!= edges.size() ) {
      std::vector<MBEntityHandle> loop_of_edges;
      // pick initial edge and point
      MBEntityHandle edge = edges.front();

      // 20091201 Update: Pinch points may not be important. If not, there is no
      // purpose detecting them. Instead assume that pinch points coincide with 
      // the endpoints of geometric curves. Also assume that the loop creation at
      // pinch points does not matter. Pinch points can result in one or more 
      // loops, depending upon the path of traversal through the point.

      // Check to make sure the beginning endpt of the first edge is not a pinch 
      // point. If it is a pinch point the loop is ambiguous. Maybe--see watertightness notes for 20091201
      {
        const MBEntityHandle *end_verts;
        int n_verts;
        result = MBI()->get_connectivity( edge, end_verts, n_verts );
        assert(MB_SUCCESS==result);
        assert( 2 == n_verts );
        // get the edges adjacent to the back vertex
        MBRange adj_edges;
        result = MBI()->get_adjacencies( &(end_verts[0]), 1, 1, false, adj_edges );
        assert(MB_SUCCESS==result);
        // keep the edges that are part of the input range
        adj_edges = intersect( adj_edges, edges );
        if(2!=adj_edges.size() && debug) {
          std::cout << "  create_loops: adj_edges.size()=" << adj_edges.size() << std::endl;
	  std::cout << "  create_loops: pinch point exists" << std::endl;
          result = MBI()->list_entity( end_verts[0] );
          assert(MB_SUCCESS == result);
          //return MB_MULTIPLE_ENTITIES_FOUND;
        }
      }

      // add it to the loop
      loop_of_edges.push_back( edge );
      if(debug) std::cout << "push_back: " << edge << std::endl;
      n_edges_out++;
      edges.erase( edge );

      // find connected edges and add to the loop
      MBEntityHandle next_edge = 0;
      while (true) {

	// get the next vertex and next edge
	result = get_next_oriented_edge( edges, edge, next_edge );
        if(MB_ENTITY_NOT_FOUND == result) {
          return result;
        } else if(MB_SUCCESS != result) {
          gen::print_arc_of_edges( loop_of_edges );
          return result;
        }
	//assert( MB_SUCCESS == result );

	// if the next edge was found
	if ( 0!=next_edge ) {
	  // add it to the loop
	  loop_of_edges.push_back( next_edge );
          if(debug) std::cout << "push_back: " << next_edge << std::endl;
          n_edges_out++;

	  // remove the edge from the possible edges
	  edges.erase( next_edge );

	  // set the new reference vertex
	  //vert = next_vert;
	  edge = next_edge;

	  // if another edge was not found
	} else {
	  break;

	}
      }

      // check to ensure the arc is closed
      MBRange first_edge;
      first_edge.insert( loop_of_edges.front() );
      result = get_next_oriented_edge( first_edge, loop_of_edges.back(), next_edge );
      assert(MB_SUCCESS == result);
      if(next_edge != first_edge.front()) {
	std::cout << "create_loops: loop is not closed" << std::endl;
	gen::print_arc_of_edges(loop_of_edges);
        return MB_FAILURE;
      }

      // add the current arc to the vector of arcs
      loops_of_edges.push_back(loop_of_edges);
    }

    // check to make sure that we have the same number of verts as we started with
    if(gen::error(n_edges_in!=n_edges_out,"edges not conserved")) return MB_FAILURE;
    assert( n_edges_in == n_edges_out );

    return MB_SUCCESS;
  }

  // return a set of ordered_verts and remaining unordered_edges
  MBErrorCode order_verts_by_edge( MBRange unordered_edges, 
                                   std::vector<MBEntityHandle> &ordered_verts ) {
    if(unordered_edges.empty()) return MB_SUCCESS;
 
    // get the endpoints of the curve. It should have 2 endpoints, unless is it a circle.
    MBRange end_verts;
    MBSkinner tool(MBI());
    MBErrorCode result;
    result = tool.find_skin( unordered_edges, 0, end_verts, false );
    if(MB_SUCCESS != result) gen::print_range_of_edges( unordered_edges );
    assert(MB_SUCCESS == result);

    // start with one endpoint
    MBEntityHandle vert, edge;
    if(2 == end_verts.size()) {
      vert = end_verts.front();
    } else if (0 == end_verts.size()) {
      result = MBI()->get_adjacencies( &unordered_edges.front(), 1, 0, false, end_verts );
      assert(MB_SUCCESS == result);
      assert(2 == end_verts.size());
      vert = end_verts.front();
    } else return MB_FAILURE;
    
    // build the ordered set of verts. It will be as large as the number
    // of edges, plus one extra endpoint.
    ordered_verts.clear();
    ordered_verts.push_back( vert );
   
    // this cannot be used if multiple loops exist
    while(!unordered_edges.empty()) {
      // get an edge of the vert
      MBRange adj_edges;
      result = MBI()->get_adjacencies( &vert, 1, 1, false, adj_edges );
      assert(MB_SUCCESS == result);
      adj_edges = intersect( adj_edges, unordered_edges );
      //assert(!adj_edges.empty());
      if(adj_edges.empty()) {
	std::cout << "    order_verts_by_edgs: adj_edges is empty" << std::endl;
        return MB_FAILURE;
      }
      edge = adj_edges.front(); 
      unordered_edges.erase( edge );

      // get the next vert
      end_verts.clear();
      result = MBI()->get_adjacencies( &edge, 1, 0, false, end_verts );
      assert(MB_SUCCESS == result);
      if(2 != end_verts.size()) {
        std::cout << "end_verts.size()=" << end_verts.size() << std::endl;
	gen::print_edge( edge );
      }
      assert(2 == end_verts.size());
      vert = end_verts.front()==vert ? end_verts.back() : end_verts.front();
      ordered_verts.push_back( vert );
    }
    return MB_SUCCESS;
  }
     
  MBErrorCode get_meshset( const MBEntityHandle set, std::vector<MBEntityHandle> &vec) {
    MBErrorCode result;
    vec.clear();
    result = MBI()->get_entities_by_handle( set, vec );
    assert(MB_SUCCESS == result);  
    return MB_SUCCESS;
  }

  MBErrorCode set_meshset( const MBEntityHandle set, const std::vector<MBEntityHandle> vec) {
    MBErrorCode result;
    result = MBI()->clear_meshset( &set, 1 );
    assert(MB_SUCCESS == result);  
    result = MBI()->add_entities( set, &vec[0], vec.size() );       
    assert(MB_SUCCESS == result);  
    return MB_SUCCESS;
  }

  MBErrorCode merge_curves( MBRange curve_sets, const double facet_tol, 
                            MBTag id_tag, MBTag merge_tag, const bool debug ) {
    // find curve endpoints to add to kd tree
    MBErrorCode result;
    const double SQR_TOL = facet_tol*facet_tol;
    MBRange endpoints;
    for(MBRange::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++) {
      std::vector<MBEntityHandle> curve;
      result = get_meshset( *i, curve );
      assert(MB_SUCCESS == result);
      //if(2 > curve.size()) continue;
      assert(1 < curve.size());
      MBEntityHandle front_endpt = curve[0];
      MBEntityHandle back_endpt  = curve[curve.size()-1];
      // ADD CODE TO HANDLE SPECIAL CASES!!
      if(front_endpt == back_endpt) { // special case
        if(0 == gen::length(curve)) { // point curve
        } else {                      // circle
        }
      } else {                        // normal curve
        endpoints.insert( front_endpt );
        endpoints.insert( back_endpt );
      }
    }

    // add endpoints to kd-tree. Tree must track ownership to know when verts are
    // merged away (deleted).  
    MBAdaptiveKDTree kdtree(MBI(), 0, MESHSET_TRACK_OWNER);                           
    MBEntityHandle root;                                                         
    MBAdaptiveKDTree::Settings settings;                
    settings.maxEntPerLeaf = 1;                                                     
    settings.candidateSplitsPerDir = 1;                           
    settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION; 
    result = kdtree.build_tree( endpoints, root, &settings);            
    assert(MB_SUCCESS == result);                                               
    MBAdaptiveKDTreeIter tree_iter;                                     
    kdtree.get_tree_iterator( root, tree_iter );     

    // search for other endpoints that match each curve's endpoints
    for(MBRange::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++) {
      std::vector<MBEntityHandle> curve_i_verts;
      result = get_meshset( *i, curve_i_verts );
      assert(MB_SUCCESS == result);
      double curve_length = gen::length( curve_i_verts );
      //if(2 > curve.size()) continue; // HANDLE SPECIAL CASES (add logic)     
      if(curve_i_verts.empty()) continue;
      MBEntityHandle endpts[2] = { curve_i_verts.front(), curve_i_verts.back() };
      MBCartVect endpt_coords;
      //if( endpts[0] == endpts[1]) continue; // special case of point curve or circle
      std::vector<MBEntityHandle> leaves;
      MBRange adj_curves[2];
      // match the front then back endpts                        
      for(unsigned int j=0; j<2; j++) {                                  
	//gen::print_vertex_coords( endpts[j] );
        result = MBI()->get_coords( &endpts[j], 1, endpt_coords.array());
        assert(MB_SUCCESS == result);                                        
        result = kdtree.leaves_within_distance( root, endpt_coords.array(), 
                                                facet_tol, leaves);
        assert(MB_SUCCESS == result);                                
        for(unsigned int k=0; k<leaves.size(); k++) {           
	  std::vector<MBEntityHandle> leaf_verts;               
          result = MBI()->get_entities_by_type( leaves[k], MBVERTEX, leaf_verts);    
          assert(MB_SUCCESS == result);             
          for(unsigned int l=0; l<leaf_verts.size(); l++) {               
            double sqr_dist;                                             
            result = gen::squared_dist_between_verts( endpts[j], leaf_verts[l], sqr_dist);        
            assert(MB_SUCCESS == result);
            /* Find parent curves. There will be no parent curves if the curve has
	       already been merged and no longer exists. */     
            if(SQR_TOL >= sqr_dist) {
              MBRange temp_curves;
              result = MBI()->get_adjacencies( &leaf_verts[l], 1, 4, false, temp_curves);
	      assert(MB_SUCCESS == result);
              // make sure the sets are curve sets before adding them to the list
              // of candidates.
              temp_curves = intersect(temp_curves, curve_sets);
              adj_curves[j].merge( temp_curves );
            }
          }
        }
      }

      // now find the curves that have matching endpts
      MBRange candidate_curves;
      candidate_curves =intersect( adj_curves[0], adj_curves[1] );
      if(candidate_curves.empty()) continue;

      // subtract the current curve
      int n_before = candidate_curves.size();
      candidate_curves.erase( *i );
      int n_after = candidate_curves.size();
      //if(gen::error(n_before!=n_after+1,"could not erase curve")) {
      //  return MB_FAILURE;
      //}

      // find matching curves
      for(MBRange::iterator j=candidate_curves.begin(); j!=candidate_curves.end(); j++) {
	std::vector<MBEntityHandle> curve_j_verts;
        result = get_meshset( *j, curve_j_verts );
        assert(MB_SUCCESS == result);
        double j_curve_length = gen::length( curve_j_verts );

        int i_id, j_id;
        result = MBI()->tag_get_data( id_tag, &(*i), 1, &i_id );                         
        assert(MB_SUCCESS == result);   
        result = MBI()->tag_get_data( id_tag, &(*j), 1, &j_id );                         
        assert(MB_SUCCESS == result);   
	if(debug) {
          std::cout << "curve i_id=" << i_id << " j_id=" << j_id 
                    << " leng0=" << curve_length << " leng1=" << j_curve_length << std::endl;
        }
        // reject curves with significantly different length (for efficiency)
        if( facet_tol < abs(curve_length - j_curve_length)) continue;

        // align j_curve to be the same as i_curve
        bool reversed;
        if(gen::dist_between_verts(curve_i_verts.front(), curve_j_verts.front()) >
           gen::dist_between_verts(curve_i_verts.front(), curve_j_verts.back())) {
	  reverse( curve_j_verts.begin(), curve_j_verts.end() );
          reversed = true;
        } else {
          reversed = false;
        }

        // Reject curves if the average distance between them is greater than 
        // merge_tol.
        double dist;
        result = gen::dist_between_arcs( debug, curve_i_verts, curve_j_verts, dist );
        assert(MB_SUCCESS == result);
        if( facet_tol < dist ) continue;

        // THE CURVE WILL BE MERGED
	std::cout << "  merging curve " << j_id << " to curve " << i_id 
                  << ", dist_between_curves=" << dist << " cm" << std::endl;

        // Merge the endpts of the curve to preserve topology. Merging (and deleting)
        // the endpoints will also remove them from the KDtree so that the merged
        // curve cannot be selected again. Update the curves when merging to avoid
        // stake info.
	if(curve_i_verts.front() != curve_j_verts.front()) { 
          result = zip::merge_verts( curve_i_verts.front(), curve_j_verts.front(), 
                                     curve_i_verts, curve_j_verts );
	  if(MB_SUCCESS != result) std::cout << result << std::endl;
	  assert(MB_SUCCESS == result);
	}
	if(curve_i_verts.back() != curve_j_verts.back()) {
          result = zip::merge_verts( curve_i_verts.back(), curve_j_verts.back(),
                                     curve_i_verts, curve_j_verts );
	  if(MB_SUCCESS != result) std::cout << result << std::endl;
	  assert(MB_SUCCESS == result);
	}

        // Tag the curve that is merged away. We do not delete it so that its 
        // parent-child links are preserved. Later, a surface's facets will be
        // deleted if all of its curves are deleted or merged with themselves.
        // Only after this can the merged away curves be deleted.
        result = MBI()->tag_set_data( merge_tag, &(*j), 1, &(*i) );
        assert(MB_SUCCESS == result);

        // clear the sets contents
        curve_j_verts.clear();
        result = set_meshset( *j, curve_j_verts ); 
        assert(MB_SUCCESS == result);

        // reverse the curve-surf senses if the curves do not share the same orientation
        if(reversed) {
  	  std::vector<MBEntityHandle> surfs;
          std::vector<int> senses;
          moab::GeomTopoTool gt(MBI(), false);
          result = gt.get_senses( *j, surfs, senses );
          if(gen::error(MB_SUCCESS!=result,"failed to get senses")) return result;
          for(unsigned k=0; k<surfs.size(); ++k) {
            if(0==senses[k])
              senses[k] = 1;
            else if(1==senses[k])
              senses[k] = 0;
            else if(-1==senses[k])
              senses[k] = -1;
            else
              if(gen::error(true,"unrecognized sense")) return MB_FAILURE;
          }   
          result = gt.set_senses( *j, surfs, senses );
          if(gen::error(MB_SUCCESS!=result,"failed to set senses")) return result;
        }

      }
    }
    return MB_SUCCESS;
  }



}
