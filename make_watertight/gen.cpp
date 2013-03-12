#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/GeomTopoTool.hpp"

#include "gen.hpp"
#include "zip.hpp"
#include "MBSkinner.hpp"
#include "DagMC.hpp"

namespace gen {

  bool error( const bool error_has_occured, const std::string message ) {
    if(error_has_occured) {
      if("" == message) {
	std::cout << "Error at " << __FILE__ << ":" << __LINE__ 
                  << std::endl;
      } else {
        std::cout << message << std::endl;
      }
      return true;
    } else {
      return false;
    }
  }

  void print_vertex_cubit( const MBEntityHandle vertex ) {

    MBErrorCode result;
    double coords[3];
    int n_precision = 20;
    result = MBI()->get_coords( &vertex, 1, coords );
    assert(MB_SUCCESS == result);
    std::cout << "  create vertex " 
	      << std::setprecision(n_precision)
	      << coords[0] << " " << coords[1] << " " << coords[2]
	      << std::endl;
  }

  void print_vertex_coords( const MBEntityHandle vertex ) {

    MBErrorCode result;
    double coords[3];
    int n_precision = 20;
    result = MBI()->get_coords( &vertex, 1, coords );
    if(MB_SUCCESS!=result) std::cout << "vert=" << vertex << std::endl;
    assert(MB_SUCCESS == result);
    std::cout << "    vertex " << vertex << " coords= (" 
      //<< std::setprecision(n_precision)
	      << coords[0] << "," << coords[1] << "," << coords[2] << ")" 
	      << std::endl;
  }

  void print_triangles( const MBRange tris ) {
    for(MBRange::const_iterator i=tris.begin(); i!=tris.end(); i++) {
      print_triangle( *i, false );
    }
  }
  // If the edges of the tri are ambiguous, do not print edges!
  void print_triangle( const MBEntityHandle tri, bool print_edges ) {
    MBErrorCode result;
    double area;
    result = triangle_area( tri, area );
    assert(MB_SUCCESS == result);
    std::cout << "    triangle " << tri << " area=" << area << std::endl;
    const MBEntityHandle *conn;
    int n_verts;
    result = MBI()->get_connectivity( tri, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    for(int i=0; i<3; i++) print_vertex_coords( conn[i] );
    
    if(print_edges) {
      MBRange edges;
      result = MBI()->get_adjacencies( &tri, 1, 1, true, edges );
      if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
      assert(MB_SUCCESS == result);
      //std::cout << "      edges: ";
      for(MBRange::iterator i=edges.begin(); i!=edges.end(); i++) {
        //std::cout << *i << " ";
        print_edge( *i );
      }
      //std::cout << std::endl;
    }
  }

  void print_edge( const MBEntityHandle edge ) {
    const MBEntityHandle *conn;
    int n_verts;
    std::cout << "    edge " << edge << std::endl;
    MBErrorCode result = MBI()->get_connectivity( edge, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(2 == n_verts);
    print_vertex_coords( conn[0] );   
    print_vertex_coords( conn[1] ); 

    MBRange tris;
    result = MBI()->get_adjacencies( &edge, 1, 2, false, tris );
    assert(MB_SUCCESS == result);
    std::cout << "     tris: ";
    for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
      std::cout << *i << " ";
    }
    std::cout << std::endl; 
  }

  void print_range( const MBRange range ) {
    std::cout << "print range:" << std::endl;
    MBRange::iterator i;
    for(i=range.begin(); i!=range.end(); i++) {
      std::cout << "    " << *i << std::endl;
    }
  }

  void print_range_of_edges( const MBRange range ) {
    std::cout << "print range:" << std::endl;
    MBRange::const_iterator i;
    for(i=range.begin(); i!=range.end(); i++) {
      print_edge( *i );
    }
  }


  void print_vertex_count(const MBEntityHandle input_meshset) {
  
    // get the range of facets of the surface meshset
    MBErrorCode result;
    MBRange vertices;
    result = MBI()->get_entities_by_type(0, MBVERTEX, vertices);
    assert( MB_SUCCESS == result );
  
    std::cout<< "    " << vertices.size() << " vertices found." << std::endl;
  }

  void print_arcs( const std::vector< std::vector<MBEntityHandle> > arcs ) {
    for(unsigned int i=0; i<arcs.size(); i++) {
      std::cout << "arc " << i << std::endl;
      print_loop( arcs[i] );
    }
  }

  void print_arc_of_edges( const std::vector<MBEntityHandle> arc_of_edges ) {
  
    MBErrorCode result;
    std::vector<MBEntityHandle>::const_iterator i;
    double dist = 0;
    for( i=arc_of_edges.begin(); i!=arc_of_edges.end(); i++ ) {
      int n_verts;
      const MBEntityHandle *conn;
      result = MBI()->get_connectivity( *i, conn, n_verts ); 
      assert(MB_SUCCESS == result);
      assert( 2 == n_verts );
      dist += dist_between_verts( conn[0], conn[1] );
      print_vertex_coords( conn[0] );
      print_vertex_coords( conn[1] );
    }
    std::cout << "  dist= " << dist << std::endl;
  }

  void print_loop( const std::vector<MBEntityHandle> loop_of_verts ) {
   
    std::cout << "  size=" << loop_of_verts.size() << std::endl;
    double dist = 0;
    //std::vector<MBEntityHandle>::iterator i;
    //for( i=loop_of_verts.begin(); i!=loop_of_verts.end(); i++ ) {
    for(unsigned int i=0; i<loop_of_verts.size(); i++) {
      print_vertex_coords( loop_of_verts[i] );
      if(i != loop_of_verts.size()-1) {
	dist += dist_between_verts( loop_of_verts[i], loop_of_verts[i+1] );
      }
    }
    std::cout << "  dist=" << dist << std::endl;
  }


// Return the closest vertex to the arc.
// For efficiency: only get_coords on the reference vertex once
//                 if specified, limit search length along curve
MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
                               const std::vector<MBEntityHandle> arc_of_verts,
                               unsigned &position,
                               const double dist_limit ) {
  MBErrorCode rval;
  const bool debug = false;
  double min_dist_sqr = std::numeric_limits<double>::max();
  MBCartVect ref_coords;
  rval = MBI()->get_coords( &reference_vert, 1, ref_coords.array() );
  if(gen::error(MB_SUCCESS!=rval,"failed to get ref coords")) return rval;
  double length = 0;
  MBCartVect prev_coords;

  for(unsigned i=0; i<arc_of_verts.size(); ++i) {
    MBCartVect coords;
    rval = MBI()->get_coords( &arc_of_verts[i], 1, coords.array() );
    if(gen::error(MB_SUCCESS!=rval,"failed to get coords")) return rval;

    // use dist_limit to exit early; avoid checking the entire arc
    if(0!=i) {
      MBCartVect temp = prev_coords - coords;
      length += temp.length();
      if(length>dist_limit && debug) 
        std::cout << "length=" << length << " dist_limit=" << dist_limit << std::endl;
      if(length > dist_limit) return MB_SUCCESS;
    }
    prev_coords = coords;

    // get distance to ref_vert
    MBCartVect temp = ref_coords - coords;
    double dist_sqr = temp.length_squared();
    if(dist_sqr < min_dist_sqr) {
      position = i;
      min_dist_sqr = dist_sqr;
      if(debug) std::cout << "min_dist_sqr=" << min_dist_sqr << std::endl;
    }
  }
 
  return MB_SUCCESS;
}



  // Return the closest vert and all within tol. This is needed because sometimes
  // the correct vert is not the closest. For example, iter_surf4010 the skin
  // loop has the same point in it twice, at two different locations (center of L).
  // This ensure that both are returned as candidates.
  MBErrorCode find_closest_vert( const double tol,
                                 const MBEntityHandle reference_vert,
				 const std::vector<MBEntityHandle> loop_of_verts,
				 std::vector<unsigned> &positions, 
				 std::vector<double> &dists) {

    MBErrorCode rval;
    positions.clear();
    dists.clear();
    const double TOL_SQR = tol*tol;
    unsigned min_pos;
    double sqr_min_dist = std::numeric_limits<double>::max();
    for(unsigned int i=0; i<loop_of_verts.size(); i++) {
	double sqr_dist = std::numeric_limits<double>::max();
	rval = squared_dist_between_verts(reference_vert, loop_of_verts[i], sqr_dist);
        if(gen::error(MB_SUCCESS!=rval,"could not get dist")) return rval;
	if(sqr_dist < sqr_min_dist) {
          if(sqr_dist >= TOL_SQR) {
            sqr_min_dist = sqr_dist;
            min_pos = i;
          } else {
            sqr_min_dist = TOL_SQR;
            positions.push_back(i);
            dists.push_back( sqrt(sqr_dist) );
          }
	}
    }
 
    if(dists.empty()) {
      dists.push_back( sqrt(sqr_min_dist) );
      positions.push_back( min_pos );
    }

            
	  //sqr_min_dist = sqr_dist;
	  //if(sqr_min_dist == 0.0) {
          //  min_dist = 0.0;
          //  return MB_SUCCESS; // can't do better than this
          //}
    //min_dist = sqrt( sqr_min_dist );
    return MB_SUCCESS;
  }
  /*  MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
				 const std::vector<std::vector<MBEntityHandle> > loops_of_verts,
				 unsigned int &loop, unsigned int &position, 
				 double &min_dist) {

    MBErrorCode result;
    min_dist = std::numeric_limits<double>::max();
    for(unsigned int i=0; i<loops_of_verts.size(); i++) {
      unsigned int p;
      double d;
      result = find_closest_vert( reference_vert, loops_of_verts[i], p, d );
      assert(MB_SUCCESS == result);  
      if(d < min_dist) {
	min_dist = d;
        loop = i;
        position = p;
        if(min_dist == 0.0) return MB_SUCCESS; // can't do better than this
      }
    }
    return MB_SUCCESS;
  }
  */
  MBErrorCode merge_vertices( MBRange verts /* in */, const double tol /* in */ ) {

    MBErrorCode result;
    const double SQR_TOL = tol*tol;
    // Clean up the created tree, and track verts so that if merged away they are
    // removed from the tree.
    MBAdaptiveKDTree kdtree(MBI(), true, 0, MESHSET_TRACK_OWNER);
    MBEntityHandle root;
    MBAdaptiveKDTree::Settings settings;
    settings.maxEntPerLeaf = 6;
    settings.maxTreeDepth  = 50;
    settings.candidateSplitsPerDir = 1;
    settings.candidatePlaneSet = MBAdaptiveKDTree::VERTEX_MEDIAN;
    result = kdtree.build_tree( verts, root, &settings);
    assert(MB_SUCCESS == result);
    MBAdaptiveKDTreeIter tree_iter;
    kdtree.get_tree_iterator( root, tree_iter );
 
    //for(unsigned int i=0; i<verts.size(); i++) {
    for(MBRange::iterator i=verts.begin(); i!=verts.end(); ++i) {
      double from_point[3];
      //MBEntityHandle vert = *i;
      result = MBI()->get_coords( &(*i), 1, from_point);
      assert(MB_SUCCESS == result);
      std::vector<MBEntityHandle> leaves_out;
      result = kdtree.leaves_within_distance( root, from_point, tol, leaves_out);
      assert(MB_SUCCESS == result);
      for(unsigned int j=0; j<leaves_out.size(); j++) {
	std::vector<MBEntityHandle> leaf_verts;
        result = MBI()->get_entities_by_type( leaves_out[j], MBVERTEX, leaf_verts);
        assert(MB_SUCCESS == result);
	if(100 < leaf_verts.size()) std::cout << "*i=" << *i << " leaf_verts.size()=" << leaf_verts.size() << std::endl;
        for(unsigned int k=0; k<leaf_verts.size(); k++) {
	  if( leaf_verts[k] == *i ) continue; 
	  double sqr_dist;
	  result = gen::squared_dist_between_verts( *i, leaf_verts[k], sqr_dist);
	  assert(MB_SUCCESS == result);

	  if(SQR_TOL >= sqr_dist) {
            //  std::cout << "merge_vertices: vert " << leaf_verts[k] << " merged to vert "
            //            << verts[i] << " dist=" << dist << " leaf=" << std::endl;

            // The delete_vert is automatically remove from the tree because it
            // uses tracking meshsets. merge_verts checks for degenerate tris.
            // Update the list of leaf verts to prevent stale handles.
	    std::vector<MBEntityHandle> temp_arc;
            MBEntityHandle keep_vert   = *i;
            MBEntityHandle delete_vert = leaf_verts[k];
	    result = zip::merge_verts( keep_vert, delete_vert, leaf_verts, temp_arc );
	    assert(MB_SUCCESS == result);
            // Erase delete_vert from verts
	    // Iterator should remain valid because delete_vert > keep_vert handle.
	    verts.erase( delete_vert );
	  }
	}
      }
    }   
    return MB_SUCCESS;
  }

  MBErrorCode squared_dist_between_verts( const MBEntityHandle v0, 
                                          const MBEntityHandle v1, 
                                          double &d) {
    MBErrorCode result;
    MBCartVect coords0, coords1;
    result = MBI()->get_coords( &v0, 1, coords0.array() );
    if(MB_SUCCESS != result) {
      std::cout << "dist_between_verts: get_coords on v0=" << v0 << " result=" 
                << result << std::endl; 
      return result;
    }
    result = MBI()->get_coords( &v1, 1, coords1.array() );
    if(MB_SUCCESS != result) {
      std::cout << "dist_between_verts: get_coords on v1=" << v1 << " result=" 
                << result << std::endl; 
      return result;
    }
    const MBCartVect diff = coords0 - coords1;
    d = diff.length_squared();
    return MB_SUCCESS;
  }

  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 ) {
    MBCartVect v2 = v0 - v1;
    return v2.length();
  }
  MBErrorCode dist_between_verts( const MBEntityHandle v0, const MBEntityHandle v1, double &d) {
    MBErrorCode result;
    MBCartVect coords0, coords1;
    result = MBI()->get_coords( &v0, 1, coords0.array() );
    if(MB_SUCCESS != result) {
      std::cout << "dist_between_verts: get_coords on v0=" << v0 << " result=" 
                << result << std::endl; 
      return result;
    }
    result = MBI()->get_coords( &v1, 1, coords1.array() );
    if(MB_SUCCESS != result) {
      std::cout << "dist_between_verts: get_coords on v1=" << v1 << " result=" 
                << result << std::endl; 
      return result;
    }
    d = dist_between_verts( coords0, coords1 );
    return MB_SUCCESS;
  }

  double dist_between_verts( double coords0[], double coords1[] ) {
    return sqrt( (coords0[0]-coords1[0])*(coords0[0]-coords1[0]) +
		 (coords0[1]-coords1[1])*(coords0[1]-coords1[1]) +
		 (coords0[2]-coords1[2])*(coords0[2]-coords1[2]) );
  }
  double dist_between_verts( MBEntityHandle vert0, MBEntityHandle vert1 ) {
    double coords0[3], coords1[3];
    MBErrorCode result;
    result = MBI()->get_coords( &vert0, 1, coords0 );
    if(MB_SUCCESS!=result) std::cout << "result=" << result << " vert=" 
                                     << vert0 << std::endl;
    assert(MB_SUCCESS == result);
    result = MBI()->get_coords( &vert1, 1, coords1 );
    if(MB_SUCCESS!=result) std::cout << "result=" << result << " vert=" 
                                     << vert1 << std::endl;
    assert(MB_SUCCESS == result);
    return dist_between_verts( coords0, coords1 );
  }

  // Return the length of the curve defined by MBEDGEs or ordered MBVERTEXs.
  double length( std::vector<MBEntityHandle> edges ) {
    if(edges.empty()) return 0;

    MBErrorCode result;
    std::vector<MBEntityHandle>::iterator i; 
    double dist = 0;
    MBEntityType type = MBI()->type_from_handle( edges[0] ); 

    // if vector has both edges and verts, only use edges
    // NOTE: The curve sets from ReadCGM do not contain duplicate endpoints for loops!
    MBEntityType end_type = MBI()->type_from_handle( edges.back() );
    if(type != end_type) {
      for(std::vector<MBEntityHandle>::iterator i=edges.begin(); i!=edges.end(); i++) {
        if(MBVERTEX == MBI()->type_from_handle( *i )) {
	  i = edges.erase(i) - 1;
        }
      }
    }

    // determine if vector defines an arc by edges of verts
    type = MBI()->type_from_handle( edges[0] ); 
    if        (MBEDGE == type) {
      if(edges.empty()) return 0.0; 
      for( i=edges.begin(); i!=edges.end(); i++ ) {
	int n_verts;
	const MBEntityHandle *conn;
	result = MBI()->get_connectivity( *i, conn, n_verts );
	if( MB_SUCCESS!=result ) std::cout << "result=" << result << std::endl; 
	assert(MB_SUCCESS == result);
	assert( 2 == n_verts );
        if(conn[0] == conn[1]) continue;
	dist += dist_between_verts( conn[0], conn[1] );
	//std::cout << "length: " << dist << std::endl;
      }
    } else if (MBVERTEX == type) {
      if(2 > edges.size()) return 0.0;
      MBEntityHandle front_vert = edges.front();
      for( i=edges.begin()+1; i!=edges.end(); i++) {
	dist += dist_between_verts( front_vert, *i );
	front_vert = *i;
      }
    } else return MB_FAILURE;

    return dist;
  }

  // Given a vertex and vector of edges, return the number of edges adjacent to the vertex.
  unsigned int n_adj_edges( MBEntityHandle vert, MBRange edges ) {
    MBErrorCode result;
    MBRange adj_edges;
    result = MBI()->get_adjacencies( &vert, 1, 1, false, adj_edges );
    assert(MB_SUCCESS == result);
    //adj_edges = adj_edges.intersect(edges);
    adj_edges = intersect( adj_edges, edges );
    return adj_edges.size();
  }



  // Return true if the edges share a vertex. Does not check for coincident edges.
  bool edges_adjacent( MBEntityHandle edge0, MBEntityHandle edge1 ) {
    MBErrorCode result;
    MBRange verts0, verts1;
    result = MBI()->get_adjacencies( &edge0, 1, 0, false, verts0 );
    assert( MB_SUCCESS == result );
    assert( 2 == verts0.size() );
    result = MBI()->get_adjacencies( &edge1, 1, 0, false, verts1 );
    assert( MB_SUCCESS == result );
    assert( 2 == verts1.size() );
    if      ( verts0.front() == verts1.front() ) return true;
    else if ( verts0.front() == verts1.back()  ) return true;
    else if ( verts0.back()  == verts1.back()  ) return true;
    else if ( verts0.back()  == verts1.front() ) return true;
    else                                         return false;
  }

  // get the direction unit vector from one vertex to another vertex
  MBErrorCode get_direction( const MBEntityHandle from_vert, const MBEntityHandle to_vert,
			     MBCartVect &dir ) {
    // double d[3];
    MBErrorCode result;
    MBCartVect coords0, coords1;
    result = MBI()->get_coords( &from_vert, 1, coords0.array() );
    assert(MB_SUCCESS==result);
    result = MBI()->get_coords( &to_vert, 1, coords1.array() );
    assert(MB_SUCCESS==result);
    dir = coords1 - coords0;
    if(0 == dir.length()) {
      MBCartVect zero_vector( 0.0 );
      dir = zero_vector;
      std::cout << "direction vector has 0 magnitude" << std::endl;
      return MB_SUCCESS;  
    }
    dir.normalize();
    return MB_SUCCESS;
  }    
 
  // from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
  double edge_point_dist( const MBCartVect a, const MBCartVect b, const MBCartVect c ) {
    MBCartVect ab, bc, ba, ac;
    ab = b - a;
    bc = c - b;
    ba = a - b;
    ac = c - a;
  
    // find the magnitude of the cross product and test the line
    MBCartVect cross_product = ab*ac;
    double dist = cross_product.length() / dist_between_verts(a,b);
 
    // test endpoint1
    if (ab%bc > 0) {
      //std::cout << "edge_point_dist=" << dist_between_verts(b,c) 
      //<< " at endpt1" << std::endl;
      return dist_between_verts(b,c);
    }

    // test endpoint0
    if (ba%ac > 0) {
      //std::cout << "edge_point_dist=" << dist_between_verts(a,c) 
      //<< " at endpt0" << std::endl;
      return dist_between_verts(a,c);
    }
 
    //std::cout << "edge_point_dist=" << fabs(dist) << " at middle" 
    //<< std::endl;      
    return fabs(dist);
  }
  double edge_point_dist( const MBEntityHandle endpt0, const MBEntityHandle endpt1, 
			  const MBEntityHandle pt ) {
    MBErrorCode result;
    MBCartVect a, b, c;
   result = MBI()->get_coords( &endpt0, 1, a.array() );
    assert(MB_SUCCESS==result);
    result = MBI()->get_coords( &endpt1, 1, b.array() );
    assert(MB_SUCCESS==result);
    result = MBI()->get_coords( &pt,     1, c.array() );
    assert(MB_SUCCESS==result);
    return edge_point_dist( a, b, c);
  }
  double edge_point_dist( const MBEntityHandle edge, const MBEntityHandle pt ) {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;     
    result = MBI()->get_connectivity( edge, conn, n_verts );
    assert(MB_SUCCESS==result);
    assert( 2 == n_verts );
    return edge_point_dist( conn[0], conn[1], pt );
  }
  /*
  MBErrorCode  point_curve_min_dist( const std::vector<MBEntityHandle> curve, // of verts 
				     const MBEntityHandle pt,
				     double &min_dist,
                                     const double max_dist_along_curve ) { 
  
    min_dist = std::numeric_limits<double>::max();
    double cumulative_dist = 0;
    bool last_edge = false;
  
    // it is a curve of verts or a curve of edges?
    MBEntityType type = MBI()->type_from_handle( curve.front() );
    std::vector<MBEntityHandle>::const_iterator i;
    if(MBVERTEX == type) {
      MBEntityHandle front_vert = curve.front();
      for( i=curve.begin()+1; i!=curve.end(); i++) {
	// if we are using verts, do not explicitly create the edge in MOAB
        cumulative_dist += gen::dist_between_verts( front_vert, *i );
	//std::cout << "  point_curve_min_dist: cumulative_dist="
	//        << cumulative_dist << std::endl;
	double d = edge_point_dist( front_vert, *i, pt );
	//std::cout << "  point_curve_min_dist: d=" << d << std::endl;        
	if(d < min_dist) {
	  min_dist = d;
	  //std::cout << "min_dist=" << min_dist << std::endl;
	  //print_vertex_coords( front_vert );
	  //print_vertex_coords( *i );
	}
        // check one edge past the point after max_dist_along_curve
        if(last_edge) return MB_SUCCESS;
        if(max_dist_along_curve<cumulative_dist) last_edge = true;
      
	front_vert = *i;
      }
    } //else if(MBEDGE == type) {        
      //for( i=curve.begin(); i!=curve.end(); i++) {
//	double d = edge_point_dist( *i, pt );
	//if(d < min_dist) min_dist = d;
      //}
     // } 
       else return MB_FAILURE;
    //std::cout << "point_curve_min_dist=" << min_dist << " curve.size()=" << curve.size() << std::endl;    
    return MB_SUCCESS;
  }

  MBErrorCode  point_curve_min_dist( const std::vector<MBEntityHandle> curve, // of verts 
				     const MBEntityHandle pt,
				     double &min_dist ) { 
    const double max_dist_along_curve = std::numeric_limits<double>::max();
    return point_curve_min_dist( curve, pt, min_dist, max_dist_along_curve ); 
  }
*/
  double triangle_area( const MBCartVect a, const MBCartVect b, 
                        const MBCartVect c) {
    MBCartVect d = c - a;
    MBCartVect e = c - b;
    MBCartVect f = d*e;
    return 0.5*f.length();
  }
  MBErrorCode triangle_area( const MBEntityHandle conn[], double &area ) {
    MBCartVect coords[3];
    MBErrorCode result = MBI()->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == result);
    area = triangle_area( coords[0], coords[1], coords[2] );
    return MB_SUCCESS;
  }
  MBErrorCode triangle_area( const MBEntityHandle tri, double &area ) {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;
    result = MBI()->get_connectivity( tri, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    
    result = triangle_area( conn, area );
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }
  double triangle_area( const MBRange tris ) {
    double a, area = 0;
    MBErrorCode result;
    for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
      result = triangle_area( *i, a);
      assert(MB_SUCCESS == result);
      area += a;
    }
    return area;
  }

  bool triangle_degenerate( const MBEntityHandle tri ) {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;
    result = MBI()->get_connectivity( tri, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    return triangle_degenerate( conn[0], conn[1], conn[2] );
  }

  bool triangle_degenerate( const MBEntityHandle v0, const MBEntityHandle v1,
			    const MBEntityHandle v2 ) { 
    if(v0==v1 || v1==v2 || v2==v0) return true;
    return false;
  }

  MBErrorCode triangle_normals( const MBRange tris, std::vector<MBCartVect> &normals ) {
    MBErrorCode result;
    normals.clear();
    for(MBRange::const_iterator i=tris.begin(); i!=tris.end(); i++) {
      MBCartVect normal;
      result = triangle_normal( *i, normal );
      assert(MB_SUCCESS==result || MB_ENTITY_NOT_FOUND==result);
      normals.push_back( normal );
    }
    return MB_SUCCESS;
  }

  MBErrorCode triangle_normal( const MBEntityHandle tri, MBCartVect &normal) {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;
    result = MBI()->get_connectivity( tri, conn, n_verts );
    if(MB_ENTITY_NOT_FOUND == result) {
      std::cout << "triangle_normal: triangle not found" << std::endl;
      MBCartVect zero_vector( 0.0 );
      normal = zero_vector;
      return result;
    }else if(MB_SUCCESS != result) {
      return result;
    } else {
      assert(3 == n_verts);
      return triangle_normal( conn[0], conn[1], conn[2], normal );
    }
  }

  MBErrorCode triangle_normal( const MBEntityHandle v0, const MBEntityHandle v1,
			       const MBEntityHandle v2, MBCartVect &normal ) {

    // if tri is degenerate return 0,0,0
    if( triangle_degenerate(v0, v1, v2) ) {
      MBCartVect zero_vector( 0.0 );
      normal = zero_vector;
      std::cout << "  normal=" << normal << std::endl;
      return MB_SUCCESS;
    }

    MBEntityHandle conn[3];
    conn[0] = v0;
    conn[1] = v1;
    conn[2] = v2;
    MBErrorCode result;
    MBCartVect coords[3];
    result = MBI()->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == result); 
    return triangle_normal( coords[0], coords[1], coords[2], normal );
  }

  MBErrorCode triangle_normal( const MBCartVect coords0, const MBCartVect coords1,
			       const MBCartVect coords2, MBCartVect &normal ) {
    MBCartVect edge0, edge1;
    edge0 = coords1-coords0;
    edge1 = coords2-coords0;
    normal = edge0*edge1;

    // do not normalize if magnitude is zero (avoid nans)
    if(0 == normal.length()) return MB_SUCCESS;

    normal.normalize();
    //if(debug) std::cout << "  normal=" << normal << std::endl;
    return MB_SUCCESS;
  }
    


  // Distance between a point and line. The line is defined by two verts.
  // We are using a line and not a line segment!
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  MBErrorCode line_point_dist( const MBEntityHandle line_pt1, const MBEntityHandle line_pt2, 
			       const MBEntityHandle pt0, double &dist ) {
    MBErrorCode result;
    MBCartVect x0, x1, x2;
    result = MBI()->get_coords( &line_pt1, 1, x1.array() );
    assert(MB_SUCCESS == result); 
    result = MBI()->get_coords( &line_pt2, 1, x2.array() );
    assert(MB_SUCCESS == result); 
    result = MBI()->get_coords( &pt0, 1, x0.array() );
    assert(MB_SUCCESS == result); 

    dist = ( ((x0-x1)*(x0-x2)).length() ) / ( (x2-x1).length() );
    return MB_SUCCESS;
  }

  // Project the point onto the line. Not the line segment!
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,
				     const MBEntityHandle line_pt2,
				     const MBEntityHandle pt0 ) {
    MBCartVect projected_coords;
    double parameter;
    MBErrorCode result = point_line_projection( line_pt1, line_pt2, 
						pt0, projected_coords,
						parameter );
    assert(MB_SUCCESS == result);
    result = MBI()->set_coords( &pt0, 1, projected_coords.array() );   
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }    
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,
				     const MBEntityHandle line_pt2,
				     const MBEntityHandle pt0,
				     MBCartVect &projected_coords,
				     double &parameter  ) {

    MBErrorCode result;
    MBCartVect coords[3];   
    result = MBI()->get_coords( &line_pt1, 1, coords[1].array() );
    assert(MB_SUCCESS == result);                                            
    result = MBI()->get_coords( &line_pt2, 1, coords[2].array() );           
    assert(MB_SUCCESS == result);                                     
    result = MBI()->get_coords( &pt0, 1, coords[0].array() );          
    assert(MB_SUCCESS == result);                                           

    // project the t_joint between the endpts                  
    // http://en.wikipedia.org/wiki/Vector_projection           
    MBCartVect a = coords[0] - coords[1];                 
    MBCartVect b = coords[2] - coords[1]; 
    parameter    = (a%b)/(b%b);
    MBCartVect c = parameter*b;                    
    projected_coords = c     + coords[1];      
    return MB_SUCCESS;
  }
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,
				     const MBEntityHandle line_pt2,
				     const MBEntityHandle pt0,
				     double &dist_along_edge  ) {

    MBErrorCode result;
    MBCartVect coords[3];   
    result = MBI()->get_coords( &line_pt1, 1, coords[1].array() );
    assert(MB_SUCCESS == result);                                            
    result = MBI()->get_coords( &line_pt2, 1, coords[2].array() );           
    assert(MB_SUCCESS == result);                                     
    result = MBI()->get_coords( &pt0, 1, coords[0].array() );          
    assert(MB_SUCCESS == result);                                           

    // project the t_joint between the endpts                  
    // http://en.wikipedia.org/wiki/Vector_projection           
    MBCartVect a = coords[0] - coords[1];                 
    MBCartVect b = coords[2] - coords[1]; 
    dist_along_edge = a%b / b.length();      
    return MB_SUCCESS;
  }
  

  double area2( const MBEntityHandle pt_a, const MBEntityHandle pt_b,
                const MBEntityHandle pt_c, const MBCartVect plane_normal ) {
    //std::cout << "area2: a=" << pt_a << " b=" << pt_b << " c=" << pt_c << std::endl;
    MBErrorCode result;
    MBCartVect a, b, c;
    result = MBI()->get_coords( &pt_a, 1, a.array() );
    assert(MB_SUCCESS == result);
    result = MBI()->get_coords( &pt_b, 1, b.array() );
    assert(MB_SUCCESS == result);
    result = MBI()->get_coords( &pt_c, 1, c.array() );
    assert(MB_SUCCESS == result);
    MBCartVect d = b - a;
    MBCartVect e = c - a;
    // project onto a plane defined by the plane's normal vector
    return (d*e)%plane_normal;
  }

  // Is point c to the left of line ab?
  bool left( const MBEntityHandle a, const MBEntityHandle b,
	     const MBEntityHandle c, const MBCartVect n ) {
    double area_2 = area2(a,b,c,n);
    //std::cout << "left: a=" << a << " b=" << b << " c=" << c
    //          << " area2=" << area_2 << std::endl;
    if(area_2 > 0) return true;
    else return false;
  } 

  // Is point c to the left of line ab or collinear?
  bool left_on( const MBEntityHandle a, const MBEntityHandle b,
		const MBEntityHandle c, const MBCartVect n ) {
    double area_2 = area2(a,b,c,n);
    //std::cout << "left_on: a=" << a << " b=" << b << " c=" << c
    //          << " area2=" << area_2 << std::endl;
    if(area_2 >= 0) return true;
    else return false;
  } 

  // Are pts a,b,c collinear?
  bool collinear( const MBEntityHandle a, const MBEntityHandle b,
		  const MBEntityHandle c, const MBCartVect n ) {
    double area_2 = area2(a,b,c,n);
    //std::cout << "collinear: a=" << a << " b=" << b << " c=" << c
    //          << " area2=" << area_2 << std::endl;
    if( area_2 ==0) return true;
    else return false;
  } 

  // Exclusive or: T iff exactly one argument is true
  bool logical_xor( const bool x, const bool y ) {
    return (x || y) && !(x && y);
  }

  bool intersect_prop( const MBEntityHandle a, const MBEntityHandle b,
                       const MBEntityHandle c, const MBEntityHandle d,
                       const MBCartVect n ) {
    if( collinear(a,b,c,n) ||
        collinear(a,b,d,n) ||
        collinear(c,d,a,n) ||
        collinear(c,d,b,n) ) {
      return false;
    } else {
      return logical_xor(left(a,b,c,n), left(a,b,d,n)) && 
	logical_xor(left(c,d,a,n), left(c,d,b,n));
    }
  }
      
  bool between( const MBEntityHandle pt_a, const MBEntityHandle pt_b, 
                const MBEntityHandle pt_c, const MBCartVect n) {
    if( !collinear(pt_a,pt_b,pt_c,n) ) return false;

    MBErrorCode result;
    MBCartVect a, b, c;
    result = MBI()->get_coords( &pt_a, 1, a.array() );
    assert(MB_SUCCESS == result);
    result = MBI()->get_coords( &pt_b, 1, b.array() );
    assert(MB_SUCCESS == result);
    result = MBI()->get_coords( &pt_c, 1, c.array() );
    assert(MB_SUCCESS == result);

    // if ab not vertical, check betweenness on x; else on y.
    if(a[0] != b[0]) {
      return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
    }else if(a[1] != b[1]) {
      return ((a[1] <= c[1]) && (c[1] <= b[1])) || ((a[1] >= c[1]) && (c[1] >= b[1]));
    } else {
      return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
    }
  }

  bool intersect( const MBEntityHandle a, const MBEntityHandle b,
		  const MBEntityHandle c, const MBEntityHandle d,
                  const MBCartVect n ) {
    if(intersect_prop(a,b,c,d,n)) return true;
    else if( between(a,b,c,n) ||
             between(a,b,d,n) ||
             between(c,d,a,n) ||
             between(c,d,b,n) ) return true;
    else return false;
  }

  // verts is an ordered polygon of verts
  bool diagonalie( const MBEntityHandle a, const MBEntityHandle b,
                   const MBCartVect n, 
		   const std::vector<MBEntityHandle> verts ) {
    for(unsigned int i=0; i<verts.size(); i++) {
      MBEntityHandle c = verts[i];
      MBEntityHandle c1;
      if(verts.size()-1 == i) c1 = verts[0];
      else                    c1 = verts[i+1];

      if( (c != a) && (c1 != a) &&
          (c != b) && (c1 != b) &&
          intersect( a, b, c, c1, n ) ) {
        //std::cout << "diagonalie a=" << a << " b=" << b << " c="
	//	  << c << " c1=" << c1 << " result=";
	//std::cout << "false" << std::endl;
        return false;
      }
    }
    //std::cout << "diagonalie a=" << a << " b=" << b << " result=";
    //std::cout << "true" << std::endl;
    return true;
  }

  // verts is an ordered polygon of verts
  bool in_cone( const MBEntityHandle a, const MBEntityHandle b,
                const MBCartVect n,
                const std::vector<MBEntityHandle> verts ) { 
    std::vector<MBEntityHandle>::const_iterator a_iter;
    a_iter = find( verts.begin(), verts.end(), a );
    MBEntityHandle a0, a1;
    // a0 is before a
    if(verts.begin() == a_iter) a0 = verts[verts.size()-1];
    else a0 = *(a_iter-1);
    // a1 is after a
    if(verts.end()-1 == a_iter) a1 = verts[0];
    else a1 = *(a_iter+1);

    //std::cout << "in_cone: a=" << a << " b=" << b << " a0="
    //	      << a0 << " a1=" << a1 << std::endl;
    // if a is a convex vertex
    if(left_on(a,a1,a0,n)) return left(a,b,a0,n) && left(b,a,a1,n);

    // else a is reflex
    else return !(left_on(a,b,a1,n) && left_on(b,a,a0,n));
  }

  bool diagonal( const MBEntityHandle a, const MBEntityHandle b,
                 const MBCartVect n,
                 const std::vector<MBEntityHandle> verts ) {
    bool result = in_cone(a,b,n,verts) && in_cone(b,a,n,verts) && diagonalie(a,b,n,verts);
    //std::cout << "diagonal a=" << a << " b=" << b << " result="
    //          << result << std::endl;
    return result;
  }


  // Determine if each vertex is an ear. Input an ordered polygon of verts.
  MBErrorCode ear_init( const std::vector<MBEntityHandle> verts,
                        const MBCartVect n, // plane normal vector
                        std::vector<bool> &is_ear ) {
    if(verts.size() != is_ear.size()) return MB_FAILURE;
    for(unsigned int i=0; i<verts.size(); i++) {
      MBEntityHandle prev, next;
      if(0 == i) prev = verts.back();
      else prev = verts[i-1];
      if(verts.size()-1 == i) next = verts[0];
      else next = verts[i+1];
      is_ear[i] = diagonal(prev,next,n,verts);
      //std::cout << "is_ear[" << i << "]=" << is_ear[i] << std::endl;
    }
    return MB_SUCCESS;
  }



  // Input an ordered polygon of verts and a normal vector of the plane
  // that the polygon is mostly in. The vector is required for orientation.
  MBErrorCode ear_clip_polygon( std::vector<MBEntityHandle> verts,
                                MBCartVect n, 
				MBRange &new_tris ) {

    // initialize the status of ears
    //std::cout << "begin ear clipping----------------------" << std::endl;
    //for(unsigned int i=0; i<verts.size(); i++) {
    //  print_vertex_cubit( verts[i] );
    //}

    //print_loop( verts );
    MBErrorCode result;
    std::vector<bool> is_ear( verts.size() );
    result = ear_init( verts, n, is_ear );
    assert(MB_SUCCESS == result);
 
    // if the polygon intersects itself the algorithm will not stop
    int counter = 0; 
    int n_initial_verts = verts.size();

    while(3 < verts.size()) {
      for(unsigned int i=0; i<verts.size(); i++) {
        if(is_ear[i]) {
          MBEntityHandle v0, v1, v2, v3, v4;
          if(0 == i)      v0 = verts[verts.size()-2];
          else if(1 == i) v0 = verts[verts.size()-1];
          else            v0 = verts[i-2]; 
          if(0 == i) v1 = verts[verts.size()-1];
          else       v1 = verts[i-1];
          v2 = verts[i];
          if(verts.size()-1 == i) v3 = verts[0];
          else                    v3 = verts[i+1];
          if(verts.size()-2 == i)       v4 = verts[0];
          else if (verts.size()-1 == i) v4 = verts[1];
          else                          v4 = verts[i+2];

	  //std::cout << "ear_clip_polygon: triangle=" << std::endl;
          //print_vertex_coords( v1 ); 
          //print_vertex_coords( v2 ); 
          //print_vertex_coords( v3 );
          MBEntityHandle new_tri; 
          MBEntityHandle conn[3] = {v1,v2,v3};
          result = MBI()->create_element( MBTRI, conn, 3, new_tri );
          assert(MB_SUCCESS == result);
          new_tris.insert( new_tri );

          // update ear status
	  if(0 == i) is_ear[verts.size()-1] = diagonal( v0, v3, n, verts );
          else       is_ear[i-1]            = diagonal( v0, v3, n, verts );
          if(verts.size()-1 == i) is_ear[0] = diagonal( v1, v4, n, verts );
	  else                    is_ear[i+1] = diagonal( v1, v4, n, verts );

          // cut off the ear at i
          verts.erase( verts.begin()+i );
          is_ear.erase( is_ear.begin()+i );
          //i--;
          break;
        }
      }

      // If the polygon has intersecting edges this loop will continue until it
      // hits this return.
      //std::cout << "counter=" << counter << " verts.size()=" << verts.size() << std::endl;
      if(counter > n_initial_verts) {
	//std::cout << "ear_clip_polygon: no ears found" << std::endl;
        result = MBI()->delete_entities( new_tris );
        assert(MB_SUCCESS == result);
        new_tris.clear();
        return MB_FAILURE;
      }
      counter++;
    } 
    //std::cout << "ear_clip_polygon: triangle=" << std::endl;
    //print_vertex_coords( verts[0] ); 
    //print_vertex_coords( verts[1] ); 
    //print_vertex_coords( verts[2] );
    MBEntityHandle new_tri; 
    MBEntityHandle conn[3] = {verts[0],verts[1],verts[2]};
    result = MBI()->create_element( MBTRI, conn, 3, new_tri );
    assert(MB_SUCCESS == result);
    new_tris.insert( new_tri );
    
    return MB_SUCCESS; 
  }

  int geom_id_by_handle( const MBEntityHandle set ) {
    MBErrorCode result;
    MBTag id_tag;
    result = MBI()->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                                MB_TYPE_INTEGER, id_tag, 0, true );           
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);                       
    int id;
    result = MBI()->tag_get_data( id_tag, &set, 1, &id );                  
    assert(MB_SUCCESS == result);                           
    return id;
  }
  
  MBErrorCode save_normals( MBRange tris, MBTag normal_tag ) {
    std::vector<MBCartVect> normals(tris.size());
    MBErrorCode result;
    result = triangle_normals( tris, normals );
    assert(MB_SUCCESS == result);

    result = MBI()->tag_set_data(normal_tag, tris, &normals[0]);
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }

  MBErrorCode flip(const MBEntityHandle tri, const MBEntityHandle vert0, 
		   const MBEntityHandle vert2, const MBEntityHandle surf_set) {

    // get the triangles in the surface. The tri and adj_tri must be in the surface.
    MBRange surf_tris;
    MBEntityHandle result = MBI()->get_entities_by_type( surf_set, MBTRI, surf_tris);
    assert(MB_SUCCESS == result);

    // get the triangle across the edge that will be flipped
    MBRange adj_tri;
    MBEntityHandle edge[2] = {vert0, vert2};
    result = MBI()->get_adjacencies( edge, 2, 2, false, adj_tri );
    assert(MB_SUCCESS == result);
    adj_tri = intersect(adj_tri, surf_tris);
    assert(2 == adj_tri.size());
    adj_tri.erase( tri );
    print_triangle( adj_tri.front(), false );
    //result = MBI()->list_entity( adj_tri.front() );
    //assert(MB_SUCCESS == result);

    // get the remaining tri vert
    MBRange tri_verts;
    result = MBI()->get_adjacencies( &tri, 1, 0, false, tri_verts );
    assert(MB_SUCCESS == result);
    assert(3 == tri_verts.size());
    tri_verts.erase(vert0);
    tri_verts.erase(vert2);
    assert(1 == tri_verts.size());
    MBEntityHandle vert1 = tri_verts.front();

    // get the remaining adj_tri vert
    MBRange adj_tri_verts;
    result = MBI()->get_adjacencies( &adj_tri.front(), 1, 0, false, adj_tri_verts );
    assert(MB_SUCCESS == result);
    assert(3 == adj_tri_verts.size());
    adj_tri_verts.erase(vert0);
    adj_tri_verts.erase(vert2);
    assert(1 == adj_tri_verts.size());
    MBEntityHandle vert3 = adj_tri_verts.front();

    // original tri_conn    = {vert0, vert1, vert2}
    // original adj_tri_conn= {vert2, vert3, vert0}
 
    // set the new connectivity
    MBEntityHandle tri_conn[3] = {vert0, vert1, vert3};
    result = MBI()->set_connectivity( tri, tri_conn, 3 );
    assert(MB_SUCCESS == result);
    MBEntityHandle adj_tri_conn[3] = {vert1, vert2, vert3};
    result = MBI()->set_connectivity( adj_tri.front(), adj_tri_conn, 3 );
    assert(MB_SUCCESS == result);
    print_triangle( tri, false );
    print_triangle( adj_tri.front(), false );
    return MB_SUCCESS;
  }

  MBErrorCode ordered_verts_from_ordered_edges( const std::vector<MBEntityHandle> ordered_edges,
                                                std::vector<MBEntityHandle> &ordered_verts ) {
    MBErrorCode result;
    ordered_verts.clear();
    ordered_verts.reserve(ordered_edges.size()+1);

    // Save the back of the previous edge to check for continuity.
    MBEntityHandle previous_back_vert;
    
    for(std::vector<MBEntityHandle>::const_iterator i=ordered_edges.begin(); 
        i!=ordered_edges.end(); i++) {
      const MBEntityHandle *conn;
      int n_verts;
      result = MBI()->get_connectivity( *i, conn, n_verts);
      assert(MB_SUCCESS == result);
      assert(2 == n_verts);
      if(ordered_edges.begin() == i) {
        ordered_verts.push_back(conn[0]);
      } else {
        assert(previous_back_vert == conn[0]);
      }
      ordered_verts.push_back(conn[1]);
      previous_back_vert = conn[1];
    }
    return MB_SUCCESS;
  }

  /* Find the distance between two arcs. Assume that their endpoints are somewhat
     close together. */
  MBErrorCode dist_between_arcs( bool debug,
                                 const std::vector<MBEntityHandle> arc0,
                                 const std::vector<MBEntityHandle> arc1,
                                 double &dist ) {
    dist = 0;

    // Special Case: arcs have no verts.
    if( arc0.empty() || arc1.empty() ) {
      std::cout << "arc has no vertices" << std::endl;
      return MB_FAILURE;
    }

    //print_loop(arc0);
    //print_loop(arc1);

    // for simplicity, put arcs into the same structure
    std::vector<MBEntityHandle> arcs[2] = {arc0, arc1};

    // Special Case: Remove duplicate vert handles
    for(unsigned int i=0; i<2; ++i) {
      if( 2>arcs[i].size() ) continue;
      for(std::vector<MBEntityHandle>::iterator j=arcs[i].begin()+1; j!=arcs[i].end(); ++j) {
        if(*j == *(j-1)) {
	  if(debug) {
            gen::print_loop( arcs[i] );
	    std::cout << "dist_between_arcs: found duplicate vert handle in arc" << std::endl;
          }
          j = arcs[i].erase(j) - 1;
        }
      }
    }
    
    // get the coords in one call per arc. For speed, do not ask MOAB again for coords.
    MBErrorCode result;
    std::vector<MBCartVect> coords[2];
    for(unsigned int i=0; i<2; i++) {
      coords[i].resize( arcs[i].size() );
      result = MBI()->get_coords( &arcs[i][0], arcs[i].size(), coords[i][0].array());
      assert(MB_SUCCESS == result);
    }

    // Special case: arc has 1 vert or a length of zero
    bool point_arc_exists = false;
    unsigned int point_arc_index;
    for(unsigned int i=0; i<2; ++i) {
      if( 1==arcs[i].size() ) {
        point_arc_exists = true;
        point_arc_index  = i;
        break;
      } else if ( 0==length(arcs[i]) ) {
        point_arc_exists = true;
        point_arc_index  = i;
        break;
      }
    }
    // If the special case exists, we can still get an average distance
    if(point_arc_exists) {
      unsigned int other_arc_index = (0==point_arc_index) ? 1 : 0;
      // Both are point arcs
      if(1==arcs[other_arc_index].size()) {
        dist = dist_between_verts( arcs[point_arc_index].front(), arcs[other_arc_index].front());
        return MB_SUCCESS;
	// The other arc has more than one point
      } else {
        double area = 0.0;
        for(unsigned int i=0; i<arcs[other_arc_index].size()-1; ++i) {
          area += fabs(triangle_area( coords[other_arc_index][i], 
                                      coords[other_arc_index][i+1], 
                                      coords[point_arc_index].front() ));
        }
        dist = area / gen::length(arcs[other_arc_index]);
        return MB_SUCCESS;
      }
    }

    // get the arc length and parametric coordinates. The parametrize the arcs by
    // length from 0 to 1.
    double arc_len[2] = {0.0};
    std::vector<double> params[2];
    for(unsigned int i=0; i<2; i++) {
      params[i].resize( arcs[i].size() );
      params[i][0] = 0.0;
      for(unsigned int j=0; j<coords[i].size()-1; j++) {
        double d = dist_between_verts( coords[i][j], coords[i][j+1] );
        arc_len[i] += d;
        params[i][1+j] = arc_len[i];
      }
      for(unsigned int j=0; j<coords[i].size(); j++) {
        params[i][j] /= arc_len[i];
	//std::cout << "params[" << i << "][" << j << "]=" << params[i][j] << std::endl;
      }
    }

    // Merge the two sets of parameters into one ordered set without duplicates.
    std::vector<double> mgd_params;
    mgd_params.reserve( params[0].size() + params[1].size() );
    unsigned int a=0, b=0;
    while(a<params[0].size() && b<params[1].size()) {
      if       (params[0][a] < params[1][b]) {
        mgd_params.push_back( params[0][a] );
        a++;
      } else if(params[0][a] > params[1][b]) {
        mgd_params.push_back( params[1][b] );
        b++;
      } else {
        mgd_params.push_back( params[0][a] );
        a++;
        b++;
      }
    }

    for(unsigned int i=0; i<mgd_params.size(); i++) {
      //std::cout << "mgd_params[" << i << "]=" << mgd_params[i] << std::endl;
    }

    // Insert new points to match the other arcs points, by parameter.
    for(unsigned int i=0; i<2; i++) {
      for(unsigned int j=0; j<mgd_params.size(); j++) {
	//std::cout << "params[" << i << "][" << j << "]=" << params[i][j] 
        //          << " mgd_params[" << j << "]=" << mgd_params[j] << std::endl;
        if(params[i][j] > mgd_params[j]) {
          double ratio = (mgd_params[j]-params[i][j-1]) / (params[i][j]-params[i][j-1]);
	  //std::cout << "j=" << j << " ratio=" << ratio << std::endl;
          MBCartVect pt = coords[i][j-1] + ratio*(coords[i][j]-coords[i][j-1]);
	  coords[i].insert( coords[i].begin()+j, pt);
          params[i].insert( params[i].begin()+j, mgd_params[j]); 
        }
      }
    }


    // Each arc should now have the same number of points
    if(coords[0].size() != coords[1].size()) {
      for(unsigned int i=0; i<2; i++) {
        for(unsigned int j=0; j<coords[i].size(); j++) {
	  std::cout << "coords[" << i << "][" << j << "]=" << coords[i][j] << std::endl;
        }
      }
    }
    assert( coords[0].size() == coords[1].size() );

    // Get the area between arcs. Use absolute value to prevent cancelling.
    double area = 0;
    for(unsigned int i=0; i<coords[0].size()-1; i++) {
      area += fabs(triangle_area( coords[0][i], coords[0][i+1], coords[1][i+1] ));
      //std::cout << "area0=" << area << std::endl;
      area += fabs(triangle_area( coords[0][i], coords[1][i+1], coords[1][i] ));
      //std::cout << "area1=" << area << std::endl;
    }

    // Divide the area by the average length to get the average distance between arcs.
    dist = fabs(2*area / (arc_len[0] + arc_len[1] ));
    //std::cout << "dist_between_arcs=" << dist << std::endl;
    return MB_SUCCESS;
  }


    // qsort struct comparision function        
    // If the first handle is the same, compare the second                                        
    int compare_edge(const void *a, const void *b) {                                 
      struct edge *ia = (struct edge *)a;           
      struct edge *ib = (struct edge *)b;
      if(ia->v0 == ib->v0) {        
        return (int)(ia->v1 - ib->v1);                   
      } else {
        return (int)(ia->v0 - ib->v0);                   
      }
      /* float comparison: returns negative if b > a           
      and positive if a > b. We multiplied result by 100.0        
      to preserve decimal fraction */                
    }  

  // WARNING: This skinner goes 10x faster by assuming that no edges already exist
  // in the MOAB instance. Otherwise checking to see if an edge exists before
  // creating a new one if very slow. This is partly the reason that MBSkinner is
  // very slow.
  MBErrorCode find_skin( MBRange tris, const int dim, 
 //                         std::vector<std::vector<MBEntityHandle> > &skin_edges, 
                         MBRange &skin_edges,                         
                         const bool temp_bool ) {    

    const bool local_debug = false;
    //MBSkinner tool(MBI());
    //MBRange skin_verts;
    //return tool.find_skin( tris, dim, skin_edges, temp_bool );
    //return tool.find_skin_vertices( tris, skin_verts, &skin_edges, true );
  
    if(1 != dim) return MB_NOT_IMPLEMENTED;
    if(MBTRI != MBI()->type_from_handle(tris.front())) return MB_NOT_IMPLEMENTED;

    skin_edges.clear();
    if(tris.empty()) return MB_ENTITY_NOT_FOUND;

    // This implementation gets some of its speed due to not checking for edges
    MBErrorCode result;
    int n_edges;
    result = MBI()->get_number_entities_by_type( 0, MBEDGE, n_edges );
    assert(MB_SUCCESS == result);
    if(0 != n_edges) {
      MBRange temp_edges;
      result = MBI()->get_entities_by_type( 0, MBEDGE, temp_edges);
      assert(MB_SUCCESS == result);
      result = MBI()->list_entities( temp_edges );
      assert(MB_SUCCESS == result);
    }      
    assert(0 == n_edges);

    // Get connectivity. Do not create edges.
    edge *edges = new edge[3*tris.size()];
    int n_verts;
    int ii = 0;
    for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
      const MBEntityHandle *conn;
      result = MBI()->get_connectivity( *i, conn, n_verts );
      assert(MB_SUCCESS == result);
      assert(3 == n_verts);
      // shouldn't be degenerate
      assert(conn[0] != conn[1]);
      assert(conn[1] != conn[2]);
      assert(conn[2] != conn[0]);
      // make edges
      edges[3*ii+0].v0 = conn[0];
      edges[3*ii+0].v1 = conn[1];
      edges[3*ii+1].v0 = conn[1];
      edges[3*ii+1].v1 = conn[2];
      edges[3*ii+2].v0 = conn[2];
      edges[3*ii+2].v1 = conn[0];
      ii++;
    }

    // Change the first handle to be lowest
    for(unsigned int i=0; i<3*tris.size(); ++i) {
      if(edges[i].v0 > edges[i].v1) {
        MBEntityHandle temp = edges[i].v0;
        edges[i].v0 = edges[i].v1;
        edges[i].v1 = temp;
      }
    }

    // Sort by first handle, then second handle. Do not sort the extra edge on the
    // back.
    qsort(edges, 3*tris.size(), sizeof(struct edge), compare_edge);    

    // Go through array, saving edges that are not paired.
    for(unsigned int i=0; i<3*tris.size(); i++) {
      // If the last edge has not been paired, create it. This avoids overrunning
      // the edges array with i+1.
      if(3*tris.size()-1 == i) {
        const MBEntityHandle conn[2] = {edges[i].v0, edges[i].v1};
        MBEntityHandle edge;
        result = MBI()->create_element( MBEDGE, conn, 2, edge );
        assert(MB_SUCCESS == result);
        skin_edges.insert(edge);
    
      // If a match exists, skip ahead
      } else if(edges[i].v0==edges[i+1].v0 && edges[i].v1==edges[i+1].v1) {
        i++;
        // test to make sure surface is manifold
        if(3*tris.size() != i+1) { // avoid overrunning array
          while( edges[i].v0==edges[i+1].v0 && edges[i].v1==edges[i+1].v1 ) {
	    if(local_debug) {
              std::cout << "find_skin WARNING: non-manifold edge" << std::endl;
              MBI()->list_entity( edges[i].v0 );
              MBI()->list_entity( edges[i].v1 );
	    }
            ++i;
          }
        }
        //assert( !(edges[i].v0==edges[i+2].v0 && edges[i].v1==edges[i+2].v1) );
	// otherwise a skin edge has been found
      } else {
        const MBEntityHandle conn[2] = {edges[i].v0, edges[i].v1};
	//std::vector<MBEntityHandle> the_edge(conn, conn+2);     
        //skin_edges.push_back(the_edge);
        // see if the edge already exists
        //MBRange the_edge;
        //result = MBI()->get_adjacencies( conn, 2, 1, true, the_edge );
        //assert(MB_SUCCESS == result);
        //assert(1 == the_edge.size());
        MBEntityHandle edge;
        result = MBI()->create_element( MBEDGE, conn, 2, edge );
        //if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
        assert(MB_SUCCESS == result);
        skin_edges.insert( edge );
      } 
    }
    delete[] edges;
    return MB_SUCCESS;
  }
  /*  MBErrorCode find_skin( MBRange tris, const int dim, MBRange &skin_edges, const bool temp ) {
    std::vector<std::vector<MBEntityHandle> > skin_edges_vctr;
    MBErrorCode result = find_skin( tris, dim, skin_edges_vctr, temp );
    assert(MB_SUCCESS == result);
    for(std::vector<std::vector<MBEntityHandle> >::const_iterator i=skin_edges_vctr.begin(); 
        i!=skin_edges_vctr.end(); i++) {
      MBEntityHandle edge;
      result = MBI()->create_element( MBEDGE, &(*i)[0], 2, edge );
      if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
      assert(MB_SUCCESS == result);
      skin_edges.insert( edge );
    }
    return MB_SUCCESS;
    }*/

// calculate volume of polyhedron
// Copied from DagMC, without index_by_handle. The dagmc function will
// segfault if build_indices is not first called. For sealing there is
// no need to build_indices.
MBErrorCode measure_volume( const MBEntityHandle volume, double& result )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle> surfaces, surf_volumes;
  result = 0.0;
  
   // don't try to calculate volume of implicit complement
  //if (volume == impl_compl_handle) {
  //  result = 1.0;
  //  return MB_SUCCESS;
  //}

    // get surfaces from volume
  rval = MBI()->get_child_meshsets( volume, surfaces );
  if (MB_SUCCESS != rval) return rval;
  
    // get surface senses
  std::vector<int> senses( surfaces.size() );
  moab::DagMC &dagmc = *moab::DagMC::instance( MBI() );
  rval = dagmc.surface_sense( volume, surfaces.size(), &surfaces[0], &senses[0] );
  if (MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }
  
  for (unsigned i = 0; i < surfaces.size(); ++i) {
      // skip non-manifold surfaces
    if (!senses[i])
      continue;
    
      // get triangles in surface
    MBRange triangles;
    rval = MBI()->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (MB_SUCCESS != rval) 
      return rval;
    if (!triangles.all_of_type(MBTRI)) {
      std::cout << "WARNING: Surface " << geom_id_by_handle(surfaces[i])
                << " contains non-triangle elements. Volume calculation may be incorrect."
                << std::endl;
      triangles.clear();
      rval = MBI()->get_entities_by_type( surfaces[i], MBTRI, triangles );
      if (MB_SUCCESS != rval) return rval;
    }
    
      // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const MBEntityHandle *conn;
    int len;
    MBCartVect coords[3];
    for (MBRange::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = MBI()->get_connectivity( *j, conn, len, true );
      if (MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = MBI()->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;
    
      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }
  
  result /= 6.0;
  return MB_SUCCESS;
}

  /* Calculate the signed volumes beneath the surface (x 6.0). Use the triangle's
     cannonical sense. Do not take sense tags into account. Code taken from 
     DagMC::measure_volume. */    
  MBErrorCode get_signed_volume( const MBEntityHandle surf_set, double &signed_volume) {
    MBErrorCode rval;                                     
    MBRange tris;                          
    rval = MBI()->get_entities_by_type( surf_set, MBTRI, tris );
    if(MB_SUCCESS != rval) return rval;       
    signed_volume = 0.0;
    const MBEntityHandle *conn;                                  
    int len;
    MBCartVect coords[3];                                              
    for (MBRange::iterator j = tris.begin(); j != tris.end(); ++j) {         
      rval = MBI()->get_connectivity( *j, conn, len, true );           
      if (MB_SUCCESS != rval) return rval;            
      assert(3 == len);                    
      rval = MBI()->get_coords( conn, 3, coords[0].array() );              
      if (MB_SUCCESS != rval) return rval;         
                            
      coords[1] -= coords[0];                              
      coords[2] -= coords[0];                     
      signed_volume += (coords[0] % (coords[1] * coords[2]));       
    }                                     
    return MB_SUCCESS;                              
  }                 

  MBErrorCode measure( const MBEntityHandle set, const MBTag geom_tag, double &size ) {
    MBErrorCode result;
    int dim;
    result = MBI()->tag_get_data( geom_tag, &set, 1, &dim );                  
    assert(MB_SUCCESS == result);                           
    if(0 == dim) {
      std::cout << "measure: cannot measure vertex" << std::endl;
      return MB_FAILURE;

    } else if(1 == dim) {
      std::vector<MBEntityHandle> vctr;
      result = arc::get_meshset( set, vctr );
      assert(MB_SUCCESS == result);
      size = length( vctr );

    } else if(2 == dim) {
      MBRange tris;
      result = MBI()->get_entities_by_type( set, MBTRI, tris );
      assert(MB_SUCCESS == result);
      size = triangle_area( tris );

    } else if(3 == dim) {
      //moab::DagMC &dagmc = *moab::DagMC::instance( MBI() );
      //result = dagmc.measure_volume( set, size );
      result = measure_volume( set, size );
      if(MB_SUCCESS != result) {
        std::cout << "result=" << result << " vol_id=" 
                  << gen::geom_id_by_handle(set) << std::endl;
      }
    } else {
      std::cout << "measure: incorrect dimension" << std::endl;
      return MB_FAILURE;
    }
    return MB_SUCCESS;
  }

  // From CGMA/builds/dbg/include/CubitDefines, CUBIT_UNKNOWN=-1, CUBIT_FORWARD=0, CUBIT_REVERSED=1
  MBErrorCode get_curve_surf_sense( const MBEntityHandle surf_set, const MBEntityHandle curve_set,
                                    int &sense ) {
    std::vector<MBEntityHandle> surfs;
    std::vector<int> senses;
    MBErrorCode rval;
    moab::GeomTopoTool gt( MBI(), false);
    rval = gt.get_senses( curve_set, surfs, senses );
    if(gen::error(MB_SUCCESS!=rval,"failed to get_senses")) return rval;
    
    unsigned counter = 0;
    for(unsigned i=0; i<surfs.size(); ++i) {
      if(surf_set==surfs[i]) {
        sense = senses[i];
        ++counter;
      }
    }

    // special case: parent surface does not exist
    if(gen::error(0==counter,"failed to find a surf in sense list")) return MB_FAILURE;

    // special case: ambiguous
    if(1<counter) sense = -1;
    
    return MB_SUCCESS;
  }

}
