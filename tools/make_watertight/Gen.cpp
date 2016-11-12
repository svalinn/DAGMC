#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/GeomTopoTool.hpp"
#include "moab/FileOptions.hpp"

#include "Gen.hpp"
#include "moab/Skinner.hpp"
#include "moab/Range.hpp"


const char GEOM_SENSE_2_TAG_NAME[] = "GEOM_SENSE_2";
const char GEOM_SENSE_N_ENTS_TAG_NAME[] = "GEOM_SENSE_N_ENTS";
const char GEOM_SENSE_N_SENSES_TAG_NAME[] = "GEOM_SENSE_N_SENSES";

void Gen::moab_printer(moab::ErrorCode error_code)
{
  if ( error_code == moab::MB_INDEX_OUT_OF_RANGE ) {
    std::cerr << "ERROR: moab::MB_INDEX_OUT_OF_RANGE" << std::endl;
  }
  if ( error_code == moab::MB_MEMORY_ALLOCATION_FAILED ) {
    std::cerr << "ERROR: moab::MB_MEMORY_ALLOCATION_FAILED" << std::endl;
  }
  if ( error_code == moab::MB_ENTITY_NOT_FOUND ) {
    std::cerr << "ERROR: moab::MB_ENTITY_NOT_FOUND" << std::endl;
  }
  if ( error_code == moab::MB_MULTIPLE_ENTITIES_FOUND ) {
    std::cerr << "ERROR: moab::MB_MULTIPLE_ENTITIES_FOUND" << std::endl;
  }
  if ( error_code == moab::MB_TAG_NOT_FOUND ) {
    std::cerr << "ERROR: moab::MB_TAG_NOT_FOUND" << std::endl;
  }
  if ( error_code == moab::MB_FILE_DOES_NOT_EXIST ) {
    std::cerr << "ERROR: moab::MB_FILE_DOES_NOT_EXIST" << std::endl;
  }
  if ( error_code == moab::MB_FILE_WRITE_ERROR ) {
    std::cerr << "ERROR: moab::MB_FILE_WRITE_ERROR" << std::endl;
  }
  if ( error_code == moab::MB_ALREADY_ALLOCATED ) {
    std::cerr << "ERROR: moab::MB_ALREADY_ALLOCATED" << std::endl;
  }
  if ( error_code == moab::MB_VARIABLE_DATA_LENGTH ) {
    std::cerr << "ERROR: moab::MB_VARIABLE_DATA_LENGTH" << std::endl;
  }
  if ( error_code == moab::MB_INVALID_SIZE ) {
    std::cerr << "ERROR: moab::MB_INVALID_SIZE" << std::endl;
  }
  if ( error_code == moab::MB_UNSUPPORTED_OPERATION ) {
    std::cerr << "ERROR: moab::MB_UNSUPPORTED_OPERATION" << std::endl;
  }
  if ( error_code == moab::MB_UNHANDLED_OPTION ) {
    std::cerr << "ERROR: moab::MB_UNHANDLED_OPTION" << std::endl;
  }
  if ( error_code == moab::MB_FAILURE ) {
    std::cerr << "ERROR: moab::MB_FAILURE" << std::endl;
  }
  return;
}





void Gen::print_vertex_cubit( const moab::EntityHandle vertex )
{

  moab::ErrorCode result;
  double coords[3];
  int n_precision = 20;
  result = MBI()->get_coords( &vertex, 1, coords );
  assert(moab::MB_SUCCESS == result);
  std::cout << "  create vertex "
            << std::setprecision(n_precision)
            << coords[0] << " " << coords[1] << " " << coords[2]
            << std::endl;
}

void Gen::print_vertex_coords( const moab::EntityHandle vertex )
{

  moab::ErrorCode result;
  double coords[3];
  result = MBI()->get_coords( &vertex, 1, coords );
  if(moab::MB_SUCCESS!=result) std::cout << "vert=" << vertex << std::endl;
  assert(moab::MB_SUCCESS == result);
  std::cout << "    vertex " << vertex << " coords= ("
            << coords[0] << "," << coords[1] << "," << coords[2] << ")"
            << std::endl;
}

void Gen::print_triangles( const moab::Range tris )
{
  for(moab::Range::const_iterator i=tris.begin(); i!=tris.end(); i++) {
    print_triangle( *i, false );
  }
}
// If the edges of the tri are ambiguous, do not print edges!
void Gen::print_triangle( const moab::EntityHandle tri, bool print_edges )
{
  moab::ErrorCode result;
  double area;
  result = triangle_area( tri, area );
  assert(moab::MB_SUCCESS == result);
  std::cout << "    triangle " << tri << " area=" << area << std::endl;
  const moab::EntityHandle *conn;
  int n_verts;
  result = MBI()->get_connectivity( tri, conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert(3 == n_verts);
  for(int i=0; i<3; i++) print_vertex_coords( conn[i] );

  if(print_edges) {
    moab::Range edges;
    result = MBI()->get_adjacencies( &tri, 1, 1, true, edges );
    if(moab::MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
    assert(moab::MB_SUCCESS == result);
    for(moab::Range::iterator i=edges.begin(); i!=edges.end(); i++) {
      print_edge( *i );
    }
  }
}

void Gen::print_edge( const moab::EntityHandle edge )
{
  const moab::EntityHandle *conn;
  int n_verts;
  std::cout << "    edge " << edge << std::endl;
  moab::ErrorCode result = MBI()->get_connectivity( edge, conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert(2 == n_verts);
  print_vertex_coords( conn[0] );
  print_vertex_coords( conn[1] );

  moab::Range tris;
  result = MBI()->get_adjacencies( &edge, 1, 2, false, tris );
  assert(moab::MB_SUCCESS == result);
  std::cout << "     tris: ";
  for(moab::Range::iterator i=tris.begin(); i!=tris.end(); i++) {
    std::cout << *i << " ";
  }
  std::cout << std::endl;
}

void Gen::print_range( const moab::Range range )
{
  std::cout << "print range:" << std::endl;
  moab::Range::iterator i;
  for(i=range.begin(); i!=range.end(); i++) {
    std::cout << "    " << *i << std::endl;
  }
}

void Gen::print_range_of_edges( const moab::Range range )
{
  std::cout << "print range:" << std::endl;
  moab::Range::const_iterator i;
  for(i=range.begin(); i!=range.end(); i++) {
    print_edge( *i );
  }
}


void Gen::print_vertex_count(const moab::EntityHandle input_meshset)
{

  // get the range of facets of the surface meshset
  moab::ErrorCode result;
  moab::Range vertices;
  result = MBI()->get_entities_by_type(0, moab::MBVERTEX, vertices);
  assert( moab::MB_SUCCESS == result );

  std::cout<< "    " << vertices.size() << " vertices found." << std::endl;
}

void Gen::print_arcs( const std::vector< std::vector<moab::EntityHandle> > arcs )
{
  for(unsigned int i=0; i<arcs.size(); i++) {
    std::cout << "arc " << i << std::endl;
    print_loop( arcs[i] );
  }
}

void Gen::print_arc_of_edges( const std::vector<moab::EntityHandle> arc_of_edges )
{

  moab::ErrorCode result;
  std::vector<moab::EntityHandle>::const_iterator i;
  double dist = 0;
  for( i=arc_of_edges.begin(); i!=arc_of_edges.end(); i++ ) {
    int n_verts;
    const moab::EntityHandle *conn;
    result = MBI()->get_connectivity( *i, conn, n_verts );
    assert(moab::MB_SUCCESS == result);
    assert( 2 == n_verts );
    dist += dist_between_verts( conn[0], conn[1] );
    print_vertex_coords( conn[0] );
    print_vertex_coords( conn[1] );
  }
  std::cout << "  dist= " << dist << std::endl;
}

void Gen::print_loop( const std::vector<moab::EntityHandle> loop_of_verts )
{

  std::cout << "  size=" << loop_of_verts.size() << std::endl;
  double dist = 0;
  for(unsigned int i=0; i<loop_of_verts.size(); i++) {
    print_vertex_coords( loop_of_verts[i] );
    if(i != loop_of_verts.size()-1) {
      dist += dist_between_verts( loop_of_verts[i], loop_of_verts[i+1] );
    }
  }
  std::cout << "  dist=" << dist << std::endl;
}


/// Return the closest vertex to the arc.
/// For efficiency: only get_coords on the reference vertex once
///                 if specified, limit search length along curve
moab::ErrorCode Gen::find_closest_vert( const moab::EntityHandle reference_vert,
                                        const std::vector<moab::EntityHandle> arc_of_verts,
                                        unsigned &position,
                                        const double dist_limit )
{
  moab::ErrorCode rval;
  const bool debug = false;
  double min_dist_sqr = std::numeric_limits<double>::max();
  moab::CartVect ref_coords;
  rval = MBI()->get_coords( &reference_vert, 1, ref_coords.array() );
  MB_CHK_SET_ERR(rval,"failed to get ref coords");
  double length = 0;
  moab::CartVect prev_coords;

  for(unsigned i=0; i<arc_of_verts.size(); ++i) {
    moab::CartVect coords;
    rval = MBI()->get_coords( &arc_of_verts[i], 1, coords.array() );
    MB_CHK_SET_ERR(rval,"failed to get coords");

    // use dist_limit to exit early; avoid checking the entire arc
    if(0!=i) {
      moab::CartVect temp = prev_coords - coords;
      length += temp.length();
      if(length>dist_limit && debug)
        std::cout << "length=" << length << " dist_limit=" << dist_limit << std::endl;
      if(length > dist_limit) return moab::MB_SUCCESS;
    }
    prev_coords = coords;

    // get distance to ref_vert
    moab::CartVect temp = ref_coords - coords;
    double dist_sqr = temp.length_squared();
    if(dist_sqr < min_dist_sqr) {
      position = i;
      min_dist_sqr = dist_sqr;
      if(debug) std::cout << "min_dist_sqr=" << min_dist_sqr << std::endl;
    }
  }

  return moab::MB_SUCCESS;
}



// Return the closest vert and all within tol. This is needed because sometimes
// the correct vert is not the closest. For example, iter_surf4010 the skin
// loop has the same point in it twice, at two different locations (center of L).
// This ensure that both are returned as candidates.
moab::ErrorCode Gen::find_closest_vert( const double tol,
                                        const moab::EntityHandle reference_vert,
                                        const std::vector<moab::EntityHandle> loop_of_verts,
                                        std::vector<unsigned> &positions,
                                        std::vector<double> &dists)
{

  moab::ErrorCode rval;
  positions.clear();
  dists.clear();
  const double TOL_SQR = tol*tol;
  unsigned min_pos;
  double sqr_min_dist = std::numeric_limits<double>::max();
  for(unsigned int i=0; i<loop_of_verts.size(); i++) {
    double sqr_dist = std::numeric_limits<double>::max();
    rval = squared_dist_between_verts(reference_vert, loop_of_verts[i], sqr_dist);
    MB_CHK_SET_ERR(rval,"could not get dist");
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


  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::merge_vertices( moab::Range verts /* in */, const double tol /* in */ )
{

  moab::ErrorCode result;
  const double SQR_TOL = tol*tol;
  // Clean up the created tree, and track verts so that if merged away they are
  // removed from the tree.
  moab::AdaptiveKDTree kdtree(MBI()); //, true, 0, moab::MESHSET_TRACK_OWNER);
  // initialize the KD Tree
  moab::EntityHandle root;
  const char settings[]="MAX_PER_LEAF=6;MAX_DEPTH=50;SPLITS_PER_DIR=1;PLANE_SET=2;MESHSET_FLAGS=0x1;TAG_NAME=0";
  moab::FileOptions fileopts(settings);
  // builds the KD Tree, making the moab::EntityHandle root the root of the tree
  result = kdtree.build_tree( verts, &root, &fileopts);
  assert(moab::MB_SUCCESS == result);
  // create tree iterator to loop over all verts in the tree
  moab::AdaptiveKDTreeIter tree_iter;
  kdtree.get_tree_iterator( root, tree_iter );
  for(moab::Range::iterator i=verts.begin(); i!=verts.end(); ++i) {
    double from_point[3];
    result = MBI()->get_coords( &(*i), 1, from_point);
    assert(moab::MB_SUCCESS == result);
    std::vector<moab::EntityHandle> leaves_out;
    result = kdtree.distance_search( from_point, tol, leaves_out, root);
    assert(moab::MB_SUCCESS == result);
    for(unsigned int j=0; j<leaves_out.size(); j++) {
      std::vector<moab::EntityHandle> leaf_verts;
      result = MBI()->get_entities_by_type( leaves_out[j], moab::MBVERTEX, leaf_verts);
      assert(moab::MB_SUCCESS == result);
      if(100 < leaf_verts.size()) std::cout << "*i=" << *i << " leaf_verts.size()=" << leaf_verts.size() << std::endl;
      for(unsigned int k=0; k<leaf_verts.size(); k++) {
        if( leaf_verts[k] == *i ) continue;
        double sqr_dist;
        result = squared_dist_between_verts( *i, leaf_verts[k], sqr_dist);
        assert(moab::MB_SUCCESS == result);

        if(SQR_TOL >= sqr_dist) {
          // The delete_vert is automatically remove from the tree because it
          // uses tracking meshsets. merge_verts checks for degenerate tris.
          // Update the list of leaf verts to prevent stale handles.
          std::vector<moab::EntityHandle> temp_arc;
          moab::EntityHandle keep_vert   = *i;
          moab::EntityHandle delete_vert = leaf_verts[k];
          result = merge_verts( keep_vert, delete_vert, leaf_verts, temp_arc );
          assert(moab::MB_SUCCESS == result);
          // Erase delete_vert from verts
          // Iterator should remain valid because delete_vert > keep_vert handle.
          verts.erase( delete_vert );
        }
      }
    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::squared_dist_between_verts( const moab::EntityHandle v0,
    const moab::EntityHandle v1,
    double &d)
{
  moab::ErrorCode result;
  moab::CartVect coords0, coords1;
  result = MBI()->get_coords( &v0, 1, coords0.array() );
  if(moab::MB_SUCCESS != result) {
    std::cout << "dist_between_verts: get_coords on v0=" << v0 << " result="
              << result << std::endl;
    return result;
  }
  result = MBI()->get_coords( &v1, 1, coords1.array() );
  if(moab::MB_SUCCESS != result) {
    std::cout << "dist_between_verts: get_coords on v1=" << v1 << " result="
              << result << std::endl;
    return result;
  }
  const moab::CartVect diff = coords0 - coords1;
  d = diff.length_squared();
  return moab::MB_SUCCESS;
}

double Gen::dist_between_verts( const moab::CartVect v0, const moab::CartVect v1 )
{
  moab::CartVect v2 = v0 - v1;
  return v2.length();
}
moab::ErrorCode Gen::dist_between_verts( const moab::EntityHandle v0, const moab::EntityHandle v1, double &d)
{
  moab::ErrorCode result;
  moab::CartVect coords0, coords1;
  result = MBI()->get_coords( &v0, 1, coords0.array() );
  if(moab::MB_SUCCESS != result) {
    std::cout << "dist_between_verts: get_coords on v0=" << v0 << " result="
              << result << std::endl;
    return result;
  }
  result = MBI()->get_coords( &v1, 1, coords1.array() );
  if(moab::MB_SUCCESS != result) {
    std::cout << "dist_between_verts: get_coords on v1=" << v1 << " result="
              << result << std::endl;
    return result;
  }
  d = dist_between_verts( coords0, coords1 );
  return moab::MB_SUCCESS;
}

double Gen::dist_between_verts( double coords0[], double coords1[] )
{
  return sqrt( (coords0[0]-coords1[0])*(coords0[0]-coords1[0]) +
               (coords0[1]-coords1[1])*(coords0[1]-coords1[1]) +
               (coords0[2]-coords1[2])*(coords0[2]-coords1[2]) );
}
double Gen::dist_between_verts( moab::EntityHandle vert0, moab::EntityHandle vert1 )
{
  double coords0[3], coords1[3];
  moab::ErrorCode result;
  result = MBI()->get_coords( &vert0, 1, coords0 );
  if(moab::MB_SUCCESS!=result) std::cout << "result=" << result << " vert="
                                           << vert0 << std::endl;
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &vert1, 1, coords1 );
  if(moab::MB_SUCCESS!=result) std::cout << "result=" << result << " vert="
                                           << vert1 << std::endl;
  assert(moab::MB_SUCCESS == result);
  return dist_between_verts( coords0, coords1 );
}

// Return the length of the curve defined by moab::MBEDGEs or ordered moab::MBVERTEXs.
double Gen::length( std::vector<moab::EntityHandle> edges )
{
  if(edges.empty()) return 0;

  moab::ErrorCode result;
  std::vector<moab::EntityHandle>::iterator i;
  double dist = 0;
  moab::EntityType type = MBI()->type_from_handle( edges[0] );

  // if vector has both edges and verts, only use edges
  // NOTE: The curve sets from ReadCGM do not contain duplicate endpoints for loops!
  moab::EntityType end_type = MBI()->type_from_handle( edges.back() );
  if(type != end_type) {
    for(std::vector<moab::EntityHandle>::iterator i=edges.begin(); i!=edges.end(); i++) {
      if(moab::MBVERTEX == MBI()->type_from_handle( *i )) {
        i = edges.erase(i) - 1;
      }
    }
  }

  // determine if vector defines an arc by edges of verts
  type = MBI()->type_from_handle( edges[0] );
  if        (moab::MBEDGE == type) {
    if(edges.empty()) return 0.0;
    for( i=edges.begin(); i!=edges.end(); i++ ) {
      int n_verts;
      const moab::EntityHandle *conn;
      result = MBI()->get_connectivity( *i, conn, n_verts );
      if( moab::MB_SUCCESS!=result ) std::cout << "result=" << result << std::endl;
      assert(moab::MB_SUCCESS == result);
      assert( 2 == n_verts );
      if(conn[0] == conn[1]) continue;
      dist += dist_between_verts( conn[0], conn[1] );
    }
  } else if (moab::MBVERTEX == type) {
    if(2 > edges.size()) return 0.0;
    moab::EntityHandle front_vert = edges.front();
    for( i=edges.begin()+1; i!=edges.end(); i++) {
      dist += dist_between_verts( front_vert, *i );
      front_vert = *i;
    }
  } else return moab::MB_FAILURE;

  return dist;
}

// Given a vertex and vector of edges, return the number of edges adjacent to the vertex.
unsigned int Gen::n_adj_edges( moab::EntityHandle vert, moab::Range edges )
{
  moab::ErrorCode result;
  moab::Range adj_edges;
  result = MBI()->get_adjacencies( &vert, 1, 1, false, adj_edges );
  assert(moab::MB_SUCCESS == result);
  adj_edges = moab::intersect( adj_edges, edges );
  return adj_edges.size();
}



// Return true if the edges share a vertex. Does not check for coincident edges.
bool Gen::edges_adjacent( moab::EntityHandle edge0, moab::EntityHandle edge1 )
{
  moab::ErrorCode result;
  moab::Range verts0, verts1;
  result = MBI()->get_adjacencies( &edge0, 1, 0, false, verts0 );
  assert( moab::MB_SUCCESS == result );
  assert( 2 == verts0.size() );
  result = MBI()->get_adjacencies( &edge1, 1, 0, false, verts1 );
  assert( moab::MB_SUCCESS == result );
  assert( 2 == verts1.size() );
  if      ( verts0.front() == verts1.front() ) return true;
  else if ( verts0.front() == verts1.back()  ) return true;
  else if ( verts0.back()  == verts1.back()  ) return true;
  else if ( verts0.back()  == verts1.front() ) return true;
  else                                         return false;
}

// get the direction unit vector from one vertex to another vertex
moab::ErrorCode Gen::get_direction( const moab::EntityHandle from_vert, const moab::EntityHandle to_vert,
                                    moab::CartVect &dir )
{
  // double d[3];
  moab::ErrorCode result;
  moab::CartVect coords0, coords1;
  result = MBI()->get_coords( &from_vert, 1, coords0.array() );
  assert(moab::MB_SUCCESS==result);
  result = MBI()->get_coords( &to_vert, 1, coords1.array() );
  assert(moab::MB_SUCCESS==result);
  dir = coords1 - coords0;
  if(0 == dir.length()) {
    moab::CartVect zero_vector( 0.0 );
    dir = zero_vector;
    std::cout << "direction vector has 0 magnitude" << std::endl;
    return moab::MB_SUCCESS;
  }
  dir.normalize();
  return moab::MB_SUCCESS;
}

// from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
double Gen::edge_point_dist( const moab::CartVect a, const moab::CartVect b, const moab::CartVect c )
{
  moab::CartVect ab, bc, ba, ac;
  ab = b - a;
  bc = c - b;
  ba = a - b;
  ac = c - a;

  // find the magnitude of the cross product and test the line
  moab::CartVect cross_product = ab*ac;
  double dist = cross_product.length() / dist_between_verts(a,b);

  // test endpoint1
  if (ab%bc > 0) {
    return dist_between_verts(b,c);
  }

  // test endpoint0
  if (ba%ac > 0) {
    return dist_between_verts(a,c);
  }
  return fabs(dist);
}
double Gen::edge_point_dist( const moab::EntityHandle endpt0, const moab::EntityHandle endpt1,
                             const moab::EntityHandle pt )
{
  moab::ErrorCode result;
  moab::CartVect a, b, c;
  result = MBI()->get_coords( &endpt0, 1, a.array() );
  assert(moab::MB_SUCCESS==result);
  result = MBI()->get_coords( &endpt1, 1, b.array() );
  assert(moab::MB_SUCCESS==result);
  result = MBI()->get_coords( &pt,     1, c.array() );
  assert(moab::MB_SUCCESS==result);
  return edge_point_dist( a, b, c);
}
double Gen::edge_point_dist( const moab::EntityHandle edge, const moab::EntityHandle pt )
{
  moab::ErrorCode result;
  const moab::EntityHandle *conn;
  int n_verts;
  result = MBI()->get_connectivity( edge, conn, n_verts );
  assert(moab::MB_SUCCESS==result);
  assert( 2 == n_verts );
  return edge_point_dist( conn[0], conn[1], pt );
}
double Gen::triangle_area( const moab::CartVect a, const moab::CartVect b,
                           const moab::CartVect c)
{
  moab::CartVect d = c - a;
  moab::CartVect e = c - b;
  moab::CartVect f = d*e;
  return 0.5*f.length();
}
moab::ErrorCode Gen::triangle_area( const moab::EntityHandle conn[], double &area )
{
  moab::CartVect coords[3];
  moab::ErrorCode result = MBI()->get_coords( conn, 3, coords[0].array() );
  assert(moab::MB_SUCCESS == result);
  area = triangle_area( coords[0], coords[1], coords[2] );
  return moab::MB_SUCCESS;
}
moab::ErrorCode Gen::triangle_area( const moab::EntityHandle tri, double &area )
{
  moab::ErrorCode result;
  const moab::EntityHandle *conn;
  int n_verts;
  result = MBI()->get_connectivity( tri, conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert(3 == n_verts);

  result = triangle_area( conn, area );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}
double Gen::triangle_area( const moab::Range tris )
{
  double a, area = 0;
  moab::ErrorCode result;
  for(moab::Range::iterator i=tris.begin(); i!=tris.end(); i++) {
    result = triangle_area( *i, a);
    assert(moab::MB_SUCCESS == result);
    area += a;
  }
  return area;
}

bool Gen::triangle_degenerate( const moab::EntityHandle tri )
{
  moab::ErrorCode result;
  const moab::EntityHandle *conn;
  int n_verts;
  result = MBI()->get_connectivity( tri, conn, n_verts );
  assert(moab::MB_SUCCESS == result);
  assert(3 == n_verts);
  return triangle_degenerate( conn[0], conn[1], conn[2] );
}

bool Gen::triangle_degenerate( const moab::EntityHandle v0, const moab::EntityHandle v1,
                               const moab::EntityHandle v2 )
{
  if(v0==v1 || v1==v2 || v2==v0) return true;
  return false;
}

moab::ErrorCode Gen::triangle_normals( const moab::Range tris, std::vector<moab::CartVect> &normals )
{
  moab::ErrorCode result;
  normals.clear();
  for(moab::Range::const_iterator i=tris.begin(); i!=tris.end(); i++) {
    moab::CartVect normal;
    result = triangle_normal( *i, normal );
    assert(moab::MB_SUCCESS==result || moab::MB_ENTITY_NOT_FOUND==result);
    normals.push_back( normal );
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::triangle_normal( const moab::EntityHandle tri, moab::CartVect &normal)
{
  moab::ErrorCode result;
  const moab::EntityHandle *conn;
  int n_verts;
  result = MBI()->get_connectivity( tri, conn, n_verts );
  if(moab::MB_ENTITY_NOT_FOUND == result) {
    std::cout << "triangle_normal: triangle not found" << std::endl;
    moab::CartVect zero_vector( 0.0 );
    normal = zero_vector;
    return result;
  } else if(moab::MB_SUCCESS != result) {
    return result;
  } else {
    assert(3 == n_verts);
    return triangle_normal( conn[0], conn[1], conn[2], normal );
  }
}

moab::ErrorCode Gen::triangle_normal( const moab::EntityHandle v0, const moab::EntityHandle v1,
                                      const moab::EntityHandle v2, moab::CartVect &normal )
{

  // if tri is degenerate return 0,0,0
  if( triangle_degenerate(v0, v1, v2) ) {
    moab::CartVect zero_vector( 0.0 );
    normal = zero_vector;
    std::cout << "  normal=" << normal << std::endl;
    return moab::MB_SUCCESS;
  }

  moab::EntityHandle conn[3];
  conn[0] = v0;
  conn[1] = v1;
  conn[2] = v2;
  moab::ErrorCode result;
  moab::CartVect coords[3];
  result = MBI()->get_coords( conn, 3, coords[0].array() );
  assert(moab::MB_SUCCESS == result);
  return triangle_normal( coords[0], coords[1], coords[2], normal );
}

moab::ErrorCode Gen::triangle_normal( const moab::CartVect coords0, const moab::CartVect coords1,
                                      const moab::CartVect coords2, moab::CartVect &normal )
{
  moab::CartVect edge0, edge1;
  edge0 = coords1-coords0;
  edge1 = coords2-coords0;
  normal = edge0*edge1;

  // do not normalize if magnitude is zero (avoid nans)
  if(0 == normal.length()) return moab::MB_SUCCESS;

  normal.normalize();
  return moab::MB_SUCCESS;
}



// Distance between a point and line. The line is defined by two verts.
// We are using a line and not a line segment!
// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
moab::ErrorCode Gen::line_point_dist( const moab::EntityHandle line_pt1, const moab::EntityHandle line_pt2,
                                      const moab::EntityHandle pt0, double &dist )
{
  moab::ErrorCode result;
  moab::CartVect x0, x1, x2;
  result = MBI()->get_coords( &line_pt1, 1, x1.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &line_pt2, 1, x2.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt0, 1, x0.array() );
  assert(moab::MB_SUCCESS == result);

  dist = ( ((x0-x1)*(x0-x2)).length() ) / ( (x2-x1).length() );
  return moab::MB_SUCCESS;
}

// Project the point onto the line. Not the line segment!
moab::ErrorCode Gen::point_line_projection( const moab::EntityHandle line_pt1,
    const moab::EntityHandle line_pt2,
    const moab::EntityHandle pt0 )
{
  moab::CartVect projected_coords;
  double parameter;
  moab::ErrorCode result = point_line_projection( line_pt1, line_pt2,
                           pt0, projected_coords,
                           parameter );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->set_coords( &pt0, 1, projected_coords.array() );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::point_line_projection( const moab::EntityHandle line_pt1,
    const moab::EntityHandle line_pt2,
    const moab::EntityHandle pt0,
    moab::CartVect &projected_coords,
    double &parameter  )
{

  moab::ErrorCode result;
  moab::CartVect coords[3];
  result = MBI()->get_coords( &line_pt1, 1, coords[1].array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &line_pt2, 1, coords[2].array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt0, 1, coords[0].array() );
  assert(moab::MB_SUCCESS == result);

  // project the t_joint between the endpts
  // http://en.wikipedia.org/wiki/Vector_projection
  moab::CartVect a = coords[0] - coords[1];
  moab::CartVect b = coords[2] - coords[1];
  parameter    = (a%b)/(b%b);
  moab::CartVect c = parameter*b;
  projected_coords = c     + coords[1];
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::point_line_projection( const moab::EntityHandle line_pt1,
    const moab::EntityHandle line_pt2,
    const moab::EntityHandle pt0,
    double &dist_along_edge  )
{

  moab::ErrorCode result;
  moab::CartVect coords[3];
  result = MBI()->get_coords( &line_pt1, 1, coords[1].array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &line_pt2, 1, coords[2].array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt0, 1, coords[0].array() );
  assert(moab::MB_SUCCESS == result);

  // project the t_joint between the endpts
  // http://en.wikipedia.org/wiki/Vector_projection
  moab::CartVect a = coords[0] - coords[1];
  moab::CartVect b = coords[2] - coords[1];
  dist_along_edge = a%b / b.length();
  return moab::MB_SUCCESS;
}

double Gen::area2( const moab::EntityHandle pt_a, const moab::EntityHandle pt_b,
                   const moab::EntityHandle pt_c, const moab::CartVect plane_normal )
{
  moab::ErrorCode result;
  moab::CartVect a, b, c;
  result = MBI()->get_coords( &pt_a, 1, a.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt_b, 1, b.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt_c, 1, c.array() );
  assert(moab::MB_SUCCESS == result);
  moab::CartVect d = b - a;
  moab::CartVect e = c - a;
  // project onto a plane defined by the plane's normal vector
  return (d*e)%plane_normal;
}

// Is point c to the left of line ab?
bool Gen::left( const moab::EntityHandle a, const moab::EntityHandle b,
                const moab::EntityHandle c, const moab::CartVect n )
{
  double area_2 = area2(a,b,c,n);
  if(area_2 > 0) return true;
  else return false;
}

// Is point c to the left of line ab or collinear?
bool Gen::left_on( const moab::EntityHandle a, const moab::EntityHandle b,
                   const moab::EntityHandle c, const moab::CartVect n )
{
  double area_2 = area2(a,b,c,n);
  if(area_2 >= 0) return true;
  else return false;
}

// Are pts a,b,c collinear?
bool Gen::collinear( const moab::EntityHandle a, const moab::EntityHandle b,
                     const moab::EntityHandle c, const moab::CartVect n )
{
  double area_2 = area2(a,b,c,n);
  if( area_2 ==0) return true;
  else return false;
}

// Exclusive or: T iff exactly one argument is true
bool logical_xor( const bool x, const bool y )
{
  return (x || y) && !(x && y);
}

bool Gen::intersect_prop( const moab::EntityHandle a, const moab::EntityHandle b,
                          const moab::EntityHandle c, const moab::EntityHandle d,
                          const moab::CartVect n )
{
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

bool Gen::between( const moab::EntityHandle pt_a, const moab::EntityHandle pt_b,
                   const moab::EntityHandle pt_c, const moab::CartVect n)
{
  if( !collinear(pt_a,pt_b,pt_c,n) ) return false;

  moab::ErrorCode result;
  moab::CartVect a, b, c;
  result = MBI()->get_coords( &pt_a, 1, a.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt_b, 1, b.array() );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->get_coords( &pt_c, 1, c.array() );
  assert(moab::MB_SUCCESS == result);

  // if ab not vertical, check betweenness on x; else on y.
  if(a[0] != b[0]) {
    return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
  } else if(a[1] != b[1]) {
    return ((a[1] <= c[1]) && (c[1] <= b[1])) || ((a[1] >= c[1]) && (c[1] >= b[1]));
  } else {
    return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
  }
}

bool Gen::intersect( const moab::EntityHandle a, const moab::EntityHandle b,
                     const moab::EntityHandle c, const moab::EntityHandle d,
                     const moab::CartVect n )
{
  if(intersect_prop(a,b,c,d,n)) return true;
  else if( between(a,b,c,n) ||
           between(a,b,d,n) ||
           between(c,d,a,n) ||
           between(c,d,b,n) ) return true;
  else return false;
}

// verts is an ordered polygon of verts
bool Gen::diagonalie( const moab::EntityHandle a, const moab::EntityHandle b,
                      const moab::CartVect n,
                      const std::vector<moab::EntityHandle> verts )
{
  for(unsigned int i=0; i<verts.size(); i++) {
    moab::EntityHandle c = verts[i];
    moab::EntityHandle c1;
    if(verts.size()-1 == i) c1 = verts[0];
    else                    c1 = verts[i+1];

    if( (c != a) && (c1 != a) &&
        (c != b) && (c1 != b) &&
        intersect( a, b, c, c1, n ) ) {
      return false;
    }
  }
  return true;
}

// verts is an ordered polygon of verts
bool Gen::in_cone( const moab::EntityHandle a, const moab::EntityHandle b,
                   const moab::CartVect n,
                   const std::vector<moab::EntityHandle> verts )
{
  std::vector<moab::EntityHandle>::const_iterator a_iter;
  a_iter = find( verts.begin(), verts.end(), a );
  moab::EntityHandle a0, a1;
  // a0 is before a
  if(verts.begin() == a_iter) a0 = verts[verts.size()-1];
  else a0 = *(a_iter-1);
  // a1 is after a
  if(verts.end()-1 == a_iter) a1 = verts[0];
  else a1 = *(a_iter+1);

  // if a is a convex vertex
  if(left_on(a,a1,a0,n)) return left(a,b,a0,n) && left(b,a,a1,n);

  // else a is reflex
  else return !(left_on(a,b,a1,n) && left_on(b,a,a0,n));
}

bool Gen::diagonal( const moab::EntityHandle a, const moab::EntityHandle b,
                    const moab::CartVect n,
                    const std::vector<moab::EntityHandle> verts )
{
  bool result = in_cone(a,b,n,verts) && in_cone(b,a,n,verts) && diagonalie(a,b,n,verts);
  return result;
}


// Determine if each vertex is an ear. Input an ordered polygon of verts.
moab::ErrorCode Gen::ear_init( const std::vector<moab::EntityHandle> verts,
                               const moab::CartVect n, // plane normal vector
                               std::vector<bool> &is_ear )
{
  if(verts.size() != is_ear.size()) return moab::MB_FAILURE;
  for(unsigned int i=0; i<verts.size(); i++) {
    moab::EntityHandle prev, next;
    if(0 == i) prev = verts.back();
    else prev = verts[i-1];
    if(verts.size()-1 == i) next = verts[0];
    else next = verts[i+1];
    is_ear[i] = diagonal(prev,next,n,verts);
  }
  return moab::MB_SUCCESS;
}

// Input an ordered polygon of verts and a normal vector of the plane
// that the polygon is mostly in. The vector is required for orientation.
moab::ErrorCode Gen::ear_clip_polygon( std::vector<moab::EntityHandle> verts,
                                       moab::CartVect n,
                                       moab::Range &new_tris )
{

  moab::ErrorCode result;
  std::vector<bool> is_ear( verts.size() );
  result = ear_init( verts, n, is_ear );
  assert(moab::MB_SUCCESS == result);

  // if the polygon intersects itself the algorithm will not stop
  int counter = 0;
  int n_initial_verts = verts.size();

  while(3 < verts.size()) {
    for(unsigned int i=0; i<verts.size(); i++) {
      if(is_ear[i]) {
        moab::EntityHandle v0, v1, v2, v3, v4;
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

        moab::EntityHandle new_tri;
        moab::EntityHandle conn[3] = {v1,v2,v3};
        result = MBI()->create_element( moab::MBTRI, conn, 3, new_tri );
        assert(moab::MB_SUCCESS == result);
        new_tris.insert( new_tri );

        // update ear status
        if(0 == i) is_ear[verts.size()-1] = diagonal( v0, v3, n, verts );
        else       is_ear[i-1]            = diagonal( v0, v3, n, verts );
        if(verts.size()-1 == i) is_ear[0] = diagonal( v1, v4, n, verts );
        else                    is_ear[i+1] = diagonal( v1, v4, n, verts );

        // cut off the ear at i
        verts.erase( verts.begin()+i );
        is_ear.erase( is_ear.begin()+i );
        break;
      }
    }

    // If the polygon has intersecting edges this loop will continue until it
    // hits this return.
    if(counter > n_initial_verts) {
      result = MBI()->delete_entities( new_tris );
      assert(moab::MB_SUCCESS == result);
      new_tris.clear();
      return moab::MB_FAILURE;
    }
    counter++;
  }
  moab::EntityHandle new_tri;
  moab::EntityHandle conn[3] = {verts[0],verts[1],verts[2]};
  result = MBI()->create_element( moab::MBTRI, conn, 3, new_tri );
  assert(moab::MB_SUCCESS == result);
  new_tris.insert( new_tri );

  return moab::MB_SUCCESS;
}

int Gen::geom_id_by_handle( const moab::EntityHandle set )
{
  moab::ErrorCode result;
  moab::Tag id_tag;
  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,id_tag,moab::MB_TAG_DENSE);
  assert(moab::MB_SUCCESS==result || moab::MB_ALREADY_ALLOCATED==result);
  int id;
  result = MBI()->tag_get_data( id_tag, &set, 1, &id );
  return id;
}

moab::ErrorCode Gen::save_normals( moab::Range tris, moab::Tag normal_tag )
{
  std::vector<moab::CartVect> normals(tris.size());
  moab::ErrorCode result;
  result = triangle_normals( tris, normals );
  assert(moab::MB_SUCCESS == result);

  result = MBI()->tag_set_data(normal_tag, tris, &normals[0]);
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::flip(const moab::EntityHandle tri, const moab::EntityHandle vert0,
                          const moab::EntityHandle vert2, const moab::EntityHandle surf_set)
{

  // get the triangles in the surface. The tri and adj_tri must be in the surface.
  moab::Range surf_tris;
  moab::EntityHandle result = MBI()->get_entities_by_type( surf_set, moab::MBTRI, surf_tris);
  assert(moab::MB_SUCCESS == result);

  // get the triangle across the edge that will be flipped
  moab::Range adj_tri;
  moab::EntityHandle edge[2] = {vert0, vert2};
  result = MBI()->get_adjacencies( edge, 2, 2, false, adj_tri );
  assert(moab::MB_SUCCESS == result);
  adj_tri = moab::intersect(adj_tri, surf_tris);
  assert(2 == adj_tri.size());
  adj_tri.erase( tri );
  print_triangle( adj_tri.front(), false );

  // get the remaining tri vert
  moab::Range tri_verts;
  result = MBI()->get_adjacencies( &tri, 1, 0, false, tri_verts );
  assert(moab::MB_SUCCESS == result);
  assert(3 == tri_verts.size());
  tri_verts.erase(vert0);
  tri_verts.erase(vert2);
  assert(1 == tri_verts.size());
  moab::EntityHandle vert1 = tri_verts.front();

  // get the remaining adj_tri vert
  moab::Range adj_tri_verts;
  result = MBI()->get_adjacencies( &adj_tri.front(), 1, 0, false, adj_tri_verts );
  assert(moab::MB_SUCCESS == result);
  assert(3 == adj_tri_verts.size());
  adj_tri_verts.erase(vert0);
  adj_tri_verts.erase(vert2);
  assert(1 == adj_tri_verts.size());
  moab::EntityHandle vert3 = adj_tri_verts.front();

  // set the new connectivity
  moab::EntityHandle tri_conn[3] = {vert0, vert1, vert3};
  result = MBI()->set_connectivity( tri, tri_conn, 3 );
  assert(moab::MB_SUCCESS == result);
  moab::EntityHandle adj_tri_conn[3] = {vert1, vert2, vert3};
  result = MBI()->set_connectivity( adj_tri.front(), adj_tri_conn, 3 );
  assert(moab::MB_SUCCESS == result);
  print_triangle( tri, false );
  print_triangle( adj_tri.front(), false );
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::ordered_verts_from_ordered_edges( const std::vector<moab::EntityHandle> ordered_edges,
    std::vector<moab::EntityHandle> &ordered_verts )
{
  moab::ErrorCode result;
  ordered_verts.clear();
  ordered_verts.reserve(ordered_edges.size()+1);

  // Save the back of the previous edge to check for continuity.
  moab::EntityHandle previous_back_vert;

  for(std::vector<moab::EntityHandle>::const_iterator i=ordered_edges.begin();
      i!=ordered_edges.end(); i++) {
    const moab::EntityHandle *conn;
    int n_verts;
    result = MBI()->get_connectivity( *i, conn, n_verts);
    assert(moab::MB_SUCCESS == result);
    assert(2 == n_verts);
    if(ordered_edges.begin() == i) {
      ordered_verts.push_back(conn[0]);
    } else {
      assert(previous_back_vert == conn[0]);
    }
    ordered_verts.push_back(conn[1]);
    previous_back_vert = conn[1];
  }
  return moab::MB_SUCCESS;
}

/* Find the distance between two arcs. Assume that their endpoints are somewhat
   close together. */
moab::ErrorCode Gen::dist_between_arcs( bool debug,
                                        const std::vector<moab::EntityHandle> arc0,
                                        const std::vector<moab::EntityHandle> arc1,
                                        double &dist )
{
  dist = 0;

  // Special Case: arcs have no verts.
  if( arc0.empty() || arc1.empty() ) {
    std::cout << "arc has no vertices" << std::endl;
    return moab::MB_FAILURE;
  }

  // for simplicity, put arcs into the same structure
  std::vector<moab::EntityHandle> arcs[2] = {arc0, arc1};

  // Special Case: Remove duplicate vert handles
  for(unsigned int i=0; i<2; ++i) {
    if( 2>arcs[i].size() ) continue;
    for(std::vector<moab::EntityHandle>::iterator j=arcs[i].begin()+1; j!=arcs[i].end(); ++j) {
      if(*j == *(j-1)) {
        if(debug) {
          print_loop( arcs[i] );
          std::cout << "dist_between_arcs: found duplicate vert handle in arc" << std::endl;
        }
        j = arcs[i].erase(j) - 1;
      }
    }
  }

  // get the coords in one call per arc. For speed, do not ask MOAB again for coords.
  moab::ErrorCode result;
  std::vector<moab::CartVect> coords[2];
  for(unsigned int i=0; i<2; i++) {
    coords[i].resize( arcs[i].size() );
    result = MBI()->get_coords( &arcs[i][0], arcs[i].size(), coords[i][0].array());
    assert(moab::MB_SUCCESS == result);
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
      return moab::MB_SUCCESS;
      // The other arc has more than one point
    } else {
      double area = 0.0;
      for(unsigned int i=0; i<arcs[other_arc_index].size()-1; ++i) {
        area += fabs(triangle_area( coords[other_arc_index][i],
                                    coords[other_arc_index][i+1],
                                    coords[point_arc_index].front() ));
      }
      dist = area / length(arcs[other_arc_index]);
      return moab::MB_SUCCESS;
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

  }

  // Insert new points to match the other arcs points, by parameter.
  for(unsigned int i=0; i<2; i++) {
    for(unsigned int j=0; j<mgd_params.size(); j++) {
      //std::cout << "params[" << i << "][" << j << "]=" << params[i][j]
      //          << " mgd_params[" << j << "]=" << mgd_params[j] << std::endl;
      if(params[i][j] > mgd_params[j]) {
        double ratio = (mgd_params[j]-params[i][j-1]) / (params[i][j]-params[i][j-1]);
        //std::cout << "j=" << j << " ratio=" << ratio << std::endl;
        moab::CartVect pt = coords[i][j-1] + ratio*(coords[i][j]-coords[i][j-1]);
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
    area += fabs(triangle_area( coords[0][i], coords[1][i+1], coords[1][i] ));
  }

  // Divide the area by the average length to get the average distance between arcs.
  dist = fabs(2*area / (arc_len[0] + arc_len[1] ));
  return moab::MB_SUCCESS;
}

// qsort struct comparision function
// If the first handle is the same, compare the second
int compare_edge(const void *a, const void *b)
{
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
moab::ErrorCode Gen::find_skin( moab::Range tris, const int dim,
                                moab::Range &skin_edges,
                                const bool temp_bool )
{

  const bool local_debug = false;

  if(1 != dim) return moab::MB_NOT_IMPLEMENTED;
  if(moab::MBTRI != MBI()->type_from_handle(tris.front())) return moab::MB_NOT_IMPLEMENTED;

  skin_edges.clear();
  if(tris.empty()) return moab::MB_ENTITY_NOT_FOUND;

  // This implementation gets some of its speed due to not checking for edges
  moab::ErrorCode result;
  int n_edges;
  result = MBI()->get_number_entities_by_type( 0, moab::MBEDGE, n_edges );
  assert(moab::MB_SUCCESS == result);
  if(0 != n_edges) {
    moab::Range temp_edges;
    result = MBI()->get_entities_by_type( 0, moab::MBEDGE, temp_edges);
    assert(moab::MB_SUCCESS == result);
    result = MBI()->list_entities( temp_edges );
    assert(moab::MB_SUCCESS == result);
  }
  assert(0 == n_edges);

  // Get connectivity. Do not create edges.
  edge *edges = new edge[3*tris.size()];
  int n_verts;
  int ii = 0;
  for(moab::Range::iterator i=tris.begin(); i!=tris.end(); i++) {
    const moab::EntityHandle *conn;
    result = MBI()->get_connectivity( *i, conn, n_verts );
    assert(moab::MB_SUCCESS == result);
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
      moab::EntityHandle temp = edges[i].v0;
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
      const moab::EntityHandle conn[2] = {edges[i].v0, edges[i].v1};
      moab::EntityHandle edge;
      result = MBI()->create_element( moab::MBEDGE, conn, 2, edge );
      assert(moab::MB_SUCCESS == result);
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
      // otherwise a skin edge has been found
    } else {
      const moab::EntityHandle conn[2] = {edges[i].v0, edges[i].v1};
      moab::EntityHandle edge;
      result = MBI()->create_element( moab::MBEDGE, conn, 2, edge );
      MB_CHK_SET_ERR(result, "could not create edge element");
      skin_edges.insert( edge );
    }
  }
  delete[] edges;
  return moab::MB_SUCCESS;
}

// calculate volume of polyhedron
// Copied from DagMC, without index_by_handle. The dagmc function will
// segfault if build_indices is not first called. For sealing there is
// no need to build_indices.
moab::ErrorCode Gen::measure_volume( const moab::EntityHandle volume, double& result, bool debug, bool verbose )
{
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> surfaces, surf_volumes;
  result = 0.0;



  // get surfaces from volume
  rval = MBI()->get_child_meshsets( volume, surfaces );
  if (moab::MB_SUCCESS != rval) {
    return rval;
  }
  if(debug) std::cout << "in measure_volume 1" << std::endl;

  // get surface senses
  std::vector<int> senses( surfaces.size() );


  if (rval != moab::MB_SUCCESS) {
    std::cout << "Could not measure volume" << std::endl;
    std::cout << "This can happen for 2 reasons, there are no volumes" << std::endl;
    std::cout << "or the pointer to the Moab instance is poor" << std::endl;
    exit(rval);
  }
  if(debug) std::cout << surfaces.size() << " " << result << std::endl;

  rval = surface_sense( volume, surfaces.size(), &surfaces[0], &senses[0] );


  if(debug) std::cout << "in measure_volume 2" << std::endl;

  if (moab::MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }


  for (unsigned i = 0; i < surfaces.size(); ++i) {
    // skip non-manifold surfaces
    if (!senses[i])
      continue;

    // get triangles in surface
    moab::Range triangles;
    rval = MBI()->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (moab::MB_SUCCESS != rval) {
      return rval;
    }
    if (!triangles.all_of_type(moab::MBTRI)) {
      std::cout << "WARNING: Surface " << geom_id_by_handle(surfaces[i])
                << " contains non-triangle elements. Volume calculation may be incorrect."
                << std::endl;
      triangles.clear();
      rval = MBI()->get_entities_by_type( surfaces[i], moab::MBTRI, triangles );
      if (moab::MB_SUCCESS != rval) {
        return rval;
      }
    }

    // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const moab::EntityHandle *conn;
    int len;
    moab::CartVect coords[3];
    for (moab::Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = MBI()->get_connectivity( *j, conn, len, true );
      if (moab::MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = MBI()->get_coords( conn, 3, coords[0].array() );
      if (moab::MB_SUCCESS != rval) return rval;

      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }

  result /= 6.0;
  if(debug)  std::cout << result << std::endl;
  return moab::MB_SUCCESS;
}

/// Calculate the signed volumes beneath the surface (x 6.0). Use the triangle's
///   cannonical sense. Do not take sense tags into account. Code taken from
///   DagMC::measure_volume.
moab::ErrorCode Gen::get_signed_volume( const moab::EntityHandle surf_set, double &signed_volume)
{
  moab::ErrorCode rval;
  moab::Range tris;
  rval = MBI()->get_entities_by_type( surf_set, moab::MBTRI, tris );
  if(moab::MB_SUCCESS != rval) return rval;
  signed_volume = 0.0;
  const moab::EntityHandle *conn;
  int len;
  moab::CartVect coords[3];
  for (moab::Range::iterator j = tris.begin(); j != tris.end(); ++j) {
    rval = MBI()->get_connectivity( *j, conn, len, true );
    if (moab::MB_SUCCESS != rval) return rval;
    assert(3 == len);
    rval = MBI()->get_coords( conn, 3, coords[0].array() );
    if (moab::MB_SUCCESS != rval) return rval;

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    signed_volume += (coords[0] % (coords[1] * coords[2]));
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::measure( const moab::EntityHandle set, const moab::Tag geom_tag, double &size, bool debug,  bool verbose )
{
  moab::ErrorCode result;
  int dim;
  result = MBI()->tag_get_data( geom_tag, &set, 1, &dim );
  assert(moab::MB_SUCCESS == result);
  if(0 == dim) {
    std::cout << "measure: cannot measure vertex" << std::endl;
    return moab::MB_FAILURE;

  } else if(1 == dim) {
    std::vector<moab::EntityHandle> vctr;
    result = get_meshset( set, vctr );
    assert(moab::MB_SUCCESS == result);
    size = length( vctr );

  } else if(2 == dim) {
    moab::Range tris;
    result = MBI()->get_entities_by_type( set, moab::MBTRI, tris );
    assert(moab::MB_SUCCESS == result);
    size = triangle_area( tris );

  } else if(3 == dim) {

    result = measure_volume( set, size, debug, verbose );

    if(moab::MB_SUCCESS != result) {
      std::cout << "result=" << result << " vol_id="
                << geom_id_by_handle(set) << std::endl;
    }
  } else {
    std::cout << "measure: incorrect dimension" << std::endl;
    return moab::MB_FAILURE;
  }
  return moab::MB_SUCCESS;
}

// From CGMA/builds/dbg/include/CubitDefines
/// gets the surface sense with respect to the curve and returns the value to sense
moab::ErrorCode Gen::get_curve_surf_sense( const moab::EntityHandle surf_set, const moab::EntityHandle curve_set,
    int &sense, bool debug )
{
  std::vector<moab::EntityHandle> surfs;
  std::vector<int> senses;
  moab::ErrorCode rval;
  moab::GeomTopoTool gt( MBI(), false);
  rval = gt.get_senses( curve_set, surfs, senses );
  MB_CHK_SET_ERR(rval,"failed to get_senses");

  unsigned counter = 0;
  int edim;
  for(unsigned i=0; i<surfs.size(); i++) {
    edim=gt.dimension(surfs[i]);
    if( edim == -1) {
      surfs.erase(surfs.begin()+i);
      senses.erase(senses.begin()+i);
      i=0;
    }
    if(surf_set==surfs[i]) {
      sense = senses[i];

      ++counter;
    }
  }

  if(debug) {
    for(unsigned int index=0; index < surfs.size() ; index++) {
      std::cout << "surf = " << geom_id_by_handle(surfs[index]) << std::endl;
      std::cout << "sense = " << senses[index] << std::endl;
    }
  }
  // special case: parent surface does not exist
  if(0==counter) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"failed to find a surf in sense list");
  }

  // special case: ambiguous
  if(1<counter) sense = moab::SENSE_BOTH;
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::surface_sense( moab::EntityHandle volume,
                                    int num_surfaces,
                                    const moab::EntityHandle* surfaces,
                                    int* senses_out )
{
  std::vector<moab::EntityHandle> surf_volumes( 2*num_surfaces );
  moab::Tag senseTag = get_tag( "GEOM_SENSE_2", 2, moab::MB_TAG_SPARSE, moab::MB_TYPE_HANDLE, NULL, false );
  moab::ErrorCode rval = MBI()->tag_get_data( senseTag , surfaces, num_surfaces, &surf_volumes[0] );
  if (moab::MB_SUCCESS != rval) {
    return rval;
  }

  const moab::EntityHandle* end = surfaces + num_surfaces;
  std::vector<moab::EntityHandle>::const_iterator surf_vols = surf_volumes.begin();
  while (surfaces != end) {
    moab::EntityHandle forward = *surf_vols;
    ++surf_vols;
    moab::EntityHandle reverse = *surf_vols;
    ++surf_vols;
    if (volume == forward) {
      *senses_out = (volume != reverse); // zero if both, otherwise 1
    } else if (volume == reverse) {
      *senses_out = moab::SENSE_BOTH;
    } else {
      return moab::MB_ENTITY_NOT_FOUND;
    }

    ++surfaces;
    ++senses_out;
  }

  return moab::MB_SUCCESS;
}

/// get sense of surface(s) wrt volume
moab::ErrorCode Gen::surface_sense( moab::EntityHandle volume,
                                    moab::EntityHandle surface,
                                    int& sense_out )
{
  // get sense of surfaces wrt volumes
  moab::EntityHandle surf_volumes[2];
  moab::Tag senseTag = get_tag( "GEOM_SENSE_2", 2, moab::MB_TAG_SPARSE, moab::MB_TYPE_HANDLE, NULL, false );
  moab::ErrorCode rval = MBI()->tag_get_data( senseTag , &surface, 1, surf_volumes );
  if (moab::MB_SUCCESS != rval) {
    return rval;
  }

  if (surf_volumes[0] == volume) {
    sense_out = (surf_volumes[1] != volume); // zero if both, otherwise 1
  } else if (surf_volumes[1] == volume) {
    sense_out = moab::SENSE_BOTH;
  } else {
    return moab::MB_ENTITY_NOT_FOUND;
  }

  return moab::MB_SUCCESS;
}

moab::Tag Gen::get_tag( const char* name, int size, moab::TagType store,
                        moab::DataType type, const void* def_value,
                        bool create_if_missing)
{
  moab::Tag retval = 0;
  unsigned flags = store|moab::MB_TAG_CREAT;
  if (!create_if_missing) {
    flags |= moab::MB_TAG_EXCL;
  }
  moab::ErrorCode result = MBI()->tag_get_handle(name, size, type, retval, flags, def_value);
  if (create_if_missing && moab::MB_SUCCESS != result) {
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  }
  return retval;
}

moab::ErrorCode Gen::delete_surface( moab::EntityHandle surf, moab::Tag geom_tag, moab::Range tris, int id, bool debug, bool verbose )
{

  moab::ErrorCode result;

  //measure area of the surface
  double area;
  result = measure( surf, geom_tag, area, debug, verbose );
  MB_CHK_SET_ERR(result,"could not measure area");
  assert(moab::MB_SUCCESS == result);

  //remove triagngles from the surface
  result = MBI()->remove_entities( surf, tris);
  MB_CHK_SET_ERR(result,"could not remove tris");
  assert(moab::MB_SUCCESS == result);
  //print information about the deleted surface if requested by the user
  if(debug) std::cout << "  deleted surface " << id
                        << ", area=" << area << " cm^2, n_facets=" << tris.size() << std::endl;

  // delete triangles from mesh data
  result = MBI()->delete_entities( tris );
  MB_CHK_SET_ERR(result,"could not delete tris");
  assert(moab::MB_SUCCESS == result);

  //remove the sense data for the surface from its child curves
  result = remove_surf_sense_data( surf, debug);
  MB_CHK_SET_ERR(result, "could not remove surface's sense data");

  //if this was the last surface in a volume, then delete the volume
  std::vector<moab::EntityHandle> parent_volumes;
  result = MBI()->get_parent_meshsets( surf, parent_volumes);
  MB_CHK_SET_ERR(result,"could not get the surface's parent meshsets");
  assert(moab::MB_SUCCESS == result);

  //check each parent volume
  std::vector<moab::EntityHandle>::iterator i;
  for( i=parent_volumes.begin(); i != parent_volumes.end(); i++) {
    moab::EntityHandle parent_vol = *i;
    std::vector<moab::EntityHandle> child_surfs;
    result = MBI()->get_child_meshsets(parent_vol, child_surfs);
    MB_CHK_SET_ERR(result, "could not get the child surfaces of the volume");
    assert(moab::MB_SUCCESS == result);

    // delete volume if it only has this surface as its child
    if( child_surfs.size() == 1 && child_surfs[0] == surf ) {
      result = delete_vol(parent_vol);
      MB_CHK_SET_ERR(result, "could not delete the parent volume");
    }
  }

  //remove the surface set itself
  result = MBI()->delete_entities( &(surf), 1);
  MB_CHK_SET_ERR(result,"could not delete surface set");
  assert(moab::MB_SUCCESS == result);


  return moab::MB_SUCCESS;
}

/// removes sense data from all curves associated with the surface given to the function
moab::ErrorCode Gen::remove_surf_sense_data(moab::EntityHandle del_surf, bool debug)
{

  moab::ErrorCode result;
  moab::GeomTopoTool gt(MBI(), false);
  int edim = gt.dimension(del_surf);

  if(edim!=2) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"could not remove sense data: entity is of the wrong dimension");
  }

// get the curves of the surface
  moab::Range del_surf_curves;
  result = MBI() -> get_child_meshsets( del_surf, del_surf_curves);
  MB_CHK_SET_ERR(result,"could not get the curves of the surface to delete");
  if (debug) std::cout << "got the curves" << std::endl;

  if (debug) {
    std::cout << "number of curves to the deleted surface = " << del_surf_curves.size() << std::endl;
    for(unsigned int index =0 ; index < del_surf_curves.size() ; index++) {
      std::cout << "deleted surface's child curve id " << index << " = " << geom_id_by_handle(del_surf_curves[index]) << std::endl;
    }
  }
  //get the sense data for each curve

  //get sense_tag handles from MOAB
  moab::Tag senseEnts, senseSenses;
  unsigned flags = moab::MB_TAG_SPARSE;

  //get tag for the entities with sense data associated with a given moab entity
  result = MBI()-> tag_get_handle(GEOM_SENSE_N_ENTS_TAG_NAME, 0, moab::MB_TYPE_HANDLE, senseEnts, flags);
  MB_CHK_SET_ERR(result, "could not get senseEnts tag");

  //get tag for the sense data associated with the senseEnts entities for a given moab entity
  result = MBI()-> tag_get_handle(GEOM_SENSE_N_SENSES_TAG_NAME, 0, moab::MB_TYPE_INTEGER, senseSenses, flags);
  MB_CHK_SET_ERR(result,"could not get senseSenses tag");

  //initialize vectors for entities and sense data
  std::vector<moab::EntityHandle> surfaces;
  std::vector<int> senses;
  for(moab::Range::iterator i=del_surf_curves.begin(); i!=del_surf_curves.end(); i++ ) {

    result = gt.get_senses(*i, surfaces, senses);
    MB_CHK_SET_ERR(result, "could not get the senses for the del_surf_curve");
    // if the surface to be deleted (del_surf) exists in the sense data (which it should), then remove it
    for(unsigned int index = 0; index < senses.size() ; index++) {
      if(surfaces[index]==del_surf) {
        surfaces.erase(surfaces.begin() + index);
        senses.erase(senses.begin() +index);
        index=-1;
      }
    }
    //remove existing sense entity data for the curve
    result= MBI()-> tag_delete_data( senseEnts, &*i, 1);
    MB_CHK_SET_ERR(result, "could not delete sense entity data");

    //remove existing sense data for the curve
    result = MBI()-> tag_delete_data(senseSenses, &*i, 1);
    MB_CHK_SET_ERR(result, "could not delete sense data");

    //reset the sense data for each curve
    result = gt.set_senses( *i, surfaces, senses);
    MB_CHK_SET_ERR(result, "could not update sense data for surface deletion");

  }

  return moab::MB_SUCCESS;
}

/// combines the senses of any curves tagged as merged in the vector curves
moab::ErrorCode Gen::combine_merged_curve_senses( std::vector<moab::EntityHandle> &curves, moab::Tag merge_tag, bool debug)
{

  moab::ErrorCode result;

  for(std::vector<moab::EntityHandle>::iterator j=curves.begin(); j!=curves.end(); j++) {

    moab::EntityHandle merged_curve;
    result = MBI() -> tag_get_data( merge_tag, &(*j), 1, &merged_curve);
    if(moab::MB_SUCCESS!=result && moab::MB_TAG_NOT_FOUND!=result) {
      MB_CHK_SET_ERR(result,"could not get the merge_tag data of the curve");
    }

    if(moab::MB_SUCCESS==result) { // we have found a merged curve pairing
      // add the senses from the curve_to_delete to curve_to keep
      // create vectors for the senses and surfaces of each curve
      std::vector<moab::EntityHandle> curve_to_keep_surfs, curve_to_delete_surfs, combined_surfs;
      std::vector<int> curve_to_keep_senses, curve_to_delete_senses, combined_senses;

      //initialize GeomTopoTool.cpp instance in MOAB
      moab::GeomTopoTool gt(MBI(), false);
      // get senses of the iterator curve and place them in the curve_to_delete vectors
      result = gt.get_senses( *j, curve_to_delete_surfs, curve_to_delete_senses);
      MB_CHK_SET_ERR(result, "could not get the surfs/senses of the curve to delete");
      // get surfaces/senses of the merged_curve and place them in the curve_to_keep vectors
      result = gt.get_senses( merged_curve, curve_to_keep_surfs, curve_to_keep_senses);
      MB_CHK_SET_ERR(result, "could not get the surfs/senses of the curve to delete");

      if(debug) {
        std::cout << "curve to keep id = " << geom_id_by_handle(merged_curve) << std::endl;
        std::cout << "curve to delete id = " << geom_id_by_handle(*j) << std::endl;
        for(unsigned int index=0; index < curve_to_keep_surfs.size(); index++) {
          std::cout << "curve_to_keep_surf " << index << " id = " << geom_id_by_handle(curve_to_keep_surfs[index]) << std::endl;
          std::cout << "curve_to_keep_sense " << index << " = " << curve_to_keep_senses[index] << std::endl;

        }
        for(unsigned int index=0; index < curve_to_keep_surfs.size(); index++) {

          std::cout << "curve_to_delete_surf " << index << " id = " << geom_id_by_handle(curve_to_delete_surfs[index]) << std::endl;
          std::cout << "curve_to_delete_sense " << index << " = " << curve_to_delete_senses[index] << std::endl;
        }
      } // end of debug st.

      // combine the surface and sense data for both curves into the same vector

      combined_surfs.insert( combined_surfs.end(), curve_to_keep_surfs.begin(),
                             curve_to_keep_surfs.end() );
      combined_surfs.insert( combined_surfs.end(), curve_to_delete_surfs.begin(),
                             curve_to_delete_surfs.end() );
      combined_senses.insert(combined_senses.end(),curve_to_keep_senses.begin(),
                             curve_to_keep_senses.end() );
      combined_senses.insert(combined_senses.end(),curve_to_delete_senses.begin(),
                             curve_to_delete_senses.end() );

      if(debug) {
        std::cout << combined_surfs.size() << std::endl;
        std::cout << combined_senses.size() << std::endl;
        for(unsigned int index=0; index < combined_senses.size(); index++) {

          std::cout << "combined_surfs{"<< index << "] = " << geom_id_by_handle(combined_surfs[index]) << std::endl;
          std::cout << "combined_sense["<< index << "] = " << combined_senses[index] << std::endl;
        }
      } // end debug st.

      result = gt.set_senses(merged_curve, combined_surfs, combined_senses);
      if(moab::MB_SUCCESS!=result && moab::MB_MULTIPLE_ENTITIES_FOUND!=result) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"failed to set senses: ");
      }



      // add the duplicate curve_to_keep to the vector of curves
      *j = merged_curve;

    } //end merge_tag result if st.


  } //end curves loop

  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::get_sealing_mesh_tags( double &facet_tol,
    double &sme_resabs_tol,
    moab::Tag &geom_tag,
    moab::Tag &id_tag,
    moab::Tag &normal_tag,
    moab::Tag &merge_tag,
    moab::Tag &faceting_tol_tag,
    moab::Tag &geometry_resabs_tag,
    moab::Tag &size_tag,
    moab::Tag &orig_curve_tag)
{

  moab::ErrorCode result;

  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
                                  moab::MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1,
                                  moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( "NORMAL", sizeof(moab::CartVect), moab::MB_TYPE_OPAQUE,
                                  normal_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( "MERGE", 1, moab::MB_TYPE_HANDLE,
                                  merge_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( "FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
                                  faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( "GEOMETRY_RESABS", 1,     moab::MB_TYPE_DOUBLE,
                                  geometry_resabs_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT  );
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  result = MBI()->tag_get_handle( "GEOM_SIZE", 1, moab::MB_TYPE_DOUBLE,
                                  size_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT  );
  assert( (moab::MB_SUCCESS == result) );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  int true_int = 1;
  result = MBI()->tag_get_handle( "ORIG_CURVE", 1,
                                  moab::MB_TYPE_INTEGER, orig_curve_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT, &true_int );
  assert( moab::MB_SUCCESS == result );
  if ( result != moab::MB_SUCCESS ) {
    moab_printer(result);
  }
  // PROBLEM: MOAB is not consistent with file_set behavior. The tag may not be
  // on the file_set.
  moab::Range file_set;
  result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &faceting_tol_tag,
           NULL, 1, file_set );

  MB_CHK_SET_ERR(result,"could not get faceting_tol_tag");

  if(file_set.empty()) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"file set not found");
  }

  if(1!=file_set.size()) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"Refacet with newer version of ReadCGM.");
  }

  result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1,
                                &facet_tol );
  assert(moab::MB_SUCCESS == result);
  result = MBI()->tag_get_data( geometry_resabs_tag, &file_set.front(), 1,
                                &sme_resabs_tol );
  if(moab::MB_SUCCESS != result) {
    std::cout <<  "absolute tolerance could not be read from file" << std::endl;
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::get_geometry_meshsets( moab::Range geometry_sets[], moab::Tag geom_tag, bool verbose)
{

  moab::ErrorCode result;

  // get all geometry sets
  for(unsigned dim=0; dim<4; dim++) {
    void *val[] = {&dim};
    result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &geom_tag,
             val, 1, geometry_sets[dim] );
    assert(moab::MB_SUCCESS == result);

    // make sure that sets TRACK membership and curves are ordered
    // moab::MESHSET_TRACK_OWNER=0x1, moab::MESHSET_SET=0x2, moab::MESHSET_ORDERED=0x4

    for(moab::Range::iterator i=geometry_sets[dim].begin(); i!=geometry_sets[dim].end(); i++) {
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

  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::check_for_geometry_sets(moab::Tag geom_tag, bool verbose)
{

  moab::ErrorCode result;
  // go get all geometry sets
  moab::Range geometry_sets[4];
  result = get_geometry_meshsets( geometry_sets, geom_tag, false);
  MB_CHK_SET_ERR(result,"could not get the geometry meshsets");

  //make sure they're there
  for(unsigned dim=2; dim<4; dim++) {

    if(geometry_sets[dim].size() == 0) return moab::MB_FAILURE;

  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::delete_vol(moab::EntityHandle volume)
{
  // Remove the volume set. This also removes parent-child relationships.
  moab::ErrorCode result;
  std::cout << "  deleting volume " << geom_id_by_handle(volume)  << std::endl;
  result = MBI()->delete_entities(&volume, 1);
  assert(moab::MB_SUCCESS == result);
  return result;
}


moab::ErrorCode Gen::get_meshset( const moab::EntityHandle set, std::vector<moab::EntityHandle> &vec)
{
  moab::ErrorCode result;
  vec.clear();
  result = MBI()->get_entities_by_handle( set, vec );
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode Gen::merge_verts( const moab::EntityHandle keep_vert,
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

moab::ErrorCode Gen::delete_degenerate_tris( moab::Range tris )
{
  moab::ErrorCode result;
  for(moab::Range::iterator i=tris.begin(); i!=tris.end(); i++) {
    result = delete_degenerate_tris( *i );
    assert(moab::MB_SUCCESS == result);
  }
  return moab::MB_SUCCESS;
}

// Delete degenerate triangles in the range.
moab::ErrorCode Gen::delete_degenerate_tris( moab::EntityHandle tri )
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
