#ifndef GEN_HPP
#define GEN_HPP

#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"

MBInterface *MBI(); 
namespace gen {
  bool error( const bool error_has_occured, const std::string message="" );
  void print_vertex_cubit( const MBEntityHandle vertex );
  void print_vertex_coords( const MBEntityHandle vertex );

  void print_triangles( const MBRange tris );
  void print_triangle( const MBEntityHandle triangle, bool print_edges );

  void print_edge( const MBEntityHandle edge );

  void print_vertex_count( const MBEntityHandle input_meshset);

  void print_range( const MBRange range );

  void print_range_of_edges( const MBRange range );

  void print_arc_of_edges( const std::vector<MBEntityHandle> arc_of_edges );

  void print_arcs( const std::vector < std::vector<MBEntityHandle> > arcs );

  void print_loop( const std::vector<MBEntityHandle> loop_of_verts );

MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
                               const std::vector<MBEntityHandle> arc_of_verts,
                               unsigned &position,
                               const double dist_limit );

  MBErrorCode find_closest_vert( const double tol,
                                 const MBEntityHandle reference_vert,
                                 const std::vector<MBEntityHandle> loop_of_verts,
                                 std::vector<unsigned> &positions, 
                                 std::vector<double>   &dists);
  /*  MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
                                  const std::vector<std::vector<MBEntityHandle> > loops_of_verts,
                                  unsigned int &loop, unsigned int &position, 
                                  double &min_dist);
  */
  // Merge the range of vertices. We do not merge by edges (more
  // stringent) because we do not want to miss corner vertices.
  MBErrorCode merge_vertices( MBRange vertices /* in */, 
			      const  double tol       /* in */);
			      //bool &merge_vertices_again /* out */);

  MBErrorCode squared_dist_between_verts( const MBEntityHandle v0, 
                                          const MBEntityHandle v1, 
                                          double &d);
  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 );
  MBErrorCode dist_between_verts( const MBEntityHandle v0, const MBEntityHandle v1,
                                  double &d );
  double dist_between_verts( double coords0[], double coords1[] );
  double dist_between_verts( MBEntityHandle vert0, MBEntityHandle vert1 );                             

  // Return the length of the curve defined by MBEDGEs or ordered MBVERTEXs.
  double length( std::vector<MBEntityHandle> curve );

  // Given a vertex and vector of edges, return the number of edges adjacent to the vertex.
  unsigned int n_adj_edges( MBEntityHandle vert, MBRange edges );

// Return true if the edges share a vertex. Does not check for coincident edges.
  bool edges_adjacent( MBEntityHandle edge0, MBEntityHandle edge1 );

// get the direction unit vector from one vertex to another vertex
  //MBErrorCode get_direction( MBEntityHandle from_vert, MBEntityHandle to_vert, double dir[] );
  MBErrorCode get_direction( const MBEntityHandle from_vert, const MBEntityHandle to_vert,
                           MBCartVect &dir ); 

// from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
  double edge_point_dist( const MBCartVect a, const MBCartVect b, const MBCartVect c );
  double edge_point_dist( const MBEntityHandle endpt0, const MBEntityHandle endpt1, 
                          const MBEntityHandle pt );
  double edge_point_dist( const MBEntityHandle edge, const MBEntityHandle pt );

  MBErrorCode point_curve_min_dist( const std::vector<MBEntityHandle> curve, 
                                    const MBEntityHandle pt, double &min_dist,
                                    const double max_dist_along_curve );
  MBErrorCode point_curve_min_dist( const std::vector<MBEntityHandle> curve, 
                                    const MBEntityHandle pt, double &min_dist );

  double triangle_area( const MBCartVect a, const MBCartVect b, const MBCartVect c);
  MBErrorCode triangle_area( const MBEntityHandle conn[], double &area );
  MBErrorCode triangle_area( const MBEntityHandle triangle, double &area );
  double triangle_area( MBRange triangles );
  
  bool triangle_degenerate( const MBEntityHandle triangle );
  bool triangle_degenerate( const MBEntityHandle v0, const MBEntityHandle v1, const MBEntityHandle v2);

  MBErrorCode triangle_normals( const MBRange triangles, std::vector<MBCartVect> &normals );
  MBErrorCode triangle_normal( const MBEntityHandle triangle, MBCartVect &normal );
  MBErrorCode triangle_normal( const MBEntityHandle v0, const MBEntityHandle v1,
                               const MBEntityHandle v2, MBCartVect &normal );
  MBErrorCode triangle_normal( const MBCartVect v0, const MBCartVect v1, 
                               const MBCartVect v2, MBCartVect &normal ); 

  // Distance between a point and line. The line is defined by two verts.
  // We are using a line and not a line segment!
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  MBErrorCode line_point_dist( const MBEntityHandle line_pt1, const MBEntityHandle line_pt2,
                               const MBEntityHandle pt0, double &dist );

  // Project the point onto the line. Not the line segment! 
  // Change the coordinates of the pt0 to the projection.
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,     
                                     const MBEntityHandle line_pt2,             
                                     const MBEntityHandle pt0 );

  // Do not change the coords of pt0. Instead return the projected coords.              
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,     
                                     const MBEntityHandle line_pt2,             
                                     const MBEntityHandle pt0,              
                                     MBCartVect &projected_coords,
                                     double &parameter );  
  // Get the distance of pt0 from line_pt1 if pt0 is projected. No coords are
  // changed in MOAB.
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,
				     const MBEntityHandle line_pt2,
				     const MBEntityHandle pt0,
				     double &dist_along_edge  );
  
  /*  double area2( const MBEntityHandle a, const MBEntityHandle b, const MBEntityHandle c);
  bool left( const MBEntityHandle a, const MBEntityHandle b, const MBEntityHandle c);
  bool left_on( const MBEntityHandle a, const MBEntityHandle b, const MBEntityHandle c);
  bool collinear( const MBEntityHandle a, const MBEntityHandle b, const MBEntityHandle c);
  bool logical_xor( const bool x, const bool y );
  bool intersect_prop( const MBEntityHandle a, const MBEntityHandle b, 
                       const MBEntityHandle c, const MBEntityHandle d);
  bool between( const MBEntityHandle a, const MBEntityHandle b, const MBEntityHandle c);
  bool intersect( const MBEntityHandle a, const MBEntityHandle b, 
                       const MBEntityHandle c, const MBEntityHandle d);
  bool diagonalie( const MBEntityHandle a, const MBEntityHandle b,
                   const std::vector<MBEntityHandle> verts ); 
  bool in_cone( const MBEntityHandle a, const MBEntityHandle b,
		const std::vector<MBEntityHandle> verts ); 
  bool diagonal( const MBEntityHandle a, const MBEntityHandle b,
		 const std::vector<MBEntityHandle> verts ); 
  MBErrorCode ear_init( std::vector<MBEntityHandle> polygon_of_verts,
                        std::vector<bool> &is_ear );
  */  MBErrorCode ear_clip_polygon( std::vector<MBEntityHandle> polygon_of_verts,
                                    const MBCartVect plane_normal_vector, MBRange &new_tris );

  int geom_id_by_handle( const MBEntityHandle set);

  MBErrorCode save_normals( MBRange tris, MBTag normal_tag );

  MBErrorCode flip(const MBEntityHandle tri, const MBEntityHandle vert0, 
                   const MBEntityHandle vert2, const MBEntityHandle surf_set);

  MBErrorCode ordered_verts_from_ordered_edges( const std::vector<MBEntityHandle> ordered_edges,
                                                std::vector<MBEntityHandle> &ordered_verts );

  MBErrorCode dist_between_arcs( bool debug,
                         const std::vector<MBEntityHandle> arc0,
                         const std::vector<MBEntityHandle> arc1,
                         double &dist );

  // skin edges are a vector of two vertex handles
    // Hold edges in an array of handles.
  struct edge {
    MBEntityHandle edge, v0, v1;
  };
  int compare_edge(const void *a, const void *b);
  MBErrorCode find_skin( MBRange tris, const int dim,                     
			 // std::vector<std::vector<MBEntityHandle> > &skin_edges,    
			 MBRange &skin_edges,                         
                         const bool );
  //MBErrorCode find_skin( const MBRange tris, const int dim, MBRange &skin_edges, const bool );
  MBErrorCode measure( const MBEntityHandle set, const MBTag geom_tag, double &size );

  // Given a curve and surface set, get the relative sense.
  // From CGMA/builds/dbg/include/CubitDefines, CUBIT_UNKNOWN=-1, CUBIT_FORWARD=0, CUBIT_REVERSED=1
  MBErrorCode get_curve_surf_sense( const MBEntityHandle surf_set, const MBEntityHandle curve_set,
                                    int &sense );


}

#endif
