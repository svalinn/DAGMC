#ifndef ARC_HPP
#define ARC_HPP

#include <vector>
#include "MBCore.hpp"
#include "gen.hpp"
#include "MBSkinner.hpp"
#include "MBTagConventions.hpp"

MBInterface *MBI(); 
namespace arc {
  MBErrorCode orient_edge_with_tri( const MBEntityHandle edge, 
                                    const MBEntityHandle tri );

  MBErrorCode remove_degenerate_edges( MBRange &edges, const bool debug );

  MBErrorCode remove_opposite_pairs_of_edges( MBRange &edges, const bool debug );
  MBErrorCode remove_opposite_pairs_of_edges_fast( MBRange &edges, const bool debug );

  MBErrorCode get_next_oriented_edge( const MBRange edges, 
                                      const MBEntityHandle edge,
				      MBEntityHandle &next_edge );

  // Given a range of edges and a vertex, find the edge the contains the
  // endpoint. Also return the opposite endpoint of the edge. This checks
  // to ensure that only one edge is found.
  MBErrorCode get_next_edge_and_vert_by_edge( const MBRange edges_in,
					      const MBEntityHandle edge_in,
					      const MBEntityHandle vertex_in,
					      MBEntityHandle &edge_out,
					      MBEntityHandle &vertex_out       );

  MBErrorCode create_loops_from_oriented_edges_fast( MBRange edges,
						std::vector< std::vector<MBEntityHandle> > &loops_of_edges,
                                                const bool debug );
  MBErrorCode create_loops_from_oriented_edges( MBRange edges,
						std::vector< std::vector<MBEntityHandle> > &loops_of_edges,
                                                const bool debug );

  MBErrorCode order_verts_by_edge( MBRange unordered_edges, std::vector<MBEntityHandle> &ordered_verts );

  MBErrorCode get_meshset( const MBEntityHandle set, std::vector<MBEntityHandle> &vec);
  MBErrorCode set_meshset( const MBEntityHandle set, const std::vector<MBEntityHandle> vec );

  MBErrorCode merge_curves(MBRange curve_sets, const double FACET_TOL,
                           MBTag idTag, MBTag merge_tag, const bool debug );
}

#endif
