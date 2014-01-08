#ifndef ARC_HPP
#define ARC_HPP

#include <vector>
#include "MBCore.hpp"
#include "gen.hpp"
#include "MBSkinner.hpp"
#include "MBTagConventions.hpp"

MBInterface *MBI(); 
namespace arc {

/// check that edge is going in the same direction as one of the edges on tri.
/// If this is not the case, the edge is reversed.
  MBErrorCode orient_edge_with_tri( const MBEntityHandle edge, 
                                    const MBEntityHandle tri );
// checks for degeneracy of edges in the MBRange edges and deletes degenerates if found
  MBErrorCode remove_degenerate_edges( MBRange &edges, const bool debug );

/// deletes any duplicate edges in MBRange edges for which one goes from vertex a to 
/// vertex b and the other from b to a
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

  MBErrorCode order_verts_by_edge( MBEntityHandle meshset, MBRange unordered_edges, std::vector<MBEntityHandle> &ordered_verts );

/// gets the moab entities in the meshset, set, and returns them to vec
  MBErrorCode get_meshset( const MBEntityHandle set, std::vector<MBEntityHandle> &vec);

/// clears the given meshset set and then adds the entities desired to the meshset
/// (apparently child_parent_relations are taken care of here? Edges are created how?)
  MBErrorCode set_meshset( const MBEntityHandle set, const std::vector<MBEntityHandle> vec );

/// goes through curve_sets and finds any curves with coincident ( dist. apart <= FACET_TOL) front and back points.
/// it then merges the curves topologically. Any merged curves aren't deleted until prepare surfaces. 
  MBErrorCode merge_curves(MBRange curve_sets, const double FACET_TOL,
                           MBTag idTag, MBTag merge_tag, const bool debug );
}

#endif
