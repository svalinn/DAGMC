#ifndef ARC_HPP
#define ARC_HPP

#include <vector>
#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "MBTagConventions.hpp"
#include "Gen.hpp"
#include "Zip.hpp"

class Arc
{
  friend class Gen;

 public:
  Arc(moab::Interface *mbInterface) : mbi(mbInterface) {
    gen = new Gen(mbInterface);
    zip = new Zip(mbInterface);
  };
  ~Arc() {};

  Zip* zip;
  Gen* gen;
  moab::Interface* mbi;
  moab::Interface* MBI() {
    return mbi;
  };

  /// check that edge is going in the same direction as one of the edges on tri.
  /// If this is not the case, the edge is reversed.
  moab::ErrorCode orient_edge_with_tri( const moab::EntityHandle edge,
                                        const moab::EntityHandle tri );
  // checks for degeneracy of edges in the moab::Range edges and deletes degenerates if found
  moab::ErrorCode remove_degenerate_edges( moab::Range &edges, const bool debug );

  /// deletes any duplicate edges in moab::Range edges for which one goes from vertex a to
  /// vertex b and the other from b to a
  moab::ErrorCode remove_opposite_pairs_of_edges( moab::Range &edges, const bool debug );
  moab::ErrorCode remove_opposite_pairs_of_edges_fast( moab::Range &edges, const bool debug );

  moab::ErrorCode get_next_oriented_edge( const moab::Range edges,
                                          const moab::EntityHandle edge,
                                          moab::EntityHandle &next_edge );

  // Given a range of edges and a vertex, find the edge the contains the
  // endpoint. Also return the opposite endpoint of the edge. This checks
  // to ensure that only one edge is found.
  moab::ErrorCode get_next_edge_and_vert_by_edge( const moab::Range edges_in,
      const moab::EntityHandle edge_in,
      const moab::EntityHandle vertex_in,
      moab::EntityHandle &edge_out,
      moab::EntityHandle &vertex_out       );

  moab::ErrorCode create_loops_from_oriented_edges_fast( moab::Range edges,
      std::vector< std::vector<moab::EntityHandle> > &loops_of_edges,
      const bool debug );
  moab::ErrorCode create_loops_from_oriented_edges( moab::Range edges,
      std::vector< std::vector<moab::EntityHandle> > &loops_of_edges,
      const bool debug );

  moab::ErrorCode order_verts_by_edge( moab::Range unordered_edges, std::vector<moab::EntityHandle> &ordered_verts );

  /// gets the moab entities in the meshset, set, and returns them to vec
  moab::ErrorCode get_meshset( const moab::EntityHandle set, std::vector<moab::EntityHandle> &vec);

  /// clears the given meshset set and then adds the entities desired to the meshset
  /// (apparently child_parent_relations are taken care of here? Edges are created how?)
  moab::ErrorCode set_meshset( const moab::EntityHandle set, const std::vector<moab::EntityHandle> vec );

  moab::ErrorCode merge_verts( const moab::EntityHandle keep_vert,
                               const moab::EntityHandle delete_vert,
                               std::vector<moab::EntityHandle> &arc0,
                               std::vector<moab::EntityHandle> &arc1 );

  /// goes through curve_sets and finds any curves with coincident ( dist. apart <= FACET_TOL) front and back points.
  /// it then merges the curves topologically. Any merged curves aren't deleted until prepare surfaces.
  moab::ErrorCode merge_curves(moab::Range curve_sets, const double FACET_TOL,
                               moab::Tag idTag, moab::Tag merge_tag, const bool debug );
};

#endif
