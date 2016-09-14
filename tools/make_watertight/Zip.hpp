#ifndef ZIP_HPP
#define ZIP_HPP

#include "moab/Core.hpp"
#include "Gen.hpp"

class Zip
{
 public:
  Zip(moab::Interface *mbInterface) : mbi(mbInterface) {
    gen = new Gen(mbInterface);
  };
  ~Zip() {};

  Gen* gen;
  moab::Interface* mbi;
  moab::Interface* MBI() {
    return mbi;
  };

  moab::ErrorCode order_verts_by_edge( moab::Range unordered_edges, std::vector<moab::EntityHandle> &ordered_verts );

  moab::ErrorCode t_joint( moab::Tag normal_tag,
                           const moab::EntityHandle vert0,
                           const moab::EntityHandle vert1,
                           const moab::EntityHandle vert2,
                           bool debug );
  /// removes the entitiy handle tri from the loaded mesh
  moab::ErrorCode delete_degenerate_tris( moab::EntityHandle tri );
  /// checks that no triangles in the moab::Range tris are degenterate. If
  /// degenerate triangles are found, they are deleted from the mesh.
  moab::ErrorCode delete_degenerate_tris( moab::Range tris );

  moab::ErrorCode delete_adj_degenerate_tris( const moab::EntityHandle adj_vert );


  /// merges two vertices by updating the entity handle of the deleted vert to the vert to
  /// keep in the correct arc. Uses MOAB function merge_entities to merge the vertices in
  /// the database. Also deletes the triangles adjacent to the merged vertices if one
  /// becomes degenerate.
  moab::ErrorCode merge_verts( const moab::EntityHandle keep_vert,
                               const moab::EntityHandle delete_vert,
                               std::vector<moab::EntityHandle> &arc0,
                               std::vector<moab::EntityHandle> &arc1 );

  /// test two normal vectors to see if they point in the same direction
  moab::ErrorCode test_normals( const std::vector<moab::CartVect> norms0,
                                const std::vector<moab::CartVect> norms1,
                                std::vector<int> &inverted_tri_indices );
  moab::ErrorCode test_normals( const             moab::CartVect  norms0,
                                const             moab::CartVect  norms1 );

  moab::ErrorCode remove_inverted_tris(moab::Tag normal_tag, moab::Range tris, const bool debug );

  /// tests the watertightness of all arcs in the vector-array of moab entity handles arcs
  moab::ErrorCode test_zipping( const double FACET_TOL,
                                const std::vector< std::vector<moab::EntityHandle> > arcs );


};


struct triangles {
  moab::EntityHandle before_tri;
  const moab::EntityHandle *before;
  moab::CartVect     before_norm;
  moab::EntityHandle after0[3];
  moab::EntityHandle after1[3];
  moab::CartVect     after0_norm;
  moab::CartVect     after1_norm;
  moab::EntityHandle surf_set;
};

#endif

