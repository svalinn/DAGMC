#ifndef CLEANUP_HPP
#define CLEANUP_HPP

// moab includes
#include "moab/Core.hpp"

#include "Gen.hpp"
#include "Arc.hpp"
#include "Zip.hpp"



class Cleanup
{
 public:
  Cleanup(moab::Interface *mbInterface) : mbi(mbInterface) {
    gen = new Gen(mbInterface);
  };
  ~Cleanup() {};

  Gen* gen;
  moab::Interface* mbi;
  moab::Interface* MBI() {
    return mbi;
  };

  /// The obbtrees are no longer valid because the triangles have been altered.
  ///  -Surface and volume sets are tagged with tags holding the obb tree
  ///   root handles.
  ///  -Surface/volume set handles are added to the root meshset.
  /// Somehow, delete the old tree without deleting the
  /// surface and volume sets, then build a new tree.
  moab::ErrorCode remove_obb_tree(bool verbose = false);

  moab::ErrorCode delete_small_edge_and_tris( const moab::EntityHandle vert0,
      moab::EntityHandle &vert1, const double tol);
  moab::ErrorCode delete_small_edges( const moab::Range &surfaces, const double MERGE_TOL);

  /// Lots of edges have been created but are no longer needed.
  /// Delete edges that are not in curves. These should be the only edges
  /// that remain. This incredibly speeds up the watertight_check tool (100x?).
  moab::ErrorCode cleanup_edges( moab::Range curve_meshsets );
};

#endif
