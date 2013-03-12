#ifndef CLEANUP_HPP
#define CLEANUP_HPP

#include "MBCore.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

MBInterface *MBI();
namespace cleanup {

  // The obbtrees are no longer valid because the triangles have been altered.
  //  -Surface and volume sets are tagged with tags holding the obb tree
  //   root handles.
  //  -Surface/volume set handles are added to the root meshset.
  // Somehow, delete the old tree without deleting the
  // surface and volume sets, then build a new tree.
  MBErrorCode remove_obb_tree();

  MBErrorCode delete_small_edge_and_tris( const MBEntityHandle vert0, 
                                          MBEntityHandle &vert1, const double tol);
  MBErrorCode delete_small_edges( const MBRange &surfaces, const double MERGE_TOL);

  // Lots of edges have been created but are no longer needed.
  // Delete edges that are not in curves. These should be the only edges
  // that remain. This incredibly speeds up the watertight_check tool (100x?).
  MBErrorCode cleanup_edges( MBRange curve_meshsets );  
}

#endif
