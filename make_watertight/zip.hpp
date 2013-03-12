#ifndef ZIP_HPP
#define ZIP_HPP

#include "MBCore.hpp"
#include "gen.hpp"
#include "arc.hpp"

MBInterface *MBI();
namespace zip {
  MBErrorCode t_joint( MBTag normal_tag, 
                       const MBEntityHandle vert0,              
                       const MBEntityHandle vert1,                         
                       const MBEntityHandle vert2 );
        
  MBErrorCode delete_degenerate_tris( MBEntityHandle tri );
  MBErrorCode delete_degenerate_tris( MBRange tris );
  MBErrorCode delete_adj_degenerate_tris( const MBEntityHandle adj_vert );

  MBErrorCode merge_verts( const MBEntityHandle keep_vert, 
                           const MBEntityHandle delete_vert,
                           std::vector<MBEntityHandle> &arc0,
                           std::vector<MBEntityHandle> &arc1 );

  // test two normal vectors to see if they point in the same direction
  MBErrorCode test_normals( const std::vector<MBCartVect> norms0, 
                            const std::vector<MBCartVect> norms1,
                            std::vector<int> &inverted_tri_indices );
  MBErrorCode test_normals( const             MBCartVect  norms0, 
                            const             MBCartVect  norms1 );

  MBErrorCode remove_inverted_tris(MBTag normal_tag, MBRange tris, const bool debug );

  MBErrorCode test_zipping( const double FACET_TOL,
                            const std::vector< std::vector<MBEntityHandle> > arcs );

}

#endif
