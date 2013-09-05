#ifndef DAGMC_TRACK_LENGTH_MESH_TALLY_H
#define DAGMC_TRACK_LENGTH_MESH_TALLY_H


#include <string>
#include <cassert>
#include <set>

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"

#include "Matrix3.hpp"
#include "MeshTally.hpp"


struct MCNPTrackParam {
  int *fmesh_index;
  double *erg;
  double *wgt;

  MCNPTrackParam(int *i, double *e, double *w): 
    fmesh_index(i), erg(e), wgt(w) {}
};

namespace moab{

/* Forward Declarations */
class AdaptiveKDTree;
class OrientedBoxTreeTool;


class TrackLengthMeshTally : public MeshTally{ 

public:
  /**
   * Public constructor interface- actual constructor is protected.
   */
  static TrackLengthMeshTally* setup( const fmesh_card& params, Interface* mbi, const int* cur_mcnp_cell );

  /**
   * Tally a single particle track segment
   */
  void add_track_segment( CartVect& start, CartVect& vec, double length, int ebin,  MCNPTrackParam* score_params );
  
  virtual void end_history ();
    
  virtual void print(  double sp_norm, double mult_fact );

  ~TrackLengthMeshTally();

 

protected:

  Range get_adjacency_info(Range input_handles);
  
  bool point_in_tet_ad(CartVect pos, EntityHandle tet);

  int tet_index(Range tets, EntityHandle tet_handle);

  EntityHandle point_in_which_tet (CartVect point);

  EntityHandle remainder( CartVect start, CartVect dir, double distance, double left_over);

  ErrorCode load_mesh( const std::string& input_filename, 
                       std::string tag_name, std::vector<std::string>& tag_values );  
  ErrorCode write_results( double sp_norm, double mult_fact, 
                           const std::string* override_output_filename = NULL );

  bool point_in_tet( const CartVect& point, const EntityHandle* tet );

  void add_score_to_mesh_cell( EntityHandle mesh_cell, double score, int ebin );
 
  void get_skin_triangle_adjacencies( EntityHandle triangle, 
                                      EntityHandle& tetrahedron, EntityHandle vertices[3] );


  EntityHandle find_next_tet_by_ray_fire( CartVect& start, CartVect& vec, double length, 
                                          EntityHandle first_tri[3], double& first_t, 
                                          EntityHandle last_crossing = 0 );

  EntityHandle get_starting_tet( CartVect& start, CartVect& vec, double length, 
                                 EntityHandle first_tri[3], double& first_t );
    

  EntityHandle get_starting_tet_conformal( CartVect& start, EntityHandle first_tri[3] );


  void set_convex_flag( bool c ){ convex = c; } 

  Interface* mb;

  std::string output_filename;

  EntityHandle tally_set;

  AdaptiveKDTree* kdtree;
  EntityHandle kdtree_root;

  OrientedBoxTreeTool* obb_tool;
  EntityHandle obbtree_root;

  EntityHandle last_visited_tet;

  std::vector< Matrix3 > tet_baryc_data;

  bool convex; // true if user asserts this mesh tally has convex geometry

  // if non-null, user has asserted mesh tally geometry
  // conforms to the cells identified in this set
  std::set<int>* conformality; 
  bool conformal_surface_source;

  const int* mcnp_current_cell; // non-null if user has asserted conformal geometry
  int last_cell;
  std::set<EntityHandle> visited_this_history; 
    
  TrackLengthMeshTally( const fmesh_card& fmesh_params, Interface* mb_p, const std::string& output_filename );

private:
  TrackLengthMeshTally& operator=( const TrackLengthMeshTally& mt ); // unimplemented
  TrackLengthMeshTally( const TrackLengthMeshTally& mt ); // unimplemented

};

} // end namespace moab


#endif /* DAGMC_TRACK_LENGTH_MESH_TALLY_H */
