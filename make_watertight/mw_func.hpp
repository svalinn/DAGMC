#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"




MBInterface *MOAB();

namespace mw_func {

void moab_printer(MBErrorCode error_code);

MBErrorCode delete_all_edges();

/// gets all triangles from the mesh set and checks them for degeneracy.
/// if any degenerate triangles are found, the program will exit
MBErrorCode find_degenerate_tris();

/// prepares all curves for the make_watertight algorithm. Merges curves that are 
/// coincident and deletes curves that are smallerthan the faceting tolerance
MBErrorCode prepare_curves(MBRange &curve_sets, 
                           MBTag geom_tag, 
                           MBTag id_tag, 
                           MBTag merge_tag, 
                           const double FACET_TOL, 
                           const bool debug,
                           bool verbose = true );

MBErrorCode create_arc_pair(  const double FACET_TOL,
                              const MBEntityHandle surf_set,
			      std::vector<MBEntityHandle> &skin_loop,
			      std::vector<MBEntityHandle> &curve_sets,
			      const MBEntityHandle front_endpt,
                              const bool debug,
			      MBEntityHandle &curve_set,
			      bool &curve_is_reversed,
			      std::vector<MBEntityHandle> &curve,
			      std::vector<MBEntityHandle> &skin_arc );

MBErrorCode seal_arc_pair( const bool debug,
                           const double FACET_TOL,
                           const MBTag normal_tag,
                           std::vector<MBEntityHandle> &edge, /* in */
                           std::vector<MBEntityHandle> &skin /* in/out */,
                           const int surf_id );

MBErrorCode seal_loop( bool debug,
                       const double FACET_TOL,
                       const MBTag normal_tag,
                       const MBTag orig_curve_tag,
                       const MBEntityHandle surf_set,
                       std::vector<MBEntityHandle> &curve_sets,
                       std::vector<MBEntityHandle> &skin_loop,
                       bool verbose = false );

MBErrorCode prepare_surfaces(MBRange &surface_sets,
                             MBTag geom_tag, 
                             MBTag id_tag, 
                             MBTag normal_tag, 
                             MBTag merge_tag,
                             MBTag orig_curve_tag,
                             const double SME_RESABS_TOL,
                             const double FACET_TOL, 
                             const bool debug,
                             bool verbose = true);

/// re-orients triangles with inverted normal vectors after being sealed
MBErrorCode fix_normals(MBRange surface_sets, 
                        MBTag id_tag, 
                        MBTag normal_tag,
                        const bool debug,
                        const bool verbose);

MBErrorCode restore_moab_curve_representation( const MBRange curve_sets );

/// prints the size of every entity in geom_sets
MBErrorCode get_geom_size_before_sealing( const MBRange geom_sets[], 
                                          const MBTag geom_tag,
                                          const MBTag size_tag,
                                          bool verbose = true );
/// prints changes in size to the mesh after sealing
MBErrorCode get_geom_size_after_sealing( const MBRange geom_sets[], 
                                         const MBTag geom_tag,
                                         const MBTag size_tag,
                                         const double FACET_TOL,
                                         bool verbose = true );

/// deletes all curves with a merge_tag and removes them from the curve sets of the mesh
MBErrorCode delete_merged_curves(MBRange &existing_curve_sets, MBTag merge_tag, bool debug = false);

/// deletes all tags created for use in sealing the model
MBErrorCode delete_sealing_tags( MBTag normal_tag, MBTag merge_tag, MBTag size_tag, MBTag orig_curve_tag);

/// gets all curve sets an returns them in curves. Places any unmerged curves in unmerged_curves. 
MBErrorCode get_unmerged_curves( MBEntityHandle surface, 
                                 std::vector<MBEntityHandle> &curves, 
                                 std::vector<MBEntityHandle> &unmerged_curves, 
                                 MBTag merge_tag, 
                                 bool verbose);

/// takes the skin_edges from the moab skinner and creates loops of vertices between the facets and geometric curves.
/// The vertex loops are returned in the vector array, skin. 
MBErrorCode create_skin_vert_loops( MBRange &skin_edges, MBRange tris, std::vector < std::vector <MBEntityHandle> > &skin, int surf_id, bool &cont, bool debug);

/// merges any skin vertices closer in proximity than the SME_RESABS_TOL.
/// It then checks the skins for any degenerate edges resultant of vertex merging.
MBErrorCode merge_skin_verts ( MBRange &skin_verts, MBRange &skin_edges, double SME_RESABS_TOL, int surf_id, bool cont, bool debug);

/// runs the make_watertight algorithm on each set of skin_loops for the surface, surf.
MBErrorCode seal_surface_loops ( MBEntityHandle surf,
                                 MBEntityHandle skin_loops[], 
                                 std::vector < std::vector<MBEntityHandle> > skin, 
                                 std::vector<MBEntityHandle> curves, 
                                 MBTag normal_tag, 
                                 MBTag orig_curve_tag, 
                                 double FACET_TOL, 
                                 int surf_id, 
                                 bool debug);

/// takes the mesh in input_set and makes it watertight
MBErrorCode make_mesh_watertight(MBEntityHandle input_set, double &facet_tol, bool verbose = true);




}
