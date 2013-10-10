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

MBErrorCode find_degenerate_tris();

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
                       std::vector<MBEntityHandle> &skin_loop );

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

MBErrorCode remove_surf_sense_data(MBEntityHandle del_surf);

MBErrorCode fix_normals(MBRange surface_sets, 
                        MBTag id_tag, 
                        MBTag normal_tag,
                        const bool debug,
                        const bool verbose);

MBErrorCode restore_moab_curve_representation( const MBRange curve_sets );

MBErrorCode get_geom_size_before_sealing( const MBRange geom_sets[], 
                                          const MBTag geom_tag,
                                          const MBTag size_tag,
                                          bool verbose = true );

MBErrorCode get_geom_size_after_sealing( const MBRange geom_sets[], 
                                         const MBTag geom_tag,
                                         const MBTag size_tag,
                                         const double FACET_TOL,
                                         bool verbose = true );

MBErrorCode make_mesh_watertight(MBEntityHandle input_set, double &facet_tol, bool verbose = true);




}
