#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "Arc.hpp"
#include "Zip.hpp"
#include "Cleanup.hpp"
#include "Gen.hpp"

class MakeWatertight
{


 public:

  MakeWatertight(moab::Interface* mbInterface): mbi(mbInterface) {
    gen = new Gen(mbInterface);
    arc = new Arc(mbInterface);
    zip = new Zip(mbInterface);

  };
  ~MakeWatertight() {};

  void moab_printer(moab::ErrorCode error_code);
  Gen* gen;
  Arc* arc;
  Zip* zip;
  moab::Interface* mbi;
  moab::Interface* MBI() {
    return mbi;
  }

  moab::ErrorCode delete_all_edges();

  /// gets all triangles from the mesh set and checks them for degeneracy.
  /// if any degenerate triangles are found, the program will exit
  moab::ErrorCode find_degenerate_tris();

  /// prepares all curves for the make_watertight algorithm. Merges curves that are
  /// coincident and deletes curves that are smallerthan the faceting tolerance
  moab::ErrorCode prepare_curves(moab::Range &curve_sets,
                                 moab::Tag geom_tag,
                                 moab::Tag id_tag,
                                 moab::Tag merge_tag,
                                 const double FACET_TOL,
                                 const bool debug,
                                 bool verbose = true );

  moab::ErrorCode create_arc_pair(  const double FACET_TOL,
                                    const moab::EntityHandle surf_set,
                                    std::vector<moab::EntityHandle> &skin_loop,
                                    std::vector<moab::EntityHandle> &curve_sets,
                                    const moab::EntityHandle front_endpt,
                                    const bool debug,
                                    moab::EntityHandle &curve_set,
                                    bool &curve_is_reversed,
                                    std::vector<moab::EntityHandle> &curve,
                                    std::vector<moab::EntityHandle> &skin_arc );

  moab::ErrorCode seal_arc_pair( const bool debug,
                                 const double FACET_TOL,
                                 const moab::Tag normal_tag,
                                 std::vector<moab::EntityHandle> &edge, /* in */
                                 std::vector<moab::EntityHandle> &skin /* in/out */,
                                 const int surf_id );

  moab::ErrorCode seal_loop( bool debug,
                             const double FACET_TOL,
                             const moab::Tag normal_tag,
                             const moab::Tag orig_curve_tag,
                             const moab::EntityHandle surf_set,
                             std::vector<moab::EntityHandle> &curve_sets,
                             std::vector<moab::EntityHandle> &skin_loop,
                             bool verbose = false );

  moab::ErrorCode prepare_surfaces(moab::Range &surface_sets,
                                   moab::Tag geom_tag,
                                   moab::Tag id_tag,
                                   moab::Tag normal_tag,
                                   moab::Tag merge_tag,
                                   moab::Tag orig_curve_tag,
                                   const double SME_RESABS_TOL,
                                   const double FACET_TOL,
                                   const bool debug,
                                   bool verbose = true);

  /// re-orients triangles with inverted normal vectors after being sealed
  moab::ErrorCode fix_normals(moab::Range surface_sets,
                              moab::Tag id_tag,
                              moab::Tag normal_tag,
                              const bool debug,
                              const bool verbose);

  moab::ErrorCode restore_moab_curve_representation( const moab::Range curve_sets );

  /// prints the size of every entity in geom_sets
  moab::ErrorCode get_geom_size_before_sealing( const moab::Range geom_sets[],
      const moab::Tag geom_tag,
      const moab::Tag size_tag,
      bool debug,
      bool verbose);
  /// prints changes in size to the mesh after sealing
  moab::ErrorCode get_geom_size_after_sealing( const moab::Range geom_sets[],
      const moab::Tag geom_tag,
      const moab::Tag size_tag,
      const double FACET_TOL,
      bool debug,
      bool verbose );

  /// deletes all curves with a merge_tag and removes them from the curve sets of the mesh
  moab::ErrorCode delete_merged_curves(moab::Range &existing_curve_sets, moab::Tag merge_tag, bool debug = false);

  /// deletes all tags created for use in sealing the model
  moab::ErrorCode delete_sealing_tags( moab::Tag normal_tag, moab::Tag merge_tag, moab::Tag size_tag, moab::Tag orig_curve_tag);

  /// gets all curve sets an returns them in curves. Places any unmerged curves in unmerged_curves.
  moab::ErrorCode get_unmerged_curves( moab::EntityHandle surface,
                                       std::vector<moab::EntityHandle> &curves,
                                       std::vector<moab::EntityHandle> &unmerged_curves,
                                       moab::Tag merge_tag,
                                       bool verbose,
                                       bool debug);

  /// takes the skin_edges from the moab skinner and creates loops of vertices between the facets and geometric curves.
  /// The vertex loops are returned in the vector array, skin.
  moab::ErrorCode create_skin_vert_loops( moab::Range &skin_edges, moab::Range tris, std::vector < std::vector <moab::EntityHandle> > &skin, int surf_id, bool &cont, bool debug);

  /// merges any skin vertices closer in proximity than the SME_RESABS_TOL.
  /// It then checks the skins for any degenerate edges resultant of vertex merging.
  moab::ErrorCode merge_skin_verts ( moab::Range &skin_verts, moab::Range &skin_edges, double SME_RESABS_TOL, int surf_id, bool cont, bool debug);

  /// runs the make_watertight algorithm on each set of skin_loops for the surface, surf.
  moab::ErrorCode seal_surface_loops ( moab::EntityHandle surf,
                                       moab::EntityHandle skin_loops[],
                                       std::vector < std::vector<moab::EntityHandle> > skin,
                                       std::vector<moab::EntityHandle> curves,
                                       moab::Tag normal_tag,
                                       moab::Tag orig_curve_tag,
                                       double FACET_TOL,
                                       int surf_id,
                                       bool debug);

  /// takes the mesh in input_set and makes it watertight
  moab::ErrorCode make_mesh_watertight(moab::EntityHandle input_set, double &facet_tol, bool verbose = true);

};
