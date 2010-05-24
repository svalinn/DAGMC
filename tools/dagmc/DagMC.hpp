#ifndef MOABMC_HPP
#define MOABMC_HPP

#include "moab/Interface.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <string>
#include <assert.h>

#include "moab/OrientedBoxTreeTool.hpp"

#ifdef CGM
class RefEntity;
#endif

namespace moab {

class CartVect;

#define DAGMC_VERSION 0.99
#define DAGMC_VERSION_STRING "0.99"
#define DAGMC_INTERFACE_REVISION "$Rev$"

class DagMC 
{
public:
  static DagMC *instance(Interface *mb_impl = NULL);
  
  ~DagMC();
  
    // Return the version of this library
  float version(std::string *version_string = NULL);

  ErrorCode ray_fire(const EntityHandle cell, const EntityHandle last_surf_hit, 
                       const int num_pts,
                       const double uuu, const double vvv, const double www,
                       const double xxx, const double yyy, const double zzz,
                       const double huge_val,
                       double &dist_traveled, EntityHandle &next_surf_hit,
                       OrientedBoxTreeTool::TrvStats* stats = NULL );

    // Test if point is inside or outside of a volume using unit sphere area method
    // Requires sense of surfaces wrt volume.
    // Return values: {0 : outside, 1: inside, -1: on boundary}
  ErrorCode point_in_volume_slow( EntityHandle volume, 
                                    double x, double y, double z,
                                    int& result );

    // Test if point is inside or outside of a volume using an OBB tree
    // Requires sense of surfaces wrt volume.
    // Return values: {0 : outside, 1: inside, -1: on boundary}
  ErrorCode point_in_volume( EntityHandle volume, 
                               double x, double y, double z,
                               int& result,
                               double u, double v, double w);

    // Find the distance to the nearest surface
  ErrorCode closest_to_location( EntityHandle volume,
                                   double* point,
                                   double& result);

    // Calculate the volume contained in a 'volume'
  ErrorCode measure_volume( EntityHandle volume, double& result );

    // Calculate sum of area of triangles
  ErrorCode measure_area( EntityHandle surface, double& result );

    // Get the sense of surfaces wrt a volume.  Sense values are:
    // {-1 -> reversed, 0 -> both, 1 -> forward}
  ErrorCode surface_sense( EntityHandle volume, 
                             int num_surfaces,
                             const EntityHandle* surfaces,
                             int* senses_out );

    // Get the sense of surfaces wrt a volume.  Sense values are:
    // {-1 -> reversed, 0 -> both, 1 -> forward}
  ErrorCode surface_sense( EntityHandle volume, EntityHandle surface, int& sense_out );

    // load mesh
  ErrorCode load_file(const char* cfile,
			const double facet_tolerance = 0);

  ErrorCode write_mesh(const char* ffile,
			 const int flen);


    // initialize data structures and OBB tree
  ErrorCode init_OBBTree();

    // map between EntityHandle, base-1 index, and GLOBAL_ID
  EntityHandle entity_by_index( int dimension, int index );
  EntityHandle entity_by_id( int dimension, int id );
  int index_by_handle( EntityHandle handle );
  int id_by_index( int dimension, int index );

    // get number of geometric sets corresponding to geometry
    // of specified dimension (i.e. number of geometric surfaces
    // if dimension == 2).
  int num_entities( int dimension );

    // read/write settings from file
  void read_settings( const char* filename );
  void write_settings( FILE* filp, bool with_description = true );
  void parse_settings();

    // pass settings from calling application
  void set_settings(int source_cell, int use_cad, int use_dist_limit,
		    double add_distance_tolerance,
		    double discard_distance_tolerance);
  void get_settings(int* source_cell, int* use_cad, int* use_dist_limit,
		    double* discard_distance_tolerance);

  char *get_spec_reflect();
  char *get_white_reflect();

  OrientedBoxTreeTool *obb_tree() {return &obbTree;}

    // Get the tag used to associate OBB trees with geometry in load_file(..).
  Tag obb_tag() {return obbTag;}

  Tag geom_tag() {return geomTag;}

  Tag id_tag() {return idTag;}

  Tag name_tag() {return nameTag;}
  
  Tag sense_tag() { return senseTag; }

  double distance_limit() {return distanceLimit;}
  void distance_limit(double this_limit) {distanceLimit = this_limit;}
  
  std::vector<EntityHandle>& surf_handles() {return entHandles[2];}
      
  std::vector<EntityHandle>& vol_handles() {return entHandles[3];}
      
  std::vector<EntityHandle>& group_handles() {return entHandles[4];}
      
  ErrorCode write_mcnp(std::string ifile, const bool overwrite = true);
  ErrorCode parse_metadata();
  void tokenize( const std::string& str,
			std::vector<std::string>& tokens,
			const char* delimiters );
  
  ErrorCode poly_solid_angle( EntityHandle face, const CartVect& point, double& area );

  ErrorCode get_angle(EntityHandle surf, 
                        double xxx, double yyy, double zzz, double *ang);


    // get the corners of the OBB for a given volume
  ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]);

    // get the center point and three vectors for the OBB of a given volume
  ErrorCode getobb(EntityHandle volume, double center[3], 
                     double axis1[3], double axis2[3], double axis3[3]);

    // get the root of the obbtree for a given entity
  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle &root);

  int get_entity_id(EntityHandle this_ent);

    // Get the instance of MOAB used by functions in this file.
  Interface* moab_instance() {return mbImpl;}

  double add_dist_tol() {return addDistTol;}

  double discard_dist_tol() {return discardDistTol;}

  double faceting_tolerance() {return facetingTolerance;}

  int source_cell() {return sourceCell;}

  bool use_dist_limit() {return useDistLimit;}

  bool use_cad() {return useCAD;}

  ErrorCode CAD_ray_intersect(const double *point, 
                                const double *dir, 
                                const double huge_val,
                                std::vector<double> &distances,
                                std::vector<EntityHandle> &surfaces, 
                                double &len);

  ErrorCode boundary_case( EntityHandle volume, int& result, 
                             double u, double v, double w,
                             EntityHandle facet,
                             EntityHandle surface);
  
    /** Get subversion revision of this file (DagMC.hpp) */
  static unsigned int interface_revision();
  
private:
  bool have_obb_tree();

  ErrorCode get_impl_compl();
  
  Tag get_tag( const char* name, int size, TagType store, DataType type,
                 const void* def_value = NULL, bool create_if_missing = true);

  bool get_group_names(EntityHandle group_set, std::vector<std::string> &grp_names);
  
  static void create_instance(Interface *mb_impl = NULL);
  
    // build internal index vectors that speed up handle-by-id, etc.
  ErrorCode build_indices(Range &surfs, Range &vols);
  
    // build obb structure
  ErrorCode build_obbs(Range &surfs, Range &vols);
  ErrorCode build_obb_impl_compl(Range &surfs);
  
  // attempt to set useCAD, first checking for availability
  void set_useCAD( bool use_cad ); 

  class Option {
  public:
    Option(){}
    Option( const char* n, const char* d, const char* v )
        : name(n), desc(d), value(v), user_set(false) {}
    std::string name, desc, value;
    bool user_set;
  };

  std::vector<int> tallyList;

  std::string itos(int ival);

  Interface *mbImpl;

  OrientedBoxTreeTool obbTree;
  Tag obbTag, geomTag, idTag, nameTag, senseTag;
  // metadata
  Tag matTag, densTag, bcTag, impTag, tallyTag;
  
  Option options[6];

  char specReflectName[NAME_TAG_SIZE];
  char whiteReflectName[NAME_TAG_SIZE];
  char implComplName[NAME_TAG_SIZE];

  std::vector<EntityHandle> entHandles[5];

  EntityHandle impl_compl_handle;

  double discardDistTol;
  double addDistTol;
  double facetingTolerance;
  int sourceCell;
  bool useDistLimit;
  bool useCAD;         /// true if user requested CAD-based ray firing
  bool have_cgm_geom;  /// true if CGM contains problem geometry; required for CAD-based ray firing.

  double distanceLimit;

    // store some lists indexed by handle
    // this is the lowest-valued handle among entity sets representing
    // surfaces & volumes
  EntityHandle setOffset;

    // list of obbTree root sets for surfaces and volumes, 
    // indexed by [surf_or_vol_handle - setOffset]
  std::vector<EntityHandle> rootSets;

    // entity index (contiguous 1-N indices) indexed like rootSets are
  std::vector<int> entIndices;

#ifdef CGM
    // corresponding geometric entities indexed like rootSets are
  std::vector<RefEntity *> geomEntities;

#endif
  
  DagMC(Interface *mb_impl);
  
  static DagMC *instance_;
  
    // temporary storage so functions don't have to reallocate vectors
  std::vector<EntityHandle> triList, surfList, facetList;
  std::vector<double> distList;
};

inline float DagMC::version(std::string *version_string) 
{
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

inline char *DagMC::get_spec_reflect() 
{
  return specReflectName;
}

inline char *DagMC::get_white_reflect() 
{
  return whiteReflectName;
}

inline DagMC *DagMC::instance(Interface *mb_impl) 
{
  if (NULL == instance_) create_instance(mb_impl);
  
  return instance_;
}

inline int DagMC::index_by_handle( EntityHandle handle )
{
  assert(handle-setOffset < entIndices.size());
  return entIndices[handle-setOffset];
}

inline int DagMC::num_entities( int dimension )
{
  assert(0 <= dimension && 3 >= dimension);
  
  return entHandles[dimension].size() - 1;
}

inline EntityHandle DagMC::entity_by_index( int dimension, int index )
{
  assert(2 <= dimension && 3 >= dimension && (unsigned) index < entHandles[dimension].size());
  return entHandles[dimension][index];
}

    // get the root of the obbtree for a given entity
inline ErrorCode DagMC::get_root(EntityHandle vol_or_surf, EntityHandle &root) 
{
  unsigned int index = vol_or_surf - setOffset;
  root = (index < rootSets.size() ? rootSets[index] : 0);
  return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

} // namespace moab

#endif

