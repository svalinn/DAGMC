#ifndef MOABMC_HPP
#define MOABMC_HPP

#include "MBInterface.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <string>
#include <assert.h>

#include "MBOrientedBoxTreeTool.hpp"

class MBCartVect;

#ifdef CGM
class RefEntity;
#endif

#define DAGMC_VERSION 0.99
#define DAGMC_VERSION_STRING "0.99"
#define DAGMC_INTERFACE_REVISION "$Rev$"

class DagMC 
{
public:
  static DagMC *instance(MBInterface *mb_impl = NULL);
  
  ~DagMC();
  
    // Return the version of this library
  float version(std::string *version_string = NULL);

  MBErrorCode ray_fire(const MBEntityHandle cell, const MBEntityHandle last_surf_hit, 
                       const int num_pts,
                       const double uuu, const double vvv, const double www,
                       const double xxx, const double yyy, const double zzz,
                       const double huge_val,
                       double &dist_traveled, MBEntityHandle &next_surf_hit,
                       MBOrientedBoxTreeTool::TrvStats* stats = NULL );

    // Test if point is inside or outside of a volume using unit sphere area method
    // Requires sense of surfaces wrt volume.
    // Return values: {0 : outside, 1: inside, -1: on boundary}
  MBErrorCode point_in_volume_slow( MBEntityHandle volume, 
                                    double x, double y, double z,
                                    int& result );

    // Test if point is inside or outside of a volume using an OBB tree
    // Requires sense of surfaces wrt volume.
    // Return values: {0 : outside, 1: inside, -1: on boundary}
  MBErrorCode point_in_volume( MBEntityHandle volume, 
                               double x, double y, double z,
                               int& result,
                               double u, double v, double w);

    // Find the distance to the nearest surface
  MBErrorCode closest_to_location( MBEntityHandle volume,
                                   double* point,
                                   double& result);

    // Calculate the volume contained in a 'volume'
  MBErrorCode measure_volume( MBEntityHandle volume, double& result );

    // Calculate sum of area of triangles
  MBErrorCode measure_area( MBEntityHandle surface, double& result );

    // Get the sense of surfaces wrt a volume.  Sense values are:
    // {-1 -> reversed, 0 -> both, 1 -> forward}
  MBErrorCode surface_sense( MBEntityHandle volume, 
                             int num_surfaces,
                             const MBEntityHandle* surfaces,
                             int* senses_out );

    // Get the sense of surfaces wrt a volume.  Sense values are:
    // {-1 -> reversed, 0 -> both, 1 -> forward}
  MBErrorCode surface_sense( MBEntityHandle volume, MBEntityHandle surface, int& sense_out );

    // load mesh
  MBErrorCode load_file(const char* cfile,
			const double facet_tolerance = 0);

  MBErrorCode write_mesh(const char* ffile,
			 const int flen);


    // initialize data structures and OBB tree
  MBErrorCode init_OBBTree();

    // map between MBEntityHandle, base-1 index, and GLOBAL_ID
  MBEntityHandle entity_by_index( int dimension, int index );
  MBEntityHandle entity_by_id( int dimension, int id );
  int index_by_handle( MBEntityHandle handle );
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

  char *get_spec_reflect();
  char *get_white_reflect();

  MBOrientedBoxTreeTool *obb_tree() {return &obbTree;}

    // Get the tag used to associate OBB trees with geometry in load_file(..).
  MBTag obb_tag() {return obbTag;}

  MBTag geom_tag() {return geomTag;}

  MBTag id_tag() {return idTag;}

  MBTag name_tag() {return nameTag;}
  
  MBTag sense_tag() { return senseTag; }

  double distance_limit() {return distanceLimit;}
  void distance_limit(double this_limit) {distanceLimit = this_limit;}
  
  std::vector<MBEntityHandle>& surf_handles() {return entHandles[2];}
      
  std::vector<MBEntityHandle>& vol_handles() {return entHandles[3];}
      
  std::vector<MBEntityHandle>& group_handles() {return entHandles[4];}
      
  MBErrorCode write_mcnp(std::string ifile, const bool overwrite = true);
  MBErrorCode parse_metadata();
  void tokenize( const std::string& str,
			std::vector<std::string>& tokens,
			const char* delimiters );
  
  MBErrorCode poly_solid_angle( MBEntityHandle face, const MBCartVect& point, double& area );

  MBErrorCode get_angle(MBEntityHandle surf, 
                        double xxx, double yyy, double zzz, double *ang);


    // get the corners of the OBB for a given volume
  MBErrorCode getobb(MBEntityHandle volume, double minPt[3], double maxPt[3]);

    // get the center point and three vectors for the OBB of a given volume
  MBErrorCode getobb(MBEntityHandle volume, double center[3], 
                     double axis1[3], double axis2[3], double axis3[3]);

    // get the root of the obbtree for a given entity
  MBErrorCode get_root(MBEntityHandle vol_or_surf, MBEntityHandle &root);

  int get_entity_id(MBEntityHandle this_ent);

    // Get the instance of MOAB used by functions in this file.
  MBInterface* moab_instance() {return mbImpl;}

  double add_dist_tol() {return addDistTol;}

  double discard_dist_tol() {return discardDistTol;}

  double faceting_tolerance() {return facetingTolerance;}

  int source_cell() {return sourceCell;}

  bool use_dist_limit() {return useDistLimit;}

  bool use_cad() {return useCAD;}

  MBErrorCode CAD_ray_intersect(const double *point, 
                                const double *dir, 
                                const double huge_val,
                                std::vector<double> &distances,
                                std::vector<MBEntityHandle> &surfaces, 
                                double &len);

  MBErrorCode boundary_case( MBEntityHandle volume, int& result, 
                             double u, double v, double w,
                             MBEntityHandle facet,
                             MBEntityHandle surface);
  
    /** Get subversion revision of this file (DagMC.hpp) */
  static unsigned int interface_revision();
  
private:
  bool have_obb_tree();

  MBErrorCode get_impl_compl();
  
  MBTag get_tag( const char* name, int size, MBTagType store, MBDataType type,
                 const void* def_value = NULL, bool create_if_missing = true);

  bool get_group_names(MBEntityHandle group_set, std::vector<std::string> &grp_names);
  
  static void create_instance(MBInterface *mb_impl = NULL);
  
    // build internal index vectors that speed up handle-by-id, etc.
  MBErrorCode build_indices(MBRange &surfs, MBRange &vols);
  
    // build obb structure
  MBErrorCode build_obbs(MBRange &surfs, MBRange &vols);
  MBErrorCode build_obb_impl_compl(MBRange &surfs);
  
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

  MBInterface *mbImpl;

  MBOrientedBoxTreeTool obbTree;
  MBTag obbTag, geomTag, idTag, nameTag, senseTag;
  // metadata
  MBTag matTag, densTag, bcTag, impTag, tallyTag;
  
  Option options[6];

  char specReflectName[NAME_TAG_SIZE];
  char whiteReflectName[NAME_TAG_SIZE];
  char implComplName[NAME_TAG_SIZE];

  std::vector<MBEntityHandle> entHandles[5];

  MBEntityHandle impl_compl_handle;

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
  MBEntityHandle setOffset;

    // list of obbTree root sets for surfaces and volumes, 
    // indexed by [surf_or_vol_handle - setOffset]
  std::vector<MBEntityHandle> rootSets;

    // entity index (contiguous 1-N indices) indexed like rootSets are
  std::vector<int> entIndices;

#ifdef CGM
    // corresponding geometric entities indexed like rootSets are
  std::vector<RefEntity *> geomEntities;

#endif
  
  DagMC(MBInterface *mb_impl);
  
  static DagMC *instance_;
  
    // temporary storage so functions don't have to reallocate vectors
  std::vector<MBEntityHandle> triList, surfList, facetList;
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

inline DagMC *DagMC::instance(MBInterface *mb_impl) 
{
  if (NULL == instance_) create_instance(mb_impl);
  
  return instance_;
}

inline int DagMC::index_by_handle( MBEntityHandle handle )
{
  assert(handle-setOffset < entIndices.size());
  return entIndices[handle-setOffset];
}

inline int DagMC::num_entities( int dimension )
{
  assert(0 <= dimension && 3 >= dimension);
  
  return entHandles[dimension].size() - 1;
}

inline MBEntityHandle DagMC::entity_by_index( int dimension, int index )
{
  assert(2 <= dimension && 3 >= dimension && (unsigned) index < entHandles[dimension].size());
  return entHandles[dimension][index];
}

    // get the root of the obbtree for a given entity
inline MBErrorCode DagMC::get_root(MBEntityHandle vol_or_surf, MBEntityHandle &root) 
{
  unsigned int index = vol_or_surf - setOffset;
  root = (index < rootSets.size() ? rootSets[index] : 0);
  return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

#endif

