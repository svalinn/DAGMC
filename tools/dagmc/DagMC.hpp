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

class DagMC 
{
public:
  static DagMC *instance(MBInterface *mb_impl = NULL);
  
  ~DagMC();
  
  MBErrorCode ray_fire(const MBEntityHandle cell, const MBEntityHandle last_surf_hit, 
                       const int num_pts,
                       const double uuu, const double vvv, const double www,
                       const double xxx, const double yyy, const double zzz,
                       const double huge,
                       double &dist_traveled, MBEntityHandle &next_surf_hit);

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
                               int& result );

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

    // load mesh and initialize
  MBErrorCode load_file_and_init(const char* cfile,
                                 const int clen,
                                 const char* ffile,
                                 const int flen);

    // map between MBEntityHandle, base-1 index, and GLOBAL_ID
  MBEntityHandle entity_by_index( int dimension, int index );
  MBEntityHandle entity_by_id( int dimension, int id );
  int index_by_handle( MBEntityHandle handle );
  int id_by_index( int dimension, int index );

    // get number of geometric sets corresponding to geometry
    // of specified dimension (i.e. number of geometric surfaces
    // if dimension == 2).
  int num_entities( int dimension );

    // read/write settings 
  void read_settings( const char* filename );
  void write_settings( FILE* filp, bool with_description = true );
  void parse_settings();

  char *get_spec_reflect();
  char *get_white_reflect();

  MBOrientedBoxTreeTool *obb_tree() {return &obbTree;}

    // Get the tag used to associate OBB trees with geometry in load_file(..).
  MBTag obb_tag() {return obbTag;}

  MBTag geom_tag() {return geomTag;}

  MBTag id_tag() {return idTag;}

  MBTag name_tag() {return nameTag;}

  double distance_limit() {return distanceLimit;}
  void distance_limit(double this_limit) {distanceLimit = this_limit;}
  
  std::vector<MBEntityHandle>& surf_handles() {return entHandles[2];}
      
  std::vector<MBEntityHandle>& vol_handles() {return entHandles[3];}
      
  std::vector<MBEntityHandle>& group_handles() {return entHandles[4];}
      
  void write_log(std::string ifile, const bool overwrite = true);
  
  MBErrorCode poly_solid_angle( MBEntityHandle face, const MBCartVect& point, double& area );

  MBErrorCode get_angle(MBEntityHandle surf, 
                        double xxx, double yyy, double zzz, double *ang);
  
  int get_entity_id(MBEntityHandle this_ent);

    // Get the instance of MOAB used by functions in this file.
  MBInterface* moab_instance() {return mbImpl;}

  double tolerance() {return moabMCNPTolerance;}
  
  int source_cell() {return moabMCNPSourceCell;}

  bool use_dist_limit() {return moabMCNPUseDistLimit;}

  bool use_cad() {return useCAD;}

  MBErrorCode CAD_ray_intersect(const double *point, 
                                const double *dir, 
                                const double huge,
                                std::vector<double> &distances,
                                std::vector<MBEntityHandle> &surfaces, 
                                double &len);

private:
  MBTag get_tag( const char* name, int size, MBTagType store, MBDataType type,
                 bool create_if_missing = true);

  bool get_group_names(MBEntityHandle group_set, std::vector<std::string> &grp_names);
  
  static void create_instance(MBInterface *mb_impl = NULL);
  
    // build internal index vectors that speed up handle-by-id, etc.
  MBErrorCode build_indices(MBRange &surfs, MBRange &vols,
                            bool is_geom);
  
    // build obb structure
  MBErrorCode build_obbs(MBRange &surfs, MBRange &vols);
  
  class Option {
  public:
    Option(){}
    Option( const char* n, const char* d, const char* v )
        : name(n), desc(d), value(v) {}
    std::string name, desc, value;
  };

  std::string itos(int ival);

  MBInterface *mbImpl;

  MBOrientedBoxTreeTool obbTree;
  MBTag obbTag, geomTag, idTag, nameTag;
  
  Option options[4];

  char specReflectName[NAME_TAG_SIZE];
  char whiteReflectName[NAME_TAG_SIZE];

  std::vector<MBEntityHandle> entHandles[5];

  double moabMCNPTolerance;
  int moabMCNPSourceCell;
  bool moabMCNPUseDistLimit;
  bool useCAD;

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
};

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

#endif

