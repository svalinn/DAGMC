#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "moab/GeomQueryTool.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <vector>
#include <iostream>
#include <math.h>
#include <limits>
#include <algorithm>
#include <stdio.h> // for remove()

#define CHKERR if (MB_SUCCESS != rval) return rval

const double ROOT2 = 1.4142135623730951;

using namespace moab;

// Create file containing geometry for 2x2x2 cube 
// centered at origin and having a convex +Z face
// (center of face at origin).
ErrorCode write_geometry( const char* output_file_name );

ErrorCode test_ray_fire( DagMC * );

ErrorCode test_point_in_volume( DagMC * );

ErrorCode test_measure_volume( DagMC * );

ErrorCode test_measure_area( DagMC * );

ErrorCode test_surface_sense( DagMC * );

ErrorCode overlap_write_geometry( const char* output_file_name );
ErrorCode overlap_test_ray_fire( DagMC * );
ErrorCode overlap_test_point_in_volume( DagMC * );
ErrorCode overlap_test_measure_volume( DagMC * );
ErrorCode overlap_test_measure_area( DagMC * );
ErrorCode overlap_test_surface_sense( DagMC * );
ErrorCode overlap_test_tracking( DagMC * );

ErrorCode write_geometry( const char* output_file_name )
{
  ErrorCode rval;

  Interface *moab = new Core();
  
  // Define a 2x2x2 cube centered at orgin
  // with concavity in +Z face.
  const double coords[] = {
    1, -1, -1, 
    1,  1, -1,
   -1,  1, -1,
   -1, -1, -1,
    1, -1,  1, 
    1,  1,  1,
   -1,  1,  1,
   -1, -1,  1,
    0,  0,  0 };
  const int connectivity[] = {
    0, 3, 1,  3, 2, 1, // -Z
    0, 1, 4,  5, 4, 1, // +X
    1, 2, 6,  6, 5, 1, // +Y
    6, 2, 3,  7, 6, 3, // -X
    0, 4, 3,  7, 3, 4, // -Y
    4, 5, 8,  5, 6, 8, // +Z
    6, 7, 8,  7, 4, 8  // +Z
  };
  const unsigned tris_per_surf[] = { 2, 2, 2, 2, 2, 4 };
  
    // Create the geometry
  const unsigned num_verts = sizeof(coords) / (3*sizeof(double));
  const unsigned num_tris = sizeof(connectivity) / (3*sizeof(int));
  const unsigned num_surfs = sizeof(tris_per_surf) / sizeof(unsigned);
  EntityHandle verts[num_verts], tris[num_tris], surfs[num_surfs];
  for (unsigned i = 0; i < num_verts; ++i) {
    rval = moab->create_vertex( coords + 3*i, verts[i] ); 
    CHKERR;
  }
  for (unsigned i = 0; i < num_tris; ++i) {
    const EntityHandle conn[] = { verts[connectivity[3*i  ]], 
                                    verts[connectivity[3*i+1]], 
                                    verts[connectivity[3*i+2]] };
    rval = moab->create_element( MBTRI, conn, 3, tris[i] );
    CHKERR;
  }
  
    // create CAD topology
  EntityHandle* tri_iter = tris;
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = moab->create_meshset( MESHSET_SET, surfs[i] );
    CHKERR;
    rval = moab->add_entities( surfs[i], tri_iter, tris_per_surf[i] );
    CHKERR;
    tri_iter += tris_per_surf[i];
  }
  
  Tag dim_tag, id_tag, sense_tag;
  rval = moab->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              dim_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  CHKERR;
  rval = moab->tag_get_handle( GLOBAL_ID_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              id_tag,
                              MB_TAG_DENSE|MB_TAG_CREAT );
  CHKERR;
  rval = moab->tag_get_handle( "GEOM_SENSE_2", 
                              2, MB_TYPE_HANDLE, 
                              sense_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  CHKERR;

  std::vector<int> dims( num_surfs, 2 );
  rval = moab->tag_set_data( dim_tag, surfs, num_surfs, &dims[0] );
  CHKERR;
  std::vector<int> ids( num_surfs );
  for (size_t i = 0; i < ids.size(); ++i) ids[i] = i+1;
  rval = moab->tag_set_data( id_tag, surfs, num_surfs, &ids[0] );
  CHKERR;

  EntityHandle volume;
  rval = moab->create_meshset( MESHSET_SET, volume );
  CHKERR;
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = moab->add_parent_child( volume, surfs[i] );
    CHKERR;
  }
  
  std::vector<EntityHandle> senses( 2*num_surfs, 0 );
  for (size_t i = 0; i < senses.size(); i += 2)
    senses[i] = volume;
  rval = moab->tag_set_data( sense_tag, surfs, num_surfs, &senses[0] );
  CHKERR;
  
  const int three = 3;
  const int one = 1;
  rval = moab->tag_set_data( dim_tag, &volume, 1, &three );
  CHKERR;
  rval = moab->tag_set_data( id_tag, &volume, 1, &one );
  CHKERR;
  
  rval = moab->write_mesh( output_file_name );
  CHKERR;
  delete moab;
  return MB_SUCCESS;
}

ErrorCode overlap_write_geometry( const char* output_file_name )
{
  ErrorCode rval;
  Interface *moab = new Core();
  
  // Define two 1x2x2 cubes that overlap from 0 <= x <= 0.01
  // cube 0 centered at (0.5,0,0)
  const double coords[] = {
    1, -1, -1, 
    1,  1, -1,
    0,  1, -1,
    0, -1, -1,
    1, -1,  1, 
    1,  1,  1,
    0,  1,  1,
    0, -1,  1,
  // cube 1 centered near (-0.5,0,0)
    0.01, -1, -1, 
    0.01,  1, -1,
   -1,     1, -1,
   -1,    -1, -1,
    0.01, -1,  1, 
    0.01,  1,  1,
   -1,     1,  1,
   -1,    -1,  1 };
  const int connectivity[] = {
    0, 3, 1,  3, 2, 1, // -Z
    0, 1, 4,  5, 4, 1, // +X
    1, 2, 6,  6, 5, 1, // +Y
    6, 2, 3,  7, 6, 3, // -X
    0, 4, 3,  7, 3, 4, // -Y
    4, 5, 6,  6, 7, 4};// +Z
  
  // Create the geometry
  const unsigned tris_per_surf = 2;
  const unsigned num_cubes     = 2;
  const unsigned num_verts = sizeof(coords) / (3*sizeof(double));
  const unsigned num_tris  = num_cubes*sizeof(connectivity) / (3*sizeof(int));
  const unsigned num_surfs = num_tris / tris_per_surf;
  EntityHandle verts[num_verts], tris[num_tris], surfs[num_surfs];

  for (unsigned i = 0; i < num_verts; ++i) {
    rval = moab->create_vertex( coords + 3*i, verts[i] ); 
    CHKERR;
  }
  // cube0
  for (unsigned i = 0; i < num_tris/2; ++i) {
    const EntityHandle conn[] = { verts[connectivity[3*i  ]], 
                                  verts[connectivity[3*i+1]], 
                                  verts[connectivity[3*i+2]] };
    rval = moab->create_element( MBTRI, conn, 3, tris[i] );
    CHKERR;
  }
  // cube1
  for (unsigned i = 0; i < num_tris/2; ++i) {
    const EntityHandle conn[] = { verts[8 + connectivity[3*i  ]], 
                                  verts[8 + connectivity[3*i+1]], 
                                  verts[8 + connectivity[3*i+2]] };
    rval = moab->create_element( MBTRI, conn, 3, tris[num_tris/2 + i] );
    CHKERR;
  }
  
  // create CAD topology
  EntityHandle* tri_iter = tris;
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = moab->create_meshset( MESHSET_SET, surfs[i] );
    CHKERR;
    rval = moab->add_entities( surfs[i], tri_iter, tris_per_surf );
    CHKERR;
    tri_iter += tris_per_surf;
  }
  
  Tag dim_tag, id_tag, sense_tag;
  rval = moab->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              dim_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  CHKERR;
  rval = moab->tag_get_handle( GLOBAL_ID_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              id_tag,
                              MB_TAG_DENSE|MB_TAG_CREAT );
  CHKERR;
  rval = moab->tag_get_handle( "GEOM_SENSE_2", 
                              2, MB_TYPE_HANDLE, 
                              sense_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  CHKERR;

  std::vector<int> dims( num_surfs, 2 );
  rval = moab->tag_set_data( dim_tag, surfs, num_surfs, &dims[0] );
  CHKERR;
  std::vector<int> ids( num_surfs );
  for (size_t i = 0; i < ids.size(); ++i) ids[i] = i+1;
  rval = moab->tag_set_data( id_tag, surfs, num_surfs, &ids[0] );
  CHKERR;

  EntityHandle volume;
  rval = moab->create_meshset( MESHSET_SET, volume );
  CHKERR;
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = moab->add_parent_child( volume, surfs[i] );
    CHKERR;
  }
  
  std::vector<EntityHandle> senses( 2*num_surfs, 0 );
  for (size_t i = 0; i < senses.size(); i += 2)
    senses[i] = volume;
  rval = moab->tag_set_data( sense_tag, surfs, num_surfs, &senses[0] );
  CHKERR;
  
  const int three = 3;
  const int one   = 1;
  rval = moab->tag_set_data( dim_tag, &volume, 1, &three );
  CHKERR;
  rval = moab->tag_set_data( id_tag, &volume, 1, &one );
  CHKERR;
  
  rval = moab->write_mesh( output_file_name );
  CHKERR;
  delete moab;
  return MB_SUCCESS;
}

static bool run_test( std::string name, int argc, char* argv[] )
{
  if (argc == 1)
    return true;
  for (int i = 1; i < argc; ++i)
    if (name == argv[i])
      return true;
  return false;
}

#define RUN_TEST(A) do { \
  if (run_test( #A, argc, argv )) { \
    std::cout << #A << "... " << std::endl; \
    if (MB_SUCCESS != A ( dagmc ) ) { \
      ++errors; \
    } \
  } \
} while(false)

int main( int argc, char* argv[] )
{
  ErrorCode rval;
  const char* filename = "test_geom.h5m";
  
#ifdef MOAB_HAVE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

  rval = write_geometry( filename );
  if (MB_SUCCESS != rval) {
    remove( filename );
    std::cerr << "Failed to create input file: " << filename << std::endl;
    return 1;
  }
  
  DagMC *dagmc = new DagMC();

  int errors = 0;
  //rval = dagmc.moab_instance()->load_file( filename );
  rval = dagmc->load_file( filename );
  remove( filename );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to load file." << std::endl;
    return 2;
  }
  rval = dagmc->init_OBBTree();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC." << std::endl;
    return 2;
  }
  
  RUN_TEST( test_ray_fire );
  RUN_TEST( test_point_in_volume );
  RUN_TEST( test_measure_volume );
  RUN_TEST( test_measure_area );
  RUN_TEST( test_surface_sense );
 
  // change settings to use overlap-tolerant mode (arbitrary thickness)
  double overlap_thickness = 0.1;
  dagmc->set_overlap_thickness( overlap_thickness );
  RUN_TEST( test_ray_fire );
  RUN_TEST( test_point_in_volume );

  // clear moab and dagmc instance
  rval = dagmc->moab_instance()->delete_mesh();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to delete mesh." << std::endl;
    return 2;
  }

  // Now load a different geometry: two cubes that slightly overlap
  rval = overlap_write_geometry( filename );
  if (MB_SUCCESS != rval) {
    remove( filename );
    std::cerr << "Failed to create input file: " << filename << std::endl;
    return 1;
  }
  
  delete dagmc;
  
  dagmc = new DagMC();
  //  rval = dagmc->moab_instance()->load_file( filename );
  rval = dagmc->load_file( filename );
  remove( filename );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to load file with overlaps." << std::endl;
    return 2;
  }
  rval = dagmc->init_OBBTree();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC with overlaps." << std::endl;
    return 2;
  }
  // change settings to use overlap-tolerant mode (with a large enough thickness)
  overlap_thickness = 3;
  dagmc->set_overlap_thickness( overlap_thickness );
  RUN_TEST( overlap_test_ray_fire );
  RUN_TEST( overlap_test_point_in_volume );
  RUN_TEST( overlap_test_measure_volume );
  RUN_TEST( overlap_test_measure_area );
  RUN_TEST( overlap_test_surface_sense );
  RUN_TEST( overlap_test_tracking );

#ifdef MOAB_HAVE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  delete dagmc;
  
  return errors;
}

ErrorCode test_surface_sense( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs, vols;
  const int two = 2, three = 3;
  const void* ptrs[] = { &two, &three };
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, ptrs, 1, surfs );
  CHKERR;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, ptrs+1, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() << std::endl;
    return MB_FAILURE;
  }
  if (surfs.size() != 6) {
    std::cerr << "ERROR: Expected 6 surfaces in input, found " << surfs.size() << std::endl;
    return MB_FAILURE;
  }
  
  for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    int sense = 0;
    rval = dagmc->surface_sense( vols.front(), 1, &*i, &sense );
    if (MB_SUCCESS != rval || sense != 1) {
      std::cerr << "ERROR: Expected 1 for surface sense, got " << sense << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}  

ErrorCode overlap_test_surface_sense( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs, vols;
  const int two = 2, three = 3;
  const void* ptrs[] = { &two, &three };
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs, 1, surfs );
  CHKERR;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs+1, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() 
              << std::endl;
    return MB_FAILURE;
  }
  if (surfs.size() != 12) {
    std::cerr << "ERROR: Expected 12 surfaces in input, found " << surfs.size() 
              << std::endl;
    return MB_FAILURE;
  }
  
  for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    int sense = 0;
    rval = dagmc->surface_sense( vols.front(), 1, &*i, &sense );
    if (MB_SUCCESS != rval || sense != 1) {
      std::cerr << "ERROR: Expected 1 for surface sense, got " << sense << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}  
ErrorCode test_measure_volume( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();
  
  Tag dim_tag = dagmc->geom_tag();
  Range vols;
  const int three = 3;
  const void* ptr = &three;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, &ptr, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() << std::endl;
    return MB_FAILURE;
  }
  
    // input file is 2x2x2 cube with convexity in +Z face that touches the origin.
    // expected volume is 8 (2x2x2) less the volume of the pyrimid concavity
  double result;
  const double vol = 2*2*2 - 1*4./3;

  rval = dagmc->measure_volume( vols.front(), result );
  CHKERR;
  if (fabs(result - vol) > 10*std::numeric_limits<double>::epsilon()) {
    std::cerr << "ERROR: Expected " << vol << " as measure of volume, got " << result << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}
ErrorCode overlap_test_measure_volume( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();
  
  Tag dim_tag = dagmc->geom_tag();
  Range vols;
  const int three = 3;
  const void* ptr = &three;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            &ptr, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() 
              << std::endl;
    return MB_FAILURE;
  }
  
  // The volume has two regions that overlap.
  double result;
  const double vol = (1+1.01)*2*2;

  rval = dagmc->measure_volume( vols.front(), result );
  CHKERR;
  if (fabs(result - vol) > 2*std::numeric_limits<double>::epsilon()) {
    std::cerr << "ERROR: Expected " << vol << " as measure of volume, got " 
              << result << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

ErrorCode test_measure_area( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs;
  const int two = 2;
  const void* ptr = &two;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, &ptr, 1, surfs );
  CHKERR;

  if (surfs.size() != 6) {
    std::cerr << "ERROR: Expected 6 surfaces in input, found " << surfs.size() 
              << std::endl;
    return MB_FAILURE;
  }
  
  int ids[6];
  rval = moab->tag_get_data( dagmc->id_tag(), surfs, ids );
  CHKERR;
  
    // expect area of 4 for all faces except face 6.
    // face 6 should have area == 4*sqrt(2)
  Range::iterator iter = surfs.begin();
  for (unsigned i = 0; i < 6; ++i, ++iter) {
    double expected = 4.0;
    if (ids[i] == 6)
      expected *= ROOT2;
    
    double result;
    
    rval = dagmc->measure_area( *iter, result );
    CHKERR;
    if (fabs(result - expected) > std::numeric_limits<double>::epsilon()) {
      std::cerr << "ERROR: Expected area of surface " << ids[i] << " to be " 
                << expected << ".  Got " << result << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode overlap_test_measure_area( DagMC * dagmc )
{
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs;
  const int two = 2;
  const void* ptr = &two;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            &ptr, 1, surfs );
  CHKERR;

  const unsigned num_surfs = 12;
  if (surfs.size() != num_surfs) {
    std::cerr << "ERROR: Expected " << num_surfs << " surfaces in input, found " 
              << surfs.size() << std::endl;
    return MB_FAILURE;
  }
  
  int ids[num_surfs];
  rval = moab->tag_get_data( dagmc->id_tag(), surfs, ids );
  CHKERR;
  
  const double x_area   = 2*2;
  const double yz_area0 = 2*1;
  const double yz_area1 = 2*1.01;
  Range::iterator iter = surfs.begin();
  for (unsigned i = 0; i < num_surfs; ++i, ++iter) {
    double expected, result;
    if      (0==i || 2==i || 4 ==i || 5 ==i) expected = yz_area0;
    else if (1==i || 3==i || 7 ==i || 9 ==i) expected = x_area;
    else if (6==i || 8==i || 10==i || 11==i) expected = yz_area1;
    
    rval = dagmc->measure_area( *iter, result );
    CHKERR;
    if (fabs(result - expected) > std::numeric_limits<double>::epsilon()) {
      std::cerr << "ERROR: Expected area of surface " << ids[i] << " to be " 
                << expected << ".  Got " << result << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

struct ray_fire {
  int prev_surf;
  double origin[3], direction[3];
  int hit_surf;
  double distance;
};

ErrorCode test_ray_fire( DagMC * dagmc )
{
  // Glancing ray-triangle intersections are not valid exit intersections. 
  // Piercing ray-triangle intersections are valid exit intersections.
  // "0" destination surface implies that it is ambiguous.
  const struct ray_fire tests[] = {
  /* src_srf origin               direction                 dest dist */
    // piercing edge
    { 1, { 0.0, 0.0, -1. }, { -1.0/ROOT2, 0.0, 1.0/ROOT2 }, 4, ROOT2 },
    // piercing edge
    { 1, { 0.0, 0.0, -1. }, {  1.0/ROOT2, 0.0, 1.0/ROOT2 }, 2, ROOT2 },
    // piercing edge
    { 1, { 0.0, 0.0, -1. }, {  0.0, 1.0/ROOT2, 1.0/ROOT2 }, 3, ROOT2 },
    // piercing edge
    { 1, { 0.5, 0.5, -1. }, {  0.0, 0.0, 1.0 },             6, 1.5   },
    // interior
    { 2, { 1.0, 0.0, 0.5 }, { -1.0, 0.0, 0.0 },             6, 0.5   },
    // glancing node then piercing edge
    { 2, { 1.0, 0.0, 0.0 }, { -1.0, 0.0, 0.0 },             4, 2.0   },
    // piercing node
    { 1, { 0.0, 0.0, -1. }, {  0.0, 0.0, 1.0 },             6, 1.0   },
    // glancing edge then interior
    { 2, { 1.0, 0.0, 0.5 }, { -1.0/ROOT2, 1.0/ROOT2, 0.0 }, 3, ROOT2 } };

  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs, vols;
  const int two = 2, three = 3;
  const void* ptrs[] = { &two, &three };
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, ptrs, 1, surfs );
  CHKERR;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, ptrs+1, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() << std::endl;
    return MB_FAILURE;
  }
  if (surfs.size() != 6) {
    std::cerr << "ERROR: Expected 6 surfaces in input, found " << surfs.size() << std::endl;
    return MB_FAILURE;
  }
  
  int ids[6];
  rval = moab->tag_get_data( dagmc->id_tag(), surfs, ids );
  CHKERR;
  EntityHandle surf[6];
  std::copy( surfs.begin(), surfs.end(), surf );
  
  const int num_test = sizeof(tests) / sizeof(tests[0]);
  for (int i = 0; i < num_test; ++i) {
    int* ptr = std::find( ids, ids+6, tests[i].prev_surf );
    int idx = ptr - ids;
    if (idx >= 6) {
      std::cerr << "Surface " << tests[i].prev_surf << " not found." << std::endl;
      return MB_FAILURE;
    }
    //const EntityHandle src_surf = surf[idx];
    
    ptr = std::find( ids, ids+6, tests[i].hit_surf );
    idx = ptr - ids;
    if (idx >= 6) {
      std::cerr << "Surface " << tests[i].hit_surf << " not found." << std::endl;
      return MB_FAILURE;
    }
    const EntityHandle hit_surf = surf[idx];

    double dist;
    EntityHandle result;
    GeomQueryTool::RayHistory history;
    rval = dagmc->ray_fire( vols.front(), 
                           tests[i].origin, tests[i].direction,
                           result, dist, &history );
    
    if (result != hit_surf || fabs(dist - tests[i].distance) > 1e-6) {
      EntityHandle *p = std::find( surf, surf+6, result );
      idx = p - surf;
      int id = idx > 5 ? 0 : ids[idx];
      
      std::cerr << "Rayfire test failed for " << std::endl
                << "\t ray from (" << tests[i].origin[0] 
                << ", " << tests[i].origin[1] << ", "
                << tests[i].origin[2] << ") going ["
                << tests[i].direction[0] << ", " 
                << tests[i].direction[1] << ", "
                << tests[i].direction[2] << "]" << std::endl
                << "\t Beginning on surface " << tests[i].prev_surf << std::endl
                << "\t Expected to hit surface " << tests[i].hit_surf << " after " 
                << tests[i].distance << " units." << std::endl
                << "\t Actually hit surface " << id << " after " << dist << " units."
                << std::endl;
      return MB_FAILURE;
    }

    CartVect loc = CartVect(tests[i].origin) + (dist * CartVect(tests[i].direction));
    
    std::vector< std::pair<int,GeomQueryTool::RayHistory*> > boundary_tests;
    boundary_tests.push_back( std::make_pair( 1, &history ) );
    boundary_tests.push_back( std::make_pair( 0, &history ) );
    boundary_tests.push_back( std::make_pair( 1, (GeomQueryTool::RayHistory*)NULL ) );
    boundary_tests.push_back( std::make_pair( 0, (GeomQueryTool::RayHistory*)NULL ) );


    for( unsigned int bt = 0; bt < boundary_tests.size(); ++bt ) {
      
      int expected = boundary_tests[bt].first;
      GeomQueryTool::RayHistory* h = boundary_tests[bt].second;
      
      // pick the direction based on expected result of test. Either reuse the ray_fire
      // vector, or reverse it to check for a vector that enters the cell
      CartVect uvw( tests[i].direction );
      if( expected == 1 )
        uvw = -uvw; 

      int boundary_result = -1;
      
      rval = dagmc->test_volume_boundary( vols.front(), result, loc.array(), 
                                         uvw.array(), boundary_result, h );
      
      
      if( boundary_result != expected ){
        std::cerr << "DagMC::test_volume_boundary failed (" << ( (expected==0)?"+":"-" )
                  << " dir," << ( (h)?"+":"-" ) << " history, i=" << i << ")" <<  std::endl;
        return MB_FAILURE;
      }
      
    }
    
 
  }
  
  return MB_SUCCESS;
}

ErrorCode overlap_test_ray_fire( DagMC * dagmc )
{
  // Glancing ray-triangle intersections are not valid exit intersections. 
  // Piercing ray-triangle intersections are valid exit intersections.
  // "0" destination surface implies that it is ambiguous.
  const struct ray_fire tests[] = {
    /*          ____________________
                |                   |
       -x <---  | region1 | overlap | region0 |   ---> +x
                          |                   |
                          ---------------------
       -x <--- -1.0      0.0      0.01        1.0 ---> +x
     surf_id:   10        4         8         2

     The zero-distance advance would occur in the implicit volume between
     these two regions---not in this volume.

     src_srf origin               direction                 dest dist */
    // numerical location on surface 1 of region0
    { 4, { 1.0,  0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 2,  0.0   },
    { 2, { 1.0,  0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 4,  1.0   },
    // numerical location inside region0
    { 4, { 0.5,  0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 2,  0.5   },
    { 2, { 0.5,  0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 4,  0.5   },
    // numerical location on surface 7 of region1
    { 4, { 0.01, 0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 8,  0.0   },
    { 2, { 0.01, 0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 4,  0.01  },
    // numerical location inside overlap
    { 10,{ 0.005,0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 8,  0.005 },
    { 2, { 0.005,0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 4,  0.005 },
    // numerical location on surface 3 of region0
    { 10,{ 0.0,  0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 8,  0.01  },
    { 8, { 0.0,  0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 4,  0.0   },
    // numerical location inside region1
    { 10,{-0.5,  0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 8,  0.51  },
    { 8, {-0.5,  0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 10, 0.5   },
    // numerical location on surface 9 of region1
    { 10,{-1.0,  0.0, 0.0 }, {  1.0, 0.0, 0.0 }            , 8,  1.01  },
    { 8, {-1.0,  0.0, 0.0 }, { -1.0, 0.0, 0.0 }            , 10, 0.0   }  };

  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range surfs, vols;
  const int two = 2, three = 3;
  const void* ptrs[] = { &two, &three };
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs, 1, surfs );
  CHKERR;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs+1, 1, vols );
  CHKERR;
  
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() 
              << std::endl;
    return MB_FAILURE;
  }
  const unsigned num_surf = 12;
  if (surfs.size() != num_surf) {
    std::cerr << "ERROR: Expected " << num_surf << " surfaces in input, found " 
              << surfs.size() << std::endl;
    return MB_FAILURE;
  }
  
  int ids[num_surf];
  rval = moab->tag_get_data( dagmc->id_tag(), surfs, ids );
  CHKERR;
  EntityHandle surf[num_surf];
  std::copy( surfs.begin(), surfs.end(), surf );
  
  const int num_test = sizeof(tests) / sizeof(tests[0]);
  for (int i = 0; i < num_test; ++i) {
    int* ptr = std::find( ids, ids+num_surf, tests[i].prev_surf );
    unsigned idx = ptr - ids;
    if (idx >= num_surf) {
      std::cerr << "Surface " << tests[i].prev_surf << " not found." << std::endl;
      return MB_FAILURE;
    }
    //const EntityHandle src_surf = surf[idx];
    
    ptr = std::find( ids, ids+num_surf, tests[i].hit_surf );
    idx = ptr - ids;
    if (idx >= num_surf) {
      std::cerr << "Surface " << tests[i].hit_surf << " not found." << std::endl;
      return MB_FAILURE;
    }
    const EntityHandle hit_surf = surf[idx];

    double dist;
    EntityHandle result;
    rval = dagmc->ray_fire( vols.front(), 
                           tests[i].origin, tests[i].direction,
                           result, dist );
    
    if (result != hit_surf || fabs(dist - tests[i].distance) > 1e-6) {
      EntityHandle *p = std::find( surf, surf+6, result );
      idx = p - surf;
      int id = idx > num_surf-1 ? 0 : ids[idx];
      
      std::cerr << "Rayfire test failed for " << std::endl
                << "\t ray from (" << tests[i].origin[0] 
                << ", " << tests[i].origin[1] << ", "
                << tests[i].origin[2] << ") going ["
                << tests[i].direction[0] << ", " 
                << tests[i].direction[1] << ", "
                << tests[i].direction[2] << "]" << std::endl
                << "\t Beginning on surface " << tests[i].prev_surf << std::endl
                << "\t Expected to hit surface " << tests[i].hit_surf << " after " 
                << tests[i].distance << " units." << std::endl
                << "\t Actually hit surface " << id << " after " << dist << " units."
                << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

struct PointInVol { double coords[3]; int result; double dir[3]; };

ErrorCode test_point_in_volume( DagMC * dagmc )
{
  const char* const NAME_ARR[] = { "Boundary", "Outside", "Inside" };
  const char* const* names = NAME_ARR + 1;
  const int INSIDE = 1, OUTSIDE = 0, BOUNDARY = -1;
  const struct PointInVol tests[] = {
    { { 0.0, 0.0, 0.5 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 0.0, 0.0,-0.5 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.7, 0.0, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { {-0.7, 0.0, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.0,-0.7, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.0,-0.7, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 1.1, 1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1, 1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1,-1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1,-1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1, 1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1, 1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1,-1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1,-1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
  // Add some directions to test special cases of edge/node intersection
    { { 0.5, 0.0, 0.0 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.5, 0.0, 0.0 }, INSIDE,  { 1.0, 0.0, 0.0} },
    { { 0.0, 0.0, 2.0 }, OUTSIDE, { 0.0, 0.0,-1.0} },
    { { 0.5, 0.0,-0.5 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.5,-0.5,-2.0 }, OUTSIDE, { 0.0, 0.0, 1.0} } };

    //    { { 1.0, 0.0, 0.0 }, BOUNDARY}, MCNP doesn't return on boundary
    //{ {-1.0, 0.0, 0.0 }, BOUNDARY},
    //{ { 0.0, 1.0, 0.0 }, BOUNDARY},
    //{ { 0.0,-1.0, 0.0 }, BOUNDARY},
    //{ { 0.0, 0.0, 0.0 }, BOUNDARY},
    //{ { 0.0, 0.0,-1.0 }, BOUNDARY} };
  const int num_test = sizeof(tests) / sizeof(tests[0]);

  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();

  Range vols;
  const int three = 3;
  const void* ptr = &three;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, &ptr, 1, vols );
  CHKERR;
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() << std::endl;
    return MB_FAILURE;
  }
  const EntityHandle vol = vols.front();

  for (int i = 0; i < num_test; ++i) {
    int result;
    rval = dagmc->point_in_volume( vol, tests[i].coords,
                                  result, tests[i].dir );
    CHKERR;
    if (result != tests[i].result) {
      std::cerr << "ERROR testing point_in_volume[" << i << "]:" << std::endl
                << "\tExpected " << names[tests[i].result] 
                << " for (" << tests[i].coords[0] << ", "
                << tests[i].coords[1] << ", " << tests[i].coords[2]
                << ").  Got " << names[result] << std::endl;
      return MB_FAILURE;
    }
    
      // point_in_volume_slow doesn't to boundary.
    if (tests[i].result == BOUNDARY)
      continue;
     
    rval = dagmc->point_in_volume_slow( vol, tests[i].coords, result );
    CHKERR;
      
    if (result != tests[i].result) {
      std::cerr << "ERROR testing point_in_volume_slow[" << i << "]:" << std::endl
                << "\tExpected " << names[tests[i].result] 
                << " for (" << tests[i].coords[0] << ", "
                << tests[i].coords[1] << ", " << tests[i].coords[2]
                << ").  Got " << names[result] << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode overlap_test_point_in_volume( DagMC * dagmc )
{
  const char* const NAME_ARR[] = { "Boundary", "Outside", "Inside" };
  const char* const* names = NAME_ARR + 1;
  const int INSIDE = 1, OUTSIDE = 0, BOUNDARY = -1;
  const struct PointInVol tests[] = {
    { { 0.5, 0.0, 0.5 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.5, 0.0,-0.5 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.5, 0.0, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { {-0.5, 0.0, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.5, 0.5, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 0.5,-0.5, 0.0 }, INSIDE , { 0.0, 0.0, 0.0} },
    { { 1.1, 1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1, 1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1,-1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1,-1.1, 1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1, 1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1, 1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { {-1.1,-1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
    { { 1.1,-1.1,-1.1 }, OUTSIDE, { 0.0, 0.0, 0.0} },
  // Add some directions to test special cases of edge/node intersection
    { { 0.5, 0.0, 0.0 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.5, 0.0, 0.0 }, INSIDE,  { 1.0, 0.0, 0.0} },
    { { 0.0, 0.0, 2.0 }, OUTSIDE, { 0.0, 0.0,-1.0} },
    { { 0.5, 0.0,-0.5 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.5,-0.5,-2.0 }, OUTSIDE, { 0.0, 0.0, 1.0} },
  // Test some points in the overlap
    { { 0.005, 0.0, 0.0 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.005, 0.0, 0.0 }, INSIDE,  { 1.0, 0.0, 0.0} },
    { { 0.005, 0.0, 2.0 }, OUTSIDE, { 0.0, 0.0,-1.0} },
    { { 0.005, 0.0,-0.5 }, INSIDE,  {-1.0, 0.0, 0.0} },
    { { 0.005,-0.5,-2.0 }, OUTSIDE, { 0.0, 0.0, 1.0} } };

  const int num_test = sizeof(tests) / sizeof(tests[0]);

  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();

  Tag dim_tag = dagmc->geom_tag();
  Range vols;
  const int three = 3;
  const void* ptr = &three;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag,
                                            &ptr, 1, vols );
  CHKERR;
  if (vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 volumes in input, found " << vols.size() 
              << std::endl;
    return MB_FAILURE;
  }
  const EntityHandle vol = vols.front();
  for (int i = 0; i < num_test; ++i) {
    int result;
    rval = dagmc->point_in_volume( vol, tests[i].coords,
                                  result, tests[i].dir );
    CHKERR;
    if (result != tests[i].result) {
      std::cerr << "ERROR testing point_in_volume[" << i << "]:" << std::endl
                << "\tExpected " << names[tests[i].result] 
                << " for (" << tests[i].coords[0] << ", "
                << tests[i].coords[1] << ", " << tests[i].coords[2]
                << ").  Got " << names[result] << std::endl;
      return MB_FAILURE;
    }
    
      // point_in_volume_slow doesn't to boundary.
    if (tests[i].result == BOUNDARY)
      continue;
     
    rval = dagmc->point_in_volume_slow( vol, tests[i].coords, result ); 
    CHKERR;
      
    if (result != tests[i].result) {
      std::cerr << "ERROR testing point_in_volume_slow[" << i << "]:" << std::endl
                << "\tExpected " << names[tests[i].result] 
                << " for (" << tests[i].coords[0] << ", "
                << tests[i].coords[1] << ", " << tests[i].coords[2]
                << ").  Got " << names[result] << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode overlap_test_tracking( DagMC * dagmc )
{
  /* Track a particle from left (-x) to right (+x) through an overlap.
                ____________________
                |                   |
       -x <---  | region1 | overlap | region0 |   ---> +x
                          |                   |
                          ---------------------
       -x <--- -1.0      0.0      0.01        1.0 ---> +x
     surf_id:   10        4         8         2                     */

  // get the surfaces and volumes
  Tag dim_tag = dagmc->geom_tag();
  Range surfs, explicit_vols;
  const int two = 2, three = 3;
  const void* ptrs[] = { &two, &three };
  ErrorCode rval;
  Interface *moab = dagmc->moab_instance();
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs, 1, surfs );
  CHKERR;
  rval = moab->get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 
                                            ptrs+1, 1, explicit_vols );
  CHKERR;
  
  if (explicit_vols.size() != 2) {
    std::cerr << "ERROR: Expected 2 explicit volumes in input, found " 
              << explicit_vols.size() << std::endl;
    return MB_FAILURE;
  }
  const unsigned num_surf = 12;
  if (surfs.size() != num_surf) {
    std::cerr << "ERROR: Expected " << num_surf << " surfaces in input, found " 
              << surfs.size() << std::endl;
    return MB_FAILURE;
  }
  const EntityHandle explicit_vol = explicit_vols.front();

  // start particle
  double point[] = { -0.9, 0, 0 };
  const double dir[] = { 1, 0, 0 };
  EntityHandle vol = explicit_vol;
  int result;
  const int INSIDE = 1; // OUTSIDE = 0, BOUNDARY = -1;
  rval = dagmc->point_in_volume( explicit_vol, point, result, dir );
  CHKERR;
  if (result != INSIDE) {
    std::cerr << "ERROR: particle not inside explicit volume" << std::endl;
    return MB_FAILURE;
  }

  // get next surface
  double dist;
  EntityHandle next_surf;
  GeomQueryTool::RayHistory history;
  rval = dagmc->ray_fire( vol, point, dir, next_surf, dist, &history );
  CHKERR;    
  if (next_surf != surfs[7] || fabs(dist - 0.91) > 1e-6) {
    std::cerr << "ERROR: failed on advance 1" << std::endl;
    return MB_FAILURE;
  }
  for(unsigned i=0; i<3; i++) point[i]+=dist*dir[i];

  // get the next volume (implicit complement)
  EntityHandle next_vol;
  rval = dagmc->next_vol( next_surf, vol, next_vol ); 
  CHKERR;

  // get the next surface (behind numerical location)
  vol       = next_vol;
  rval = dagmc->ray_fire( vol, point, dir, next_surf, dist, &history );
  CHKERR;    
  if (next_surf != surfs[3] || fabs(dist - 0.0) > 1e-6) {
    std::cerr << "ERROR: failed on advance 2" << std::endl;
    return MB_FAILURE;
  }
  for(unsigned i=0; i<3; i++) point[i]+=dist*dir[i];

  // get the next volume (the explicit volume)
  rval = dagmc->next_vol( next_surf, vol, next_vol );
  CHKERR;

  // get the next surface
  vol       = next_vol;
  rval = dagmc->ray_fire( vol, point, dir, next_surf, dist, &history );
  CHKERR;    
  if (next_surf != surfs[1] || fabs(dist - 0.99) > 1e-6) {
    std::cerr << "ERROR: failed on advance 3" << std::endl;
    return MB_FAILURE;
  }
  for(unsigned i=0; i<3; i++) point[i]+=dist*dir[i];

  return MB_SUCCESS;
}


