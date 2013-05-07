#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"

MBInterface *MOAB(); 
namespace geometry 
{
  MBErrorCode measure( const MBEntityHandle set, const MBTag geom_tag, double &size );
  
  double length( std::vector<MBEntityHandle> edges );

  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 );
  MBErrorCode dist_between_verts( const MBEntityHandle v0, const MBEntityHandle v1,double &d);
  double dist_between_verts( double coords0[], double coords1[] );
  double dist_between_verts( MBEntityHandle vert0, MBEntityHandle vert1 ); 
  
  double triangle_area( const MBCartVect a, const MBCartVect b, const MBCartVect c);
  MBErrorCode triangle_area( const MBEntityHandle conn[], double &area );
  MBErrorCode triangle_area( const MBEntityHandle triangle, double &area );
  double triangle_area( MBRange triangles );

  MBErrorCode measure_volume( MBEntityHandle volume, double& result );

  MBErrorCode surface_sense( MBEntityHandle volume, int num_surfaces,const MBEntityHandle* surfaces,int* senses_out );
  MBErrorCode surface_sense( MBEntityHandle volume, MBEntityHandle surface, int& sense_out );

  MBTag get_tag( const char* name, int size, MBTagType store,MBDataType type, const void* def_value, bool create_if_missing);

  bool triangle_degenerate( const MBEntityHandle tri );
  bool triangle_degenerate( const MBEntityHandle v0, const MBEntityHandle v1, const MBEntityHandle v2 );
}
