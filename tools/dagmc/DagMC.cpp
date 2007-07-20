#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "MBCartVect.hpp"
#include "MBRange.hpp"
#include "MBCore.hpp"
#include "MBGeomUtil.hpp"

#ifdef CGM
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "CubitVector.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "cgm2moab.hpp"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define CATEGORY_TAG_LENGTH 32

DagMC *DagMC::instance_ = NULL;

const bool debug = false;

DagMC::DagMC(MBInterface *mb_impl) 
    : mbImpl(mb_impl), obbTree(mb_impl), 
      moabMCNPTolerance(1e-6), moabMCNPSourceCell(0), moabMCNPUseDistLimit(false)
{
  options[0] = Option( "source_cell",        "source cell ID, or zero if unknown", "0" );
  options[1] = Option( "distance_tolerance", "positive real value", "0.001" );
  options[2] = Option( "use_distance_limit", "one to enable distance limit optimization, zero otherwise", "0" );
  options[3] = Option( "use_cad", "one to ray-trace to cad, zero for just facets", "0" );
  
  memset( specReflectName, 0, NAME_TAG_SIZE );
  strcpy( specReflectName, "spec_reflect" );
  memset( whiteReflectName, 0, NAME_TAG_SIZE );
  strcpy( whiteReflectName, "white_reflect" );

  distanceLimit = std::numeric_limits<double>::max();

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, false);
  
  idTag = get_tag( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER );
  
  geomTag = get_tag( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER );

  obbTag = get_tag( MB_OBB_TREE_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE );
}

MBErrorCode DagMC::ray_fire(const MBEntityHandle vol, const MBEntityHandle last_surf_hit, 
                             const int num_pts,
                             const double uuu, const double vvv, const double www,
                             const double xxx, const double yyy, const double zzz,
                             const double huge,
                             double &dist_traveled, MBEntityHandle &next_surf_hit) 
{
  if (debug) {
    std::cout << "Vol " << id_by_index(3, index_by_handle(vol)) << ", xyz = " 
              << xxx << " " << yyy << " " << zzz 
              << ", uvw = " 
              << uuu << " " << vvv << " " << www << std::endl;
  }
  
  assert(vol - setOffset < rootSets.size());
  
  MBEntityHandle root = rootSets[vol - setOffset];
  MBErrorCode rval;
    // delcare some stuff static so we don't need to re-created
    // it for every call
  static std::vector<double> distances;
  static std::vector<MBEntityHandle> surfaces;
  
  
    // do ray fire
  const double point[] = {xxx, yyy, zzz};
  const double dir[] = {uuu, vvv, www};
  distances.clear();
  surfaces.clear();
  double len = use_dist_limit() ? distance_limit() : huge;

  rval = obbTree.ray_intersect_sets( distances,
                                     surfaces, 
                                     root, 
                                     tolerance(),
                                     2,
                                     point, dir,
                                     &len);
  assert( MB_SUCCESS == rval );

#ifdef CGM
  if (useCAD) {
#ifdef CUBIT_CGM
    std::cout << "USE_CAD = 1 not supported with this build of CGM/DagMC;"
              << std:: endl
              << "need an ACIS-based install of CGM." << std::endl;
    return MB_NOT_IMPLEMENTED;
#else    
    rval = CAD_ray_intersect(point, dir, huge,
                             distances, surfaces, len);
    if (MB_SUCCESS != rval) return rval;
#endif
  }
#endif
  
    // Find smallest intersection
  if (distances.empty()) {
    next_surf_hit = 0;
    dist_traveled = (use_dist_limit() ? len*10.0 : huge);
    if (debug) {
      std::cout << "next surf hit = " << 0 << ", dist = (huge)" << std::endl;
    }
    return MB_SUCCESS;
  }
  int smallest = std::min_element( distances.begin(), distances.end() ) - distances.begin();
  
    // If intersected previous surface near start of ray, reject it
  if (surfaces[smallest] == last_surf_hit && distances[smallest] < tolerance()) {
    distances.erase( distances.begin() + smallest );
    surfaces.erase( surfaces.begin() + smallest );
  
      // Find smallest intersection
    if (distances.empty()) {
      next_surf_hit = 0;
      dist_traveled = (use_dist_limit() ? distance_limit()*10.0 : huge);
      if (debug) {
        std::cout << "next surf hit = " << 0 << ", dist = (huge)" << std::endl;
      }
      return MB_SUCCESS;
    }
    smallest = std::min_element( distances.begin(), distances.end() ) - distances.begin();
  }
  
    // return results
  dist_traveled = distances[smallest];
  next_surf_hit = surfaces[smallest];
  
  if (debug) {
    std::cout << "next surf hit = " <<  id_by_index(2, index_by_handle(next_surf_hit)) 
              << ", dist = " << dist_traveled << std::endl;
  }

  return MB_SUCCESS;
}

void DagMC::create_instance(MBInterface *mb_impl) 
{
  if (NULL == mb_impl) mb_impl = new MBCore();
  instance_ = new DagMC(mb_impl);
}

void DagMC::write_settings( FILE* fptr, bool desc ) {
  int num_opt = sizeof(options) / sizeof(options[0]);
  if (desc) for (int i = 0; i < num_opt; ++i) 
    fprintf( fptr, "%s = %s  # %s\n", options[i].name.c_str(), options[i].value.c_str(), options[i].desc.c_str() );
  else for (int i = 0; i < num_opt; ++i) 
    fprintf( fptr, "%s = %s\n", options[i].name.c_str(), options[i].value.c_str() );
}

void DagMC::parse_settings() {
  moabMCNPSourceCell = atoi( options[0].value.c_str() );
  if (moabMCNPSourceCell < 0) {
    std::cerr << "Invalid source_cell = " << moabMCNPSourceCell << std::endl;
    exit(2);
  }
  moabMCNPTolerance = strtod( options[1].value.c_str(), 0 );
  if (moabMCNPTolerance <= 0 || moabMCNPTolerance > 1) {
    std::cerr << "Invalid distance_tolerance = " << moabMCNPTolerance << std::endl;
    exit(2);
  }
  moabMCNPUseDistLimit = !!atoi( options[2].value.c_str() );
}

void DagMC::read_settings( const char* filename )
{
  int num_opt = sizeof(options) / sizeof(options[0]);
  FILE* file;
  if (filename && (file = fopen( filename, "r" ))) {
    int line = 0;
    char buffer[256];
    bool havenl = true;
    while ( fgets( buffer, sizeof(buffer), file ) ) {
      if (havenl) {
        havenl = false;
        ++line;
      }
      int len = strlen(buffer);
      if (len && buffer[len-1] == '\n')
        havenl = true;
        
        // chop off any thing after the '#' comman indicator
      char* c = strchr( buffer, '#' );
      if (c) *c = '\0';
        // skip leading white space
      char* p = buffer;
      while (isspace(*p)) ++p;
        // if empty line, done
      if (!*p)
        continue;
        // find '='
      c = strchr( p, '=' );
      if (!c) {
        fprintf(stderr, "Invalid option at line %d of config file '%s'\n", line, filename );
        exit(2);
      }
        // get rid of white space around '='
      char* v = c+1;
      *c = ' ';
      while (c > p && isspace(*c)) {
        *c = '\0';
        --c;
      }
      while (isspace(*v))
        ++v;
      
        // search for option
      bool found = false;
      for (int i = 0; i < num_opt; ++i) {
        if (options[i].name == p) {
          found = true;
          options[i].value = v;
        }
      }
      if (!found) 
        fprintf(stderr,"Warning: unknown option at line %d of '%s': %s\n", line, filename, p );
    } // while( fgets() )
    
    fclose(file);
  } // file(file)
  
  // always do this, even if no file was found.
  // it will either read the default values, or re-parse
  // the values previously read from some other file.
  parse_settings();
}

// point_in_volume_slow, including poly_solid_angle helper subroutine
// are adapted from "Point in Polyhedron Testing Using Spherical Polygons", Paulo Cezar 
// Pinto Carvalho and Paulo Roma Cavalcanti, _Graphics Gems V_, pg. 42.  Original algorithm
// was described in "An Efficient Point In Polyhedron Algortihm", Jeff Lane, Bob Magedson, 
// and Mike Rarick, _Computer Visoin, Graphics, and Image Processing 26_, pg. 118-225, 1984.

// helper function for point_in_volume_slow.  calculate area of a polygon 
// projected into a unit-sphere space
MBErrorCode DagMC::poly_solid_angle( MBEntityHandle face, const MBCartVect& point, double& area )
{
  MBErrorCode rval;
  
    // Get connectivity
  const MBEntityHandle* conn;
  int len;
  rval = moab_instance()->get_connectivity( face, conn, len, true );
  if (MB_SUCCESS != rval)
    return rval;
  
    // Allocate space to store vertices
  static MBCartVect coords_static[4];
  static std::vector<MBCartVect> coords_dynamic;
  MBCartVect* coords = coords_static;
  if ((unsigned)len > (sizeof(coords_static)/sizeof(coords_static[0]))) {
    coords_dynamic.resize(len);
    coords = &coords_dynamic[0];
  }
  
    // get coordinates
  rval = moab_instance()->get_coords( conn, len, coords->array() );
  if (MB_SUCCESS != rval)
    return rval;
  
    // calculate normal
  MBCartVect norm(0.0), v1, v0 = coords[1] - coords[0];
  for (int i = 2; i < len; ++i) {
    v1 = coords[i] - coords[0];
    norm += v0 * v1;
    v0 = v1;
  }
  
    // calculate area
  double s, ang;
  area = 0.0;
  MBCartVect r, n1, n2, b, a = coords[len-1] - coords[0];
  for (int i = 0; i < len; ++i) {
    r = coords[i] - point;
    b = coords[(i+1)%len] - coords[i];
    n1 = a * r;
    n2 = r * b;
    s = (n1 % n2) / (n1.length() * n2.length());
    ang = s <= -1.0 ? M_PI : s >= 1.0 ? 0.0 : acos(s);
    s = (b * a) % norm;
    area += s > 0.0 ? M_PI - ang : M_PI + ang;
    a = -b;
  }
  
  area -= M_PI * (len - 2);
  if ((norm % r) > 0)
    area = -area;
  return MB_SUCCESS;
}
  
  
// use spherical area test to determine inside/outside of a polyhedron.
MBErrorCode DagMC::point_in_volume_slow( MBEntityHandle volume,
                                  double x, double y, double z,
                                  int& result )
{
  MBErrorCode rval;
  MBRange faces;
  std::vector<MBEntityHandle> surfs;
  std::vector<int> senses;
  double sum = 0.0;
  const MBCartVect point(x,y,z);
  
  rval = moab_instance()->get_child_meshsets( volume, surfs );
  if (MB_SUCCESS != rval)
    return rval;
  
  senses.resize( surfs.size() );
  rval = surface_sense( volume, surfs.size(), &surfs[0], &senses[0] );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (unsigned i = 0; i < surfs.size(); ++i) {
    if (!senses[i])  // skip non-manifold surfaces
      continue;
     
    double surf_area = 0.0, face_area;
    faces.clear();
    rval = moab_instance()->get_entities_by_dimension( surfs[i], 2, faces );
    if (MB_SUCCESS != rval)
      return rval;
    
    for (MBRange::iterator j = faces.begin(); j != faces.end(); ++j) {
      rval = poly_solid_angle( *j, point, face_area );
      if (MB_SUCCESS != rval)
        return rval;
      
      surf_area += face_area;
    }
    
    sum += senses[i] * surf_area;
  }
  
  result = fabs(sum) > 2.0*M_PI;
  return MB_SUCCESS;
}


// Fast point_in_volume code.  This function assumes that there is an 
// OBB tree containing only manifold surfaces in the volume, such that
// the closest facet to the input position found using the OBB tree is
// guaranteed to be on the boundary of the volume.  
//
// This function will find the closest point on the volume boundary.
// If that point is closest to the interior of a facet, the direction
// of the surface normal of that facet can be used to determine if the
// point is inside or outside of the volume.  If the point is closest to
// an edge joining two facets, then the point is inside the volume if 
// that edge is at a concavity wrt the volume.  If the point is closest
// to a vertex in the facetting, it relies on the fact that the closest
// vertex cannot be a saddle: it must be a strict concavity or convexity.
MBErrorCode DagMC::point_in_volume( MBEntityHandle volume, 
                             double x, double y, double z,
                             int& result )
{
  const double epsilon = moabMCNPTolerance;
  
    // Get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  MBEntityHandle root = rootSets[volume - setOffset];

    // Get closest triangles in volume
  static std::vector<MBEntityHandle> tris, surfs;
  tris.clear();
  surfs.clear();
  const MBCartVect point(x,y,z);
  MBErrorCode rval = obbTree.closest_to_location( point.array(), root, epsilon, tris, &surfs );
  if (MB_SUCCESS != rval) return rval;
  
    // Get sense of each surface wrt volume
  std::vector<int> senses( surfs.size() );
  rval = surface_sense( volume, surfs.size(), &surfs[0], &senses[0] );
  if (MB_SUCCESS != rval) return rval;
  
    // Get closest point on each triangle, and group triangles
    // by intersection type
  std::vector<int> edge_tris, vertex_tris;
  
    // Get corrected normal, closest point and topolgocal closest for each triangle
  std::vector<MBCartVect> normals( tris.size() ), closest( tris.size() );
  std::vector<int> topo( tris.size() );
  const MBEntityHandle* conn;
  int len;
  MBCartVect coords[3];
  for (unsigned i = 0; i < tris.size(); ++i) {
      // volume tree shouldn't contain non-manifold surfaces
    assert(senses[i]);
    
    rval = moab_instance()->get_connectivity( tris[i], conn, len );
    if (MB_SUCCESS != rval || len != 3) 
      return MB_FAILURE;
    rval = moab_instance()->get_coords( conn, len, coords[0].array() );
    if (MB_SUCCESS != rval)
      return rval;
    
    MBGeomUtil::closest_location_on_tri( point, coords, epsilon, closest[i], topo[i] );
    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normals[i] = senses[i] * coords[1] * coords[2];
  }
  
    // If within tolerance of triangle, then return that point is on boundary
  if ((point - closest[0]).length() < epsilon) {
    result = -1;
    return MB_SUCCESS;
  }
  
    // if only one, then return 
  if (tris.size() == 1) {
      // Ignore the possibility of being in the plane of the 
      // triangle.  If parallel to one triangle, should have
      // at least one other triangle.
    result = (normals[0] % (closest[0] - point) > 0.0);
    return MB_SUCCESS;
  }
  
    // If we found any triangle for which the closest point was in
    // the interior, just use that one
  std::vector<int>::iterator ii = std::find( topo.begin(), topo.end(), 6 );
  if (ii != topo.end()) {
    int idx = ii - topo.begin();
    result = (normals[idx] % (closest[idx] - point) > 0.0);
    return MB_SUCCESS;
  }
  
    // Otherwise if we found any triangle for which the
    // closest point was at an edge, find it's neighbor and use the pair
  for (ii = topo.begin(); ii != topo.end() && *ii < 3; ++ii);
  if (ii != topo.end()) // found one
  {
    unsigned idx1 = ii - topo.end();
    int idx2 = -1;
    MBEntityHandle edge[2];
    moab_instance()->get_connectivity( tris[idx1], conn, len, true );
    edge[0] = conn[*ii - 3];
    edge[1] = conn[(*ii - 2) % 3];
    for (unsigned i = 0; i < tris.size(); ++i) {
      if (i == idx1)
        continue;
      moab_instance()->get_connectivity( tris[i], conn, len, true );
      int j;
      for (j = 0; j < 3; ++j)
        if (conn[j] == edge[0])
          if (conn[(j+1)%3] == edge[1] || conn[(j+2)%3] == edge[1])
            break;
      if (j < 3) {
        idx2 = i;
        break;
      }
    }
    
    if (idx2 < 0) { // if no adjacent triangle
      assert(false);
      result = (normals[idx1] % (closest[idx1] - point) > 0.0);
      return MB_SUCCESS;
    }
    
      // Check if inside or outside of both facets first.
      // This avoids potential rounding/accuracy issues 
      // when testing for convexity/concavity if the facets
      // are near co-planer, and handles the special case 
      // where both facets are exactly coplanar.
    const double dot1 = normals[idx1] % (closest[idx1] - point);
    const double dot2 = normals[idx2] % (closest[idx2] - point);
    const bool inside1 = dot1 > -std::numeric_limits<double>::epsilon();
    const bool inside2 = dot2 > -std::numeric_limits<double>::epsilon();
    const bool outside1 = dot1 < std::numeric_limits<double>::epsilon();
    const bool outside2 = dot2 < std::numeric_limits<double>::epsilon();
    if (inside1 && inside2)
      result = 1;
    else if (outside1 && outside2)
      result = 0;
    else {
        // Edge is at a concavity if the projection of the normal of 
        // triangle A into the plane of triangle B points into triangle
        // B.  It is a convexity if it points out of triangle B (from the
        // shared edge).  This test can be reduced to checking the sign of
        // dot product of the normal of triangle A with any vector pointing
        // into triangle B from the shared edge.  
        
        // Get triangle coordinates
      moab_instance()->get_connectivity( tris[idx2], conn, len, true );
      moab_instance()->get_coords( conn, len, coords[0].array() );
        // Get a vector from edge end point to opposite vertex
      const MBCartVect vect2 = coords[(topo[idx2] - 1) % 3] - coords[topo[idx2]];
      result = (vect2 % normals[idx1] >= 0.0); // true if a concavity
        result = inside1 || inside2;
    }
    return MB_SUCCESS;
  }
  
    // Otherwise pick a vertex that the point was closest to remove
    // any triangles from lists that are not adjacent to that vertex.
  assert( topo[0] < 3 ); // if here, must be closest to vertex
  rval = moab_instance()->get_connectivity( tris[0], conn, len );
  if (MB_SUCCESS != rval || len != 3) 
    return MB_FAILURE;
  const MBEntityHandle closest_vert = conn[topo[0]];
  unsigned w = 1;
  for (unsigned r = 1; r < tris.size(); ++r) {
    assert( topo[r] < 3 ); // if here, must be closest to vertex
    moab_instance()->get_connectivity( tris[r], conn, len );
    if (conn[topo[r]] == closest_vert) {
      tris[w] = tris[r];
      surfs[w] = surfs[r];
      normals[w] = normals[r];
      ++w;
    }
  }
  tris.resize(w);
  surfs.resize(w);
  normals.resize(w);
  
    // use algorithm from:
    // (referance)
    // to determine inside vs. outside.
  
    // first find single triangle for which all other triangles are to
    // one side.
  bool some_above, some_below;
  MBCartVect vertcoords, closestcoords;
  moab_instance()->get_coords( &closest_vert, 1, closestcoords.array() );
  for (unsigned i = 0; i < tris.size(); ++i) {
    some_above = some_below = false;
    for (unsigned j = 0; j < tris.size(); ++j) {
      if (j == i)
        continue;
      
      moab_instance()->get_connectivity( tris[j], conn, len );
      for (int k = 0; k < len; ++k) {
        moab_instance()->get_coords( conn+k, 1, vertcoords.array() );
        const double dot = normals[i] % (closestcoords - vertcoords);
        if (dot * dot / (normals[i] % normals[i]) > epsilon * epsilon) {
          if (dot < 0.0)
            some_above = true;
          else
            some_below = true;
        }
      }
    }
    
      // All triangles are roughly co-planar: if we're inside of
      // one then we're inside of all.
    if (!some_above && !some_below) {
      result = (normals[i] % (closestcoords - point)) > 0.0;
      return MB_SUCCESS;
    }
      // All other triangles to one side of this triangle:
    else if (some_above != some_below) {
        // If all other triangles are above this one (some_above == true),
        // then the vertex is a concavity so the input point is inside.
        // If all other triangles are below this one (some_above == false),
        // then the vertex is a convexity and the input point must be outside.
      result = some_above;
      return MB_SUCCESS;
    }
  }
  
    // If we got here, then something's wrong somewhere.
    // We appear to be closest to a "saddle" vertex in the
    // triangulation.  That shoudn't be possible (must be 
    // closer to at least one of the adjacent triagles than
    // to the saddle vertex.)
  assert(false /*shouldn't be here*/);
  return point_in_volume_slow( volume, x, y, z, result );
}

// detemine distance to nearest surface
MBErrorCode DagMC::closest_to_location( MBEntityHandle volume, double* coords, double& result)
{

  const double epsilon = moabMCNPTolerance;
  
    // Get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  MBEntityHandle root = rootSets[volume - setOffset];

    // Get closest triangles in volume
  const MBCartVect point(coords);
  MBCartVect nearest;
  MBEntityHandle facet_out;
  MBErrorCode rval = obbTree.closest_to_location( point.array(), root, epsilon, nearest.array(), facet_out );
  if (MB_SUCCESS != rval) return rval;

  // calculate distance between point and nearest facet
  result = (point-nearest).length();
  
  return MB_SUCCESS;

}




// calculate volume of polyhedron
MBErrorCode DagMC::measure_volume( MBEntityHandle volume, double& result )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle> surfaces, surf_volumes;
  result = 0.0;
  
    // get surfaces from volume
  rval = moab_instance()->get_child_meshsets( volume, surfaces );
  if (MB_SUCCESS != rval) return rval;
  
    // get surface senses
  std::vector<int> senses( surfaces.size() );
  rval = surface_sense( volume, surfaces.size(), &surfaces[0], &senses[0] );
  if (MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }
  
  for (unsigned i = 0; i < surfaces.size(); ++i) {
      // skip non-manifold surfaces
    if (!senses[i])
      continue;
    
      // get triangles in surface
    MBRange triangles;
    rval = moab_instance()->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (MB_SUCCESS != rval) 
      return rval;
    if (!triangles.all_of_type(MBTRI)) {
      std::cerr << "ERROR: Surface contains non-triangle elements.  Cannot calculate volume." << std::endl;
      return MB_FAILURE;
    }
    
      // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const MBEntityHandle *conn;
    int len;
    MBCartVect coords[3];
    for (MBRange::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = moab_instance()->get_connectivity( *j, conn, len, true );
      if (MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = moab_instance()->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;
    
      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }
  
  result /= 6.0;
  return MB_SUCCESS;
}

// sum area of elements in surface
MBErrorCode DagMC::measure_area( MBEntityHandle surface, double& result )
{
    // get triangles in surface
  MBRange triangles;
  MBErrorCode rval = moab_instance()->get_entities_by_dimension( surface, 2, triangles );
  if (MB_SUCCESS != rval) 
    return rval;
  if (!triangles.all_of_type(MBTRI)) {
    std::cerr << "ERROR: Surface contains non-triangle elements.  Cannot calculate area." << std::endl;
    return MB_FAILURE;
  }

    // calculate sum of area of triangles
  result = 0.0;
  const MBEntityHandle *conn;
  int len;
  MBCartVect coords[3];
  for (MBRange::iterator j = triangles.begin(); j != triangles.end(); ++j) {
    rval = moab_instance()->get_connectivity( *j, conn, len, true );
    if (MB_SUCCESS != rval) return rval;
    assert(3 == len);
    rval = moab_instance()->get_coords( conn, 3, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    coords[0] = coords[1] * coords[2];
    result += coords[0].length();
  }
  result *= 0.5;
  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
MBErrorCode DagMC::surface_sense( MBEntityHandle volume, 
                           int num_surfaces,
                           const MBEntityHandle* surfaces,
                           int* senses_out )
{
  
    // get sense of surfaces wrt volumes
  static const MBTag tag = get_tag( "GEOM_SENSE_2", 2*sizeof(MBEntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE );

  std::vector<MBEntityHandle> surf_volumes( 2*num_surfaces );
  MBErrorCode rval = moab_instance()->tag_get_data( tag, surfaces, num_surfaces, &surf_volumes[0] );
  if (MB_SUCCESS != rval)  return rval;
  
  const MBEntityHandle* end = surfaces + num_surfaces;
  std::vector<MBEntityHandle>::const_iterator surf_vols = surf_volumes.begin();
  while (surfaces != end) {
    MBEntityHandle forward = *surf_vols; ++surf_vols;
    MBEntityHandle reverse = *surf_vols; ++surf_vols;
    if (volume == forward) 
      *senses_out = (volume != reverse); // zero if both, otherwise 1
    else if (volume == reverse)
      *senses_out = -1;
    else 
      return MB_ENTITY_NOT_FOUND;
    
    ++surfaces;
    ++senses_out;
  }
  
  return MB_SUCCESS;
}

MBErrorCode DagMC::load_file_and_init(const char* cfile,
                                       const int clen,
                                       const char* ffile,
                                       const int flen)
{
    // Always do this first to make sure we have at least the
    // default values if noone ever calls read_settings(). If
    // read_settings() has already been called, it will just
    // re-parse the same values.
  read_settings(0);

  MBErrorCode rval;
  
  std::string scfile = cfile;                                           

    // save whether it's geometry or mesh, building obb later depends
    // on this setting
  bool is_geom = false;
  
  if (scfile.find(".sat",0) > (unsigned int) 0 &&
      scfile.find(".sat",0) <= scfile.length()) {

    is_geom = true;
    
#ifdef CGM
      // initialize cgm
    InitCGMA::initialize_cgma();
    InitCGMA::initialize_engine("ACIS");
    CGMApp::instance()->attrib_manager()->set_all_auto_read_flags(true);
    CGMApp::instance()->attrib_manager()->set_all_auto_actuate_flags(true);
    CubitStatus s = GeometryQueryTool::instance()->import_solid_model( scfile.c_str(), "ACIS_SAT");
    if (CUBIT_SUCCESS != s) {
      std::cerr << "Failed to read '" << scfile << "' of type 'ACIS_SAT'" 
                << std::endl;
      exit(2);
    }
    if (!cgm2moab(moab_instance())) {
      std::cerr << "Internal error copying geometry" << std::endl;
      exit(5);
    }

    if (!useCAD) {
      // ok, we're done with CGM
      CGMApp::instance()->shutdown();
    }
#else
    std::cerr << "CGM not supported in this version." << std::endl;
#endif
  }
  else if (scfile.find(".cub",0) > (unsigned int) 0 &&
           scfile.find(".cub",0) <= scfile.length()) {
    std::cerr << "Error: can't read .cub files yet." << std::endl;
    exit(7);
  }
  else {
    rval = moab_instance()->load_mesh(scfile.c_str());
    if (MB_SUCCESS != rval) {
      std::cerr << "Couldn't read file " << scfile << std::endl;
      return rval;
    }
    else if (useCAD) {
      std::cerr << "useCAD parameter turned on but geometry read from "
                << "a mesh file;" << std::endl
                << " turning that parameter off for you." << std::endl;
      useCAD = 0;
    }
  }

  MBRange surfs, vols;
  const int three = 3;
  const void* const three_val[] = {&three};
  rval = 
    moab_instance()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, 
                                                   three_val, 1, vols );
  if (MB_SUCCESS != rval)
    return rval;

  const int two = 2;
  const void* const two_val[] = {&two};
  rval = moab_instance()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, two_val, 1, surfs );
  if (MB_SUCCESS != rval)
    return rval;

  // Build OBB trees for everything, but only if we only read geometry
  // Changed to build obb tree if tree does not already exist. -- JK
  if (!have_obb_tree()) {
    rval = build_obbs(surfs, vols);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to build obb." << std::endl;
      return rval;
    }
  }

    // build the various index vectors used for efficiency
  rval = build_indices(surfs, vols, is_geom);
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to build surface/volume indices." << std::endl;
    return rval;
  }
  
    // finally, if user entered a facet file name, export facets to
    // that file
  if (ffile && 0 < flen) {
    rval = moab_instance()->write_mesh(ffile);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return rval;
    }
  }
  
  return MB_SUCCESS;
}

bool DagMC::have_obb_tree()
{
  MBRange entities;
  MBErrorCode rval = mbImpl->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &obbTag, 0, 1,
                                                           entities );
  return MB_SUCCESS == rval && !entities.empty();
}                                                    

MBEntityHandle DagMC::entity_by_id( int dimension, int id )
{
  assert(0 <= dimension && 3 >= dimension);
  const MBTag tags[] = { idTag, geomTag };
  const void* const vals[] = { &id, &dimension };
  MBErrorCode rval;
  
  MBRange results;
  rval = moab_instance()->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, results );
  if (MB_SUCCESS != rval || results.empty())
    return 0;
  
  return results.front();
}

int DagMC::id_by_index( int dimension, int index )
{
  MBEntityHandle h = entity_by_index( dimension, index );
  if (!h)
    return 0;
  
  int result = 0;
  moab_instance()->tag_get_data( idTag, &h, 1, &result );
  return result;
}

void DagMC::write_log(std::string ifile, const bool overwrite) 
{
  if (!overwrite) {
      // if not overwriting, test for file, and if it exists, just return;
      // this replaces the old inquir function
    std::ifstream testfile;
    testfile.open(ifile.c_str());

    if (!testfile.fail()) return;
  }
  
  const char *cifile = ifile.c_str();
  std::ofstream cgmfile;
  cgmfile.open(cifile);

  int ncells;
  int nsurfs;
  int ngroups;
  int maxital;
// Get the number of surfaces and cells for dynamic memory                     
  ncells = vol_handles().size()-1;
  nsurfs = surf_handles().size()-1;
  ngroups = group_handles().size()-1;
  MBRange vols_save, surfs_save;
  std::copy(vol_handles().begin()+1, vol_handles().end(), 
            mb_range_inserter(vols_save));
  std::copy(surf_handles().begin()+1, surf_handles().end(), 
            mb_range_inserter(surfs_save));
  
  std::cout << "ncells: " << ncells << std::endl;
  std::cout << "nsurfs: " << nsurfs << std::endl;
  std::cout << "ngroups: " << ngroups << std::endl;
  std::vector<int> mats(ncells);
  std::vector<double> dens(ncells);
  std::vector<int> imp(ncells);
  std::vector<int> surflist(nsurfs);
  std::vector<int> surfbcs(nsurfs);
  std::vector<std::string> gruplist;
// Zero out material and density info                                          
  for (int ii=1; ii <=ncells; ii++) {
    mats[ii-1] = 0;
    dens[ii-1] = 0.0;
    imp[ii-1] = 1;
  }
  maxital = -1;
  for (int ii=1; ii <=nsurfs; ii++)
    surfbcs[ii-1] = 0;
// Read in each group name, check if a material specification                  
// or surface boundary condition                                               
    // get all the groups
  std::vector<std::string> grp_names;
  for (int ii=1; ii <=ngroups; ii++) {
      // get group & names
    MBEntityHandle group = group_handles()[ii];
    grp_names.clear();
    bool success = get_group_names(group, grp_names);
    assert(success);
    if (grp_names.empty()) continue;

      // get sets in the group
    MBRange grp_sets;
    MBErrorCode result = moab_instance()->get_entities_by_type(group, MBENTITYSET, grp_sets);
    if (MB_SUCCESS != result) continue;

    if (grp_names[0].find("mat",0) == 0 && (grp_names[0].find("rho",0) > 0)) {
        // Extract material ID and density                                         
      int cpos = grp_names[0].find("mat",0) + 4;
      std::string matid_str;
      while ((int(grp_names[0][cpos]) >= 48) && (int(grp_names[0][cpos]) <= 57)) {
        matid_str += grp_names[0][cpos];
        cpos++;
      }
      cpos = grp_names[0].find("rho",0) + 4;
      std::string density_str;
      while (((int(grp_names[0][cpos]) >= 48) && (int(grp_names[0][cpos]) <= 57)) ||
             (int(grp_names[0][cpos]) == 46) || (int(grp_names[0][cpos]) == 45) ||
             (int(grp_names[0][cpos]) == 69) || (int(grp_names[0][cpos]) == 101) ||
             (int(grp_names[0][cpos]) == 68) || (int(grp_names[0][cpos]) == 100)) {
        density_str += grp_names[0][cpos];
        cpos++;
      }
      const char *cmatid = matid_str.c_str();
      const char *cdensity = density_str.c_str();
      int matid = atoi(cmatid);
      double density = atof(cdensity);
        // Get the volumes in the current group                                    
      MBRange grp_vols = grp_sets.intersect(vols_save);
      for (MBRange::iterator rit = grp_vols.begin(); rit != grp_vols.end(); rit++) {
        MBEntityHandle vol = *rit;
        int cell_num = index_by_handle(vol);
        mats[cell_num-1]=matid;
        dens[cell_num-1]=density;
      }
    }
    else if ((grp_names[0].find("spec_reflect",0)==0) ||
             (grp_names[0].find("white_reflect",0)==0)) {
        // Check whether group is reflecing or white                               
      int bc_id;
      if (grp_names[0].find("spec_reflect",0)==0)
        bc_id = 1;
      else
        bc_id = 2;
        // Get the surfaces in the current group                                   
      MBRange grp_faces = grp_sets.intersect(surfs_save);
      for (MBRange::iterator rit = grp_faces.begin(); rit != grp_faces.end(); rit++) {
        MBEntityHandle surf = *rit;
        int surf_num = index_by_handle(surf);
        surfbcs[surf_num-1]=bc_id;
      }
    }
    else if (grp_names[0].find("periodic",0)==0) {
        // Get the periodic surface id                                             
      int bc_id;
      std::string bc_str;
      int cpos = 9;
      bc_str.clear();
      while ((int(grp_names[0][cpos]) >= 48) && (int(grp_names[0][cpos]) <= 57)) {
        bc_str += grp_names[0][cpos];
        cpos++;
      }
      const char *cbc_id = bc_str.c_str();
      bc_id = -1*atoi(cbc_id);
        // Get the surfaces in the current group                                   
      MBRange grp_faces = grp_sets.intersect(surfs_save);
      for (MBRange::iterator rit = grp_faces.begin(); rit != grp_faces.end(); rit++) {
        MBEntityHandle surf = *rit;
        int surf_num = index_by_handle(surf);
        surfbcs[surf_num-1]=bc_id;
      }
    }
// Set graveyard to zero importance                                            
    else if ((grp_names[0].find("graveyard",0)==0)||(grp_names[0].find("outside_world",0)==0)
             ||(grp_names[0].find("rest_of_world",0)==0)) {
        // Get the volumes in the current group                                    
      MBRange grp_vols = grp_sets.intersect(vols_save);
      for (MBRange::iterator rit = grp_vols.begin(); rit != grp_vols.end(); rit++) {
        MBEntityHandle vol = *rit;
        int vol_num = index_by_handle(vol);
        imp[vol_num-1]=0;
      }
    }
// Check for tally groups, assign maximum CUBIT index                          
    else if (grp_names[0].find("tally_",0)==0) {
      std::string stal;
      int cpos = 6;
      while ((int(grp_names[0][cpos]) >= 48) && (int(grp_names[0][cpos]) <= 57)) {
        stal += grp_names[0][cpos];
        cpos++;
      }
      const char *ctal = stal.c_str();
      if (atoi(ctal) > maxital)
        maxital = atoi(ctal);
    }
  }
// Generate list of surfaces                                                   
  for (int ii=1; ii <= nsurfs; ii++) {
    MBEntityHandle surf = surf_handles()[ii];
    surflist[ii-1] = get_entity_id(surf);
//    cout << "surf#: " << ii << " = " << surflist[ii-1] << endl;              
  }

    // write in cell information                                                 
  for (int ii=1; ii <= ncells; ii++) {
    if (mats[ii-1] == 0) {
      cgmfile << id_by_index(3, ii) << " " << mats[ii-1] << " " <<
        "imp:n=" << imp[ii-1] << std::endl;
    }
    else {
      cgmfile << id_by_index(3, ii) << " " << mats[ii-1] << " " <<
        dens[ii-1] << " " << "imp:n=" << imp[ii-1] << std::endl;
    }
  }
    // skip a line                                                               
  cgmfile << std::endl;
    // write out surface information (need boundary conditions)                  
  for (int ii=1; ii <= nsurfs; ii++) {
    if (surfbcs[ii-1] == 1)
      cgmfile << "*" << surflist[ii-1] << std::endl;
    else if (surfbcs[ii-1] == 2)
      cgmfile << "+" << surflist[ii-1] << std::endl;
    else if (surfbcs[ii-1] < 0)
      cgmfile << surflist[ii-1] << " " << surfbcs[ii-1] << std::endl;
    else
      cgmfile << surflist[ii-1] << std::endl;
  }
    // add a final blank line                                                    
  cgmfile << " ";
    // if tallies were specified, write those out at end of .cgm file            
  if (maxital >= 0) {
    cgmfile << std::endl;
      // loop through each possible CUBIT tally index, write information           
    for (int ii=0; ii <= maxital; ii++) {
        // loop through each group name, check to see if the desired tally type    
      for (int jj=1; jj <= ngroups; jj++) {
        MBEntityHandle group = group_handles()[jj];
        grp_names.clear();
        bool success = get_group_names(group, grp_names);
        assert(success);
        if (!success || grp_names.empty()) continue;
        std::string talstr;
        talstr = "tally_" + itos(ii) + "_";
        const char *ctalstr = talstr.c_str();
        if (grp_names[0].find(ctalstr,0)==0) {
          unsigned int cpos = talstr.length();
          int etal = 0;
            // check if next character is "e" or "q" denoting energy or charge   
	  // need to also check if next character conforms to acceptable argument     
          if ((int(grp_names[0][cpos]) == 101) && 
              ((int(grp_names[0][cpos+1]) == 115) ||
               (int(grp_names[0][cpos+1]) == 99) || (int(grp_names[0][cpos+1]) == 112))) {
            etal = 1;
            cpos++;
          }
          else if ((int(grp_names[0][cpos]) == 113) && (int(grp_names[0][cpos+1]) == 112)) {
            etal = -1;
            cpos++;
          }
            // check to obtain tally type                                        
          int ital = 0;
          int tidx = 0;
          int tklen = 0;
          if (grp_names[0].find("surf_current",cpos)==cpos) {
            tidx = 1;
            ital = 10*ii + tidx;
            tklen = 13;
          }
          else if (grp_names[0].find("surf_flux",cpos)==cpos) {
            tidx = 2;
            ital = 10*ii + tidx;
            tklen = 10;
          }
          else if (grp_names[0].find("cell_flux",cpos)==cpos) {
            tidx = 4;
            ital = 10*ii + tidx;
            tklen = 10;
          }
          else if (grp_names[0].find("cell_heating",cpos)==cpos) {
            tidx = 6;
            ital = 10*ii + tidx;
            tklen = 13;
          }
          else if (grp_names[0].find("cell_fission",cpos)==cpos) {
            tidx = 7;
            ital = 10*ii + tidx;
            tklen = 13;
          }
          else if (grp_names[0].find("pulse_height",cpos)==cpos) {
            tidx = 8;
            ital = 10*ii + tidx;
            tklen = 13;
          }
          if (tidx > 0) {
            int ipos = 1;
              // place the tally header in the .cgm file                         
            if (etal == 1) {
              cgmfile << "*";
              ipos++;
            }
            else if ((etal == -1) && (tidx == 8)) {
              cgmfile << "+";
              ipos++;
            }
            cgmfile << "f" << ital << ":";
            std::stringstream spos;
            spos << ital;
            ipos += spos.str().length() + 2;
	    // get the particles the tally covers (default is n)                      
            cpos += tklen;
            std::string talptls;
            while((cpos < grp_names[0].length()) && (int(grp_names[0][cpos]) >= 65) &&
                  (int(grp_names[0][cpos]) <= 122) && (int(grp_names[0][cpos]) != 95)) {
              talptls += grp_names[0][cpos];
              cpos++;
            }
	    // read the tally particle string and place in .cgm file                  
            if (talptls.empty()) {
              cgmfile << "n" << " ";
              ipos += 2;
            }
            else {
              for(unsigned int kk=1; kk <= talptls.length(); kk++) {
                cgmfile << talptls[kk-1];
                if (kk < talptls.length())
                  cgmfile << ",";
                else
                  cgmfile << " ";
                ipos += 2*talptls.length();
              }
            }
	    // begin read of surfaces in current group                                
            MBRange grp_sets;
            MBErrorCode result = moab_instance()->get_entities_by_type(group, MBENTITYSET, grp_sets);
            if (MB_SUCCESS != result) continue;
            if ((tidx == 1) || (tidx == 2)) {
              MBRange grp_faces = grp_sets.intersect(surfs_save);
              for (MBRange::iterator rit = grp_faces.begin(); rit != grp_faces.end(); rit++) {
                int jtal = get_entity_id(*rit);
                std::stringstream spos1;
                spos1 << jtal;
                ipos += spos1.str().length() + 1;
                if (ipos > 78) {
                  cgmfile << "&" << std::endl << "     ";
                  ipos = spos1.str().length() + 7;
                }
                cgmfile << jtal << " ";
              }
            }
            else if ((tidx == 4) || ((tidx >= 6) && (tidx <= 8))) {
              MBRange grp_vols = grp_sets.intersect(vols_save);
              for (MBRange::iterator rit = grp_vols.begin(); rit != grp_vols.end(); rit++) {
                int jtal = get_entity_id(*rit);
                std::stringstream spos1;
                spos1 << jtal;
                ipos += spos1.str().length() + 1;
                if (ipos > 78) {
                  cgmfile << "&" << std::endl << "     ";
                  ipos = spos1.str().length() + 7;
                }
                cgmfile << jtal << " ";
              }
            }
	    // set default behavior to provide total over all tally regions
	    if (ipos+2 > 78) {
	      cgmfile << "&" << std::endl << "     ";
	    }
            cgmfile << " T" << std::endl;
          }
        }
      }
    }
  }
  cgmfile.close();
}

std::string DagMC::itos(int ival) {
  std::stringstream s;
  s << ival;
  return s.str();
}

bool DagMC::get_group_names(MBEntityHandle group_set, 
                             std::vector<std::string> &grp_names) 
{
    // get names
  char name0[NAME_TAG_SIZE];
  std::fill(name0, name0+NAME_TAG_SIZE, '\0');
  MBErrorCode result = moab_instance()->tag_get_data(name_tag(), &group_set, 1,
                                                     &name0);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return false;

  if (MB_TAG_NOT_FOUND != result) grp_names.push_back(std::string(name0));

  int extra_num = 0;
  while (true) {
    sprintf(name0, "%s%s%d", "EXTRA_", NAME_TAG_NAME, extra_num);
    extra_num++;
    MBTag this_tag = get_tag(name0, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, 
                             false);
    if (0 == this_tag) break;
    std::fill(name0, name0+NAME_TAG_SIZE, '\0');
    result = moab_instance()->tag_get_data(this_tag, &group_set, 1, &name0);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return false;
    if (MB_TAG_NOT_FOUND == result) break;
    else grp_names.push_back(std::string(name0));
  }
  
  return true;
}

MBTag DagMC::get_tag( const char* name, int size, MBTagType store, MBDataType type,
                       bool create_if_missing) 
{
  MBTag retval = 0;
  MBErrorCode result = moab_instance()->tag_create(name, size, store, type,
                                                   retval, NULL, create_if_missing);
  if (create_if_missing && MB_SUCCESS != result) 
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  
  return retval;
}

int DagMC::get_entity_id(MBEntityHandle this_ent) 
{
  int id = 0;
  MBErrorCode result = moab_instance()->tag_get_data(idTag, &this_ent, 1, &id);
  if (MB_TAG_NOT_FOUND == result)
    id = moab_instance()->id_from_handle(this_ent);
    
  return id;
}

MBErrorCode DagMC::get_angle(MBEntityHandle surf, 
                              double xxx, double yyy, double zzz, double *ang)
{
  MBEntityHandle root = rootSets[surf - setOffset];
  
  const double in_pt[] = { xxx, yyy, zzz };
  static std::vector<MBEntityHandle> facets;
  facets.clear();
  MBErrorCode rval = obbTree.closest_to_location( in_pt, root, tolerance(), facets );
  assert(MB_SUCCESS == rval);
  
  MBCartVect coords[3], normal(0.0);
  const MBEntityHandle *conn;
  int len;
  for (unsigned i = 0; i < facets.size(); ++i) {
    rval = mbImpl->get_connectivity( facets[i], conn, len );
    assert( MB_SUCCESS == rval );
    assert( 3 == len );
  
    rval = mbImpl->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == rval);
    
    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal += coords[1] * coords[2];
  }
  
  normal.normalize();
  normal.get( ang );

  return MB_SUCCESS;
}

MBErrorCode DagMC::CAD_ray_intersect(const double *point, 
                                      const double *dir, 
                                      const double huge,
                                      std::vector<double> &distances,
                                      std::vector<MBEntityHandle> &surfaces, 
                                      double &len) 
{
#ifdef CGM
#ifdef CUBIT_CGM
  return MB_NOT_IMPLEMENTED;
#else
  std::vector<double>::iterator dit = distances.begin();
  std::vector<MBEntityHandle>::iterator sit = surfaces.begin();
  static DLIList<double> ray_params;
  
  for (; dit != distances.end(); dit++, sit++) {
      // get the RefFace
    RefEntity *this_face = geomEntities[*sit - setOffset];
      // get the ray distance to this face
    ray_params.clean_out();
    int result = GeometryQueryTool::instance()->
      fire_ray(dynamic_cast<RefFace*>(this_face), CubitVector(point),
               CubitVector(dir), ray_params);
    assert(CUBIT_SUCCESS == result);
    if (ray_params.size() != 0) {
      ray_params.reset();
      *dit = ray_params.get();
    }
    else *dit = huge;
  }
  
    // now bubble sort list
  bool done = false;
  while (!done) {
    dit = distances.begin();
    sit = surfaces.begin();
    done = true;
    for (; dit != distances.end(); dit++, sit++) {
      if (dit+1 != distances.end() && *dit > *(dit+1)) {
        double tmp_dist = *dit;
        *dit = *(dit+1);
        *(dit+1) = tmp_dist;
        MBEntityHandle tmp_hand = *sit;
        *sit = *(sit+1);
        *(sit+1) = tmp_hand;
        done = false;
      }
    }
  }

  if (!distances.empty()) len = distances[0];
  
  return MB_SUCCESS;
#endif
#else
  return MB_FAILURE;
#endif
}

MBErrorCode DagMC::build_indices(MBRange &surfs, MBRange &vols,
                                  bool is_geom) 
{
  MBErrorCode rval = MB_SUCCESS;
  
    // surf/vol offsets are just first handles
  setOffset = (*surfs.begin() < *vols.begin() ? *surfs.begin() : *vols.begin());
  unsigned int tmp_offset = (surfs.back() > vols.back() ? surfs.back() : vols.back())
    - setOffset + 1;
  rootSets.resize(tmp_offset);
  entIndices.resize(rootSets.size());

    // store surf/vol handles lists (surf/vol by index) and
    // index by handle lists
  surf_handles().resize( surfs.size() + 1 );
  std::vector<MBEntityHandle>::iterator iter = entHandles[2].begin();
  *(iter++) = 0;
  std::copy( surfs.begin(), surfs.end(), iter );
  int idx = 1;
  for (MBRange::iterator rit = surfs.begin(); rit != surfs.end(); rit++)
    entIndices[*rit-setOffset] = idx++;
  
  entHandles[3].resize( vols.size() + 1 );
  iter = entHandles[3].begin();
  *(iter++) = 0;
  std::copy( vols.begin(), vols.end(), iter );
  idx = 1;
  for (MBRange::iterator rit = vols.begin(); rit != vols.end(); rit++)
    entIndices[*rit-setOffset] = idx++;

#ifdef CGM
  if (is_geom) {
    geomEntities.resize(rootSets.size());
      // get geometry entities by id and cache in this vector
    std::vector<int> ids;
    ids.resize(surfs.size());
    rval = moab_instance()->tag_get_data(id_tag(), surfs, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    int i = 0;
    MBRange::iterator rit = surfs.begin();
    for (; rit != surfs.end(); rit++, i++) {
      RefEntity *this_surf = GeometryQueryTool::instance()->
        get_ref_face(ids[i]);
      assert(NULL != this_surf);
      geomEntities[*rit - setOffset] = this_surf;
    }
    ids.resize(vols.size());
    rval = moab_instance()->tag_get_data(id_tag(), vols, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    i = 0;
    rit = vols.begin();
    for (; rit != vols.end(); rit++, i++) {
      RefEntity *this_vol = GeometryQueryTool::instance()->
        get_ref_volume(ids[i]);
      assert(NULL != this_vol);
      geomEntities[*rit - setOffset] = this_vol;
    }
  }
#endif  

    // get group handles
  MBTag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_LENGTH, 
                               MB_TAG_SPARSE, MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category+CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};
  MBRange groups;
  rval = moab_instance()->get_entities_by_type_and_tag(0, MBENTITYSET, &category_tag, 
                                                       group_val, 1, groups);
  if (MB_SUCCESS != rval)
    return rval;
  entHandles[4].resize(groups.size()+1);
  entHandles[4][0] = 0;
  std::copy(groups.begin(), groups.end(), &entHandles[4][1]);

    // populate root sets vector
  std::vector<MBEntityHandle> rsets;
  rsets.resize(surfs.size());
  rval = moab_instance()->tag_get_data(obb_tag(), surfs, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  MBRange::iterator rit;
  int i;
  for (i = 0, rit = surfs.begin(); rit != surfs.end(); rit++, i++)
    rootSets[*rit-setOffset] = rsets[i];

  rsets.resize(vols.size());
  rval = moab_instance()->tag_get_data(obb_tag(), vols, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  for (i = 0, rit = vols.begin(); rit != vols.end(); rit++, i++)
    rootSets[*rit-setOffset] = rsets[i];

  return MB_SUCCESS;
}

MBErrorCode DagMC::build_obbs(MBRange &surfs, MBRange &vols) 
{
  MBErrorCode rval = MB_SUCCESS;
  
  for (MBRange::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    MBEntityHandle root;
    MBRange tris;
    rval = moab_instance()->get_entities_by_dimension( *i, 2, tris );
    if (MB_SUCCESS != rval) 
      return rval;
    if (tris.empty()) 
      std::cerr << "WARNING: Surface " << get_entity_id(*i) << " has no facets." << std::endl;
    rval = obbTree.build( tris, root );
    if (MB_SUCCESS != rval) 
      return rval;
    rval = moab_instance()->add_entities( root, &*i, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    rval = moab_instance()->tag_set_data( obbTag, &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  for (MBRange::iterator i = vols.begin(); i != vols.end(); ++i) {
      // get all surfaces in volume
    MBRange tmp_surfs;
    rval = moab_instance()->get_child_meshsets( *i, tmp_surfs );
    if (MB_SUCCESS != rval)
      return rval;
    
      // get OBB trees for each surface
    MBEntityHandle root;
    MBRange trees;
    for (MBRange::iterator j = tmp_surfs.begin();  j != tmp_surfs.end(); ++j) {
        // skip any surfaces that are non-manifold in the volume
        // because point containment code will get confused by them
      int sense;
      rval = surface_sense( *i, 1, &*j, &sense );
      if (MB_SUCCESS != rval) {
        std::cerr << "Surface/Volume sense data missing." << std::endl;
        return rval;
      }
      if (!sense)
        continue;
      
      rval = moab_instance()->tag_get_data( obbTag, &*j, 1, &root );
      if (MB_SUCCESS != rval || !root)  
        return MB_FAILURE;
      trees.insert( root );
    }
    
      // build OBB tree for volume
    rval = obbTree.join_trees( trees, root );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = moab_instance()->tag_set_data( obbTag, &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;

  }

  return MB_SUCCESS;
}
