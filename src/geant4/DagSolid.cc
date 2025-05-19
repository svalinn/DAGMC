//
// Date:                27/03/2014
// Author:              A. Davis M. C. Han, C. H. Kim, J. H. Jeong, Y. S. Yeom,
// S.
//                      Kim, Paul. P. H. Wilson, J. Apostolakis
// Organisation:        The University of Wisconsin-Madison, USA & Hanyang
// Univ., KR
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 19th May 2025, A Davis - version 2.0
// 27 March 2014, A Davis, UW - Updated description text
// 31 October 2010, J. H. Jeong, Hanyang Univ., KR
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "DagSolid.hh"

#include <iostream>

#include "DagMC.hpp"
#include "DagSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "G4VisExtent.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

// #define DAGDEBUG 1

// Constructor for empty DagSolid
DagSolid::DagSolid() : G4VSolid("dummy"), cubicVolume(0.), surfaceArea(0.) {
  geometryType = "DagSolid";

  xMinExtent = kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent = kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent = kInfinity;
  zMaxExtent = -kInfinity;
}

// Main Constructor
DagSolid::DagSolid(const G4String& name, moab::DagMC* dagmc, int volID)
    : G4VSolid(name), cubicVolume(0.), surfaceArea(0.) {
  geometryType = "DagSolid";

  Myname = name;   // set the name
  fdagmc = dagmc;  // set the dagmc instance
  // get the moab instance from the DAGMC pointer
  moab = fdagmc->moab_instance();

  fvolID = volID;
  fvolEntity = fdagmc->entity_by_index(3, volID);

  // get extents from the bounding box
  double min[3], max[3];
  fdagmc->getobb(fvolEntity, min, max);

  xMinExtent = min[0];
  xMaxExtent = max[0];
  yMinExtent = min[1];
  yMaxExtent = max[1];
  zMinExtent = min[2];
  zMaxExtent = max[2];

  // if surfaces have not been instantiated yet, do so
  if (fSurfaces.size() == 0) {
    moab->get_child_meshsets(fvolEntity, fSurfaces, 1);
  }

  // cache entity handles of triangles for later
  for (moab::EntityHandle surf : fSurfaces) {
    G4cout << fvolEntity << " " << surf << G4endl;

    std::vector<moab::EntityHandle> triangles;
    moab->get_entities_by_type(surf, moab::MBTRI, triangles);
    fTriangles.insert(fTriangles.end(), triangles.begin(), triangles.end());
  }
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
DagSolid::DagSolid(__void__& a)
    : G4VSolid(a),
      geometryType("DagSolid"),
      cubicVolume(0.),
      surfaceArea(0.),
      xMinExtent(0.),
      xMaxExtent(0.),
      yMinExtent(0.),
      yMaxExtent(0.),
      zMinExtent(0.),
      zMaxExtent(0.) {}

// Destructor
DagSolid::~DagSolid() {
  if (fPolyhedron) delete fPolyhedron;
}

// Determine if the point p is inside the solid or not
EInside DagSolid::Inside(const G4ThreeVector& p) const {
  G4double point[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm

  G4double minDist = 0.0;

  int inside;  // test to see if the point is inside the volume

  // are we inside
  fdagmc->point_in_volume(fvolEntity, point, inside);
  // also need to know how far from the surface
  fdagmc->closest_to_location(fvolEntity, point, minDist);

  // if on surface
  EInside result;
  if (minDist * cm <= kCarToleranceHalf) {
    result = kSurface;
  } else {
    if (inside == 0) {
      result = kOutside;
    } else {
      result = kInside;
    }
  }

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "Inside(p)" << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "pos: " << p << G4endl;
  G4cout << "return: " << result << G4endl;
  G4cout << ">>>>>" << G4endl;
#endif
  return result;
}

// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.
G4ThreeVector DagSolid::SurfaceNormal(const G4ThreeVector& p) const {
  G4double ang[3] = {0, 0, 1};
  G4double position[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm
  G4double distance;
  moab::EntityHandle surface;

  // find out the nearest surface
  fdagmc->closest_to_location(fvolEntity, position, distance, &surface);

  // now figure out the normal
  // TODO check expectation about sign flips for shared surfaces
  moab::ErrorCode rval = fdagmc->get_angle(surface, position, ang);

  // backup estimate of normal - a la G4TessellatedSolid
  if (rval != moab::MB_SUCCESS) {
    return (p.z() > 0 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(0, 0, -1));
  } else {
    // return the normal
    G4ThreeVector normal = G4ThreeVector(ang[0], ang[1], ang[2]);

    G4int sense;  // sense is either +1 or -1
    fdagmc->surface_sense(fvolEntity, surface, sense);
    if (sense == -1) normal = -normal;
    return normal;
  }

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "SurfaceNormal(p)" << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "normal: " << normal << G4endl;
  G4cout << ">>>>>" << G4endl;
#endif
}

// Given a point, p and a vector, v, determine the distance
// from the point oustide the volume until we enter
G4double DagSolid::DistanceToIn(const G4ThreeVector& p,
                                const G4ThreeVector& v) const {
  G4double position[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm
  G4ThreeVector vec = v.unit();
  G4double dir[3] = {vec.x(), vec.y(), vec.z()};
  moab::EntityHandle next_surf = 0;
  G4double distance;

  // look for an entering intesection, i.e. opposing the ray direction
  fdagmc->ray_fire(fvolEntity, position, dir, next_surf, distance, NULL, 0, -1);

  // if we are close the surface set surfaace distance 0
  if (distance * cm <= kCarToleranceHalf) {
    distance = 0.;
  }

  // no hits
  if (next_surf == 0) return kInfinity;

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "DistanceToIn(p,v) " << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "pos: " << p << G4endl;
  G4cout << "direction: " << v << G4endl;
  G4cout << "distance: " << distance * cm << G4endl;
  G4cout << ">>>>>" << G4endl;
#endif

  return distance * cm;
}

// Calculate distance to nearest surface of shape from an outside point p.
// The distance can be an underestimate
G4double DagSolid::DistanceToIn(const G4ThreeVector& p) const {
  G4double minDist = kInfinity;
  G4double point[3] = {p.x() / cm, p.y() / cm,
                       p.z() / cm};  // convert position to cm

  fdagmc->closest_to_location(fvolEntity, point, minDist);

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "DistanceToIn(p)" << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "pos: " << p << G4endl;
  G4cout << "return: " << minDist * cm << G4endl;
  G4cout << ">>>>>" << G4endl;
#endif

  return minDist * cm;
}

// Return distance along the normalised vector v to the shape, from a
// point at an offset p inside or on the surface of the shape.
// Intersections with surfaces, when the point is not greater
// than kCarTolerance/2 from a surface, must be ignored.
//     If calcNorm is true, then it must also set validNorm to either
//     * true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//     * false, if the solid does not lie entirely behind or on the
//       exiting surface.
// If calcNorm is false, then validNorm and n are unused.
// distqance to out - when inside - if the distaance is within halfkcartolerance
// then return 0, if aasked to calculate the normal then
// set validnorm true and set the normal
G4double DagSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                 const G4bool calcNorm, G4bool* validNorm,
                                 G4ThreeVector* n) const {
  double position[3] = {p.x() / cm, p.y() / cm,
                        p.z() / cm};  // convert position to cm
  G4ThreeVector vec = v.unit();
  double dir[3] = {vec.x(), vec.y(), vec.z()};

  moab::EntityHandle surface;
  G4double distance;
  moab::DagMC::RayHistory history;

  // fire a ray
  fdagmc->ray_fire(fvolEntity, position, dir, surface, distance, &history, 0,
                   1);

  // if we hit nothing
  if (!surface) distance = kInfinity;

  // if we are asked to calculate the normal
  if (calcNorm) {
    *validNorm = true;
    // get the direction vector;
    double normal[3];
    double hit[3];

    hit[0] = position[0] + (distance * dir[0]);
    hit[1] = position[1] + (distance * dir[1]);
    hit[2] = position[2] + (distance * dir[2]);

    // the the angle
    moab::ErrorCode rval = fdagmc->get_angle(surface, hit, normal, &history);
    if (rval != moab::MB_SUCCESS) {
      *n = (p.z() > 0 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(0, 0, -1));
    } else {
      *n = G4ThreeVector(normal[0], normal[1], normal[2]);
    }

    G4int sense;
    // get the surface sense
    fdagmc->surface_sense(fvolEntity, surface, sense);
    if (sense == -1) *n = -(*n);
  }

  // if hit is too close
  if (distance * cm <= kCarToleranceHalf) {
    distance = 0.;
  }

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "DistanceToOut(p,v) " << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "point: " << p << G4endl;
  G4cout << "direction: " << v << G4endl;
  G4cout << "distance: " << distance * cm << G4endl;
  G4cout << "calcNorm: " << calcNorm << G4endl;
  if (calcNorm) {
    G4cout << "validNorm: " << *validNorm << G4endl;
    G4cout << "normal: " << n->x() << " " << n->y() << " " << n->z() << G4endl;
  }
  G4cout << ">>>>>" << G4endl;
#endif

  // return the hit distance
  return distance * cm;
}

// Calculate distance to nearest surface of shape from an inside point.
// The distance can be an underestimate.
G4double DagSolid::DistanceToOut(const G4ThreeVector& p) const {
  // same as DistanceToIn - calculate safety
  G4double minDist = kInfinity;
  G4double point[3] = {p.x() / cm, p.y() / cm,
                       p.z() / cm};  // convert position to cm

  fdagmc->closest_to_location(fvolEntity, point, minDist);

#ifdef DAGDEBUG
  G4cout << "<<<<<" << G4endl;
  G4cout << "DistanceToOut(p)" << G4endl;
  G4cout << "Name: " << Myname << G4endl;
  G4cout << "pos: " << p << G4endl;
  G4cout << "return: " << minDist * cm << G4endl;
  G4cout << ">>>>>" << G4endl;
#endif

  return minDist * cm;
  // return DistanceToIn(p);
}

// Return the entity type
G4GeometryType DagSolid::GetEntityType() const { return geometryType; }

// for visualiation adds in our case the polyhedron
void DagSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
  scene.AddSolid(*this);
}

// string based description of the solid
std::ostream& DagSolid::StreamInfo(std::ostream& os) const {
  os << G4endl;
  os << "Geometry Type    = " << geometryType << G4endl;
  return os;
}

// Based on correction provided by Stan Seibert, University of Texas.
G4bool DagSolid::CalculateExtent(const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const {
  // Calculate the minimum and maximum extent of the solid,
  // when under the specified transform, and within the specified limits.
  // If the solid is not intersected by the region, return false, else return
  // true.

  G4ThreeVector minExtent(xMinExtent, yMinExtent, zMinExtent);
  G4ThreeVector maxExtent(xMaxExtent, yMaxExtent, zMaxExtent);

  // Check for containment and clamp to voxel boundaries
  for (G4int axis = G4ThreeVector::X; axis < G4ThreeVector::SIZE; axis++) {
    EAxis geomAxis = kXAxis;  // G4 geom classes use different index type
    switch (axis) {
      case G4ThreeVector::X:
        geomAxis = kXAxis;
        break;
      case G4ThreeVector::Y:
        geomAxis = kYAxis;
        break;
      case G4ThreeVector::Z:
        geomAxis = kZAxis;
        break;
    }
    G4bool isLimited = pVoxelLimit.IsLimited(geomAxis);
    G4double voxelMinExtent = pVoxelLimit.GetMinExtent(geomAxis);
    G4double voxelMaxExtent = pVoxelLimit.GetMaxExtent(geomAxis);

    if (isLimited) {
      if (minExtent[axis] > voxelMaxExtent + kCarTolerance ||
          maxExtent[axis] < voxelMinExtent - kCarTolerance) {
        return false;
      } else {
        if (minExtent[axis] < voxelMinExtent) {
          minExtent[axis] = voxelMinExtent;
        }
        if (maxExtent[axis] > voxelMaxExtent) {
          maxExtent[axis] = voxelMaxExtent;
        }
      }
    }
  }

  // Convert pAxis into G4ThreeVector index
  G4int vecAxis = 0;
  switch (pAxis) {
    case kXAxis:
      vecAxis = G4ThreeVector::X;
      break;
    case kYAxis:
      vecAxis = G4ThreeVector::Y;
      break;
    case kZAxis:
      vecAxis = G4ThreeVector::Z;
      break;
    default:
      break;
  }

  pMin = minExtent[vecAxis] - kCarTolerance;
  pMax = maxExtent[vecAxis] + kCarTolerance;

  return true;
}

// sampling the solid, get a random point on the surface
// used during visualisation if GetPolyhedron fails
G4ThreeVector DagSolid::GetPointOnSurface() const {
  // number of triangles
  G4int num_tri = fTriangles.size();
  // sampled triangle
  G4int tri_idx = (G4int)G4RandFlat::shoot(0., num_tri);

  // get the triangle connectivity
  const moab::EntityHandle* tri_conn;
  G4int n_verts;
  moab->get_connectivity(fTriangles[tri_idx], tri_conn, n_verts);
  // get the spatial coordinates
  std::vector<moab::CartVect> coords(3);
  moab->get_coords(tri_conn, n_verts, coords[0].array());

  // get the random point on the triangle
  G4double a = G4RandFlat::shoot(0., 1.);
  G4double b = G4RandFlat::shoot(0., 1.);

  //
  G4ThreeVector x1, x2, x3;

  // vertices
  x1 = G4ThreeVector(coords[0][0] * cm, coords[0][1] * cm, coords[0][2] * cm);
  x2 = G4ThreeVector(coords[1][0] * cm, coords[1][1] * cm, coords[1][2] * cm);
  x3 = G4ThreeVector(coords[2][0] * cm, coords[2][1] * cm, coords[2][2] * cm);

  // the sample point
  G4ThreeVector samplePoint;

  samplePoint = (1 - std::sqrt(a)) * x1 + (std::sqrt(a) * (1 - b)) * x2 +
                (std::sqrt(a) * b) * x3;

  // return the point
  return samplePoint;
}

// Create a polyhedron for visualising the the solid in the
// geant4 viewer
G4Polyhedron* DagSolid::CreatePolyhedron() const {
  // Create a polyhedron from the facets and vertices of the solid.
  // get the moab instance

  std::unordered_set<moab::EntityHandle> facetSet;
  // loop over the surfaces - we need to get all the triangles
  for (unsigned i = 0; i < fSurfaces.size(); i++) {
    // get the triangle entities
    std::vector<moab::EntityHandle> facets;
    // get the triangles
    moab->get_entities_by_type(fSurfaces[i], moab::MBTRI, facets);
    // insert them for later
    // facetCollection.insert(facetCollection.end(), facets.begin(),
    // facets.end());
    facetSet.insert(facets.begin(), facets.end());
  }
  // make the collection removing duplicates
  std::vector<moab::EntityHandle> facetCollection(facetSet.begin(),
                                                  facetSet.end());

  // get the vertices used by the facets
  std::vector<moab::EntityHandle> vertexCollection;
  moab->get_connectivity(&facetCollection[0], facetCollection.size(),
                         vertexCollection);

  // remove duplicates
  std::unordered_set<moab::EntityHandle> vertexSet(vertexCollection.begin(),
                                                   vertexCollection.end());

  // clearout the vector
  vertexCollection.clear();
  // insert the set
  vertexCollection.assign(vertexSet.begin(), vertexSet.end());

  // instanciate the polyhedron
  G4int nFacets = facetCollection.size();
  G4int nVertices = vertexCollection.size();

  #ifdef GEANT4_GT_11
  // make a new polyhedron container
  auto polyhedron = new G4Polyhedron(nVertices, nFacets);
  // now loop over the surfaces creating the facets
  // using the threevector versions

  // create index maps
  std::map<moab::EntityHandle, int> g4VertexMap;

  // create a map of ThreeVectors for the polyhedron
  for (int i = 0; i < nVertices; i++) {
    moab::EntityHandle vertex = vertexCollection[i];
    // get the coordinates
    G4double coords[3];
    // get the coordinates of the vertex
    moab->get_coords(&vertex, 1, coords);

    // create the G4 vertex
    G4ThreeVector g4vertex =
        G4ThreeVector(coords[0] * cm, coords[1] * cm, coords[2] * cm);

    polyhedron->SetVertex(i + 1, g4vertex);
    g4VertexMap[vertex] = i + 1;
  }

  G4int facet_index = 0;
  // loop over the surfaces
  for (unsigned i = 0; i < fSurfaces.size(); i++) {
    // get the triangle entities
    std::vector<moab::EntityHandle> facets;

    // get the triangles on the surface
    moab->get_entities_by_type(fSurfaces[i], moab::MBTRI, facets);
    G4int surfaceSense;
    fdagmc->surface_sense(fvolEntity, fSurfaces[i], surfaceSense);

    // for each triangle
    for (moab::EntityHandle facet : facets) {
      // get the connectivity of the triangle
      const moab::EntityHandle* tri_conn;
      G4int n_verts;
      moab->get_connectivity(facet, tri_conn, n_verts);

      // get each vertex
      G4int vertex[3];
      for (int j = 0; j < n_verts; j++) {
        vertex[j] = g4VertexMap[tri_conn[j]];
      }
      // get the sense in case we need to flip the triangle
      if (surfaceSense > 0) {
        polyhedron->SetFacet(facet_index + 1, vertex[0], vertex[1], vertex[2],
                             0);
      } else {
        polyhedron->SetFacet(facet_index + 1, vertex[2], vertex[1], vertex[0],
                             0);
      }
      facet_index++;
    }
  }
  //  finalise the polyhedron
  # else
  // make a HepPolyhedron first, then make a Polyhedron from it

  //std::vector<std::array<G4double,3>> xyz_arr;
  //xyz_arr.reserve(nVertices);
  G4double xyz_arr[nVertices][3];

  std::map<moab::EntityHandle,G4int> vertex_lookup;
  // create a map of ThreeVectors for the polyhedron
  for (int i = 0; i < nVertices; i++) {
    moab::EntityHandle vertex = vertexCollection[i];
    // get the coordinates
    G4double coords[3];
    // get the coordinates of the vertex
    moab->get_coords(&vertex, 1, coords);

    // create the G4 vertex
    G4double xyz[3] = {coords[0]*cm, coords[1]*cm, coords[2]*cm};
    //xyz_arr[i] = xyz;
    xyz_arr[i][0] = xyz[0];
    xyz_arr[i][1] = xyz[1];
    xyz_arr[i][2] = xyz[2];
    
    vertex_lookup[vertex] = i+1; // note uses a 1 based indexing
  }

  //std::vector<std::array<G4int,4>> faces_arr;
  //faces_arr.reserve(nFacets);
  G4int faces_arr[nFacets][4];
  // loop over the surfaces
  G4int facet_idx = 0;
  for (unsigned i = 0; i < fSurfaces.size(); i++) {
    // get the triangle entities
    std::vector<moab::EntityHandle> facets;

    // get the triangles on the surface
    moab->get_entities_by_type(fSurfaces[i], moab::MBTRI, facets);
    G4int surfaceSense;
    fdagmc->surface_sense(fvolEntity, fSurfaces[i], surfaceSense);

    // for each triangle
    for (moab::EntityHandle facet : facets) {
      // get the connectivity of the triangle
      const moab::EntityHandle* tri_conn;
      G4int n_verts;
      moab->get_connectivity(facet, tri_conn, n_verts);
      
      // get each vertex
      G4int vertex[3];
      for (int j = 0; j < n_verts; j++) {
        if(surfaceSense > 0) {
          faces_arr[facet_idx][j] = vertex_lookup[tri_conn[j]];
        } else {
          G4int i = 2-j;
          faces_arr[facet_idx][j] = vertex_lookup[tri_conn[i]];
        }
      }
      faces_arr[facet_idx][3] = 0;
      //faces_arr[i] = {vertex[0],vertex[1],vertex[2],0};
      facet_idx++;
    }
  }
  HepPolyhedron *poly = new HepPolyhedron();
  poly->createPolyhedron(nVertices, nFacets,
			 xyz_arr, faces_arr);
  G4Polyhedron* polyhedron = new G4Polyhedron(*poly);
  #endif
  return polyhedron;
}

// get the Polyhedron if it doesnt exist create it
G4Polyhedron* DagSolid::GetPolyhedron() const {
  if (fPolyhedron == nullptr) fPolyhedron = CreatePolyhedron();
  return fPolyhedron;
}

// get the Cartesian extent of the solid
G4VisExtent DagSolid::GetExtent() const {
  // Define the sides of the box into which the G4Tubs instance would fit.
  return G4VisExtent(GetMinXExtent(), GetMaxXExtent(), GetMinYExtent(),
                     GetMaxYExtent(), GetMinZExtent(), GetMaxZExtent());
}

// get the bounding limits of the solid
void DagSolid::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const {
  pMin = G4ThreeVector(GetMinXExtent(), GetMinYExtent(), GetMinZExtent());
  pMax = G4ThreeVector(GetMaxXExtent(), GetMaxYExtent(), GetMaxZExtent());
}

// get the minimum extent in the x direction
G4double DagSolid::GetMinXExtent() const { return xMinExtent * cm; }

// get the maximum extent in the x direction
G4double DagSolid::GetMaxXExtent() const { return xMaxExtent * cm; }

// get the minimum extent in the y direction
G4double DagSolid::GetMinYExtent() const { return yMinExtent * cm; }

// get the maximum extent in the y direction
G4double DagSolid::GetMaxYExtent() const { return yMaxExtent * cm; }

// get the minimum extent in the z direction
G4double DagSolid::GetMinZExtent() const { return zMinExtent * cm; }

// get the maximum extent in the z direction
G4double DagSolid::GetMaxZExtent() const { return zMaxExtent * cm; }

// Return the volume of the volume
// note DAGMC is always in the units of cm
G4double DagSolid::GetCubicVolume() {
  G4double result;
  fdagmc->measure_volume(fvolEntity, result);
  return result * cm * cm * cm;
}

// Return the surface area of the volume -
// note DAGMC is always in the units of cm
G4double DagSolid::GetSurfaceArea() {
  // get the moab instance
  std::vector<moab::EntityHandle> surfaces;
  moab::ErrorCode rval = moab->get_child_meshsets(fvolEntity, surfaces);

  G4double area = 0.;
  for (auto surface : surfaces) {
    G4double surf_area;
    fdagmc->measure_area(surface, surf_area);
    area += surf_area;
  }
  return area * cm * cm;
}
