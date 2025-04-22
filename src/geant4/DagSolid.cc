//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DagSolid.hh,v 1.10 2010/12/10 16:30:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-05 $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              DagSolid.hh
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
// 27 March 2014, A Davis, UW - Updated description text
// 31 October 2010, J. H. Jeong, Hanyang Univ., KR
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "DagSolid.hh"

#include <iostream>

#include "DagMC.hpp"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
using namespace moab;

#include "DagSolid.hh"

#define plot true
#define debug false

//#define DAGDEBUG 

// #define G4SPECSDEBUG 1
///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no facets.
//
DagSolid::DagSolid()
    : G4TessellatedSolid("dummy"), cubicVolume(0.), surfaceArea(0.) {
  geometryType = "DagSolid";

  xMinExtent = kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent = kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent = kInfinity;
  zMaxExtent = -kInfinity;
}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no facets
// to define.
//
DagSolid::DagSolid(const G4String& name, DagMC* dagmc, int volID)
    : G4TessellatedSolid(name), cubicVolume(0.), surfaceArea(0.) {
  geometryType = "DagSolid";
  Myname = name;

  fdagmc = dagmc;
  fvolID = volID;
  fvolEntity = fdagmc->entity_by_index(3, volID);

  double min[3], max[3];
  fdagmc->getobb(fvolEntity, min, max);

  xMinExtent = min[0];
  xMaxExtent = max[0];
  yMinExtent = min[1];
  yMaxExtent = max[1];
  zMinExtent = min[2];
  zMaxExtent = max[2];

  int num_entities;
  std::vector<EntityHandle> surfs;
  std::vector<EntityHandle> tris;
  const EntityHandle* tri_conn;

  std::vector<CartVect> coords(3);
  G4ThreeVector vertex[3];
  int n_verts = 0;

  Interface* moab = dagmc->moab_instance();
  moab->get_child_meshsets(fvolEntity, surfs, 1);
  int sense = 0;
  if (plot) {
    for (unsigned i = 0; i < surfs.size(); i++) {
      My_sulf_hit = surfs[i];
      moab->get_number_entities_by_type(surfs[i], MBTRI, num_entities);
      moab->get_entities_by_type(surfs[i], MBTRI, tris);
      dagmc->surface_sense(fvolEntity,surfs[i],sense);
      for (unsigned j = 0; j < tris.size(); j++) {
        moab->get_connectivity(tris[j], tri_conn, n_verts);
        moab->get_coords(tri_conn, n_verts, coords[0].array());

	vertex[0] = G4ThreeVector(coords[0][0] * cm, coords[0][1] * cm,
				  coords[0][2] * cm);
	vertex[1] = G4ThreeVector(coords[1][0] * cm, coords[1][1] * cm,
				  coords[1][2] * cm);
	vertex[2] = G4ThreeVector(coords[2][0] * cm, coords[2][1] * cm,
				  coords[2][2] * cm);

	G4TriangularFacet* facet = NULL;
	if ( sense > 0 ) {
	  facet = new G4TriangularFacet(vertex[0], vertex[1], vertex[2], ABSOLUTE);
	} else {
	  facet = new G4TriangularFacet(vertex[2], vertex[1], vertex[0], ABSOLUTE);
	}
        AddFacet((G4VFacet*)facet);

        for (G4int k = 0; k < 3; k++) {
          if (vertex[k].x() < xMinExtent) xMinExtent = vertex[k].x();
          if (vertex[k].x() > xMaxExtent) xMaxExtent = vertex[k].x();
          if (vertex[k].y() < yMinExtent) yMinExtent = vertex[k].y();
          if (vertex[k].y() > yMaxExtent) yMaxExtent = vertex[k].y();
          if (vertex[k].z() < zMinExtent) zMinExtent = vertex[k].z();
          if (vertex[k].z() > zMaxExtent) zMaxExtent = vertex[k].z();
        }
      }
      tris.clear();
    }
  }

  SetSolidClosed(true);
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
DagSolid::DagSolid(__void__& a)
    : G4TessellatedSolid(a),
      geometryType("DagSolid"),
      cubicVolume(0.),
      surfaceArea(0.),
      xMinExtent(0.),
      xMaxExtent(0.),
      yMinExtent(0.),
      yMaxExtent(0.),
      zMinExtent(0.),
      zMaxExtent(0.) {}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
DagSolid::~DagSolid() {}

///////////////////////////////////////////////////////////////////////////////
//
// EInside DagSolid::Inside (const G4ThreeVector &p) const
//
EInside DagSolid::Inside(const G4ThreeVector& p) const {
  G4double point[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm

  G4double minDist = 0.0;

  int inside;
  ErrorCode ec;

  // are we inside
  fdagmc->point_in_volume(fvolEntity, point, inside);
  // also need to know how far from the surface
  fdagmc->closest_to_location(fvolEntity, point, minDist);

  // if on surface
  EInside result;
  if (minDist*cm <= kCarToleranceHalf) {
    result = kSurface;
  } else {
    if (inside == 0) {
      result = kOutside;
    } else {
      result = kInside;
    }
  }

  #ifdef DAGDEBUGA
    G4cout << "<<<<<" << G4endl;
    G4cout << "Inside(p)" << G4endl;
    G4cout << "Name: " << Myname << G4endl;
    G4cout << "pos: " << p << G4endl;
    G4cout << "return: " << result << G4endl;
    G4cout << ">>>>>" << G4endl;
  #endif
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4ThreeVector DagSolid::SurfaceNormal (const G4ThreeVector &p) const
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.

G4ThreeVector DagSolid::SurfaceNormal(const G4ThreeVector& p) const {
  G4double ang[3] = {0, 0, 1};
  G4double position[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm
  G4double distance;
  EntityHandle surface;
  fdagmc->closest_to_location(fvolEntity, position, distance, &surface);

  // currently get warnings from RTI.cpp in double down since the
  // point may not be on the surface 
  fdagmc->get_angle(surface, position, ang);

  G4ThreeVector normal = G4ThreeVector(ang[0], ang[1], ang[2]);

  #ifdef DAGDEBUG
    G4cout << "<<<<<" << G4endl;
    G4cout << "SurfaceNormal(p)" << G4endl;
    G4cout << "Name: " << Myname << G4endl;
    G4cout << "normal: " << normal << G4endl;
    G4cout << ">>>>>" << G4endl;
  #endif

  return normal;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)
//
//
G4double DagSolid::DistanceToIn(const G4ThreeVector& p,
                                const G4ThreeVector& v) const {
  G4double minDist = kInfinity;
  G4double position[3] = {p.x() / cm, p.y() / cm, p.z() / cm};  // convert to cm
  G4ThreeVector vec = v.unit();
  G4double dir[3] = {vec.x(), vec.y(), vec.z()};
  EntityHandle next_surf = 0;
  G4double distance;

  // look for an entering intesection, i.e. opposing the ray direction
  fdagmc->ray_fire(fvolEntity, position, dir, next_surf, distance, NULL, 0, -1);
  
  // if we are close the surface set surfaace distance 0
  if (distance*cm <= kCarToleranceHalf) {
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
    G4cout << "distance: " << distance*cm << G4endl;
    G4cout << ">>>>>" << G4endl;
  #endif

  return distance*cm;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an outside point p.
// The distance can be an underestimate
//
/////////////////////////////////////////////////////////////////////////
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
    G4cout << "return: " << minDist*cm << G4endl;
    G4cout << ">>>>>" << G4endl;
  #endif

  return minDist*cm;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
//                        const G4bool calcNorm=false,
//                        G4bool *validNorm=0, G4ThreeVector *n=0);
//
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

  double position[3] = {p.x()/cm, p.y()/cm, p.z()/cm};  // convert position to cm
  G4ThreeVector vec = v.unit();
  double dir[3] = {vec.x(), vec.y(), vec.z()};

  EntityHandle surface;
  G4double distance;
  moab::DagMC::RayHistory history;
    
  // fire a ray
  fdagmc->ray_fire(fvolEntity,position,dir,surface,distance,&history,0,1);

  // if we are asked to calculate the normal
  if (calcNorm) {
    *validNorm = true;
    // get the direction vector;
    double normal[3];
    double hit[3];
    hit[0] = position[0] + (distance*dir[0]);
    hit[1] = position[1] + (distance*dir[1]);
    hit[2] = position[2] + (distance*dir[2]);

    fdagmc->get_angle(surface, hit, normal, &history);
    *n = G4ThreeVector(normal[0],normal[1],normal[2]);
    G4int sense;
    fdagmc->surface_sense(fvolEntity, surface, sense);
    if (sense == -1) *n = -(*n);
  
  }

  // if hit is too close
   if ( distance*cm <= kCarToleranceHalf ) {
    distance = 0.;
  }

  #ifdef DAGDEBUG
    G4cout << "<<<<<" << G4endl;
    G4cout << "DistanceToOut(p,v) " << G4endl;
    G4cout << "Name: " << Myname << G4endl;
    G4cout << "point: " << p << G4endl;
    G4cout << "direction: " << v << G4endl;
    G4cout << "distance: " << distance*cm << G4endl;
    G4cout << "calcNorm: " << calcNorm << G4endl;
    if (calcNorm) {
      G4cout << "validNorm: " << *validNorm << G4endl;
      G4cout << "normal: " << n->x() << " " << n->y() << " " << n->z() << G4endl;
    }
    G4cout << ">>>>>" << G4endl;
  #endif 

  // return the hit distance
  return distance*cm;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
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
    G4cout << "return: " << minDist*cm << G4endl;
    G4cout << ">>>>>" << G4endl;
  #endif
  
  return minDist*cm;
  // return DistanceToIn(p);
}

///////////////////////////////////////////////////////////////////////////////
//
// G4GeometryType GetEntityType() const;
//
// Provide identification of the class of an object (required for persistency
// and STEP interface).
//
G4GeometryType DagSolid::GetEntityType() const { return geometryType; }

///////////////////////////////////////////////////////////////////////////////
//
void DagSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
  scene.AddSolid(*this);
}

std::ostream& DagSolid::StreamInfo(std::ostream& os) const {
  os << G4endl;
  os << "Geometry Type    = " << geometryType << G4endl;
  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
// CalculateExtent
//
// Based on correction provided by Stan Seibert, University of Texas.
//
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

G4double DagSolid::GetMinXExtent() const { return xMinExtent / cm; }

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxXExtent() const { return xMaxExtent / cm; }

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinYExtent() const { return yMinExtent / cm; }

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxYExtent() const { return yMaxExtent / cm; }

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinZExtent() const { return zMinExtent / cm; }

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxZExtent() const { return zMaxExtent / cm; }

/*
 * Return the volume of the volume - note DAGMC is always in the units of cm
 */
G4double DagSolid::GetCubicVolume() {
  G4double result;
  fdagmc->measure_volume(fvolEntity, result);
  return result * cm * cm * cm;
}

/*
 * Return the surface area of the volume - 
 * note DAGMC is always in the units of cm
 */
G4double DagSolid::GetSurfaceArea() {
  G4double result;
  // get the moab instance
  moab::Interface *moab = fdagmc->moab_instance();
  std::vector<moab::EntityHandle> surfaces;
  moab::ErrorCode rval = moab->get_child_meshsets(fvolEntity,surfaces);

  G4double area = 0.;
  for (auto surface : surfaces) {
    G4double surf_area;
    fdagmc->measure_area(surface,surf_area);
    area += surf_area;
  } 
  return area * cm * cm;
}
