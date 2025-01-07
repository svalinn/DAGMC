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

  int result;
  ErrorCode ec;

  if (debug) {
    std::cout << "Inside: " << std::endl;
  }
  
  // are we inside
  fdagmc->point_in_volume(fvolEntity, point, result);
  // also need to know how far from the surface
  fdagmc->closest_to_location(fvolEntity, point, minDist);

  // if on surface
  if (minDist <= kCarToleranceHalf) {
    if (debug) {
      std::cout << "dist: " << minDist << std::endl;
      std::cout << "inside: " << result << std::endl;
      std::cout << "return: (onsurf) " << kSurface << std::endl;
    }
    return kSurface;
  } else {
    // result == 0 is outside 
    if (result == 0 || minDist > kCarToleranceHalf) {
      if (debug) {
	std::cout << "dist: " << minDist << std::endl;
	std::cout << "inside: " << result << std::endl;
	std::cout << "return: (outside) " << kOutside << std::endl;
      }
      return kOutside;
    } else {
      if (debug) {
	std::cout << "dist: " << minDist << std::endl;
	std::cout << "inside: " << result << std::endl;
	std::cout << "return: (inside) " << kInside << std::endl;
      }
      return kInside;
    }
  }
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
  
  // currently get warnings from RTI.cpp in double down since the
  // point may not be on the surface 
  fdagmc->get_angle(My_sulf_hit, position, ang);

  G4ThreeVector normal = G4ThreeVector(ang[0], ang[1], ang[2]);

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
  G4double forwardDistance,reverseDistance;

  // use the safety
  G4double safety = DistanceToIn(p);
  
  // if we aren't close enough to the surface
  if (safety > kCarToleranceHalf) {
    fdagmc->ray_fire(fvolEntity, position, dir, next_surf, reverseDistance, NULL, 0, -1);
    distance = reverseDistance*cm;  // convert back to mm
    if (next_surf == 0) return kInfinity;
  } else {
    // otherwise use the safety
    distance = safety*cm;
  }

  if(debug) {
    std::cout << "DistanceToIn(raytrace) " << std::endl;
    std::cout << "Name: " << Myname << std::endl;
    std::cout << "pos: " << p.x() << " " << p.y() << "  " << p.z() << std::endl;
    std::cout << "direction: " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    std::cout << "distance: " << distance << std::endl;
  }
  return distance;

  /*   
  // perform the ray fire with modified dag call
  fdagmc->ray_fire(fvolEntity, position, dir, next_surf, forwardDistance, NULL, 0,  1);
  fdagmc->ray_fire(fvolEntity, position, dir, next_surf, reverseDistance, NULL, 0, -1);
  distance = reverseDistance;
  distance *= cm;  // convert back to mm

  // sometimes when on a surface, the entering intersection is missed
  // but the exiting intersection is found

  if(debug) {
    std::cout << "DistanceToIn(raytrace) " << std::endl;
    std::cout << "Name: " << Myname << std::endl;
    std::cout << "pos: " << p.x() << " " << p.y() << "  " << p.z() << std::endl;
    std::cout << "direction: " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    std::cout << "distance: " << distance << std::endl;
    std::cout << "forwraddistance: " << forwardDistance << std::endl;

    G4ThreeVector newpos = p - v*kCarToleranceHalf;
    G4double npos[3] = {newpos.x() / cm, newpos.y() / cm, newpos.z() / cm};
    fdagmc->ray_fire(fvolEntity, npos, dir, next_surf, reverseDistance, NULL, 0, -1);
    std::cout << "bumpeddistance: " << reverseDistance << std::endl;

  }  
  if (next_surf == 0) { // no intersection
    if(debug) {
      std::cout << "hit:nothing" << std::endl;
      std::cout << "return: " << kInfinity << std::endl;
    }
    return kInfinity;
  } else {
    if(debug) {
      std::cout << "hit:surface" << std::endl;
      std::cout << "return: " << distance << std::endl;
    }
    return distance;
  }
  */
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
  minDist *= cm;  // convert back to mm

  if(debug) {
    std::cout << "DistanceToIn(point)" << std::endl;
    std::cout << "pos: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    std::cout << "return: " << minDist << std::endl;
  }
  
  //if (minDist <= kCarToleranceHalf)
  //  return 0.0;
  //else
  return minDist;
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

G4double DagSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                 const G4bool calcNorm, G4bool* validNorm,
                                 G4ThreeVector* n) const {

  double position[3] = {p.x()/cm, p.y()/cm, p.z()/cm};  // convert position to cm
  G4ThreeVector vec = v.unit();
  double dir[3] = {vec.x(), vec.y(), vec.z()};

  EntityHandle next_surf,next_surf_b;
  G4double distance;
  G4double forwardDistance, backwardDistance;

  fdagmc->ray_fire(fvolEntity,position,dir,next_surf,forwardDistance,NULL,0,1);
  fdagmc->ray_fire(fvolEntity,position,dir,next_surf_b,backwardDistance,NULL,0,-1);

  distance = forwardDistance;

  if(debug) {
    std::cout << "DistanceToOut(trace) " << std::endl;
    std::cout << "pos: " << p.x() << " " << p.y() << "  " << p.z() << std::endl;
    std::cout << "distance: " << distance << std::endl;
    std::cout << "direction: " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    std::cout << "forwraddistance: " << forwardDistance << std::endl;
  }  

  // no surfaces to hit
  if (next_surf == 0) {
    if(debug) {
      std::cout << "return: kInfinity" << std::endl;
    }
    return kInfinity;
  }
  
  // calculate the normal
  if (calcNorm) {
    *n = SurfaceNormal(p + v*distance);
    DagMC::RayHistory history;
    // fire twice
    fdagmc->ray_fire(fvolEntity,position,dir,next_surf,forwardDistance,&history,0,1);
    fdagmc->ray_fire(fvolEntity,position,dir,next_surf,forwardDistance,&history,0,1);
    if(next_surf != 0 ) 
      *validNorm = true;
    else
      *validNorm = false;
  }

  // must convert here since previous needs it in cm
  distance *= cm;  // convert back to mm
  
  // we are on a surface
  if ( distance > 0.0 && distance <= kCarToleranceHalf ) {
    if(debug) {
      std::cout << "on a surface" << std::endl;
      std::cout << "distance: " << distance << std::endl;
    }
    return 0.0;
  }

  // distance greater than kCarTolerance
  if ( distance > 0.0 ) {
    if (debug) {
      std::cout << "normal intersection" << std::endl;
      std::cout << "distance: " << distance << std::endl;      
    }
    return distance;
  }
  std::cout << "help" << std::endl;
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an inside point.
// The distance can be an underestimate.
G4double DagSolid::DistanceToOut(const G4ThreeVector& p) const {
  // same as DistanceToIn - calculate safety
  return DistanceToIn(p);
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

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetCubicVolume() {
  G4double result;
  fdagmc->measure_volume(fvolEntity, result);
  return result * cm * cm * cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetSurfaceArea() {
  G4double result;
  fdagmc->measure_area(fvolEntity, result);
  return result * cm * cm;
}
