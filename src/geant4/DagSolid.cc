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
// Author:              A. Davis M. C. Han, C. H. Kim, J. H. Jeong, Y. S. Yeom, S.
//                      Kim, Paul. P. H. Wilson, J. Apostolakis
// Organisation:        The University of Wisconsin-Madison, USA & Hanyang Univ., KR
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

//#include "G4TessellatedSolid.hh"
#include "DagSolid.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "G4TessellatedSolid.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>

#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "DagMC.hpp"
using namespace moab;

#include "DagSolid.hh"

#define plot true

//
// Standard contructor has blank name and defines no facets.
//
DagSolid::DagSolid ()
  : G4TessellatedSolid("dummy"), cubicVolume(0.), surfaceArea(0.)
{

  geometryType = "DagSolid";

  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;

  // geometric precision for surfaces
  delta = 0.5*kCarTolerance;
  last_surf_hit = 0;
}


///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no facets
// to detine.
//
DagSolid::DagSolid (const G4String &name, DagMC* dagmc, int volID)
  : G4TessellatedSolid(name), cubicVolume(0.), surfaceArea(0.)
{
  geometryType = "DagSolid";
  Myname=name;

  fdagmc=dagmc;
  fvolID=volID;
  fvolEntity = fdagmc->entity_by_index(3, volID);

  delta = 0.5*kCarTolerance;
  last_surf_hit = 0;
  
  double min[3],max[3];
  fdagmc->getobb(fvolEntity,min,max);

  xMinExtent =  min[0]/mm;
  xMaxExtent =  max[0]/mm;
  yMinExtent =  min[1]/mm;
  yMaxExtent =  max[1]/mm;
  zMinExtent =  min[2]/mm;
  zMaxExtent =  max[2]/mm;

  int num_entities;
  std::vector<EntityHandle> surfs;
  std::vector<EntityHandle> tris;
  const EntityHandle *tri_conn;

  std::vector<CartVect> coords(3);
  G4ThreeVector vertex[3];
  int n_verts=0;

  Interface* moab = dagmc->moab_instance();
  moab->get_child_meshsets(fvolEntity, surfs, 1 );

  if ( plot ) {
    //  G4cout<<"please wait for visualization... "<<G4endl;
    for(unsigned i=0 ; i<surfs.size() ; i++) {
      moab->get_number_entities_by_type( surfs[i], MBTRI, num_entities);
      //G4cout<<"Number of triangles = "<<num_entities<<" in surface index: "<<fdagmc->index_by_handle(surfs[i])<<G4endl;
      //G4cout<<"please wait for visualization... "<<G4endl;

      moab->get_entities_by_type( surfs[i], MBTRI, tris);

      for (unsigned j=0 ; j <tris.size() ; j++) {
        moab->get_connectivity( tris[j], tri_conn, n_verts );
        moab->get_coords( tri_conn, n_verts, coords[0].array() );

        //	  G4cout<<"add facet for vis = "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<G4endl;

        vertex[0] = G4ThreeVector( coords[0][0]*cm, coords[0][1]*cm, coords[0][2]*cm );
        vertex[1] = G4ThreeVector( coords[1][0]*cm, coords[1][1]*cm, coords[1][2]*cm );
        vertex[2] = G4ThreeVector( coords[2][0]*cm, coords[2][1]*cm, coords[2][2]*cm );

        G4TriangularFacet *facet = new G4TriangularFacet (vertex[0], vertex[1], vertex[2], ABSOLUTE);
        //  G4cout << vertex[0] << " " << vertex[1] << " " << vertex[2] << G4endl;
        AddFacet((G4VFacet*)facet);

      }
      tris.clear();
    }
  }

//SetRandomVectorSet();

  SetSolidClosed(true);

//  G4cout<<"Number Of Facets = "<<GetNumberOfFacets() <<G4endl;
//  G4cout <<"maximum point = "<< xMinExtent <<" "<< yMinExtent <<" "<< zMinExtent << G4endl
//         <<"minimum point = "<< xMaxExtent <<" "<< yMaxExtent <<" "<< zMaxExtent << G4endl;

}



///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
DagSolid::DagSolid( __void__& a )
  : G4TessellatedSolid(a),
  geometryType("DagSolid"), cubicVolume(0.), surfaceArea(0.),
  xMinExtent(0.), xMaxExtent(0.),
  yMinExtent(0.), yMaxExtent(0.),
  zMinExtent(0.), zMaxExtent(0.)
{
  //SetRandomVectorSet();
}



///////////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
DagSolid::~DagSolid ()
{}


//
// 
// EInside DagSolid::Inside (const G4ThreeVector &p) const
//
// This method must return:
//  • kOutside if the point at offset p is outside the shape boundaries plus Tolerance/2,
//  • kSurface if the point is <= Tolerance/2 from a surface, or
//  • kInside otherwise.
//
EInside DagSolid::Inside (const G4ThreeVector &p) const
{
  G4double point[3]= {p.x()/cm, p.y()/cm, p.z()/cm}; //convert to cm

  double u = rand();
  double v = rand();
  double w = rand();

  const double magnitude = sqrt( u*u + v*v + w*w );
  u /= magnitude;
  v /= magnitude;
  w /= magnitude;

  G4double direction[3]= {u,v,w};

  G4double minDist = 0.0;

  int result;
  ErrorCode rval;
  rval = fdagmc->point_in_volume(fvolEntity, point, result,
				 direction);  // if uvw is not given, this function generate uvw ran
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from point_in_volume" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Proposed Result = " << result << G4endl;
    G4Exception("DagSolid::Inside(p)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }
  // note: would be nice to early exit from this function
  // but we to be more careful than that, if point is wihtin
  // kCartolerance/2.0 of surface we are actually on the surface
  rval = fdagmc->closest_to_location(fvolEntity,point,minDist);
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from clost_to_location" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << minDist/mm << " mm" << G4endl;
    G4Exception("DagSolid::Inside(p)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }
  // convert to mm
  minDist/=mm;
  // if on surface
  if (minDist <= delta) {
    return kSurface;
  } else {
    // either outside or  inside
    if ( result == 0 )
      return kOutside;
    else
      return kInside;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4ThreeVector DagSolid::SurfaceNormal (const G4ThreeVector &p) const
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.

G4ThreeVector DagSolid::SurfaceNormal (const G4ThreeVector &p) const
{
  G4double position[3]= {p.x()/cm,p.y()/cm,p.z()/cm}; //convert to cm

  // right now we are assuming that this is called only after a succesful rayfire
  // hence last_surf_hit is used. we should modify closest_to_location to allow
  // returning of the surface as well as distance
  G4double angle[3] = {0.,0.,0.};
  moab::ErrorCode rval = fdagmc->get_angle(last_surf_hit,position,angle);

  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from closest_to_location" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl;
    G4Exception("DagSolid::SurfaceNormal(p)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }
  
  G4ThreeVector normal = G4ThreeVector(angle[0],angle[1],angle[2]);

  return normal;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)
//
// Return the distance along the normalised vector v to the shape,
// from the point at offset p. If there is no intersection,
//  return kInfinity. The first intersection resulting from `leaving' a
// surface/volume is discarded. Hence, this is tolerant
// of points on surface of shape.
// 
///////////////////////////////////////////////////////////////////////////////
G4double DagSolid::DistanceToIn (const G4ThreeVector &p, 
                                 const G4ThreeVector &v) const
{

  G4double position[3] = {p.x()/cm,p.y()/cm,p.z()/cm}; //convert to cm
  G4double dir[3] = {v.x(),v.y(),v.z()};
  EntityHandle next_surf;
  G4double distance = kInfinity;

  DagMC::RayHistory history;

  // perform the ray fire with modified dag call
  moab::ErrorCode rval = fdagmc->ray_fire(fvolEntity,position,dir,next_surf,distance,&history,0,-1);
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from ray_fire" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Direction:" << G4endl << G4endl
            << "v.x() = "   << v.x() << G4endl
            << "v.y() = "   << v.y() << G4endl
            << "v.z() = "   << v.z() << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << distance/mm << " mm" << G4endl;
    G4Exception("DagSolid::DistanceToIn(p,v)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }
  
  last_surf_hit = next_surf;
  history.reset();
  distance *= cm; // convert back to mm

  if ( next_surf == 0 ) { // no intersection
    return kInfinity;
  } else if (distance <= delta) {    
    return 0.0;
  } else {
    return distance;
  }
  // if less than delta, return 0
}



///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an outside point p.
// The distance can be an underestimate.
//
//////////////
G4double DagSolid::DistanceToIn (const G4ThreeVector &p) const
{
  G4double safety = 0.0;
  G4double point[3]= {p.x()/cm, p.y()/cm, p.z()/cm}; // convert position to cm
  moab::ErrorCode rval;
  rval = fdagmc->closest_to_location(fvolEntity, point, safety);
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from clostest_to_location" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << safety/mm << " mm" << G4endl;
    G4Exception("DagSolid::DistanceToIn(p)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }
  safety *= cm; // convert back to mm

  return safety;
}


///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
//                        const G4bool calcNorm=false,
//                        G4bool *validNorm=0, G4ThreeVector *n=0);
//
// Return distance along the normalised vector v to the shape,
// from a point at an offset p inside or on the surface of
// the shape. Intersections with surfaces, when the point is not
// greater than kCarTolerance/2 from a surface, must
// be ignored.
// 
// If calcNorm is true, then it must also set validNorm to either
//    • true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//    • false, if the solid does not lie entirely behind or on the
//         exiting surface.
// If calcNorm is false, then validNorm and n are unused.
//
G4double DagSolid::DistanceToOut (const G4ThreeVector &p,
                                  const G4ThreeVector &v, const G4bool calcNorm,
                                  G4bool *validNorm, G4ThreeVector *n) const
{
  G4double minDist = kInfinity;
  double position[3]= {p.x()/cm,p.y()/cm,p.z()/cm}; //convert position to cm
  double dir[3]= {v.x(),v.y(),v.z()};

  EntityHandle next_surf;
  double next_dist;
  DagMC::RayHistory history;

  moab::ErrorCode rval = fdagmc->ray_fire(fvolEntity,position,dir,next_surf,
					    next_dist,&history,0,1);
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from ray_fire" << G4endl
	    << "ErrorCode " << rval << G4endl
	    << "Position:"  << G4endl << G4endl
	    << "p.x() = "   << p.x()/mm << " mm" << G4endl
	    << "p.y() = "   << p.y()/mm << " mm" << G4endl
	    << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
	    << "Direction:" << G4endl << G4endl
	    << "v.x() = "   << v.x() << G4endl
	    << "v.y() = "   << v.y() << G4endl
	    << "v.z() = "   << v.z() << G4endl << G4endl
	    << "Proposed distance :" << G4endl << G4endl
	    << "snxt = "    << next_dist/mm << " mm" << G4endl;
    G4Exception("DagSolid::DistanceToOut(p,v,...)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }

  last_surf_hit = next_surf;
  history.reset();
  next_dist *= cm; // convert back to mm

  // no more surfaces
  if(next_surf == 0 )
    return kInfinity;

  // true if solid lies entirely behind or on the surface
  // otherwise false
  if (calcNorm) {
    *n         = SurfaceNormal(p+minDist*v);
    *validNorm = false;
  }

  if (next_dist < minDist )
    minDist = next_dist;

  // particle considered to be on surface
  if ( minDist <= delta ) {
    *validNorm = true;
    return 0.0;
  } else if ( minDist > delta && minDist < kInfinity) {
    if(calcNorm) // if calc norm true, s should set validNorm true
      *validNorm = true;
    return minDist;
  } else {
    G4cout << minDist << G4endl;
    return kInfinity; //was kinfinity
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an inside point.

G4double DagSolid::DistanceToOut (const G4ThreeVector &p) const
{
  G4double safety = 0.0;
  G4double point[3]= {p.x()/cm, p.y()/cm, p.z()/cm}; // convert to cm

  moab::ErrorCode rval = fdagmc->closest_to_location(fvolEntity, point, safety);
  if(rval != moab::MB_SUCCESS) {
    std::ostringstream message;
    message << "Failure from clostest_to_location" << G4endl
	    << "ErrorCode " << rval << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << safety/mm << " mm" << G4endl;
    G4Exception("DagSolid::DistanceToOut(p)","GeomSolidsFailure",
		EventMustBeAborted,message);
  }

  safety *= cm; // convert back to mm
  if ( safety < delta )
    return 0.0;
  else
    return safety;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4GeometryType GetEntityType() const;
//
// Provide identification of the class of an object (required for persistency and STEP interface).
//
G4GeometryType DagSolid::GetEntityType () const
{
  return geometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
void DagSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}


std::ostream &DagSolid::StreamInfo(std::ostream &os) const
{
  os << G4endl;
  os << "Geometry Type    = " << geometryType  << G4endl;
  return os;
}


///////////////////////////////////////////////////////////////////////////////
//
// CalculateExtent
//
// Based on correction provided by Stan Seibert, University of Texas.
//
G4bool
DagSolid::CalculateExtent(const EAxis pAxis,
                          const G4VoxelLimits& pVoxelLimit,
                          const G4AffineTransform& pTransform,
                          G4double& pMin, G4double& pMax) const
{


// Calculate the minimum and maximum extent of the solid,
// when under the specified transform, and within the specified limits.
// If the solid is not intersected by the region, return false, else return true.


  G4ThreeVector minExtent(xMinExtent, yMinExtent, zMinExtent);
  G4ThreeVector maxExtent(xMaxExtent, yMaxExtent, zMaxExtent);


  // Check for containment and clamp to voxel boundaries
  for (G4int axis=G4ThreeVector::X; axis < G4ThreeVector::SIZE; axis++) {
    EAxis geomAxis = kXAxis; // G4 geom classes use different index type
    switch(axis) {
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
      if ( minExtent[axis] > voxelMaxExtent+kCarTolerance ||
           maxExtent[axis] < voxelMinExtent-kCarTolerance    ) {
        return false ;
      } else {
        if (minExtent[axis] < voxelMinExtent) {
          minExtent[axis] = voxelMinExtent ;
        }
        if (maxExtent[axis] > voxelMaxExtent) {
          maxExtent[axis] = voxelMaxExtent;
        }
      }
    }
  }

  // Convert pAxis into G4ThreeVector index
  G4int vecAxis=0;
  switch(pAxis) {
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

G4double DagSolid::GetMinXExtent () const
{
  return xMinExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxXExtent () const
{
  return xMaxExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinYExtent () const
{
  return yMinExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxYExtent () const
{
  return yMaxExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinZExtent () const
{
  return zMinExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxZExtent () const
{
  return zMaxExtent/cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetCubicVolume()
{
  G4double result;
  fdagmc->measure_volume(fvolEntity, result);
  return result*cm*cm*cm;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetSurfaceArea()
{
  G4double result;
  fdagmc->measure_area(fvolEntity, result);
  return result*cm*cm;
}

