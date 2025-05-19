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
// Date:                20/12/2010
// Author:              M. C. Han, C. H. Kim, J. H. Jeong, Y. S. Yeom, S. Kim,
//                      Paul. P. H. Wilson, J. Apostolakis
// Organisation:        Hanyang Univ., KR
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2010, J. H. Jeong, Hanyang Univ., KR
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

///////////////////////////////////////////////////////////////////////////////
#ifndef DagSolid_hh
#define DagSolid_hh 1

#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include "DagMC.hpp"
#include "G4AffineTransform.hh"
//#include "G4TessellatedSolid.hh"
#include "G4VSolid.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "globals.hh"

class DagSolid : public G4VSolid {
 public:  // with description
  DagSolid();
  DagSolid(const G4String& name, moab::DagMC* dagmc, int volID);
  virtual ~DagSolid();

  // Mandatory Functions

  virtual EInside Inside(const G4ThreeVector& p) const;

  virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
  virtual G4double DistanceToIn(const G4ThreeVector& p,
                                const G4ThreeVector& v) const;
  virtual G4double DistanceToIn(const G4ThreeVector& p) const;
  virtual G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                 const G4bool calcNorm = false,
                                 G4bool* validNorm = 0,
                                 G4ThreeVector* n = 0) const;
  virtual G4double DistanceToOut(const G4ThreeVector& p) const;
  virtual G4GeometryType GetEntityType() const;

  virtual G4bool CalculateExtent(const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;
                                 
  virtual std::ostream& StreamInfo(std::ostream& os) const;

  G4ThreeVector GetPointOnSurface() const;

  G4Polyhedron* CreatePolyhedron() const;

  G4Polyhedron* GetPolyhedron() const;

  virtual G4double GetCubicVolume();
  virtual G4double GetSurfaceArea();
  
  G4VisExtent GetExtent() const;

  void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

  G4double GetMinXExtent() const;
  G4double GetMaxXExtent() const;
  G4double GetMinYExtent() const;
  G4double GetMaxYExtent() const;
  G4double GetMinZExtent() const;
  G4double GetMaxZExtent() const;
  // Functions for visualization

  virtual void DescribeYourselfTo(G4VGraphicsScene& scene) const;

 public:  // without description
  DagSolid(__void__&);
  // Fake default constructor for usage restricted to direct object
  // persistency for clients requiring preallocation of memory for
  // persistifiable objects.

 protected:  // with description
  void DeleteObjects();
  void CopyObjects(const DagSolid& s);

 private:

  G4GeometryType geometryType; // the geant4 solid type
  G4double cubicVolume; // the volumes volume
  G4double surfaceArea; // the surface area of the volume
  G4double xMinExtent; // lowest x extent
  G4double xMaxExtent; // highest x extent
  G4double yMinExtent; // lowest y extent
  G4double yMaxExtent; // highest y extent
  G4double zMinExtent; // lowest z extent
  G4double zMaxExtent; // highest z extent

  G4String Myname; // the name of the solid

  G4double kCarToleranceHalf = kCarTolerance/2.; // half the tracking tolerance
  moab::DagMC* fdagmc; // the DAGMC pointer
  moab::Interface* moab; // the MOAB pointer associated with DAGMC
  G4int fvolID; // the ID of the volume 
  moab::EntityHandle fvolEntity; // its MOAB entity handle

  mutable G4Polyhedron *fPolyhedron = nullptr; // a pointer 
  
  std::vector<moab::EntityHandle> fSurfaces; // surfaces of the volume
  std::vector<moab::EntityHandle> fTriangles; // triangles of the volume

};

#endif
