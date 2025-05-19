#ifndef DagSolid_hh
#define DagSolid_hh 1

#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include "DagMC.hpp"
#include "G4Polyhedron.hh"
#include "G4AffineTransform.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VoxelLimits.hh"
#include "globals.hh"

class DagSolid : public G4VSolid {
 public:  // with description
  // empty constructor
  DagSolid();
  // main constructor
  DagSolid(const G4String& name, moab::DagMC* dagmc, int volID);
  // destructor
  virtual ~DagSolid();

  // given a point p determine if the point is inside the volume
  virtual EInside Inside(const G4ThreeVector& p) const;

  // return the surface normal given the point p
  virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;

  // given the point p and the vector v, determine the distance until
  // the ray enters the solid
  virtual G4double DistanceToIn(const G4ThreeVector& p,
                                const G4ThreeVector& v) const;

  // given the point p, determine an estimate until the ray
  // would enter the solid
  virtual G4double DistanceToIn(const G4ThreeVector& p) const;

  // given the point p and the vector v, determine the distance until
  // the ray leaves the solid, could be asked to provide a valid normal
  virtual G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                 const G4bool calcNorm = false,
                                 G4bool* validNorm = 0,
                                 G4ThreeVector* n = 0) const;

  // given the point p, provide an estimate of the distance until
  // a ray would leave
  virtual G4double DistanceToOut(const G4ThreeVector& p) const;

  // return the GeometryType
  virtual G4GeometryType GetEntityType() const;

  // calculate the extent of the solid given an axis and rotation
  virtual G4bool CalculateExtent(const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

  // text dump of the DagSolid object
  virtual std::ostream& StreamInfo(std::ostream& os) const;

  // return a valid point on the surface of the solid
  G4ThreeVector GetPointOnSurface() const;

  // return the volume of the solid
  virtual G4double GetCubicVolume();

  // returns the total surface areas of the
  // solid
  virtual G4double GetSurfaceArea();

  // return the bounding limits of the solid
  void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

  // return the minimum extent of the solid in the x direction
  G4double GetMinXExtent() const;

  // return the maximum extent of the solid in the x direction
  G4double GetMaxXExtent() const;

  // return the minimum extent of the solid in the y direction
  G4double GetMinYExtent() const;

  // return the maxmimum extent of the solid in the y direction
  G4double GetMaxYExtent() const;

  // return the minimum extent of the solid in the z direction
  G4double GetMinZExtent() const;

  // return the maximum extent of the solid in the z direction
  G4double GetMaxZExtent() const;

  // Functions for visualization

  // add the solid to a G4 scene
  virtual void DescribeYourselfTo(G4VGraphicsScene& scene) const;

  // create a polyhedron for visualisation
  G4Polyhedron* CreatePolyhedron() const;

  // returns the polyhedron if it exists, if not
  // calls CreatePolyhedron
  G4Polyhedron* GetPolyhedron() const;

  // gets an extent for visualisation
  G4VisExtent GetExtent() const;

 public:  // without description
  DagSolid(__void__&);
  // Fake default constructor for usage restricted to direct object
  // persistency for clients requiring preallocation of memory for
  // persistifiable objects.

 protected:  // with description
  void DeleteObjects();
  void CopyObjects(const DagSolid& s);

 private:
  G4GeometryType geometryType;  // the geant4 solid type
  G4double cubicVolume;         // the volumes volume
  G4double surfaceArea;         // the surface area of the volume
  G4double xMinExtent;          // lowest x extent
  G4double xMaxExtent;          // highest x extent
  G4double yMinExtent;          // lowest y extent
  G4double yMaxExtent;          // highest y extent
  G4double zMinExtent;          // lowest z extent
  G4double zMaxExtent;          // highest z extent

  G4String Myname;  // the name of the solid

  G4double kCarToleranceHalf =
      kCarTolerance / 2.;         // half the tracking tolerance
  moab::DagMC* fdagmc;            // the DAGMC pointer
  moab::Interface* moab;          // the MOAB pointer associated with DAGMC
  G4int fvolID;                   // the ID of the volume
  moab::EntityHandle fvolEntity;  // its MOAB entity handle

  mutable G4Polyhedron* fPolyhedron = nullptr;  // a pointer

  std::vector<moab::EntityHandle> fSurfaces;   // surfaces of the volume
  std::vector<moab::EntityHandle> fTriangles;  // triangles of the volume
};

#endif
