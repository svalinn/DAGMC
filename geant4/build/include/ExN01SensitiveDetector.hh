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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: ExN01SensitiveDetector.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file ExN01SensitiveDetector.hh
/// \brief Definition of the ExN01SensitiveDetector class

#ifndef HistogramManager_h
#include "HistogramManager.hh"
#define HistogramManager_h 1
#endif

#ifndef ExN01SensitiveDetector_h
#define ExN01SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

#include "ExN01DetectorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// B2Tracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero
/// energy deposit.

class ExN01SensitiveDetector : public G4VSensitiveDetector
{
 public:
  ExN01SensitiveDetector(const G4String& name,
                         const G4String& collectionName,
                         const G4int     detectorIndex,
                         const G4double  detectorVolume,
                         HistogramManager* HM);
  virtual ~ExN01SensitiveDetector();

  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

 private:
  ExN01DetectorHitsCollection* fHitsCollection;
  G4int collectionID;
  G4String DetectorName;
  G4int DetectorIndex;
  G4double DetectorVolume;
  G4int* SensitiveParticles;
  G4int  NumberOfParticles;
  std::map<G4int,G4int> hist_part_map;
  G4int hist_index;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
