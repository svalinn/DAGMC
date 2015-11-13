/// \file ExN01DetectorHit.hh
/// \brief Definition of the ExN01DetectorHit class

#ifndef ExN01DetectorHit_h
#define ExN01DetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class ExN01DetectorHit : public G4VHit
{
 public:
  ExN01DetectorHit();
  ExN01DetectorHit(const ExN01DetectorHit&);
  virtual ~ExN01DetectorHit();

  // operators
  const ExN01DetectorHit& operator=(const ExN01DetectorHit&);
  G4int operator==(const ExN01DetectorHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  // methods from base class
  virtual void Draw();
  virtual void Print();

  // Set methods
  void SetTrackID  (G4int track)      {
    fTrackID = track;
  };
  /*
  void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
  void SetEdep     (G4double de)      { fEdep = de; };
  void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
  */
  void SetParticleEnergy       (G4double ke)      {
    fKe  = ke;
  };
  void SetTrackLength  (G4double tl)  {
    fTl = tl;
  };
  void SetWeight  (G4double weight)  {
    fWeight = weight;
  };
  void SetParticleName (G4String name) {
    fName = name;
  };
  void SetParticlePDG (G4int PID)      {
    fPID = PID;
  };

  // Get methods
  G4int GetTrackID() const     {
    return fTrackID;
  };
  /*
  G4int GetChamberNb() const   { return fChamberNb; };
  G4double GetEdep() const     { return fEdep; };
  G4ThreeVector GetPos() const { return fPos; };
  */
  G4double GetKE() const       {
    return fKe;
  };
  G4double GetTrackLength() const {
    return fTl;
  };
  G4double GetWeight()      const {
    return fWeight;
  };
  G4String GetParticleName() const {
    return fName;
  };
  G4int GetParticlePDG() const {
    return fPID;
  };

 private:

  G4int         fTrackID;
//      G4int         fChamberNb;
//      G4double      fEdep;
//      G4ThreeVector fPos;
  G4double      fKe;
  G4double      fTl;
  G4double      fWeight;
  G4String      fName;
  G4int         fPID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<ExN01DetectorHit> ExN01DetectorHitsCollection;

extern G4ThreadLocal G4Allocator<ExN01DetectorHit>* ExN01DetectorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* ExN01DetectorHit::operator new(size_t)
{
  if(!ExN01DetectorHitAllocator)
    ExN01DetectorHitAllocator = new G4Allocator<ExN01DetectorHit>;
  return (void *) ExN01DetectorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void ExN01DetectorHit::operator delete(void *hit)
{
  ExN01DetectorHitAllocator->FreeSingle((ExN01DetectorHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
