/// ile ExN01EventAction.hh
/// rief Definition of the ExN01EventAction class

#ifndef ExN01EventAction_h
#define ExN01EventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"
#include "globals.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class ExN01EventAction : public G4UserEventAction
{
public:
  ExN01EventAction();
  virtual ~ExN01EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

private:
  // methods
  G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
                                          const G4Event* event) const;
  G4double GetSum(G4THitsMap<G4double>* hitsMap) const;

  //  void PrintEventStatistics(std::vector<G4Double> tracklengths) const;

  // data members                   
  std::vector<G4double> tracklengths;
  std::vector<G4int> vol_tl_ids;
  
};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
