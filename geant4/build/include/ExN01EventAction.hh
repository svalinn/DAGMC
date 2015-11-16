/// ile ExN01EventAction.hh
/// brief Definition of the ExN01EventAction class

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
  std::vector<int> GetTallyID(); // get all the tally ids

  // data members
  std::vector<G4int> tally_ids;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
