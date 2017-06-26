#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4LogicalVolume.hh"


class ExN01Run : public G4Run {
   public:
     // constructor
     ExN01Run();
     // destructor
     virtual ~ExN01Run();
     //  
     void RecordEvent(const G4Event *evt);

     void Merge(const G4Run *aRun);

     G4THitsMap<G4double> GetTotal(G4int id);

     G4THitsMap<G4double> GetTotal(G4LogicalVolume*, G4String score_name);

   private:
     G4int nEvent;
     std::map<G4String,G4int> detectors;
     std::map<G4int,G4THitsMap<G4double> > total_score;
};