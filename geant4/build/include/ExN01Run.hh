#ifndef ExN01Run_h
#define ExN01Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4LogicalVolume.hh"

class G4Event;

class ExN01Run : public G4Run {
   public:
     // constructor
     ExN01Run();
     // destructor
     virtual ~ExN01Run();
     //  
     virtual void RecordEvent(const G4Event *evt);

     virtual void Merge(const G4Run *);

     G4double GetTotal(G4int id);

     G4double GetTotal(G4LogicalVolume*, G4String score_name);

   private:
     G4THitsMap<G4double> GetScore(G4int id);

     G4THitsMap<G4double> GetScore(G4LogicalVolume*, G4String score_name);

   private:
     std::map<G4String,G4int> detectors;
     std::map<G4int,G4THitsMap<G4double> > total_score;
};

#endif