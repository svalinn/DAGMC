#ifndef DAG_SOLID_STEPPING_ACTION_HH
#define DAG_SOLID_STEPPING_ACTION_HH 1

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

class TimeCutSteppingAction : public G4UserSteppingAction {
    public:
        TimeCutSteppingAction();
       ~TimeCutSteppingAction();

        virtual void UserSteppingAction(const G4Step* step);
            
};
#endif