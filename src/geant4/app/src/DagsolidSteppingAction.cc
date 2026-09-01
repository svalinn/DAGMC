#include "DagSolidSteppingAction.hh"
#include "G4Electron.hh"

TimeCutSteppingAction::TimeCutSteppingAction() {
}

TimeCutSteppingAction::~TimeCutSteppingAction() {
}

void TimeCutSteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (track->GetDefinition() == G4Electron::Definition() && track->GetLocalTime() >= 10.0 * us) {
        track->SetTrackStatus(fStopAndKill);
        G4cout << "Track killed at time: " << track->GetLocalTime() / ms << " ms";
        G4cout << " due to time cut" << G4endl;
    }
    return;
}    
