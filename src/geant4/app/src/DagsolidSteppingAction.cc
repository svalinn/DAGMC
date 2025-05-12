#include "DagSolidSteppingAction.hh"

TimeCutSteppingAction::TimeCutSteppingAction() {
}

TimeCutSteppingAction::~TimeCutSteppingAction() {
}

void TimeCutSteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (track->GetGlobalTime() >= 1.0 * us) {
        track->SetTrackStatus(fStopAndKill);
        G4cout << "Track killed at time: " << track->GetGlobalTime() / ms << " ms";
        G4cout << " due to time cut" << G4endl;
    }
    return;
}    
