//
//
// $Id$
//

#include "ExN01PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
//
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"


ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  fParticleGun = new G4GeneralParticleSource();
}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction() {
  delete fParticleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  fParticleGun->GeneratePrimaryVertex(anEvent) ;
}


