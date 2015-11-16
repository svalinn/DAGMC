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
  :G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
  fParticleGun = new G4GeneralParticleSource();
  /*
  G4int n_particle = 1;
  particleGun = new G4GeneralParticleSource();
  //  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  //  particleGun->SetParticleEnergy(14.0*MeV);
  */
}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent) ;


  //G4double y = 1.0*(G4UniformRand()-0.5);
  //  G4double z = 1.0*(G4UniformRand()-0.5);
  //  particleGun->SetParticlePosition(G4ThreeVector(1.24*cm,1.24*cm,100.0*cm));
  // particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  //  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.0,0.0,0.0));

  //  particleGun->SetParticlePosition(G4ThreeVector(-10.0*cm,y*cm,z*cm));
  //  G4int i = anEvent->GetEventID() % 3;

  /*
  G4double cosTheta = 2*G4UniformRand() - 1.;
  G4double phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
    uy = sinTheta*std::sin(phi),
    uz = cosTheta;
  //particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  */

}


