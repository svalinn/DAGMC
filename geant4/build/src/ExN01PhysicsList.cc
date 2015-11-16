//
//
// $Id$
//
//

#include "ExN01PhysicsList.hh"
#include "G4ParticleTypes.hh"


ExN01PhysicsList::ExN01PhysicsList()
{
  ;
}

ExN01PhysicsList::~ExN01PhysicsList()
{
  ;
}

void ExN01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4Geantino::GeantinoDefinition();
}

void ExN01PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
}

void ExN01PhysicsList::SetCuts()
{
  // uppress error messages even in case e/gamma/proton do not exist
  G4int temp = GetVerboseLevel();
  SetVerboseLevel(0);
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();

  // Retrieve verbose level
  SetVerboseLevel(temp);
}

