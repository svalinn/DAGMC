
#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4SDParticleFilter.hh"
#include "DagSolidTally.hh"
#include <map>
#include <string>

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
 public:

  ExN01DetectorConstruction(UWUW *uwuw_workflow_data);
  ~ExN01DetectorConstruction();

 public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  // the tally library
  std::map<std::string,pyne::Tally> tally_library;
  // dag_volumes collection mapped by id number
  std::map<int, G4LogicalVolume*> dag_logical_volumes;
  // particle filters for tallies
  std::map<std::string,G4SDParticleFilter*> particle_filters;

  void BuildParticleFilter(std::string particle_name);
  void build_histogram();
  void add_histogram_description(std::string tally_name);
  void end_histogram();
 private:
  std::string _to_string(int var);

 private:

  // Logical volumes
  //
  std::string uwuw_filename;

  G4LogicalVolume* world_volume_log;

  UWUW *workflow_data;

  // DAG Logical volumes
  // G4LogicalVolume* dag_vol_log;

  // Physical volumes
  //
  // G4VPhysicalVolume* world_volume_phys;
  // DAG Physical volumes


};

#endif
