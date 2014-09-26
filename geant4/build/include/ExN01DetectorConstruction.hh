
#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "DagSolidTally.hh"
#include <map>
#include <string>

class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  ExN01DetectorConstruction(std::string uwuw_file);
    ~ExN01DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

   std::map<std::string,pyne::Tally> tally_library;

  private:
    
    // Logical volumes
    //
    std::string uwuw_filename;

    G4LogicalVolume* world_volume_log;
    // DAG Logical volumes
    G4LogicalVolume* dag_vol_log;

    // Physical volumes
    //
    G4VPhysicalVolume* world_volume_phys;
    // DAG Physical volumes

};

#endif

