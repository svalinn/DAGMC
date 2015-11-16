#ifndef ExN01UserScoreWriter_h
#define ExN01UserScoreWriter_h 1

//#include "MBCore.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "globals.hh"
#include "G4VScoreWriter.hh"
#include <vector>

class ExN01UserScoreWriter : public G4VScoreWriter
{

 public:
  ExN01UserScoreWriter();
  virtual ~ExN01UserScoreWriter();

 public:
  // generates numbers from -end_coord to end_coord with num_bounds bin
  std::vector<double> generate_bin_bounds(int num_bounds, double end_coord);
  moab::ErrorCode generate_moab_mesh(std::vector<double> x_bins,
                                     std::vector<double> y_bins,
                                     std::vector<double> z_bins,
                                     std::vector<moab::EntityHandle> &mesh_elements);

 public:
  virtual void DumpAllQuantitiesToFile(const G4String& fileName,
                                       const G4String& option);


};

moab::Interface *MBI();
#endif
