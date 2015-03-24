#ifndef ExN01UserScoreWriter_h
#define ExN01UserScoreWriter_h 1

#include "globals.hh"
#include "G4VScoreWriter.hh"
#include <vector>

class ExN01UserScoreWriter : public G4VScoreWriter {

public:
  ExN01UserScoreWriter();
  virtual ~ExN01UserScoreWriter();

  // generates numbers from -end_coord to end_coord with num_bounds bins
  std::vector<double> generate_bin_bounds(int num_bounds, double end_coord);

public:
  virtual void DumpAllQuantitiesToFile(const G4String& fileName,
				       const G4String& option);
};

#endif
