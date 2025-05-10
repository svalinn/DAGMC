#ifndef ExN01ActionInitialization_h
#define ExN01ActionInitialization_h 1

#include <string>

#include "G4AnalysisManager.hh"
#include "G4TScoreHistFiller.hh"
#include "G4Types.hh"

#include "G4VUserActionInitialization.hh"

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

/// Action initialization class.
///

class ExN01ActionInitialization : public G4VUserActionInitialization {
 public:
  ExN01ActionInitialization(UWUW* uwuw_workflow_data);
  virtual ~ExN01ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
  UWUW* workflow_data;
private:
  G4TScoreHistFiller<G4AnalysisManager>* filler;

};

#endif
