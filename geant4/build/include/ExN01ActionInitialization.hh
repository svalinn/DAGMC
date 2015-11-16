#ifndef ExN01ActionInitialization_h
#define ExN01ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include <string>

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

/// Action initialization class.
///

class ExN01ActionInitialization : public G4VUserActionInitialization
{
 public:
  ExN01ActionInitialization(UWUW *uwuw_workflow_data);
  virtual ~ExN01ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
  UWUW *workflow_data;
};

#endif
