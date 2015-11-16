// $Id: B4RunAction.hh 74265 2013-10-02 14:41:20Z gcosmo $
//
/// ile B4RunAction.hh
/// rief Definition of the B4RunAction class

#ifndef ExN01RunAction_h
#define ExN01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

class G4Run;

/// Run action class
///
/// It accumulates statistic
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in B4Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed
/// dispersion is printed.
///

class ExN01RunAction : public G4UserRunAction
{
 public:
  ExN01RunAction(UWUW *uwuw_workflow_data);
  virtual ~ExN01RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

 private:
  UWUW *workflow_data;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
