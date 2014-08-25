#ifndef ExN01ActionInitialization_h 
#define ExN01ActionInitialization_h 1 
  
#include "G4VUserActionInitialization.hh"
  
/// Action initialization class. 
/// 
  
class ExN01ActionInitialization : public G4VUserActionInitialization 
{ 
  public: 
    ExN01ActionInitialization(); 
    virtual ~ExN01ActionInitialization(); 
  
    virtual void BuildForMaster() const; 
    virtual void Build() const; 
}; 
  
#endif 
