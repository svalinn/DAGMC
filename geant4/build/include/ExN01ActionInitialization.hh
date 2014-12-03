#ifndef ExN01ActionInitialization_h 
#define ExN01ActionInitialization_h 1 
  
#include "G4VUserActionInitialization.hh"
#include <string>
/// Action initialization class. 
/// 
  
class ExN01ActionInitialization : public G4VUserActionInitialization 
{ 
  public: 
    ExN01ActionInitialization(std::string filename); 
    virtual ~ExN01ActionInitialization(); 
  
    virtual void BuildForMaster() const; 
    virtual void Build() const; 
}; 
  
#endif 
