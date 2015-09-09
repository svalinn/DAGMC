#include <iostream>
#include <string>
#include <set>

class name_concatenator {
  public:
    // constructor
    name_concatenator();

    // destructor
   ~name_concatenator();

    // returns the a unique name
    std::string make_name_8bytes(std::string name);

  private: 
    // shifts the string should it be required
    std::string shift_and_increment(std::string name);

    // extracts the chars [A-Z] and numbers [0-9] only 
    std::string extract_alpha_num(std::string name);

  private:
    std::set<std::string> used_b8_names; // the collection of names used so far
    
};

