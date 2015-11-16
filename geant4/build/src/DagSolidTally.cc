#include "DagSolidTally.hh"
#include "uwuw.hpp"

// loads the tallies from the file
std::map<std::string,pyne::Tally> load_uwuw_tallies(std::string filepath)
{

  bool end = false;
  std::map<std::string,pyne::Tally> tally_library;
  int i = 0;

  std::cout << filepath << std::endl;
  while( !end ) {
    pyne::Tally tally; // from file
    std::cout << i << std::endl;
    tally.from_hdf5(filepath,"/tally",i++);
    if ( tally_library.count(tally.tally_name) ) {
      end = true;
    } else {
      tally_library[tally.tally_name]=tally;
    }
  }

  for(std::map<std::string,pyne::Tally>::const_iterator it = tally_library.begin() ; it != tally_library.end() ; ++it ) {
    std::cout << it->first <<  std::endl;
  }

  return tally_library;
}
