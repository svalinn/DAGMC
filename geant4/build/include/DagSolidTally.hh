#include <map>
#include "../pyne/pyne.h"

/*
 * General note about UWUW-DagSolid tallies, the pure PyNE form is not currently
 * suited to use in Geant4 directly, however, we can load the tallies from file
 * turn the id number in the DagTally into an index, so in indirect correlation.
 */

/*
 * Loads all UWUW tallies from the hdf5 file
 */
std::map<std::string,pyne::Tally> load_uwuw_tallies(std::string filepath);
