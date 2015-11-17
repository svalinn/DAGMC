#include <map>
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "../pyne/pyne.h"

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

/*
 * Load the full material collection from PyNE library, return G4Materials
 */
std::map<std::string,G4Material*> load_uwuw_materials(UWUW *filename);

/*
 * Load the material library from the h5m file into a standard map of PyNE Materials
 */
std::map<std::string,pyne::Material> load_materials(std::string filepath);

/*
 * From the map of PyNE material objects, populate a list of G4 Isotopes
 */
std::map<int,G4Isotope*> get_g4isotopes(std::map<std::string, pyne::Material> material_library);

/*
 * From the Map of G4Isotopes create an element for each Isotope
 */
std::map<int,G4Element*> get_g4elements(std::map<int,G4Isotope*> isotope_map);

/*
 * Using the MaterialLiibrary and the element map, make the G4 Materials
 */
std::map<std::string,G4Material*> get_g4materials(std::map<int,G4Element*> element_map,
    std::map<std::string, pyne::Material> material_library);
