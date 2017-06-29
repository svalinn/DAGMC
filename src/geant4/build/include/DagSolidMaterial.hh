#ifndef DAGSOLIDMATERIAL_HH
#define DAGSOLIDMATERIAL_HH

#include <map>
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "../pyne/pyne.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "uwuw.hpp"

#define frand() static_cast <float> (rand()) / static_cast <float> (RAND_MAX)

// struct to store hue, saturation and lightness
  struct hsl{
    float h;
    float s;
    float l;
  };

  // struct to store red, green and blue values
  struct rgb {
    float r;
    float g;
    float b;
  };

rgb hsl_to_rgbp(float c, float x, hsl hsl_c);
rgb rgbp_to_rgb(rgb rgbp, float m );
rgb hsl_to_rgb(hsl color);
std::vector<rgb> generate_colors(int n);

/*
 * Load the full material collection from PyNE library, return G4Materials
 */
std::map<std::string, G4Material*> load_uwuw_materials(UWUW* filename);

/*
 * Load the material library from the h5m file into a standard map of PyNE Materials
 */
std::map<std::string, pyne::Material> load_materials(std::string filepath);

/*
 * From the map of PyNE material objects, populate a list of G4 Isotopes
 */
std::map<int, G4Isotope*> get_g4isotopes(std::map<std::string, pyne::Material> material_library);

/*
 * From the Map of G4Isotopes create an element for each Isotope
 */
std::map<int, G4Element*> get_g4elements(std::map<int, G4Isotope*> isotope_map);

/*
 * Using the MaterialLiibrary and the element map, make the G4 Materials
 */
std::map<std::string, G4Material*> get_g4materials(std::map<int, G4Element*> element_map,
                                                   std::map<std::string, pyne::Material> material_library);
/*
 * Given the map of G4Materials, returns a map of colours by material name so that the colours
 * are unique 
 */                                                   
std::map<std::string, G4VisAttributes*> get_visualisation_attributes(std::map<std::string,G4Material*> materials);

#endif /* DAGSOLIDMATERIAL_HH */