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
