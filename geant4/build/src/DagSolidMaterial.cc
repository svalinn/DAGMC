#include "DagSolidMaterial.hh"

// this function takes a UWUW library, makes a collection of g4istopes
// to make a collection of g4elements, ultimately returning a map
// of g4materials
std::map<std::string,G4Material*> load_uwuw_materials(UWUW *workflow_data)
{
  // new material library
  std::map<std::string,pyne::Material> material_library;
  material_library = workflow_data->material_library;

  // make sure to expand_elements
  std::map<std::string,pyne::Material>::iterator it;
  for ( it = material_library.begin() ; it != material_library.end() ; ++it ) {
    pyne::Material new_mat = it->second;
    new_mat = new_mat.expand_elements();
    material_library[it->first] = new_mat;
  }


  std::map<std::string,G4Material*> g4_mat_empty;
  if (material_library.size() == 0 )
    return g4_mat_empty;

  // get all the nuclides in the problem
  std::map<int,G4Isotope*> g4_isotopes;
  g4_isotopes = get_g4isotopes(material_library);

  // generate an element for each isotopes, this is so we need only create
  // one element for each nuclide, this makes making G4 materials much easier
  std::map<int,G4Element*> g4_elements;
  g4_elements = get_g4elements(g4_isotopes);

  // from the element map, now generate G4Materials
  std::map<std::string,G4Material*> g4_materials;
  g4_materials = get_g4materials( g4_elements, material_library );

  return g4_materials;
}


// Given a pyne material library, get hold of the materials and return the collection
// of nuclides in the problem *note* using the name istope is not correct but follows
// the Geant4 style
std::map<int,G4Isotope*> get_g4isotopes(std::map<std::string, pyne::Material> material_library)
{
  // map to store the isotopes
  std::map<int,G4Isotope*> g4isotopes;
  pyne::Material sum_mat;

  // parameters for new isotope
  std::string name; // name of the isotope
  int iz, n; // atomic number, nucleon number
  double a; // atomic mass

  // loop over the materials present and produce a single new material which contains all the isotopes
  for(std::map<std::string,pyne::Material>::const_iterator it = material_library.begin() ; it != material_library.end() ; ++it ) {
    sum_mat = sum_mat + (it->second);
  }

  // loop through the comp map of this material
  pyne::comp_iter mat_it;
  // loop over the nuclides in the composition, making a new g4isotope for each one
  for ( mat_it = sum_mat.comp.begin() ; mat_it != sum_mat.comp.end() ; ++mat_it ) {
    int nuc_id = mat_it->first;
    G4Isotope* new_iso = new G4Isotope(name=pyne::nucname::name(nuc_id),iz=pyne::nucname::znum(nuc_id), n=pyne::nucname::anum(nuc_id),
                                       a=(pyne::atomic_mass(nuc_id)*g/mole));
    g4isotopes[nuc_id]=new_iso;
  }

  return g4isotopes;
}

// this may seem odd to a geant4 user, but since we may care about
// specific nuclides (g4isotopes) we need to make sure that our materials
// have as much data regarding them as possible, since g4materials can only
// be made of elements, we make an element for each g4isotope
std::map<int,G4Element*> get_g4elements(std::map<int,G4Isotope*> isotope_map)
{
  std::map<int,G4Element*> element_map;
  std::map<int,G4Isotope*>::iterator it;

  std::string name, symbol; // element name and symbol
  int ncomponents ; // number of constituent isotopes
  double abundance; //abundance of each isotope
  for ( it = isotope_map.begin() ; it != isotope_map.end() ; ++it ) {
    G4Element* tmp = new G4Element(name = pyne::nucname::name(it->first),
                                   symbol = pyne::nucname::name(it->first),
                                   ncomponents = 1);
    tmp->AddIsotope(isotope_map[it->first], abundance = 100.0*perCent);
    element_map[it->first]=tmp;
  }

  return element_map;
}

// using each g4element make the geant4 version of each pyne material
std::map<std::string,G4Material*> get_g4materials(std::map<int,G4Element*> element_map,
    std::map<std::string, pyne::Material> material_library)
{
  std::map<std::string, G4Material*> material_map;
  std::map<std::string, pyne::Material>::iterator it;
  pyne::comp_iter mat_it;

  std::string name;
  int ncomponents;

  // loop over the PyNE material library instanciating as needed
  for ( it = material_library.begin() ; it != material_library.end() ; ++it ) {
    // get the material object
    pyne::Material mat = it->second;

    int num_nucs = 0;
    for ( mat_it = mat.comp.begin() ; mat_it != mat.comp.end() ; ++mat_it ) {
      if( (mat_it->second) > 0.0 ) {
        num_nucs++;
      }
    }


    // create the g4 material
    std::cout << mat.metadata["name"].asString() << std::endl;
    G4Material * g4mat = new G4Material(name = mat.metadata["name"].asString(),
                                        mat.density*g/(cm*cm*cm),
                                        ncomponents = num_nucs);
    // now iterate over the composiiton
    for ( mat_it = mat.comp.begin() ; mat_it != mat.comp.end() ; ++mat_it ) {
      if( (mat_it->second) > 0.0 ) { // strip out abundance 0 nuclides
        g4mat->AddElement(element_map[mat_it->first], mat_it->second);
      }
    }
    material_map[mat.metadata["name"].asString()]=g4mat;
  }

  // Add vacuum
  double z,a,density;
  G4Material* Vacuum =
      new G4Material("mat:Vacuum", z=1., a=1.0*g/mole, density=1.0e-20*mg/cm3);
  // add to lib
  material_map["mat:Vacuum"]=Vacuum;

  return material_map;
}
