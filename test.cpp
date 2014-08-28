#include <iostream>
#include "DagSolidMaterial.hh"

int main(int argc, char* argv[])
{
  std::string filename = "atic_uwuw_zip.h5m";
  std::map<std::string,G4Material*> mat_mat = load_uwuw_materials(filename);
  return 0;
}
