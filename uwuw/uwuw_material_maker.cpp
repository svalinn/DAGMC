#include "uwuw.hpp"
#include "moab/ProgOptions.hpp"
#include "hdf5.h"
#include "pyne.h"

bool check_file_exists(std::string filename)
{
  std::ifstream infile(filename.c_str());
  return infile.good();
}


void print_mat(std::ostream& os, pyne::Material mat)
{
  //print the Mass Stream to stdout
  os << "\tMass: " << mat.mass << "\n";
  os << "\tDensity: " << mat.density << "\n";
  os << "\t---------\n";
  for(pyne::comp_iter i = mat.comp.begin(); i != mat.comp.end(); i++) {
    os << "\t" << pyne::nucname::name( i->first ) << "\t" << i->second << "\n";
  }
}

int main(int argc, char* argv[])
{

  ProgOptions po("uwuw_material_maker: a tool for creating materials and adding them to a UWUW Material Library");

  bool append = false;
  bool overwrite = false;
  bool print = false;
  bool expand = false;
  bool dump = false;

  std::string lib_file;
  std::string material_file;
  std::string out_file = "";

  po.addOpt<void>( "dump,d", "Dump material library contents to screen",&dump);
  po.addOpt<void>( "append,a", "Append to existing library", &append);
  po.addOpt<void>( "overwite,w", "Overwrite existing materials with the same name", &overwrite);
  po.addOpt<void>( "print,p", "Print the full material compositon", &print);
  po.addOpt<void>( "expand,e", "Expand elements of the material", &expand);

  po.addRequiredArg<std::string>("material_file", "Path to DAGMC file to proccess", &material_file);
  po.addOpt<std::string>("output,o", "Specify the output filename (default "")", &out_file);

  po.addOptionHelpHeading("Options for loading files");
  po.parseCommandLine(argc, argv);
 
  if(dump) {
    UWUW *uwuw = new UWUW(material_file);
    std::map<std::string,pyne::Material> mat_lib = uwuw->material_library;
    std::map<std::string,pyne::Material>::iterator it;
    std::ostringstream ostr;
    for ( it = mat_lib.begin() ; it != mat_lib.end() ; ++it ) {
      print_mat(ostr,it->second);
    }
    std::cout << ostr.str() << std::endl;
    exit(EXIT_SUCCESS);
  }
 
  if(out_file == "") {
    std::cout << "ouput file must be set" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(!append && check_file_exists(out_file)){
    std::cout << "output file already exists and you are not appending" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(append && !check_file_exists(out_file)){
    std::cout << "output doesnt already exist and you are appending" << std::endl;
    exit(EXIT_FAILURE);
  }


  
  // make a new material 
  pyne::Material material = pyne::Material();
  // retreve it
  material.from_text(material_file);
  //
  std::string mat_name = material.metadata["name"].asString();
  if(mat_name == "") {
    std::cout << "Material name not set" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // remove whitespace
    mat_name.erase(std::remove(mat_name.begin(), mat_name.end(), ' '), mat_name.end() );
    material.metadata["name"] = mat_name;
  }
  if (expand) material = material.expand_elements();
  std::ostringstream ostr;
  if (print) {
    std::cout << "Name: '"<< material.metadata["name"].asString() << "'" << std::endl; 
    print_mat(ostr,material);
    std::cout << ostr.str() << std::endl;
  }

  // need to key on Name 
  std::string new_mat_name = material.metadata["name"].asString();

    // write it
  if(!append) {
    material.write_hdf5(out_file,"/materials");
  } else {
    // if appending 
    // new uwuw class
    UWUW *uwuw = new UWUW(out_file);
    // get the library
    std::map<std::string,pyne::Material> mat_lib = uwuw->material_library;
    // check for a name clash
    if(mat_lib.count(new_mat_name) != 0 && !overwrite) {
      std::cout << "A material already exists in the library with that name" << std::endl;
      std::cout << "please rename it, if you wish to overwrite it, use the -overwrite option" << std::endl;
      exit(EXIT_FAILURE);
    }    
    // add the new material
    mat_lib[new_mat_name] = material;
    // orphan the /materials dataspace in the file.

    //Set file access properties so it closes cleanly                                                                            
    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);

    // open the file 
    hid_t id = H5Fopen(out_file.c_str(), H5F_ACC_RDWR, fapl);

    // datapaths to delete
    char* mat_path = "/materials";
    char* nuc_path = "/nucid";
    char* meta_path = "/materials_metadata";

    hid_t lapl = 0;
    // delete the links
    herr_t error = 0;
    error = H5Ldelete(id,mat_path,lapl);
    error = H5Ldelete(id,nuc_path,lapl);
    error = H5Ldelete(id,meta_path,lapl);

    // close the hdf5 file
    H5Fclose(id);

    // write the new materials
    std::map<std::string,pyne::Material>::iterator it;
    for ( it = mat_lib.begin() ; it != mat_lib.end() ; ++it ) {
      it->second.write_hdf5(out_file,"/materials");
    }
    // all done
  } 
  // thats all
  return 0;
}
