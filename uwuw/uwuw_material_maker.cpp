#include "uwuw.hpp"
#include "moab/ProgOptions.hpp"
#include "hdf5.h"

int main(int argc, char* argv[])
{

  ProgOptions po("uwuw_material_maker: a tool for creating materials and adding them to a UWUW Material Library");

  bool append = false;
  bool overwrite = false;

  std::string lib_file;
  std::string material_file;
  std::string out_file = "";

  po.addOpt<void>( "append,a", "Append to existing library", &append);
  po.addOpt<void>( "overwite,w", "Overwrite existing materials with the same name", &overwrite);

  po.addRequiredArg<std::string>("material_file", "Path to DAGMC file to proccess", &material_file);
  po.addOpt<std::string>("output,o", "Specify the output filename (default "")", &out_file);

  po.addOptionHelpHeading("Options for loading files");
  po.parseCommandLine(argc, argv);
  
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
  pyne::Material material = Material();
  // retreve it
  material.from_text(material_file);
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
    if(material_library.count(new_mat_name) != 0 && !overwrite) {
      std::cout << "A material already exists in the library with that name" << std::endl;
      std::cout << "please rename it, if you wish to overwrite it, use the -overwrite option" << std::end;
      exit(EXIT_FAILURE);
    }    
    // add the new material
    mat_lib[new_mat_name] = material;
    // orphan the /materials dataspace in the file.
    char* output_cstr = out_file.c_str();

    //Set file access properties so it closes cleanly                                                                            
    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);

    // open the file 
    hid_t id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);

    // datapaths to delete
    char* mat_path = "/materials";
    char* nuc_path = "/nucid";
    char* meta_path = "/metadata";

    hid_t lapl = 0;
    // delete the links
    herr_t error = H5Ldelete(id,mat_path,lapl);
    herr_t error = H5Ldelete(id,nuc_path,lapl);
    herr_t error = H5Ldelete(id,meta_path,lapl);

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
