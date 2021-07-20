#include "uwuw.hpp"

#ifndef _WIN32
#include <unistd.h>
#else
#include <filesystem>
#endif

#include <iostream>

// Empty Constructor
UWUW::UWUW() {
  num_tallies = 0;
  num_materials = 0;
};

// Default constructor
UWUW::UWUW(char* file) {
  std::string filename(file);

  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);

  if (!check_file_exists(full_filepath)) {
    std::cerr << "The file " << full_filepath
              << " does not exist or is read protected" << std::endl;
    exit(1);
  }

  // load materials
  material_library = load_pyne_materials(full_filepath);

  // load tallies
  tally_library = load_pyne_tallies(full_filepath);
};

// Default constructor
UWUW::UWUW(std::string filename) {
  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);

  // check for file existence
  if (!check_file_exists(full_filepath)) {
    std::cerr << "The file " << full_filepath
              << " does not exist or is read protected" << std::endl;
    exit(1);
  }

  // load materials
  material_library = load_pyne_materials(full_filepath);

  // load tallies
  tally_library = load_pyne_tallies(full_filepath);
};

// Destructor
UWUW::~UWUW(){};

// convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(char* filename) {
  std::string file(filename);
  return UWUW::get_full_filepath(file);
}

// convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(std::string filename) {
  // remove all extra whitespace
  filename.erase(std::remove(filename.begin(), filename.end(), ' '),
                 filename.end());
  // use stdlib call
#ifndef _WIN32
  const std::string full_filepath = realpath(filename.c_str(), NULL);
#else
  const std::string full_filepath =
      std::filesystem::canonical(filename.c_str()).string();
#endif
  return full_filepath;
}

// see if file exists
bool UWUW::check_file_exists(const std::string& filename) {
  // from http://stackoverflow.com/questions/12774207/
  // fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  std::ifstream infile(filename.c_str());
  return infile.good();
}

// loads all materials into map

pyne::MaterialLibrary UWUW::load_pyne_materials(std::string filename,
                                                std::string datapath) {
  pyne::MaterialLibrary library(filename, datapath);  // material library
  return library;
}

// loads all tallies into map
std::map<std::string, pyne::Tally> UWUW::load_pyne_tallies(
    std::string filename, std::string datapath) {
  std::map<std::string, pyne::Tally> library;  // material library

  if (!hdf5_path_exists(filename, datapath)) return library;

  num_tallies = get_length_of_table(filename, datapath);

  for (int i = 0; i < num_tallies; i++) {
    pyne::Tally tally;  // from file
    tally.from_hdf5(filename, datapath, i);
    library[tally.tally_name] = tally;
  }
  return library;
}

// see if path exists before we go on
bool UWUW::hdf5_path_exists(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath.c_str());

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return datapath_exists;
}

int UWUW::get_length_of_table(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);

  hsize_t arr_dims[1];
  H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return arr_dims[0];
}
