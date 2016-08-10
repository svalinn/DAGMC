// This file is composed of the following original files:

//   license.txt
//   src/utils.h
//   src/extra_types.h
//   src/h5wrap.h
//   src/nucname.h
//   src/rxname.h
//   src/particle.h
//   src/data.h
//   src/json-forwards.h
//   src/json.h
//   src/jsoncustomwriter.h
//   src/material.h
//   src/tally.h
//   src/atomic_data.cpp~
//   src/atomic_data.h
//   src/atomic_data.h~

// PyNE amalgated header http://pyne.io/
#ifndef PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY
#define PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY

#define PYNE_IS_AMALGAMATED

//
// start of license.txt
//
// Copyright 2011-2015, the PyNE Development Team. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//    1. Redistributions of source code must retain the above copyright notice, this list of
//       conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright notice, this list
//       of conditions and the following disclaimer in the documentation and/or other materials
//       provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE PYNE DEVELOPMENT TEAM ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those of the
// authors and should not be interpreted as representing official policies, either expressed
// or implied, of the stakeholders of the PyNE project or the employers of PyNE developers.
//
// -------------------------------------------------------------------------------
// The files cpp/measure.cpp and cpp/measure.hpp are covered by:
//
// Copyright 2004 Sandia Corporation.  Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
// retains certain rights in this software.
//
// http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB
// //
// end of license.txt
//


//
// start of src/utils.h
//
/// \brief This is the base PyNE library.
///
/// It contains a lot of utility functions and constants that are globaly useful
/// through out the rest of the PyNE infrastructure.
///

// Header for general library file.

#ifndef PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA
#define PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <exception>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <vector>
#include <algorithm>

#if (__GNUC__ >= 4)
#include <cmath>
#define isnan(x) std::isnan(x)
#else
#include <math.h>
#define isnan(x) __isnand((double)x)
#endif

#ifdef __WIN_MSVC__
#define isnan(x) ((x) != (x))
#endif

#ifndef JSON_IS_AMALGAMATION
#define JSON_IS_AMALGAMATION
#endif

/// The 'pyne' namespace all PyNE functionality is included in.
namespace pyne
{

void pyne_start (); ///< Initializes PyNE based on environment.

/// Path to the directory containing the PyNE data.
extern std::string PYNE_DATA;
extern std::string NUC_DATA_PATH; ///< Path to the nuc_data.h5 file.

// String Transformations
/// string of digit characters
static std::string digits = "0123456789";
/// uppercase alphabetical characters
static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
/// string of all valid word characters for variable names in programing languages.
static std::string words = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_";

/// \name String Conversion Functions
/// \{
/// Converts the variables of various types to their C++ string representation.
std::string to_str(int t);
std::string to_str(unsigned int t);
std::string to_str(double t);
std::string to_str(bool t);
/// \}

int to_int(std::string s);  ///< Converts a string of digits to an int using atoi().

double to_dbl(std::string s);  ///< Converts a valid string to a float using atof().

/// Converts a string from ENDF format to a float. Only handles E-less format
/// but is roughly 5 times faster than endftod.
double endftod_cpp(char * s);
double endftod_f(char * s); ///< Converts a string from ENDF format to a float.
extern  double (*endftod)(char * s); ///< endftod function pointer. defaults to fortran

void use_fast_endftod();/// switches endftod to fast cpp version

/// Returns an all upper case copy of the string.
std::string to_upper(std::string s);

/// Returns an all lower case copy of the string.
std::string to_lower(std::string s);

/// Returns a capitalized copy of the string.
std::string capitalize(std::string s);

/// Finds and returns the first white-space delimited token of a line.
/// \param line a character array to take the first token from.
/// \param max_l an upper bound to the length of the token.  Must be 11 or less.
/// \returns a the flag as a string
std::string get_flag(char line[], int max_l);

/// Creates a copy of \a s with all instances of \a substr taken out.
std::string remove_substring(std::string s, std::string substr);

/// Removes all characters in the string \a chars from \a s.
std::string remove_characters(std::string s, std::string chars);

/// Replaces all instance of \a substr in \a s with \a repstr.
std::string replace_all_substrings(std::string s, std::string substr,
                                   std::string repstr);

/// Returns the last character in a string.
std::string last_char(std::string s);

/// Returns the slice of a string \a s using the negative index \a n and the
/// length of the slice \a l.
std::string slice_from_end(std::string s, int n=-1, int l=1);

/// Returns true if \a a <= \a b <= \a c and flase otherwise.
bool ternary_ge(int a, int b, int c);

/// Returns true if \a substr is in \a s.
bool contains_substring(std::string s, std::string substr);

/// Calculates a version of the string \a name that is also a valid variable name.
/// That is to say that the return value uses only word characters.
std::string natural_naming(std::string name);

/// Finds the slope of a line from the points (\a x1, \a y1) and (\a x2, \a y2).
double slope (double x2, double y2, double x1, double y1);

/// Solves the equation for the line y = mx + b, given \a x and the points that
/// form the line: (\a x1, \a y1) and (\a x2, \a y2).
double solve_line (double x, double x2, double y2, double x1, double y1);

double tanh(double x);  ///< The hyperbolic tangent function.
double coth(double x);  ///< The hyperbolic cotangent function.


// File Helpers
/// Returns true if the file can be found.
bool file_exists(std::string strfilename);

// Message Helpers
extern bool USE_WARNINGS;
/// Toggles warnings on and off
bool toggle_warnings();

/// Prints a warning message.
void warning(std::string s);

/// Custom exception to be thrown in the event that a required file is not able to
/// be found.
class FileNotFound : public std::exception
{
 public:

  /// default constructor
  FileNotFound () {};

  /// default destructor
  ~FileNotFound () throw () {};

  /// constructor with the filename \a fname.
  FileNotFound(std::string fname) {
    filename = fname;
  };

  /// Creates a helpful error message.
  virtual const char* what() const throw() {
    std::string FNFstr ("File not found: ");
    if (!filename.empty())
      FNFstr += filename;

    return (const char *) FNFstr.c_str();
  };

 private:
  std::string filename; ///< unfindable filename.
};


// End PyNE namespace
}

#endif  // PYNE_KMMHYNANYFF5BFMEYIP7TUNLHA
//
// end of src/utils.h
//


//
// start of src/extra_types.h
//
/// \file extra_types.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// Provides some extra types that may be generally useful

#if !defined(_XDRESS_EXTRA_TYPES_)
#define _XDRESS_EXTRA_TYPES_

#if defined(__cplusplus)
namespace extra_types
{
/// complex type struct, matching PyTables definition
//  typedef struct {
//    double re;  ///< real part
//    double im;  ///< imaginary part
//  } complex_t;

/// Chivalrously handles C++ memory issues that Cython does
/// not yet have a syntax for.  This is a template class,
/// rather than three template functions, because Cython does
/// not yet support template function wrapping.
template <class T>
class MemoryKnight
{
 public:
  MemoryKnight() {};  ///< Default constructor
  ~MemoryKnight() {}; ///< Default Destructor

  /// Creates a new instance of type T on the heap using
  /// its default constructor.
  /// \return T *
  T * defnew() {
    return new T();
  };

  /// Creates a new instance of type T, using T's default
  /// constructor, at a given location.
  /// \param void * ptr, location to create T instance
  /// \return value of ptr recast as T *
  T * renew(void * ptr) {
    return new (ptr) T();
  };

  /// Deallocates a location in memory using delete.
  /// \param T * ptr, location to remove
  void deall(T * ptr) {
    delete ptr;
  };
};

// End namespace extra_types
}

#elif defined(__STDC__)

// de nada

#endif


/// complex type struct, matching PyTables definition
typedef struct {
  double re;  ///< real part
  double im;  ///< imaginary part
} xd_complex_t;

#endif

//
// end of src/extra_types.h
//


//
// start of src/h5wrap.h
//
/// \file h5wrap.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Provides some HDF5 helper functionality in its own namespace

#ifndef PYNE_MRNAFG5GNZDNPCRPX3UCBZ5MFE
#define PYNE_MRNAFG5GNZDNPCRPX3UCBZ5MFE

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <exception>

#include "hdf5.h"

#ifndef PYNE_IS_AMALGAMATED
#include "extra_types.h"
#endif

//! Wrapper for standard HDF5 operations
namespace h5wrap
{
/// Custom exception for HDF5 indexing errors.
class HDF5BoundsError: public std::exception
{
  /// returns error message.
  virtual const char* what() const throw() {
    return "Index of point is out of bounds.  Cannot handle in HDF5 file.";
  };
};


/// Custom exception for when an existing file is not in a valid HDF5 format.
class FileNotHDF5: public std::exception
{
 public:

  /// default constructor
  FileNotHDF5() {};

  /// default destructor
  ~FileNotHDF5() throw () {};

  /// constructor with the filename
  FileNotHDF5(std::string fname) {
    filename = fname;
  };

  /// helpful error message that includes the filename
  virtual const char* what() const throw() {
    std::string FNH5str ("Not a valid HDF5 file: ");
    if (!filename.empty())
      FNH5str += filename;

    return (const char *) FNH5str.c_str();
  };

 private:
  std::string filename; ///< the file which is not in HDF5 format.
};


/// Custom exception for when a group cannot be found in an HDF5 file.
class GroupNotFound: public std::exception
{
 public:

  /// default constructor
  GroupNotFound() {};

  /// default destructor
  ~GroupNotFound() throw () {};

  /// constructor with the filename and the groupname
  GroupNotFound(std::string fname, std::string gname) {
    filename = fname;
  };

  /// helpful error message that includes the filename and the groupname
  virtual const char* what() const throw() {
    std::string msg ("the group ");
    msg += groupname;
    msg += " not found in the file ";
    msg += filename;
    return (const char *) msg.c_str();
  };

 private:
  std::string filename;   ///< the HDF5 file
  std::string groupname;  ///< the group in the hierarchy
};

/// Custom exception for when a path is not found in an HDF5 file
class PathNotFound: public std::exception
{
 public:

  /// default constructor
  PathNotFound() {};

  /// default destructor
  ~PathNotFound() throw () {};

  /// constructor with the filename and the pathname
  PathNotFound(std::string fname, std::string pname) {
    filename = fname;
    path = pname;
  };

  /// helpful error message that includes the filename and the pathname
  virtual const char* what() const throw() {
    std::string msg ("the path ");
    msg += path;
    msg += " was not found in the HDF5 file ";
    msg += filename;
    return (const char *) msg.c_str();
  };

 private:
  std::string filename; ///< the HDF5 file
  std::string path;     ///< the path in the file
};



// Read-in Functions

/// Retrieves the \a nth index out of the dataset \a dset (which has an HDF5
/// datatype \a dtype).  The value is returned as the C/C++ type given by \a T.
template <typename T>
T get_array_index(hid_t dset, int n, hid_t dtype=H5T_NATIVE_DOUBLE)
{
  hsize_t count  [1] = {1};
  hsize_t offset [1] = {static_cast<hsize_t>(n)};

  hid_t dspace = H5Dget_space(dset);
  hsize_t npoints = H5Sget_simple_extent_npoints(dspace);

  //Handle negative indices
  if (n < 0)
    offset[0] = offset[0] + npoints;

  //If still out of range we have a problem
  if (npoints <= offset[0])
    throw HDF5BoundsError();

  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  //Set memmory hyperspace
  hsize_t dimsm[1] = {1};
  hid_t memspace = H5Screate_simple(1, dimsm, NULL);

  hsize_t count_out  [1] = {1};
  hsize_t offset_out [1] = {0};

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
                      count_out, NULL);

  T data_out [1];
  H5Dread(dset, dtype, memspace, dspace, H5P_DEFAULT, data_out);

  return data_out[0];
}


// Conversion functions

/// Reads in data from an HDF5 file as a C++ set.  \a T should roughly match
/// \a dtype.
/// \param h5file HDF5 file id for an open file.
/// \param data_path path to the data in the open file.
/// \param dtype HDF5 data type for the data set at \a data_path.
/// \return an in memory set of type \a T.
template <typename T>
std::set<T> h5_array_to_cpp_set(hid_t h5file, std::string data_path, hid_t dtype=H5T_NATIVE_DOUBLE)
{
  std::set<T> cpp_set = std::set<T>();
  hsize_t arr_len[1];
  hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(dset);
  int arr_dim = H5Sget_simple_extent_dims(arr_space, arr_len, NULL);

  // Read in data from file to memory
  T * mem_arr = new T [arr_len[0]];
  H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

  // Load new values into the set
  cpp_set.insert(&mem_arr[0], &mem_arr[arr_len[0]]);

  H5Dclose(dset);

  delete[] mem_arr;
  return cpp_set;
}


/// Reads in data from an HDF5 file as a 1 dimiensional vector.  \a T should roughly
/// match \a dtype.
/// \param h5file HDF5 file id for an open file.
/// \param data_path path to the data in the open file.
/// \param dtype HDF5 data type for the data set at \a data_path.
/// \return an in memory 1D vector of type \a T.
template <typename T>
std::vector<T> h5_array_to_cpp_vector_1d(hid_t h5file, std::string data_path,
    hid_t dtype=H5T_NATIVE_DOUBLE)
{
  std::vector<T> cpp_vec;
  hsize_t arr_dims [1];
  hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(dset);
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  // Read in data from file to memory
  T mem_arr [arr_dims[0]];
  H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

  // Load new values into the vector
  cpp_vec.assign(mem_arr, mem_arr+arr_dims[0]);

  H5Dclose(dset);
  return cpp_vec;
}


/// Reads in data from an HDF5 file as a 2 dimiensional vector.  \a T should roughly
/// match \a dtype.
/// \param h5file HDF5 file id for an open file.
/// \param data_path path to the data in the open file.
/// \param dtype HDF5 data type for the data set at \a data_path.
/// \return an in memory 2D vector of type \a T.
template <typename T>
std::vector< std::vector<T> > h5_array_to_cpp_vector_2d(hid_t h5file, std::string data_path,
    hid_t dtype=H5T_NATIVE_DOUBLE)
{
  hsize_t arr_dims [2];
  hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(dset);
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  // Read in data from file to memory
  // Have to read in as 1D array to get HDF5 and new keyword
  // to play nice with each other
  T mem_arr [arr_dims[0] * arr_dims[1]];
  H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

  // Load new values into the vector of vectors, using some indexing tricks
  std::vector< std::vector<T> > cpp_vec (arr_dims[0], std::vector<T>(arr_dims[1]));
  for(int i = 0; i < arr_dims[0]; i++) {
    cpp_vec[i].assign(mem_arr+(i*arr_dims[1]), mem_arr+((i+1)*arr_dims[1]));
  };

  H5Dclose(dset);
  return cpp_vec;
}


/// Reads in data from an HDF5 file as a 3 dimiensional vector.  \a T should roughly
/// match \a dtype.
/// \param h5file HDF5 file id for an open file.
/// \param data_path path to the data in the open file.
/// \param dtype HDF5 data type for the data set at \a data_path.
/// \return an in memory 3D vector of type \a T.
template <typename T>
std::vector< std::vector< std::vector<T> > > h5_array_to_cpp_vector_3d(hid_t h5file,
    std::string data_path,
    hid_t dtype=H5T_NATIVE_DOUBLE)
{
  hsize_t arr_dims [3];
  hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(dset);
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  // Read in data from file to memory
  // Have to read in as 1D array to get HDF5 and new keyword
  // to play nice with each other
  T mem_arr [arr_dims[0] * arr_dims[1] * arr_dims[2]];
  H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

  // Load new values into the vector of vectors of vectors, using some indexing tricks
  std::vector< std::vector< std::vector<T> > > cpp_vec (arr_dims[0], std::vector< std::vector<T> >(arr_dims[1], std::vector<T>(arr_dims[2])));
  for(int i = 0; i < arr_dims[0]; i++) {
    for(int j = 0; j < arr_dims[1]; j++) {
      cpp_vec[i][j].assign(mem_arr+((i*arr_dims[1]*arr_dims[2]) + (j*arr_dims[2])), mem_arr+((i*arr_dims[1]*arr_dims[2]) + ((j+1)*arr_dims[2])));
    };
  };

  H5Dclose(dset);
  return cpp_vec;
}



// Classes
/// A class representing a high-level table contruct whose columns all have the same
/// type \a T in C/C++ (and the analogous type in HDF5).
template <typename T>
class HomogenousTypeTable
{
 public:

  /// default constructor
  HomogenousTypeTable() {};

  /// default destructor
  ~HomogenousTypeTable() {};

  /// Constructor to load in data upon initialization.  \a T should roughly
  /// match \a dtype.
  /// \param h5file HDF5 file id for an open file.
  /// \param data_path path to the data in the open file.
  /// \param dtype HDF5 data type for the data set at \a data_path.
  HomogenousTypeTable(hid_t h5file, std::string data_path, hid_t dtype=H5T_NATIVE_DOUBLE) {
    hid_t h5_set = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);
    hid_t h5_space = H5Dget_space(h5_set);
    hid_t h5_type = H5Dget_type(h5_set);

    // set path
    path = data_path;

    // set shape
    shape[0] = H5Sget_simple_extent_npoints(h5_space);
    shape[1] = H5Tget_nmembers(h5_type);

    // set cols
    std::string * cols_buf = new std::string [shape[1]];
    for(int n = 0; n < shape[1]; n++)
      cols_buf[n] = H5Tget_member_name(h5_type, n);
    cols.assign(cols_buf, cols_buf+shape[1]);

    // set data
    hid_t col_type;
    T * col_buf = new T [shape[0]];

    data.clear();
    for(int n = 0; n < shape[1]; n++) {
      // Make a compound data type of just this column
      col_type = H5Tcreate(H5T_COMPOUND, sizeof(T));
      H5Tinsert(col_type, cols[n].c_str(), 0, dtype);

      // Read in this column
      H5Dread(h5_set, col_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, col_buf);

      // save this column as a vector in out data map
      data[cols[n]] = std::vector<T>(col_buf, col_buf+shape[0]);
    };
    delete[] col_buf;
  };

  // Metadata attributes
  std::string path; ///< path in file to the data
  int shape [2];    ///< table shape, rows x columns.
  std::vector<std::string> cols;  ///< column names
  /// mapping from column names to column data
  std::map<std::string, std::vector<T> > data;

  //
  // operator overloads
  //
  /// index into the table by column name (string)
  std::vector<T> operator[] (std::string col_name) {
    return data[col_name];
  };

  /// index into the table by row
  std::map<std::string, T> operator[] (int m) {
    std::map<std::string, T> row = std::map<std::string, T>();

    for(int n = 0; n < shape[1]; n++)
      row[cols[n]] = data[cols[n]][m];

    return row;
  };
};


/// Create an HDF5 data type for complex 128 bit data, which happens to match the
/// complex data type that is used by PyTables ^_~.
inline hid_t _get_PYTABLES_COMPLEX128()
{
  hid_t ct = H5Tcreate(H5T_COMPOUND, sizeof(xd_complex_t));
  H5Tinsert(ct, "r", HOFFSET(xd_complex_t, re), H5T_NATIVE_DOUBLE);
  H5Tinsert(ct, "i", HOFFSET(xd_complex_t, im), H5T_NATIVE_DOUBLE);
  return ct;
}

/// The HDF5 id for a complex data type compatible with PyTables generated data.
static hid_t PYTABLES_COMPLEX128 = _get_PYTABLES_COMPLEX128();


/// Determines if a path exists in an hdf5 file.
/// \param h5file HDF5 file id for an open file.
/// \param path path to the data in the open file.
/// \return true or false
inline bool path_exists(hid_t h5file, std::string path)
{
  bool rtn = false;
  hid_t ds = H5Dopen2(h5file, path.c_str(), H5P_DEFAULT);
  if (0 <= ds) {
    rtn = true;
    H5Dclose(ds);
  } else {
    hid_t grp = H5Gopen2(h5file, path.c_str(), H5P_DEFAULT);
    if (0 <= grp) {
      rtn = true;
      H5Gclose(grp);
    }
  }
  return rtn;
}


// End namespace h5wrap
}



#endif
//
// end of src/h5wrap.h
//


//
// start of src/nucname.h
//
/// \brief Converts between naming conventions for nuclides.

#ifndef PYNE_D35WIXV5DZAA5LLOWBY2BL2DPA
#define PYNE_D35WIXV5DZAA5LLOWBY2BL2DPA
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#endif

namespace pyne
{
//! Nuclide naming conventions
namespace nucname
{
typedef std::string name_t; ///< name type
typedef int zz_t;           ///< Z number type

typedef std::map<name_t, zz_t> name_zz_t; ///< name and Z num map type
typedef name_zz_t::iterator name_zz_iter; ///< name and Z num iter type
name_zz_t get_name_zz();  ///< Creates standard name to Z number mapping.
extern name_zz_t name_zz; ///< name to Z num map

typedef std::map<zz_t, name_t> zzname_t;  ///< Z num to name map type
typedef zzname_t::iterator zzname_iter;   ///< Z num to name iter type
zzname_t get_zz_name();   ///< Creates standard Z number to name mapping.
extern zzname_t zz_name;  ///< Z num to name map

name_zz_t get_fluka_zz();  ///< Creates standard fluka-name to nucid mapping.
extern name_zz_t fluka_zz; ///< fluka-name to nucid map
zzname_t get_zz_fluka();   ///< Creates standard nucid to fluka-name mapping.
extern zzname_t zz_fluka;  ///< nucid to fluka-name map
/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

/// name grouping type (for testing containment)
typedef std::set<name_t> name_group;
typedef name_group::iterator name_group_iter; ///< name grouping iter type

/// Z number grouping type (for testing containment)
typedef std::set<zz_t> zz_group;
typedef zz_group::iterator zz_group_iter; ///< Z number grouping iter

/// Converts a name group to a Z number group.
/// \param eg a grouping of nuclides by name
/// \return a Z numbered group
zz_group name_to_zz_group (name_group eg);

extern name_t LAN_array[15];  ///< array of lanthanide names
extern name_group LAN;        ///< lanthanide name group
extern zz_group lan;          ///< lanthanide Z number group

extern name_t ACT_array[15];  ///< array of actinide names
extern name_group ACT;        ///< actinide name group
extern zz_group act;          ///< actinide Z number group

extern name_t TRU_array[22];  ///< array of transuranic names
extern name_group TRU;        ///< transuranic name group
extern zz_group tru;          ///< transuranic Z number group

extern name_t MA_array[10];   ///< array of minor actinide names
extern name_group MA;         ///< minor actinide name group
extern zz_group ma;           ///< minor actinide Z number group

extern name_t FP_array[88];   ///< array of fission product names
extern name_group FP;         ///< fission product name group
extern zz_group fp;           ///< fission product Z number group


/******************/
/*** Exceptions ***/
/******************/

/// Custom expection for declaring that a value does not follow a recognizable
/// nuclide naming convention.
class NotANuclide : public std::exception
{
 public:
  /// default constructor
  NotANuclide () {};

  /// default destructor
  ~NotANuclide () throw () {};

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  NotANuclide(std::string wasptr, std::string nowptr) {
    nucwas = wasptr;
    nucnow = nowptr;
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  NotANuclide(std::string wasptr, int nowptr) {
    nucwas = wasptr;
    nucnow = pyne::to_str(nowptr);
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  NotANuclide(int wasptr, std::string nowptr) {
    nucwas = pyne::to_str(wasptr);
    nucnow = nowptr;
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  NotANuclide(int wasptr, int nowptr) {
    nucwas = pyne::to_str(wasptr);
    nucnow = pyne::to_str(nowptr);
  };

  /// Generates an informational message for the exception
  /// \return The error string
  virtual const char* what() const throw() {
    std::string NaNEstr ("Not a Nuclide! ");
    if (!nucwas.empty())
      NaNEstr += nucwas;

    if (!nucnow.empty()) {
      NaNEstr += " --> ";
      NaNEstr += nucnow;
    }
    return (const char *) NaNEstr.c_str();
  };

 private:
  std::string nucwas; ///< previous nuclide state
  std::string nucnow; ///< current nuclide state
};

/// Custom expection for declaring that a value represents one or more nuclides
/// in one or more namig conventions
class IndeterminateNuclideForm : public std::exception
{
 public:
  /// default constructor
  IndeterminateNuclideForm () {};

  /// default destuctor
  ~IndeterminateNuclideForm () throw () {};

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  IndeterminateNuclideForm(std::string wasptr, std::string nowptr) {
    nucwas = wasptr;
    nucnow = nowptr;
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  IndeterminateNuclideForm(std::string wasptr, int nowptr) {
    nucwas = wasptr;
    nucnow = pyne::to_str(nowptr);
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  IndeterminateNuclideForm(int wasptr, std::string nowptr) {
    nucwas = pyne::to_str(wasptr);
    nucnow = nowptr;
  };

  /// Constructor given previous and current state of nulide name
  /// \param wasptr Previous state, typically user input.
  /// \param nowptr Current state, as far as PyNE could get.
  IndeterminateNuclideForm(int wasptr, int nowptr) {
    nucwas = pyne::to_str(wasptr);
    nucnow = pyne::to_str(nowptr);
  };

  /// Generates an informational message for the exception
  /// \return The error string
  virtual const char* what() const throw() {
    std::string INFEstr ("Indeterminate nuclide form: ");
    if (!nucwas.empty())
      INFEstr += nucwas;

    if (!nucnow.empty()) {
      INFEstr += " --> ";
      INFEstr += nucnow;
    }
    return (const char *) INFEstr.c_str();
  }

 private:
  std::string nucwas; ///< previous nuclide state
  std::string nucnow; ///< current nuclide state
};

/// \name isnuclide functions
/// \{
/// These functions test if an input \a nuc is a valid nuclide.
/// \param nuc a possible nuclide
/// \return a bool
bool isnuclide(std::string nuc);
bool isnuclide(const char * nuc);
bool isnuclide(int nuc);
/// \}

/// \name iselement functions
/// \{
/// These functions test if an input \a nuc is a valid element.
/// \param nuc a possible element
/// \return a bool
bool iselement(std::string nuc);
bool iselement(const char * nuc);
bool iselement(int nuc);

/// \}
/// \name Identifier Form Functions
/// \{
/// The 'id' nuclide naming convention is the canonical form for representing
/// nuclides in PyNE. This is termed a ZAS, or ZZZAAASSSS, representation because
/// It stores 3 Z-number digits, 3 A-number digits, followed by 4 S-number digits
/// which the nucleus excitation state.
///
/// The id() function will always return an nuclide in id form, if successful.
/// If the input nuclide is in id form already, then this is function does no
/// work. For all other formats, the id() function provides a best-guess based
/// on a heirarchy of other formats that is used to resolve ambiguities between
/// naming conventions. For integer input the form resolution order is:
///   - id
///   - zz (elemental z-num only given)
///   - zzaaam
///   - cinder (aaazzzm)
///   - mcnp
///   - zzaaa
/// For string (or char *) input the form resolution order is as follows:
///   - ZZ-LL-AAAM
///   - Integer form in a string representation, uses interger resolution
///   - NIST
///   - name form
///   - Serpent
///   - LL (element symbol)
/// For well-defined situations where you know ahead of time what format the
/// nuclide is in, you should use the various form_to_id() functions, rather
/// than the id() function which is meant to resolve possibly ambiquous cases.
/// \param nuc a nuclide
/// \return nucid 32-bit integer identifier
int id(int nuc);
int id(const char * nuc);
int id(std::string nuc);
/// \}

/// \name Name Form Functions
/// \{
/// The 'name' nuclide naming convention is the more common, human readable
/// notation. The chemical symbol (one or two characters long) is first, followed
/// by the nucleon number. Lastly if the nuclide is metastable, the letter M is
/// concatenated to the end. For example, ‘H-1’ and ‘Am242M’ are both valid.
/// Note that nucname will always return name form with dashes removed, the
/// chemical symbol used for letter casing (ie 'Pu'), and a trailing upercase 'M'
/// for a metastable flag. The name() function first converts functions to id form
/// using the id() function. Thus the form order resolution for id() also applies
/// here.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
std::string name(int nuc);
std::string name(const char * nuc);
std::string name(std::string nuc);
/// \}

/// \name Z-Number Functions
/// \{
/// The Z-number, or charge number, represents the number of protons in a
/// nuclide.  This function returns that number.
/// \param nuc a nuclide
/// \return an integer Z-number.
int znum(int nuc);
int znum(const char * nuc);
int znum(std::string nuc);
/// \}

/// \name A-Number Functions
/// \{
/// The A-number, or nucleon number, represents the number of protons and
/// neutrons in a nuclide.  This function returns that number.
/// \param nuc a nuclide
/// \return an integer A-number.
int anum(int nuc);
int anum(const char * nuc);
int anum(std::string nuc);
/// \}

/// \name S-Number Functions
/// \{
/// The S-number, or excitation state number, represents the excitation
/// level of a nuclide.  Normally, this is zero.  This function returns
/// that number.
/// \param nuc a nuclide
/// \return an integer A-number.
int snum(int nuc);
int snum(const char * nuc);
int snum(std::string nuc);
/// \}

/// \name ZZAAAM Form Functions
/// \{
/// The ZZAAAM nuclide naming convention is the former canonical form for
/// nuclides in PyNE. This places the charge of the nucleus out front, then has
/// three digits for the atomic mass number, and ends with a metastable flag
/// (0 = ground, 1 = first excited state, 2 = second excited state, etc).
/// Uranium-235 here would be expressed as ‘922350’.
/// \param nuc a nuclide
/// \return an integer nuclide identifier.
int zzaaam(int nuc);
int zzaaam(const char * nuc);
int zzaaam(std::string nuc);
/// \}

/// \name ZZAAAM Form to Identifier Form Functions
/// \{
/// This converts from the ZZAAAM nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in ZZAAAM form.
/// \return an integer id nuclide identifier.
int zzaaam_to_id(int nuc);
int zzaaam_to_id(const char * nuc);
int zzaaam_to_id(std::string nuc);
/// \}


/// \name ZZZAAA Form Functions
/// \{
/// The ZZZAAA nuclide naming convention is a form in which the nuclides three
///digit ZZZ number is followed by the 3 digit AAA number.  If the ZZZ number
///is 2 digits, the preceding zeros are not included.
/// Uranium-235 here would be expressed as ‘92235’.
/// \param nuc a nuclide
/// \return an integer nuclide identifier.
int zzzaaa(int nuc);
int zzzaaa(const char * nuc);
int zzzaaa(std::string nuc);
/// \}


/// \name ZZZAAA Form to Identifier Form Functions
/// \{
/// This converts from the ZZZAAA nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in ZZZAAA form.
/// \return an integer id nuclide identifier.
int zzzaaa_to_id(int nuc);
int zzzaaa_to_id(const char * nuc);
int zzzaaa_to_id(std::string nuc);
/// \}


/// \name ZZLLAAAM Form Functions
/// \{
/// The ZZLLAAAM nuclide naming convention is a form in which the nuclides
/// AA number is followed by the redundant two LL characters, followed by
/// the nuclides ZZZ number.  Can also be followed with a metastable flag.
/// Uranium-235 here would be expressed as ‘92-U-235’.
/// \param nuc a nuclide
/// \return an integer nuclide identifier.
std::string zzllaaam(int nuc);
std::string zzllaaam(const char * nuc);
std::string zzllaaam(std::string nuc);
/// \}


/// \name ZZLLAAAM Form to Identifier Form Functions
/// \{
/// This converts from the ZZLLAAAM nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in ZZLLAAAM form.
/// \return an integer id nuclide identifier.
//int zzllaaam_to_id(int nuc);
int zzllaaam_to_id(const char * nuc);
int zzllaaam_to_id(std::string nuc);
/// \}


/// \name MCNP Form Functions
/// \{
/// This is the naming convention used by the MCNP suite of codes.
/// The MCNP format for entering nuclides is unfortunately non-standard.
/// In most ways it is similar to zzaaam form, except that it lacks the metastable
/// flag. For information on how metastable isotopes are named, please consult the
/// MCNP documentation for more information.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
int mcnp(int nuc);
int mcnp(const char * nuc);
int mcnp(std::string nuc);
/// \}

/// \name MCNP Form to Identifier Form Functions
/// \{
/// This converts from the MCNP nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in MCNP form.
/// \return an integer id nuclide identifier.
int mcnp_to_id(int nuc);
int mcnp_to_id(const char * nuc);
int mcnp_to_id(std::string nuc);
/// \}

/// \name FLUKA Form Functions
/// \{
/// This is the naming convention used by the FLUKA suite of codes.
/// The FLUKA format for entering nuclides requires some knowledge of FLUKA
/// The nuclide in must cases should be the atomic # times 10000000.
/// The exceptions are for FLUKA's named isotopes
/// See the FLUKA Manual for more information.
/// \param nuc a nuclide
/// \return the received FLUKA name
std::string fluka(int nuc);
/// \}

/// \name FLUKA Form to Identifier Form Functions
/// \{
/// This converts from the FLUKA name to the
/// id canonical form  for nuclides in PyNE.
/// \param name a fluka name
/// \return an integer id nuclide identifier.
int fluka_to_id(std::string name);
int fluka_to_id(char * name);
/// \}

/// \name Serpent Form Functions
/// \{
/// This is the string-based naming convention used by the Serpent suite of codes.
/// The serpent naming convention is similar to name form. However, only the first
/// letter in the chemical symbol is uppercase, the dash is always present, and the
/// the meta-stable flag is lowercase. For instance, ‘Am-242m’ is the valid serpent
/// notation for this nuclide.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
std::string serpent(int nuc);
std::string serpent(const char * nuc);
std::string serpent(std::string nuc);
/// \}

/// \name Serpent Form to Identifier Form Functions
/// \{
/// This converts from the Serpent nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in Serpent form.
/// \return an integer id nuclide identifier.
//int serpent_to_id(int nuc);  Should be ZAID
int serpent_to_id(const char * nuc);
int serpent_to_id(std::string nuc);
/// \}

/// \name NIST Form Functions
/// \{
/// This is the string-based naming convention used by NIST.
/// The NIST naming convention is also similar to the Serpent form. However, this
/// convention contains no metastable information. Moreover, the A-number comes
/// before the element symbol. For example, ‘242Am’ is the valid NIST notation.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
std::string nist(int nuc);
std::string nist(const char * nuc);
std::string nist(std::string nuc);
/// \}

/// \name NIST Form to Identifier Form Functions
/// \{
/// This converts from the NIST nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in NIST form.
/// \return an integer id nuclide identifier.
//int serpent_to_id(int nuc);  NON-EXISTANT
int nist_to_id(const char * nuc);
int nist_to_id(std::string nuc);
/// \}

/// \name CINDER Form Functions
/// \{
/// This is the naming convention used by the CINDER burnup library.
/// The CINDER format is similar to zzaaam form except that the placement of the
/// Z- and A-numbers are swapped. Therefore, this format is effectively aaazzzm.
/// For example, ‘2420951’ is the valid cinder notation for ‘AM242M’.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
int cinder(int nuc);
int cinder(const char * nuc);
int cinder(std::string nuc);
/// \}

/// \name Cinder Form to Identifier Form Functions
/// \{
/// This converts from the Cinder nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in Cinder form.
/// \return an integer id nuclide identifier.
int cinder_to_id(int nuc);
int cinder_to_id(const char * nuc);
int cinder_to_id(std::string nuc);
/// \}

/// \name ALARA Form Functions
/// \{
/// This is the format used in the ALARA activation code elements library.
/// For elements, the form is "ll" where ll is the atomic symbol. For isotopes
/// the form is "ll:AAA". No metastable isotope flag is used.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
std::string alara(int nuc);
std::string alara(const char * nuc);
std::string alara(std::string nuc);
/// \}

/// \name ALARA Form to Identifier Form Functions
/// \{
/// This converts from the ALARA nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in ALARA form.
/// \return an integer id nuclide identifier.
//int alara_to_id(int nuc); NOT POSSIBLE
int alara_to_id(const char * nuc);
int alara_to_id(std::string nuc);
/// \}

/// \name SZA Form Functions
/// \{
/// This is the new format for ACE data tables in the form SSSZZZAAA.
/// The first three digits represent the excited state (000 = ground,
/// 001 = first excited state, 002 = second excited state, etc).
/// The second three digits are the atomic number and the last three
/// digits are the atomic mass. Prepending zeros can be omitted, making
/// the SZA form equal to the MCNP form for non-excited nuclides.
/// \param nuc a nuclide
/// \return a string nuclide identifier.
int sza(int nuc);
int sza(const char * nuc);
int sza(std::string nuc);
/// \}

/// \name SZA Form to Identifier Form Functions
/// \{
/// This converts from the SZA nuclide naming convention
/// to the id canonical form  for nuclides in PyNE.
/// \param nuc a nuclide in SZA form.
/// \return an integer id nuclide identifier.
int sza_to_id(int nuc);
int sza_to_id(const char * nuc);
int sza_to_id(std::string nuc);
/// \}

/// \name Ground State Form Functions
/// \{
/// This form stores the nuclide in id form, but removes
/// the state information about the nuclide.  I is in the same
/// form as ID, but the four last digits are all zeros.
/// \param nuc a nuclide
/// \return a integer groundstate id
inline int groundstate(int nuc)
{
  return (id(nuc) / 10000 ) * 10000;
}
inline int groundstate(std::string nuc)
{
  return groundstate(id(nuc));
}
inline int groundstate(const char * nuc)
{
  return groundstate(std::string(nuc));
}
/// \}

/// \name State Map functions
/// \{
/// These convert from/to decay state ids (used in decay data)
/// to metastable ids (the PyNE default)
void _load_state_map();
int state_id_to_id(int state);
int id_to_state_id(int nuc_id);
extern std::map<int, int> state_id_map;
/// \}

/// \name ENSDF Form Functions
/// \{
/// This converts id's stored using standard ensdf syntax to nuc_id's
/// \param ensdf nuc string
/// \return PyNE nuc_id
int ensdf_to_id(const char * nuc);
int ensdf_to_id(std::string nuc);
/// \}

}
}

#endif  // PYNE_D35WIXV5DZAA5LLOWBY2BL2DPA
//
// end of src/nucname.h
//


//
// start of src/rxname.h
//
/// \brief Converts between naming conventions for reaction channels.

#ifndef PYNE_7DOEB2PKSBEFFIA3Q2NARI3KFY
#define PYNE_7DOEB2PKSBEFFIA3Q2NARI3KFY
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#include "nucname.h"
#endif

/// Number of reactions supported by default.
#define NUM_RX_NAMES 572

namespace pyne
{
//! Converts between naming conventions for reaction channels.
namespace rxname
{
extern std::string _names[NUM_RX_NAMES];  ///< Raw array of reaction names
/// Set of reaction names, must be valid variable names.
extern std::set<std::string> names;
/// Mapping from reaction ids to reaction names.
extern std::map<unsigned int, std::string> id_name;
/// Mapping from reaction names to reaction ids.
extern std::map<std::string, unsigned int> name_id;
/// Mapping between alternative names for reactions and the reaction id.
extern std::map<std::string, unsigned int> altnames;
/// Mapping from reaction ids to MT numbers.
extern std::map<unsigned int, unsigned int> id_mt;
/// Mapping from MT numbers to reaction names.
extern std::map<unsigned int, unsigned int> mt_id;
/// Mapping from reaction ids to labels (short descriptions).
extern std::map<unsigned int, std::string> labels;
/// Mapping from reaction ids to documentation strings (long descriptions).
extern std::map<unsigned int, std::string> docs;
/// Mapping from particle type and offset pairs to reaction ids.
/// Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
extern std::map<std::pair<std::string, int>, unsigned int> offset_id;
/// Mapping from particle type and reaction ids to offsets.
/// Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
extern std::map<std::pair<std::string, unsigned int>, int> id_offset;

/// A helper function to set the contents of the variables in this library.
void * _fill_maps();
extern void * _;  ///< A dummy variable used when calling #_fill_maps().

/// A helper function to compute nuclide id offsets from z-, a-, and s- deltas
inline int offset(int dz, int da, int ds=0)
{
  return dz*10000000 + da*10000 + ds;
}

/// \name Hash Functions
/// \{
/// Custom hash function for reaction name to reaction ids.
/// This functions will not return a value less than 1000, effectively reserving
/// space for the MT numbers.
unsigned int hash(std::string s);
unsigned int hash(const char * s);
/// \}

/// \name Name Functions
/// \{
/// Returns the canonical name of a reaction channel.
/// \param n Integer input of possible reaction, nominally an id or MT number.
/// \param s String input of possible reaction, often a reaction or alternate name.
/// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
///                 an integer it must be in id form.
/// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
///               an integer it must be in id form.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
std::string name(int n);
std::string name(unsigned int n);
std::string name(char * s);
std::string name(std::string s);
std::string name(int from_nuc, int to_nuc, std::string z="n");
std::string name(int from_nuc, std::string to_nuc, std::string z="n");
std::string name(std::string from_nuc, int to_nuc, std::string z="n");
std::string name(std::string from_nuc, std::string to_nuc, std::string z="n");
/// \}

/// \name ID Functions
/// \{
/// Returns the recation id of a reaction channel.  This id has been precomputed
/// from the hash of the name.
/// \param x Input reaction specification, may be a reaction name, alternate name,
///          an id, or an MT number.
/// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
///                 an integer it must be in id form.
/// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
///               an integer it must be in id form.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
unsigned int id(int x);
unsigned int id(unsigned int x);
unsigned int id(const char * x);
unsigned int id(std::string x);
unsigned int id(int from_nuc, int to_nuc, std::string z="n");
unsigned int id(int from_nuc, std::string to_nuc, std::string z="n");
unsigned int id(std::string from_nuc, int to_nuc, std::string z="n");
unsigned int id(std::string from_nuc, std::string to_nuc, std::string z="n");
/// \}

/// \name MT Number Functions
/// \{
/// Returns the MT number of a reaction channel, if available.
/// \param x Input reaction specification, may be a reaction name, alternate name,
///          an id, or an MT number.
/// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
///                 an integer it must be in id form.
/// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
///               an integer it must be in id form.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
unsigned int mt(int x);
unsigned int mt(unsigned int x);
unsigned int mt(char * x);
unsigned int mt(std::string x);
unsigned int mt(int from_nuc, int to_nuc, std::string z="n");
unsigned int mt(int from_nuc, std::string to_nuc, std::string z="n");
unsigned int mt(std::string from_nuc, int to_nuc, std::string z="n");
unsigned int mt(std::string from_nuc, std::string to_nuc, std::string z="n");
/// \}

//// \name Label Functions
/// \{
/// Returns a short description of a reaction channel.
/// \param x Input reaction specification, may be a reaction name, alternate name,
///          an id, or an MT number.
/// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
///                 an integer it must be in id form.
/// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
///               an integer it must be in id form.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
std::string label(int x);
std::string label(unsigned int x);
std::string label(char * x);
std::string label(std::string x);
std::string label(int from_nuc, int to_nuc, std::string z="n");
std::string label(int from_nuc, std::string to_nuc, std::string z="n");
std::string label(std::string from_nuc, int to_nuc, std::string z="n");
std::string label(std::string from_nuc, std::string to_nuc, std::string z="n");
/// \}

/// \name Documentation Functions
/// \{
/// Returns a short description of a reaction channel.
/// \param x Input reaction specification, may be a reaction name, alternate name,
///          an id, or an MT number.
/// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
///                 an integer it must be in id form.
/// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
///               an integer it must be in id form.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
std::string doc(int x);
std::string doc(unsigned int x);
std::string doc(char * x);
std::string doc(std::string x);
std::string doc(int from_nuc, int to_nuc, std::string z="n");
std::string doc(int from_nuc, std::string to_nuc, std::string z="n");
std::string doc(std::string from_nuc, int to_nuc, std::string z="n");
std::string doc(std::string from_nuc, std::string to_nuc, std::string z="n");
/// \}

/// \name Child Functions
/// \{
/// Returns the child nuclide comming from a parent for a reaction channel.
/// \param nuc Nuclide after reaction occurs.  When \a nuc is
///               an integer it must be in id form.
/// \param rx Input reaction specification, may be a reaction name, alternate name,
///           an id, or an MT number.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
int parent(int nuc, unsigned int rx, std::string z="n");
int parent(int nuc, std::string rx, std::string z="n");
int parent(std::string nuc, unsigned int rx, std::string z="n");
int parent(std::string nuc, std::string rx, std::string z="n");
/// \}

/// \name Parent Functions
/// \{
/// Returns the parent nuclide comming for a child and a given reaction channel.
/// \param nuc Initial target nuclide prior to reaction.  When \a nuc is
///            an integer it must be in id form.
/// \param rx Input reaction specification, may be a reaction name, alternate name,
///           an id, or an MT number.
/// \param z Flag for incident particle type.
///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
int child(int nuc, unsigned int rx, std::string z="n");
int child(int nuc, std::string rx, std::string z="n");
int child(std::string nuc, unsigned int rx, std::string z="n");
int child(std::string nuc, std::string rx, std::string z="n");
/// \}

/// Custom exception for declaring a value not to be a valid reaction.
class NotAReaction : public std::exception
{
 public:

  /// default constructor
  NotAReaction () {};

  /// default destructor
  ~NotAReaction () throw () {};

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(std::string wasptr, std::string nowptr) {
    rxwas = wasptr;
    rxnow = nowptr;
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(std::string wasptr, int nowptr) {
    rxwas = wasptr;
    rxnow = pyne::to_str(nowptr);
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(int wasptr, std::string nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = nowptr;
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(int wasptr, int nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = pyne::to_str(nowptr);
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(std::string wasptr, unsigned int nowptr) {
    rxwas = wasptr;
    rxnow = pyne::to_str(nowptr);
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(unsigned int wasptr, std::string nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = nowptr;
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  NotAReaction(unsigned int wasptr, unsigned int nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = pyne::to_str(nowptr);
  };

  /// Returns a helpful error message containing prior and current reaction state.
  virtual const char* what() const throw() {
    std::string narxstr ("Not a reaction! ");
    if (!rxwas.empty())
      narxstr += rxwas;

    if (!rxnow.empty()) {
      narxstr += " --> ";
      narxstr += rxnow;
    }
    return (const char *) narxstr.c_str();
  };

 private:
  std::string rxwas;  ///< previous reaction state
  std::string rxnow;  ///< current reaction state
};



/// Custom exception for declaring a value not to be of ambiquous reaction form.
class IndeterminateReactionForm : public std::exception
{
 public:

  /// default constructor
  IndeterminateReactionForm () {};

  /// default destructor
  ~IndeterminateReactionForm () throw () {};

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  IndeterminateReactionForm(std::string wasptr, std::string nowptr) {
    rxwas = wasptr;
    rxnow = nowptr;
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  IndeterminateReactionForm(std::string wasptr, int nowptr) {
    rxwas = wasptr;
    rxnow = pyne::to_str(nowptr);
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  IndeterminateReactionForm(int wasptr, std::string nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = nowptr;
  };

  /// Constructor using original reaction (\a wasptr) and the eventual state
  /// that PyNE calculated (\a nowptr).
  IndeterminateReactionForm(int wasptr, int nowptr) {
    rxwas = pyne::to_str(wasptr);
    rxnow = pyne::to_str(nowptr);
  };

  /// Returns a helpful error message containing prior and current reaction state.
  virtual const char* what() const throw() {
    std::string INFEstr ("Indeterminate reaction form: ");
    if (!rxwas.empty())
      INFEstr += rxwas;

    if (!rxnow.empty()) {
      INFEstr += " --> ";
      INFEstr += rxnow;
    }
    return (const char *) INFEstr.c_str();
  }

 private:
  std::string rxwas;  ///< previous reaction state
  std::string rxnow;  ///< current reaction state
};
}
}

#endif  // PYNE_7DOEB2PKSBEFFIA3Q2NARI3KFY
//
// end of src/rxname.h
//


//
// start of src/particle.h
//
/// \brief Converts betweeen naming/numbering conventions for particle types

// defines the primary particle types that are allowed by most monte carlo codes
// some monte carlo codes allow us to score so called "heavy ions", in fact we
// define heavy ions to be particles with more than one neutron or proton

#ifndef PYNE_UWPOP4EE6BEB5CZK4BQQHYTCEI
#define PYNE_UWPOP4EE6BEB5CZK4BQQHYTCEI

#include <string>
#include <map>
#include <set>

#ifndef PYNE_IS_AMALGAMATED
#include "nucname.h"
#endif

/// Number of pure particle types currently supported
#define NUM_PARTICLES 32

namespace pyne
{
//! Converts betweeen naming/numbering conventions for particle types
namespace particle
{
extern int _pdcids[NUM_PARTICLES]; ///
/// set of Particle Data Centre integer id numbers
extern std::string _docs[NUM_PARTICLES];
/// set of doc strings that describe the particle types
extern std::string _names[NUM_PARTICLES];
/// set of name strings that are the particle types
extern std::set<std::string> names;
/// set of valid names
extern std::set<int> pdc_nums;
/// set of valid pdc numbers
extern std::map<std::string,int> name_id;
/// map of name to pdc number
extern std::map<int,std::string> id_name;
/// map of pdc number to name
extern std::map<std::string,std::string> docs;
/// map of name to doc string
extern std::map<std::string,int> altnames;
/// map of alternative name to pdc number
extern std::map<std::string,std::string> part_to_mcnp;
/// map of name to mcnp string
extern std::map<std::string,std::string> part_to_mcnp6;
/// map of name to mcnp6 string
extern std::map<std::string,std::string> part_to_fluka;
/// map of name to fluka string
extern std::map<std::string,std::string> part_to_geant4;
/// map of name to geant4 string


/// \name is_hydrogen functions
/// \{
/// Returns whether or not the given particle is hydrogen or not, for
/// example, Protons (Hydrogen) are both valid nucids and fundamental
/// pdc types, all the following identify as hydrogen, Proton, Hydrogen,
/// Protium, "H1", "1H", 100001000, PDC(2212)
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
bool is_hydrogen(int n);
bool is_hydrogen(char *s);
bool is_hydrogen(std::string s);
/// \}

/// \name is_heavy_ion functions
/// \{
/// Returns whether or not the given particle is a heavy ion
/// or not. Heavy ions are not covered by the PDC scheme, therefore
/// the pyne::nucname class is used.
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
bool is_heavy_ion(int s);
bool is_heavy_ion(char *s);
bool is_heavy_ion(std::string s);
///\}

/// \name is_valid functions
/// \{
/// Returns whether or not the the given particle is a valid particle
/// in the PyNE particle class. All PDC numbers, names, altnames nucnames
/// are valid particle types
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
bool is_valid(int n);
bool is_valid(char *s);
bool is_valid(std::string s);
///\}

/// \name pdc_number functions
/// \{
/// Returns the PDC number of the particle given, if a valid pdc particle,
/// will return the number, for heavy ions will return 0.
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
int id(int n);
int id(char *s);
int id(std::string s);
///\}

/// \name name functions
/// \{
/// Returns the pyne::particle name of the particle given,
/// if a valid pdc particle number, name or nucname. Raises
/// exception if not a valid name
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
std::string name(int s);
std::string name(char *s);
std::string name(std::string s);
///\}

/// \name mcnp functions
/// \{
/// Returns the mcnp string of a valid pyne::particle name
/// \param s int, char*, String valid particle name, altname or nucid
std::string mcnp(int s);
std::string mcnp(char *s);
std::string mcnp(std::string s);
///\}

/// \name mcnp6 functions
/// \{
/// Returns the mcnp6 string of a valid pyne::particle name
/// \param s int, char*, String valid particle name, altname or nucid
std::string mcnp6(int s);
std::string mcnp6(char *s);
std::string mcnp6(std::string s);
///\}

/// \name fluka functions
/// \{
/// Returns the Fluka string of a valid pyne::particle name, or heavy ion
/// \param s int, char*, String valid particle name, altname or nucid
std::string fluka(int s);
std::string fluka(char *s);
std::string fluka(std::string s);
///\}

/// \name geant4 functions
/// \{
/// Returns the Geant4 string of a valid pyne::particle name, or heavy ion
/// \param s int, char*, String valid particle name, altname or nucid
std::string geant4(int s);
std::string geant4(char *s);
std::string geant4(std::string s);
///\}


/// \name describe functions
/// \{
/// Returns a long string that describes the particle, if given
/// a valid particle name, otherwise raises exception
/// \param n Integer PDC number or nucid
/// \param s String valid particle name, altname or nucid
std::string describe(int s);
std::string describe(char *s);
std::string describe(std::string s);

/// A helper function to set the contents of the variables in this library.
void * _fill_maps();
extern void * filler;  ///< A dummy variable used when calling #_fill_maps().


/// Custom excpeption for failed particle types
class NotAParticle : public std::exception
{
 public:
  /// Default constructor
  NotAParticle () {};

  /// Default destructor
  ~NotAParticle () throw () {};

  /// Constructor for raising the exception
  /// Spits out the particle name as input
  NotAParticle(std::string particle_name) {
    part_name = particle_name;
  }

  /// raises error message
  virtual const char* what() const throw() {
    std::string pname ("Not a valid particle name ");
    if(!part_name.empty())
      pname += part_name;
    return (const char *) pname.c_str();
  }

 private:
  std::string part_name;  /// the particle name


};
}
}

#endif
//
// end of src/particle.h
//


//
// start of src/data.h
//
/// \file data.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Implements basic nuclear data functions.

#ifndef PYNE_TEWK4A7VOFFLHDDXD5ZZ7KPXEQ
#define PYNE_TEWK4A7VOFFLHDDXD5ZZ7KPXEQ
#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <set>
#include <limits>
#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include "hdf5.h"
#include "hdf5_hl.h"

#ifndef PYNE_IS_AMALGAMATED
#include "h5wrap.h"
#include "extra_types.h"
#include "utils.h"
#include "nucname.h"
#include "rxname.h"
#endif

namespace pyne
{
/// \name Mathematical and Physical Constants
/// \{
extern const double pi;   ///< pi = 3.14159265359
extern const double N_A;  ///< Avogadro's Number
extern const double barns_per_cm2;  ///< barns per cm^2
extern const double cm2_per_barn;   ///< cm^2 per barn
extern const double sec_per_day;    ///< seconds per day
extern const double MeV_per_K;    ///< MeV per Kelvin
extern const double MeV_per_MJ;  ///< MeV per MJ
extern const double Bq_per_Ci;   ///< Becquerel per Curie
extern const double Ci_per_Bq;   ///< Curies per Becquerel
/// \}

extern std::string NUC_DATA_PATH; ///< Path to the nuc_data.h5 file.

/// Mapping from nodes in nuc_data.h5 to hashes of nodes
extern std::map<std::string, std::string> data_checksums;

/// \name Atomic Mass Data
/// \{

/// Mapping from nuclides in id form to their atomic masses.
extern std::map<int, double> atomic_mass_map;

/// a struct matching the atomic_mass table in nuc_data.h5.
typedef struct atomic_mass_data {
  int nuc;      ///< nuclide in id form
  double mass;  ///< nuclide atomic mass [amu]
  double error; ///< error in atomic mass [amu]
  double abund; ///< natural abundance of nuclide [atom fraction]
} atomic_mass_data;

// Loads preset dataset hashes into memory.
std::map<std::string, std::string> get_data_checksums();

/// Loads the atomic mass and natural abundance data from the nuc_data.h5 file
/// into memory.
void _load_atomic_mass_map();

/// \brief Returns the atomic mass of a nuclide \a nuc.
///
/// This function will first try to find the atomic mass data in the atomic_mass_map.
/// If this map is empty, it will load the data from disk.  If the nuclide is in an
/// excited state and not found in the map, it will give the value for the ground
/// state nuclide.  If the nuclide simply cannot be found, the A number is returned.
double atomic_mass(int nuc);
/// Returns the atomic mass of a nuclide \a nuc.
double atomic_mass(char * nuc);
/// Returns the atomic mass of a nuclide \a nuc.
double atomic_mass(std::string nuc);
/// \}


/// \name Natural Abundance Data
/// \{

/// Mapping from nuclides in id form to their natural abundances.
extern std::map<int, double> natural_abund_map;

/// \brief Returns the natural abundance of a nuclide \a nuc.
///
/// This follows the same the basic rules for finding or computing the natural
/// abundances as the atomic_mass() functions do.  However, if the nuclide cannot
/// be found, the default value returned is 0.0.
double natural_abund(int nuc);
/// Returns the natural abundance of a nuclide \a nuc.
double natural_abund(char * nuc);
/// Returns the natural abundance of a nuclide \a nuc.
double natural_abund(std::string nuc);
/// \}



/// \name Q_value Data
/// \{

/// Mapping from nuclides in id form to their q_values and
/// the fraction of Q that comes from gammas.
extern std::map<int, double> q_val_map;
extern std::map<int, double> gamma_frac_map;

/// a struct matching the q_value table in nuc_data.h5.
typedef struct q_val_data {
  int nuc;          ///< nuclide in id form
  double q_val;      ///< nuclide q_value [MeV/fission]
  double gamma_frac; ///< fraction of q that comes from gammas
} q_val_data;

/// Loads the q_value data from the nuc_data.h5 file into memory.
void _load_q_val_map();

/// \brief Returns the q_value of a nuclide \a nuc.
///
/// This function will first try to find the q_value data in the q_val_map.
/// If this map is empty, it will load the data from disk. If the nuclide simply
/// cannot be found, the default value returned is 0.0.
double q_val(int nuc);
double q_val(const char * nuc);
double q_val(std::string nuc);
double gamma_frac(int nuc);
double gamma_frac(const char * nuc);
double gamma_frac(std::string nuc);
/// \}


/// \name Dose Factor Data
/// \{

/// A struct matching the dose factor table in nuc_data.h5.
typedef struct dose {
  int nuc;              ///< nuclide in id form
  double ext_air_dose;  ///< nuclide ext_air dose factor [mrem/h per Ci/m^3]
  double ratio;         ///< ratio of external air dose factor to dose factor due to inhalation
  double ext_soil_dose; ///< nuclide ext_soil dose factor [mrem/h per Ci/m^2]
  double ingest_dose;   ///< nuclide dose factor due to ingestion [mrem/pCi]
  double fluid_frac;    ///< fraction of activity abosorbed in body fluids
  double inhale_dose;   ///< nuclide dose factor due to inhalation [mrem/pCi]
  char lung_mod;        ///< model of lung used (time of biological half life-- D, W, or Y)
} dose;

/// Mapping from int to dose for 3 sources
extern std::map<int, dose> epa_dose_map;
extern std::map<int, dose> doe_dose_map;
extern std::map<int, dose> genii_dose_map;

/// Loads the dose factor data from the nuc_data.h5 file into memory
/// according to the user-input source.
void _load_dose_map(std::map<int, dose>& dm, std::string source_path);

/// \brief Returns the dose factors of a nuclide.
///
/// These functions will first try to find the dose factor data in the dose_maps.
/// If the maps are empty, it will load the data from disk. If the nuclide simply
/// cannot be found, the default value returned is -1.
double ext_air_dose(int nuc, int source);
double ext_air_dose(const char * nuc, int source);
double ext_air_dose(std::string nuc, int source);
double ext_soil_dose(int nuc, int source);
double ext_soil_dose(const char * nuc, int source);
double ext_soil_dose(std::string nuc, int source);
double ingest_dose(int nuc, int source);
double ingest_dose(const char * nuc, int source);
double ingest_dose(std::string nuc, int source);
double inhale_dose(int nuc, int source);
double inhale_dose(const char * nuc, int source);
double inhale_dose(std::string nuc, int source);
double dose_ratio(int nuc, int source);
double dose_ratio(const char * nuc, int source);
double dose_ratio(std::string nuc, int source);
double dose_fluid_frac(int nuc, int source);
double dose_fluid_frac(const char * nuc, int source);
double dose_fluid_frac(std::string nuc, int source);
std::string dose_lung_model(int nuc, int source);
std::string dose_lung_model(const char * nuc, int source);
std::string dose_lung_model(std::string nuc, int source);
/// \}



/// \name Scattering Length Data
/// \{

/// Mapping from nuclides in id form to their coherent scattering length.
extern std::map<int, xd_complex_t> b_coherent_map;
/// Mapping from nuclides in id form to their incoherent scattering length.
extern std::map<int, xd_complex_t> b_incoherent_map;
/// Mapping from nuclides in id form to their scattering length.
extern std::map<int, double> b_map;

/// a struct matching the '/neutron/scattering_lengths' table in nuc_data.h5.
typedef struct scattering_lengths {
  int nuc;  ///< nuclide in id form
  xd_complex_t b_coherent;  ///< coherent scattering length [cm]
  xd_complex_t b_incoherent;  ///< incoherent scattering length [cm]
  double xs_coherent;   ///< coherent scattering cross section
  double xs_incoherent; ///< incoherent scattering cross section
  double xs;            ///< scattering cross section
} scattering_lengths;

/// Loads the scattering length data from the nuc_data.h5 file into memory.
void _load_scattering_lengths();

/// \brief Finds the coherent scattering length [cm] for a nuclide \a nuc.
///
/// This function works by first checking the b_coherent_map.  If the map is empty
/// then the data is read in from disk.  If no data is found than a value from a
/// nuclide with the same A number is returned instead.  If none of these exist,
/// then the value of a nuclide with the same Z number is used.  If none of these
/// work then 0.0 is returned.
xd_complex_t b_coherent(int nuc);
/// Finds the coherent scattering length [cm] for a nuclide \a nuc.
xd_complex_t b_coherent(char * nuc);
/// Finds the coherent scattering length [cm] for a nuclide \a nuc.
xd_complex_t b_coherent(std::string nuc);

/// \brief Finds the incoherent scattering length [cm] for a nuclide \a nuc.
///
/// This function works in the same way that b_coherent() does.
xd_complex_t b_incoherent(int nuc);
/// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
xd_complex_t b_incoherent(char * nuc);
/// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
xd_complex_t b_incoherent(std::string nuc);

/// Computes the scattering length [cm] from the coherent and incoherent components.
double b(int nuc);
/// Computes the scattering length [cm] from the coherent and incoherent components.
double b(char * nuc);
/// Computes the scattering length [cm] from the coherent and incoherent components.
double b(std::string nuc);
/// \}


/// \name Fission Product Yield Data
/// \{

/// Mapping from nuclides in id form to their scattering length.
extern std::map<std::pair<int, int>, double> wimsdfpy_data;

/// a struct matching the '/neutron/wimsd_fission_product' table in nuc_data.h5.
typedef struct wimsdfpy {
  int from_nuc;  ///< from nuclide in id form
  int to_nuc;  ///< from nuclide in id form
  double yields; ///< fission product yield, fraction [unitless]
} wimsdfpy;

/// Loads the WIMSD fission product yield data from the nuc_data.h5 file into memory.
void _load_wimsdfpy();

/// a struct matching the '/neutron/nds_fission_product' table in nuc_data.h5
typedef struct ndsfpy {
  int from_nuc; ///< id of fissioning nuclide
  int to_nuc; ///< id of fission product
  double yield_thermal; ///< thermal yield [fraction]
  double yield_thermal_err; ///< thermal yield error [fraction]
  double yield_fast; ///< fast yield [fraction]
  double yield_fast_err; ///< fast yield error [fraction]
  double yield_14MeV; ///< 14 MeV yield [fraction]
  double yield_14MeV_err; ///< 14 MeV yield error [fraction]
} ndsfpy;

/// a struct for the nds data for fpyield
typedef struct ndsfpysub {
  double yield_thermal; ///< thermal yield [fraction]
  double yield_thermal_err; ///< thermal yield error [fraction]
  double yield_fast; ///< fast yield [fraction]
  double yield_fast_err; ///< fast yield error [fraction]
  double yield_14MeV; ///< 14 MeV yield [fraction]
  double yield_14MeV_err; ///< 14 MeV yield error [fraction]
} ndsfpysub;


extern std::map<std::pair<int, int>, ndsfpysub> ndsfpy_data;

/// Loads the NDS fission product yield data from the nuc_data.h5 file into memory.
void _load_ndsfpy();

/// \brief Returns the fission product yield for a parent/child nuclide pair
///
/// This function works by first checking the fission yield data.  If this is
/// empty it loads the data from disk. If the parent/child nuclide pair
/// is still not found, then the process is assumed to be impossible
/// and 0.0 is returned. The data source is determined by the type value
/// as follows: 0 WIMS, 1 thermal NDS, 2 fast NDS, 3 14 MeV NDS.
/// negative type values return error for that data type.
double fpyield(std::pair<int, int> from_to, int source, bool get_error);
/// Returns the fission product yield for a parent/child nuclide pair
double fpyield(int from_nuc, int to_nuc, int source, bool get_error);
/// Returns the fission product yield for a parent/child nuclide pair
double fpyield(char * from_nuc, char * to_nuc, int source, bool get_error);
/// Returns the fission product yield for a parent/child nuclide pair
double fpyield(std::string from_nuc, std::string to_nuc, int source, bool get_error);

/// \}


/// \name Decay Data
/// \{

/// Data access functions

/// simple class to swap the order in which a pair is compared
class swapmapcompare
{
 public:
  /// This operator compares the second item in a pair first
  bool operator()(const std::pair<int, double>& lhs,
                  const std::pair<int, double>& rhs) const;
};

/// Access data in a std::map<std::pair<int, double> for a range of
/// values of the second member of the pair. Returns a vector of all
/// values at valoffset of class U of type T f
template<typename T, typename U> std::vector<T> data_access(double emin,
    double emax, size_t valoffset, std::map<std::pair<int, double>, U>  &data);
/// Access data in a std::map<std::pair<int, double> for a given
/// value of the first member of the pair. Returns a vector of all
/// values at valoffset of class U of type T
template<typename T, typename U> std::vector<T> data_access(int parent,
    double min, double max, size_t valoffset,
    std::map<std::pair<int, double>, U>  &data);
/// Access data in a std::map<std::pair<int, int> for a given
/// matching pair. Returns the value at valoffset of
/// class U of type T
template<typename T, typename U> T data_access(std::pair<int, int> from_to,
    size_t valoffset, std::map<std::pair<int, int>, U> &data);
/// Access data in a std::map<std::pair<int, int> for a given
/// value of the first member of the pair. Returns an array of the values
/// at valoffset of class U of type T
template<typename T, typename U> std::vector<T> data_access(int parent,
    size_t valoffset, std::map<std::pair<int, int>, U> &data);
template<typename T, typename U> std::vector<T> data_access(int parent,
    size_t valoffset, std::map<std::pair<int, unsigned int>, U> &data);

/// Access data in a std::map<int, data> format for a given first member
/// of the pair. Returns the value at valoffset of the matching datapoint.
template<typename U> double data_access(int parent,
                                        size_t valoffset, std::map<int, U> &data);

/// Structure for atomic data
typedef struct atomic {
  int z; ///< number of protons [int]
  double k_shell_fluor; ///< K-shell fluorescence [fraction]
  double k_shell_fluor_error; ///< K-shell fluorescence error [fraction]
  double l_shell_fluor; ///< L-shell fluorescence [fraction]
  double l_shell_fluor_error; ///< L-shell fluorescence error [fraction]
  double prob; ///< probability K shell hole is filled by L shell [fraction]
  double k_shell_be; ///< K-shell binding energy  [fraction]
  double k_shell_be_err; ///< K-shell binding energy error [fraction]
  double li_shell_be; ///< L-shell binding energy  [fraction]
  double li_shell_be_err; ///< L-shell binding energy error [fraction]
  double mi_shell_be; ///< M-shell binding energy  [fraction]
  double mi_shell_be_err; ///< M-shell binding energy error [fraction]
  double ni_shell_be; ///< N-shell binding energy  [fraction]
  double ni_shell_be_err; ///< N-shell binding energy error [fraction]
  double kb_to_ka; ///< ratio of Kb to Ka fluorescence [fraction]
  double kb_to_ka_err; ///< error in ratio of Kb to Ka fluorescence [fraction]
  double ka2_to_ka1; ///< Ka2 to Ka1 fluorescence ratio [fraction]
  double ka2_to_ka1_err; ///< Ka2 to Ka1 fluorescence error [fraction]
  double k_auger; ///< Auger electrons from k shell holes [fraction]
  double l_auger; ///< Auger electrons from l shell holes [fraction]
  double ka1_x_ray_en; ///< Ka1 X-ray energy [keV]
  double ka1_x_ray_en_err; ///< Ka1 X-ray energy error [keV]
  double ka2_x_ray_en; ///< Ka2 X-ray energy [keV]
  double ka2_x_ray_en_err; ///< Ka2 X-ray energy error [keV]
  double kb_x_ray_en; ///< Kb X-ray energy [keV]
  double l_x_ray_en; ///< L X-ray energy [keV]
} atomic;

// map of Z to atomic data
extern std::map<int, atomic> atomic_data_map;

template<typename T> void _load_data();
template<> void _load_data<atomic>();

// compute X-ray data
std::vector<std::pair<double, double> >
calculate_xray_data(int z, double k_conv, double l_conv);


/// a struct matching the '/decay/level_list' table in nuc_data.h5.
typedef struct level_data {
  int nuc_id; ///< state id of nuclide
  unsigned int rx_id; ///< rx id of reaction, 0 for basic level data
  double half_life; ///< half life [seconds]
  double level; ///< level energy [keV]
  double branch_ratio; ///< branch ratio [fraction]
  int metastable; ///< metastable level [int]
  char special; ///< special high-spin state [character]
} level_data;

/// Mapping from nuclides in id form to a struct containing data associated
/// with that level.
extern std::map<std::pair<int,double>, level_data> level_data_lvl_map;
extern std::map<std::pair<int,unsigned int>, level_data> level_data_rx_map;

template<> void _load_data<level_data>();

/// \brief Returns the nuc_id of an energy level
///
/// This function looks for the level that best matches the input level
/// within 1 keV of the input nuclide
int id_from_level(int nuc, double level);
int id_from_level(int nuc, double level, std::string special);
/// \brief Returns the nuc_id of a metastable state
///
/// This function looks through the level map for a given input nuc_id to find the
/// nuc_id corresponding to the level
int metastable_id(int nuc, int m);
/// Assumes the first metastable state is the desired one
int metastable_id(int nuc);

/// \brief Returns the half life for a nuclide \a nuc.
///
/// This function works by first checking the half_life_map.  If this is empty it
/// loads the data from disk.  If the nuclide is still not found, then the species
/// is assumed to be stable and infinity is returned.
double half_life(int nuc);
/// Returns the half life for a nuclide \a nuc.
double half_life(char * nuc);
/// Returns the half life for a nuclide \a nuc.
double half_life(std::string nuc);

/// \brief Returns the decay constant for a nuclide \a nuc.
///
/// This function works by first checking the decay_const_map.  If this is empty it
/// loads the data from disk.  If the nuclide is still not found, then the species
/// is assumed to be stable and 0.0 is returned.
double decay_const(int nuc);
/// Returns the decay constant for a nuclide \a nuc.
double decay_const(char * nuc);
/// Returns the decay constant for a nuclide \a nuc.
double decay_const(std::string nuc);

/// \brief Returns the branch ratio for a parent/child nuclide pair.
///
/// This function works by first checking the branch_ratio_map.  If this is empty it
/// loads the data from disk.  If the parent/child nuclide pair is still not found,
/// then the decay is assumed to be impossible and 0.0 is returned.
double branch_ratio(std::pair<int, int> from_to);
/// Returns the branch ratio for a parent/child nuclide pair.
double branch_ratio(int from_nuc, int to_nuc);
/// Returns the branch ratio for a parent/child nuclide pair.
double branch_ratio(char * from_nuc, char * to_nuc);
/// Returns the branch ratio for a parent/child nuclide pair.
double branch_ratio(std::string from_nuc, std::string to_nuc);

/// \brief Returns the excitation energy [MeV] of a \a nuc in a given state.
///
/// This function works by first checking the state_energy_map.  If this is empty it
/// loads the data from disk.  If the nuclide is still not found, then the species
/// is assumed to be in a ground state and 0.0 is returned.
double state_energy(int nuc);
/// Returns the excitation energy [MeV] of a \a nuc in a given state.
double state_energy(char * nuc);
/// Returns the excitation energy [MeV] of a \a nuc in a given state.
double state_energy(std::string nuc);

/// \brief Returns a set of decay children of a \a nuc.
///
/// This function works by first checking decay_chidlren_map.  If this is empty it
/// loads the data from disk.  If the nuclide is still not found, then the species
/// is assumed to be stable and an empty set is returned.
std::set<int> decay_children(int nuc);
/// Returns the decay constant for a nuclide \a nuc.
std::set<int> decay_children(char * nuc);
/// Returns the decay constant for a nuclide \a nuc.
std::set<int> decay_children(std::string nuc);

/// a struct matching the '/decay/decays' table in nuc_data.h5.
typedef struct decay {
  int parent; ///< state id of decay parent
  int child; ///< state id of decay child
  unsigned int decay; ///< rx id of decay
  double half_life; ///< half life of the decay [s]
  double half_life_error; ///< half life error of the decay [s]
  double branch_ratio; ///< branching ratio of this decay [fraction]
  double branch_ratio_error; ///< branching ratio of this decay [fraction]
  /// photon branching ratio of this decay [fraction]
  double photon_branch_ratio;
  /// photon branching ratio error of this decay [fraction]
  double photon_branch_ratio_error;
  /// beta branching ratio of this decay [fraction]
  double beta_branch_ratio;
  /// beta branching ratio error of this decay [fraction]
  double beta_branch_ratio_error;
} decay;

/// Loads the decay data from the nuc_data.h5 file into memory.
template<> void _load_data<decay>();
/// Mapping from a pair of nuclides in id form to a struct containing data
/// associated with the decay from the first to the second
extern std::map<std::pair<int, int>, decay> decay_data;

//
//
std::vector<int> decay_data_children(int parent);
std::pair<double, double> decay_half_life(std::pair<int,int>);
std::vector<std::pair<double, double> > decay_half_lifes(int);
std::pair<double, double> decay_branch_ratio(std::pair<int,int>);
std::vector<double> decay_branch_ratios(int parent);
std::pair<double, double> decay_photon_branch_ratio(std::pair<int,int>);
std::vector<std::pair<double, double> >decay_photon_branch_ratios(int parent);
std::pair<double, double> decay_beta_branch_ratio(std::pair<int,int>);
std::vector<std::pair<double, double> >decay_beta_branch_ratios(int parent);


/// a struct matching the '/decay/gammas' table in nuc_data.h5.
typedef struct gamma {
  int from_nuc; ///< state id of starting level
  int to_nuc; ///< state id of final level
  int parent_nuc; ///< state id of the primary decaying nucleus
  int child_nuc; ///< stateless id of the child nucleus
  double energy; ///< energy of the photon [keV]
  double energy_err; ///< energy error of the photon [keV]
  double photon_intensity; ///< photon intensity
  double photon_intensity_err; ///< photon intensity error
  double conv_intensity; ///< conversion intensity
  double conv_intensity_err; ///< conversion intensity error
  double total_intensity; ///< total decay intensity
  double total_intensity_err; ///< total decay intensity error
  double k_conv_e; ///< k conversion electron fraction
  double l_conv_e; ///< l conversion electron fraction
  double m_conv_e; ///< m conversion electron fraction
} gamma;

/// Loads the gamma ray data from the nuc_data.h5 file into memory.
template<> void _load_data<gamma>();

extern std::map<std::pair<int, double>, gamma> gamma_data;

//returns a list of gamma decay energies from input parent nuclide
std::vector<std::pair<double, double> > gamma_energy(int parent);
std::vector<std::pair<double, double> > gamma_energy(double energy,
    double error);
//returns a list of gamma photon intensities from input parent nuclide
std::vector<std::pair<double, double> > gamma_photon_intensity(int parent);
std::vector<std::pair<double, double> > gamma_photon_intensity(double energy,
    double error);
//returns a list of gamma conversion intensities from input parent nuclide
std::vector<std::pair<double, double> > gamma_conversion_intensity(int parent);
//returns a list of gamma total intensities from input parent nuclide
std::vector<std::pair<double, double> > gamma_total_intensity(int parent);
//returns a list of pairs of excited state transitions from an input parent nuclide
std::vector<std::pair<int, int> > gamma_from_to(int parent);
//returns a list of pairs of excited state transitions from an decay energy
std::vector<std::pair<int, int> > gamma_from_to(double energy, double error);
//returns a list of parent/child pairs associated with an input decay energy
std::vector<std::pair<int, int> > gamma_parent_child(double energy, double error);
//returns a list of parent nuclides associated with an input decay energy
std::vector<int> gamma_parent(double energy, double error);
// returns a list of child state_id's based on a gamma-ray energy
std::vector<int> gamma_child(double energy, double error);
// returns a list of child state_id's based on a parent state_id
std::vector<int> gamma_child(int parent);
//returns an array of arrays of X-ray energies and intesities for a
//given parent
std::vector<std::pair<double, double> > gamma_xrays(int parent);

/// Returns a list of energies and intensities normalized to branching ratios
std::vector<std::pair<double, double> > gammas(int parent_state_id);
std::vector<std::pair<double, double> > alphas(int parent_state_id);
std::vector<std::pair<double, double> > betas(int parent_state_id);
std::vector<std::pair<double, double> > xrays(int parent);

/// a struct matching the '/decay/alphas' table in nuc_data.h5.
typedef struct alpha {
  int from_nuc; ///< state id of parent nuclide
  int to_nuc; ///< state id of child nuclide
  double energy; ///< energy of alpha
  double intensity; ///< intensity of alpha decay
} alpha;

/// Loads the alpha decay data from the nuc_data.h5 file into memory.
template<> void _load_data<alpha>();

/// A vector of structs containing alpha data for access in memory
extern std::map<std::pair<int, double>, alpha> alpha_data;

//returns a list of alpha decay energies from input parent nuclide
std::vector<double > alpha_energy(int parent);
//returns a list of alpha decay intensities from input parent nuclide
std::vector<double> alpha_intensity(int parent);
//returns a list of alpha decay parents from input decay energy range
std::vector<int> alpha_parent(double energy, double error);
//returns a list of alpha decay children from input decay energy range
std::vector<int> alpha_child(double energy, double error);
//returns a list of alpha decay children from input parent nuclide
std::vector<int> alpha_child(int parent);

/// a struct matching the '/decay/betas' table in nuc_data.h5.
typedef struct beta {
  int from_nuc; ///< state id of parent nuclide
  int to_nuc; ///< state id of child nuclide
  double endpoint_energy; ///< beta decay endpoint energy
  double avg_energy; ///< beta decay average energy
  double intensity; ///< beta intensity
} beta;

/// Loads the beta decay data from the nuc_data.h5 file into memory.
template<> void _load_data<beta>();

/// A vector of structs containing beta data for access in memory
extern std::map<std::pair<int, double>, beta> beta_data;
//returns a list of beta decay endpoint energies from input parent nuclide
std::vector<double > beta_endpoint_energy(int parent);
//returns a list of beta decay average energies from input parent nuclide
std::vector<double > beta_average_energy(int parent);
//returns a list of beta decay intensities from input parent nuclide
std::vector<double> beta_intensity(int parent);
//returns a list of beta decay parents from input decay energy range
std::vector<int> beta_parent(double energy, double error);
//returns a list of beta decay children from input decay energy range
std::vector<int> beta_child(double energy, double error);
//returns a list of beta decay children from input parent nuclide
std::vector<int> beta_child(int parent);

/// A struct matching the '/decay/ecbp' table in nuc_data.h5.
typedef struct ecbp {
  int from_nuc;  ///< state id of parent nuclide
  int to_nuc; ///< state id of child nuclide
  double endpoint_energy; ///< beta decay endpoint energy
  double avg_energy; ///< beta decay average energy
  double beta_plus_intensity; ///< intensity of beta plus decay
  double ec_intensity; ///< intensity of electron capture
  double k_conv_e; ///< k conversion electron fraction
  double l_conv_e; ///< l conversion electron fraction
  double m_conv_e; ///< m conversion electron fraction
} ecbp;

/// A vector of structs containing ecbp data for access in memory
extern std::map<std::pair<int, double>, ecbp> ecbp_data;

/// Loads the electron capture and beta plus decay data from the
/// nuc_data.h5 file into memory.
template<> void _load_data<ecbp>();
///returns a list of electron capture/ beta plus decay endpoint energies from
///input parent nuclide
std::vector<double > ecbp_endpoint_energy(int parent);
//returns a list of electron capture/ beta plus decay average energies from
//input parent nuclide
std::vector<double > ecbp_average_energy(int parent);
//returns a list of electron capture decay intensities from input parent
//nuclide
std::vector<double> ec_intensity(int parent);
//returns a list of beta plus decay intensities from input parent nuclide
std::vector<double> bp_intensity(int parent);
//returns a list of electron capture /beta plus decay parents from input
//decay energy range
std::vector<int> ecbp_parent(double energy, double error);
//returns a list of electron capture /beta plus decay children from input
//decay energy range
std::vector<int> ecbp_child(double energy, double error);
//returns a list of electron capture /beta plus decay children from input
//parent nuclide
std::vector<int> ecbp_child(int parent);
//returns an array of arrays of X-ray energies and intesities for a
//given parent
std::vector<std::pair<double, double> > ecbp_xrays(int parent);
/// \}

/// map<energy, map<nuclide, map<rx, xs> > >
extern std::map<std::string, std::map<int, std::map<int, double> > >
simple_xs_map;

/// returns the microscopic cross section in barns for the specified
/// nuclide, reaction, and energy group.  energy must be one of: "thermal",
/// "thermal_maxwell_ave", "resonance_integral", "fourteen_MeV",
/// "fission_spectrum_ave".
double simple_xs(int nuc, int rx, std::string energy);
/// returns the microscopic cross section in barns for the specified
/// nuclide, reaction, and energy group.  energy must be one of: "thermal",
/// "thermal_maxwell_ave", "resonance_integral", "fourteen_MeV",
/// "fission_spectrum_ave".
double simple_xs(int nuc, std::string rx, std::string energy);
/// returns the microscopic cross section in barns for the specified
/// nuclide, reaction, and energy group.  energy must be one of: "thermal",
/// "thermal_maxwell_ave", "resonance_integral", "fourteen_MeV",
/// "fission_spectrum_ave".
double simple_xs(std::string nuc, int rx, std::string energy);
/// returns the microscopic cross section in barns for the specified
/// nuclide, reaction, and energy group.  energy must be one of: "thermal",
/// "thermal_maxwell_ave", "resonance_integral", "fourteen_MeV",
/// "fission_spectrum_ave".
double simple_xs(std::string nuc, std::string rx, std::string energy);

/// Custom exception for declaring a simple_xs request invalid
class InvalidSimpleXS : public std::exception
{
 public:
  InvalidSimpleXS () {};
  ~InvalidSimpleXS () throw () {};
  /// Exception thrown if energy group or rxname are invalid
  InvalidSimpleXS(std::string msg) : msg_(msg) {};
  /// Exception returns the string passed when thrown.
  virtual const char* what() const throw() {
    return msg_.c_str();
  };

 private:
  std::string msg_;
};

} // namespace pyne

#endif
//
// end of src/data.h
//


//
// start of src/json-forwards.h
//
/// Json-cpp amalgated forward header (http://jsoncpp.sourceforge.net/).
/// It is intented to be used with #include <json/json-forwards.h>
/// This header provides forward declaration for all JsonCpp types.

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////

/*
The JsonCpp library's source code, including accompanying documentation,
tests and demonstration applications, are licensed under the following
conditions...

The author (Baptiste Lepilleur) explicitly disclaims copyright in all
jurisdictions which recognize such a disclaimer. In such jurisdictions,
this software is released into the Public Domain.

In jurisdictions which do not recognize Public Domain property (e.g. Germany as of
2010), this software is Copyright (c) 2007-2010 by Baptiste Lepilleur, and is
released under the terms of the MIT License (see below).

In jurisdictions which recognize Public Domain property, the user of this
software may choose to accept it either as 1) Public Domain, 2) under the
conditions of the MIT License (see below), or 3) under the terms of dual
Public Domain/MIT License conditions described here, as they choose.

The MIT License is about as close to Public Domain as a license can get, and is
described in clear, concise terms at:

   http://en.wikipedia.org/wiki/MIT_License

The full text of the MIT License follows:

========================================================================
Copyright (c) 2007-2010 Baptiste Lepilleur

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
========================================================================
(END LICENSE TEXT)

The MIT license is compatible with both the GPL and commercial
software, affording one all of the rights of Public Domain with the
minor nuisance of being required to keep the above copyright notice
and license text in the source code. Note also that by accepting the
Public Domain "license" you can re-license your copy using whatever
license you like.

*/

// //////////////////////////////////////////////////////////////////////
// End of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////


#ifdef PYNE_IS_AMALGAMATED
#if !defined(JSON_IS_AMALGAMATION)
#define JSON_IS_AMALGAMATION
#endif
#endif


#ifndef JSON_FORWARD_AMALGATED_H_INCLUDED
# define JSON_FORWARD_AMALGATED_H_INCLUDED
/// If defined, indicates that the source file is amalgated
/// to prevent private header inclusion.
#define JSON_IS_AMALGATED

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/config.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_CONFIG_H_INCLUDED
# define JSON_CONFIG_H_INCLUDED

/// If defined, indicates that json library is embedded in CppTL library.
//# define JSON_IN_CPPTL 1

/// If defined, indicates that json may leverage CppTL library
//#  define JSON_USE_CPPTL 1
/// If defined, indicates that cpptl vector based map should be used instead of std::map
/// as Value container.
//#  define JSON_USE_CPPTL_SMALLMAP 1
/// If defined, indicates that Json specific container should be used
/// (hash table & simple deque container with customizable allocator).
/// THIS FEATURE IS STILL EXPERIMENTAL! There is know bugs: See #3177332
//#  define JSON_VALUE_USE_INTERNAL_MAP 1
/// Force usage of standard new/malloc based allocator instead of memory pool based allocator.
/// The memory pools allocator used optimization (initializing Value and ValueInternalLink
/// as if it was a POD) that may cause some validation tool to report errors.
/// Only has effects if JSON_VALUE_USE_INTERNAL_MAP is defined.
//#  define JSON_USE_SIMPLE_INTERNAL_ALLOCATOR 1

/// If defined, indicates that Json use exception to report invalid type manipulation
/// instead of C assert macro.
# define JSON_USE_EXCEPTION 1

/// If defined, indicates that the source file is amalgated
/// to prevent private header inclusion.
/// Remarks: it is automatically defined in the generated amalgated header.
// #define JSON_IS_AMALGAMATION


# ifdef JSON_IN_CPPTL
#  include <cpptl/config.h>
#  ifndef JSON_USE_CPPTL
#   define JSON_USE_CPPTL 1
#  endif
# endif

# ifdef JSON_IN_CPPTL
#  define JSON_API CPPTL_API
# elif defined(JSON_DLL_BUILD)
#  define JSON_API __declspec(dllexport)
# elif defined(JSON_DLL)
#  define JSON_API __declspec(dllimport)
# else
#  define JSON_API
# endif

// If JSON_NO_INT64 is defined, then Json only support C++ "int" type for integer
// Storages, and 64 bits integer support is disabled.
// #define JSON_NO_INT64 1

#if defined(_MSC_VER)  &&  _MSC_VER <= 1200 // MSVC 6
// Microsoft Visual Studio 6 only support conversion from __int64 to double
// (no conversion from unsigned __int64).
#define JSON_USE_INT64_DOUBLE_CONVERSION 1
#endif // if defined(_MSC_VER)  &&  _MSC_VER < 1200 // MSVC 6

#if defined(_MSC_VER)  &&  _MSC_VER >= 1500 // MSVC 2008
/// Indicates that the following function is deprecated.
# define JSONCPP_DEPRECATED(message) __declspec(deprecated(message))
#endif

#if !defined(JSONCPP_DEPRECATED)
# define JSONCPP_DEPRECATED(message)
#endif // if !defined(JSONCPP_DEPRECATED)

namespace Json
{
typedef int Int;
typedef unsigned int UInt;
# if defined(JSON_NO_INT64)
typedef int LargestInt;
typedef unsigned int LargestUInt;
#  undef JSON_HAS_INT64
# else // if defined(JSON_NO_INT64)
// For Microsoft Visual use specific types as long long is not supported
#  if defined(_MSC_VER) // Microsoft Visual Studio
typedef __int64 Int64;
typedef unsigned __int64 UInt64;
#  else // if defined(_MSC_VER) // Other platforms, use long long
typedef long long int Int64;
typedef unsigned long long int UInt64;
#  endif // if defined(_MSC_VER)
typedef Int64 LargestInt;
typedef UInt64 LargestUInt;
#  define JSON_HAS_INT64
# endif // if defined(JSON_NO_INT64)
} // end namespace Json


#endif // JSON_CONFIG_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/config.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/forwards.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_FORWARDS_H_INCLUDED
# define JSON_FORWARDS_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "config.h"
#endif // if !defined(JSON_IS_AMALGAMATION)

namespace Json
{

// writer.h
class FastWriter;
class StyledWriter;

// reader.h
class Reader;

// features.h
class Features;

// value.h
typedef unsigned int ArrayIndex;
class StaticString;
class Path;
class PathArgument;
class Value;
class ValueIteratorBase;
class ValueIterator;
class ValueConstIterator;
#ifdef JSON_VALUE_USE_INTERNAL_MAP
class ValueMapAllocator;
class ValueInternalLink;
class ValueInternalArray;
class ValueInternalMap;
#endif // #ifdef JSON_VALUE_USE_INTERNAL_MAP

} // namespace Json


#endif // JSON_FORWARDS_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/forwards.h
// //////////////////////////////////////////////////////////////////////





#endif //ifndef JSON_FORWARD_AMALGATED_H_INCLUDED
//
// end of src/json-forwards.h
//


//
// start of src/json.h
//
/// Json-cpp amalgated header (http://jsoncpp.sourceforge.net/).
/// It is intented to be used with #include <json.h>

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////

/*
The JsonCpp library's source code, including accompanying documentation,
tests and demonstration applications, are licensed under the following
conditions...

The author (Baptiste Lepilleur) explicitly disclaims copyright in all
jurisdictions which recognize such a disclaimer. In such jurisdictions,
this software is released into the Public Domain.

In jurisdictions which do not recognize Public Domain property (e.g. Germany as of
2010), this software is Copyright (c) 2007-2010 by Baptiste Lepilleur, and is
released under the terms of the MIT License (see below).

In jurisdictions which recognize Public Domain property, the user of this
software may choose to accept it either as 1) Public Domain, 2) under the
conditions of the MIT License (see below), or 3) under the terms of dual
Public Domain/MIT License conditions described here, as they choose.

The MIT License is about as close to Public Domain as a license can get, and is
described in clear, concise terms at:

   http://en.wikipedia.org/wiki/MIT_License

The full text of the MIT License follows:

========================================================================
Copyright (c) 2007-2010 Baptiste Lepilleur

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
========================================================================
(END LICENSE TEXT)

The MIT license is compatible with both the GPL and commercial
software, affording one all of the rights of Public Domain with the
minor nuisance of being required to keep the above copyright notice
and license text in the source code. Note also that by accepting the
Public Domain "license" you can re-license your copy using whatever
license you like.

*/

// //////////////////////////////////////////////////////////////////////
// End of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////

#ifdef PYNE_IS_AMALGAMATED
#if !defined(JSON_IS_AMALGAMATION)
#define JSON_IS_AMALGAMATION
#endif
#endif



#ifndef JSON_AMALGATED_H_INCLUDED
# define JSON_AMALGATED_H_INCLUDED
/// If defined, indicates that the source file is amalgated
/// to prevent private header inclusion.
#define JSON_IS_AMALGATED

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/config.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_CONFIG_H_INCLUDED
# define JSON_CONFIG_H_INCLUDED

/// If defined, indicates that json library is embedded in CppTL library.
//# define JSON_IN_CPPTL 1

/// If defined, indicates that json may leverage CppTL library
//#  define JSON_USE_CPPTL 1
/// If defined, indicates that cpptl vector based map should be used instead of std::map
/// as Value container.
//#  define JSON_USE_CPPTL_SMALLMAP 1
/// If defined, indicates that Json specific container should be used
/// (hash table & simple deque container with customizable allocator).
/// THIS FEATURE IS STILL EXPERIMENTAL! There is know bugs: See #3177332
//#  define JSON_VALUE_USE_INTERNAL_MAP 1
/// Force usage of standard new/malloc based allocator instead of memory pool based allocator.
/// The memory pools allocator used optimization (initializing Value and ValueInternalLink
/// as if it was a POD) that may cause some validation tool to report errors.
/// Only has effects if JSON_VALUE_USE_INTERNAL_MAP is defined.
//#  define JSON_USE_SIMPLE_INTERNAL_ALLOCATOR 1

/// If defined, indicates that Json use exception to report invalid type manipulation
/// instead of C assert macro.
# define JSON_USE_EXCEPTION 1

/// If defined, indicates that the source file is amalgated
/// to prevent private header inclusion.
/// Remarks: it is automatically defined in the generated amalgated header.
// #define JSON_IS_AMALGAMATION


# ifdef JSON_IN_CPPTL
#  include <cpptl/config.h>
#  ifndef JSON_USE_CPPTL
#   define JSON_USE_CPPTL 1
#  endif
# endif

# ifdef JSON_IN_CPPTL
#  define JSON_API CPPTL_API
# elif defined(JSON_DLL_BUILD)
#  define JSON_API __declspec(dllexport)
# elif defined(JSON_DLL)
#  define JSON_API __declspec(dllimport)
# else
#  define JSON_API
# endif

// If JSON_NO_INT64 is defined, then Json only support C++ "int" type for integer
// Storages, and 64 bits integer support is disabled.
// #define JSON_NO_INT64 1

#if defined(_MSC_VER)  &&  _MSC_VER <= 1200 // MSVC 6
// Microsoft Visual Studio 6 only support conversion from __int64 to double
// (no conversion from unsigned __int64).
#define JSON_USE_INT64_DOUBLE_CONVERSION 1
#endif // if defined(_MSC_VER)  &&  _MSC_VER < 1200 // MSVC 6

#if defined(_MSC_VER)  &&  _MSC_VER >= 1500 // MSVC 2008
/// Indicates that the following function is deprecated.
# define JSONCPP_DEPRECATED(message) __declspec(deprecated(message))
#endif

#if !defined(JSONCPP_DEPRECATED)
# define JSONCPP_DEPRECATED(message)
#endif // if !defined(JSONCPP_DEPRECATED)

namespace Json
{
typedef int Int;
typedef unsigned int UInt;
# if defined(JSON_NO_INT64)
typedef int LargestInt;
typedef unsigned int LargestUInt;
#  undef JSON_HAS_INT64
# else // if defined(JSON_NO_INT64)
// For Microsoft Visual use specific types as long long is not supported
#  if defined(_MSC_VER) // Microsoft Visual Studio
typedef __int64 Int64;
typedef unsigned __int64 UInt64;
#  else // if defined(_MSC_VER) // Other platforms, use long long
typedef long long int Int64;
typedef unsigned long long int UInt64;
#  endif // if defined(_MSC_VER)
typedef Int64 LargestInt;
typedef UInt64 LargestUInt;
#  define JSON_HAS_INT64
# endif // if defined(JSON_NO_INT64)
} // end namespace Json


#endif // JSON_CONFIG_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/config.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/forwards.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_FORWARDS_H_INCLUDED
# define JSON_FORWARDS_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "config.h"
#endif // if !defined(JSON_IS_AMALGAMATION)

namespace Json
{

// writer.h
class FastWriter;
class StyledWriter;

// reader.h
class Reader;

// features.h
class Features;

// value.h
typedef unsigned int ArrayIndex;
class StaticString;
class Path;
class PathArgument;
class Value;
class ValueIteratorBase;
class ValueIterator;
class ValueConstIterator;
#ifdef JSON_VALUE_USE_INTERNAL_MAP
class ValueMapAllocator;
class ValueInternalLink;
class ValueInternalArray;
class ValueInternalMap;
#endif // #ifdef JSON_VALUE_USE_INTERNAL_MAP

} // namespace Json


#endif // JSON_FORWARDS_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/forwards.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/features.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef CPPTL_JSON_FEATURES_H_INCLUDED
# define CPPTL_JSON_FEATURES_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "forwards.h"
#endif // if !defined(JSON_IS_AMALGAMATION)

namespace Json
{

/** \brief Configuration passed to reader and writer.
 * This configuration object can be used to force the Reader or Writer
 * to behave in a standard conforming way.
 */
class JSON_API Features
{
 public:
  /** \brief A configuration that allows all features and assumes all strings are UTF-8.
   * - C & C++ comments are allowed
   * - Root object can be any JSON value
   * - Assumes Value strings are encoded in UTF-8
   */
  static Features all();

  /** \brief A configuration that is strictly compatible with the JSON specification.
   * - Comments are forbidden.
   * - Root object must be either an array or an object value.
   * - Assumes Value strings are encoded in UTF-8
   */
  static Features strictMode();

  /** \brief Initialize the configuration like JsonConfig::allFeatures;
   */
  Features();

  /// \c true if comments are allowed. Default: \c true.
  bool allowComments_;

  /// \c true if root must be either an array or an object value. Default: \c false.
  bool strictRoot_;
};

} // namespace Json

#endif // CPPTL_JSON_FEATURES_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/features.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/value.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef CPPTL_JSON_H_INCLUDED
# define CPPTL_JSON_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "forwards.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
# include <string>
# include <vector>

# ifndef JSON_USE_CPPTL_SMALLMAP
#  include <map>
# else
#  include <cpptl/smallmap.h>
# endif
# ifdef JSON_USE_CPPTL
#  include <cpptl/forwards.h>
# endif

/** \brief JSON (JavaScript Object Notation).
 */
namespace Json
{

/** \brief Type of the value held by a Value object.
 */
enum ValueType {
    nullValue = 0, ///< 'null' value
    intValue,      ///< signed integer value
    uintValue,     ///< unsigned integer value
    realValue,     ///< double value
    stringValue,   ///< UTF-8 string value
    booleanValue,  ///< bool value
    arrayValue,    ///< array value (ordered list)
    objectValue    ///< object value (collection of name/value pairs).
};

enum CommentPlacement {
    commentBefore = 0,        ///< a comment placed on the line before a value
    commentAfterOnSameLine,   ///< a comment just after a value on the same line
    commentAfter,             ///< a comment on the line after a value (only make sense for root value)
    numberOfCommentPlacement
};

//# ifdef JSON_USE_CPPTL
//   typedef CppTL::AnyEnumerator<const char *> EnumMemberNames;
//   typedef CppTL::AnyEnumerator<const Value &> EnumValues;
//# endif

/** \brief Lightweight wrapper to tag static string.
 *
 * Value constructor and objectValue member assignement takes advantage of the
 * StaticString and avoid the cost of string duplication when storing the
 * string or the member name.
 *
 * Example of usage:
 * \code
 * Json::Value aValue( StaticString("some text") );
 * Json::Value object;
 * static const StaticString code("code");
 * object[code] = 1234;
 * \endcode
 */
class JSON_API StaticString
{
 public:
  explicit StaticString( const char *czstring )
    : str_( czstring ) {
  }

  operator const char *() const {
    return str_;
  }

  const char *c_str() const {
    return str_;
  }

 private:
  const char *str_;
};

/** \brief Represents a <a HREF="http://www.json.org">JSON</a> value.
 *
 * This class is a discriminated union wrapper that can represents a:
 * - signed integer [range: Value::minInt - Value::maxInt]
 * - unsigned integer (range: 0 - Value::maxUInt)
 * - double
 * - UTF-8 string
 * - boolean
 * - 'null'
 * - an ordered list of Value
 * - collection of name/value pairs (javascript object)
 *
 * The type of the held value is represented by a #ValueType and
 * can be obtained using type().
 *
 * values of an #objectValue or #arrayValue can be accessed using operator[]() methods.
 * Non const methods will automatically create the a #nullValue element
 * if it does not exist.
 * The sequence of an #arrayValue will be automatically resize and initialized
 * with #nullValue. resize() can be used to enlarge or truncate an #arrayValue.
 *
 * The get() methods can be used to obtanis default value in the case the required element
 * does not exist.
 *
 * It is possible to iterate over the list of a #objectValue values using
 * the getMemberNames() method.
 */
class JSON_API Value
{
  friend class ValueIteratorBase;
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  friend class ValueInternalLink;
  friend class ValueInternalMap;
# endif
 public:
  typedef std::vector<std::string> Members;
  typedef ValueIterator iterator;
  typedef ValueConstIterator const_iterator;
  typedef Json::UInt UInt;
  typedef Json::Int Int;
# if defined(JSON_HAS_INT64)
  typedef Json::UInt64 UInt64;
  typedef Json::Int64 Int64;
#endif // defined(JSON_HAS_INT64)
  typedef Json::LargestInt LargestInt;
  typedef Json::LargestUInt LargestUInt;
  typedef Json::ArrayIndex ArrayIndex;

  static const Value null;
  /// Minimum signed integer value that can be stored in a Json::Value.
  static const LargestInt minLargestInt;
  /// Maximum signed integer value that can be stored in a Json::Value.
  static const LargestInt maxLargestInt;
  /// Maximum unsigned integer value that can be stored in a Json::Value.
  static const LargestUInt maxLargestUInt;

  /// Minimum signed int value that can be stored in a Json::Value.
  static const Int minInt;
  /// Maximum signed int value that can be stored in a Json::Value.
  static const Int maxInt;
  /// Maximum unsigned int value that can be stored in a Json::Value.
  static const UInt maxUInt;

  /// Minimum signed 64 bits int value that can be stored in a Json::Value.
  static const Int64 minInt64;
  /// Maximum signed 64 bits int value that can be stored in a Json::Value.
  static const Int64 maxInt64;
  /// Maximum unsigned 64 bits int value that can be stored in a Json::Value.
  static const UInt64 maxUInt64;

 private:
#ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION
# ifndef JSON_VALUE_USE_INTERNAL_MAP
  class CZString
{
   public:
    enum DuplicationPolicy {
        noDuplication = 0,
        duplicate,
        duplicateOnCopy
    };
    CZString( ArrayIndex index );
    CZString( const char *cstr, DuplicationPolicy allocate );
    CZString( const CZString &other );
    ~CZString();
    CZString &operator =( const CZString &other );
    bool operator<( const CZString &other ) const;
    bool operator==( const CZString &other ) const;
    ArrayIndex index() const;
    const char *c_str() const;
    bool isStaticString() const;
   private:
    void swap( CZString &other );
    const char *cstr_;
    ArrayIndex index_;
  };

 public:
#  ifndef JSON_USE_CPPTL_SMALLMAP
  typedef std::map<CZString, Value> ObjectValues;
#  else
  typedef CppTL::SmallMap<CZString, Value> ObjectValues;
#  endif // ifndef JSON_USE_CPPTL_SMALLMAP
# endif // ifndef JSON_VALUE_USE_INTERNAL_MAP
#endif // ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION

 public:
  /** \brief Create a default Value of the given type.

    This is a very useful constructor.
    To create an empty array, pass arrayValue.
    To create an empty object, pass objectValue.
    Another Value can then be set to this one by assignment.
  This is useful since clear() and resize() will not alter types.

    Examples:
  \code
  Json::Value null_value; // null
  Json::Value arr_value(Json::arrayValue); // []
  Json::Value obj_value(Json::objectValue); // {}
  \endcode
  */
  Value( ValueType type = nullValue );
  Value( Int value );
  Value( UInt value );
#if defined(JSON_HAS_INT64)
  Value( Int64 value );
  Value( UInt64 value );
#endif // if defined(JSON_HAS_INT64)
  Value( double value );
  Value( const char *value );
  Value( const char *beginValue, const char *endValue );
  /** \brief Constructs a value from a static string.

   * Like other value string constructor but do not duplicate the string for
   * internal storage. The given string must remain alive after the call to this
   * constructor.
   * Example of usage:
   * \code
   * Json::Value aValue( StaticString("some text") );
   * \endcode
   */
  Value( const StaticString &value );
  Value( const std::string &value );
# ifdef JSON_USE_CPPTL
  Value( const CppTL::ConstString &value );
# endif
  Value( bool value );
  Value( const Value &other );
  ~Value();

  Value &operator=( const Value &other );
  /// Swap values.
  /// \note Currently, comments are intentionally not swapped, for
  /// both logic and efficiency.
  void swap( Value &other );

  ValueType type() const;

  bool operator <( const Value &other ) const;
  bool operator <=( const Value &other ) const;
  bool operator >=( const Value &other ) const;
  bool operator >( const Value &other ) const;

  bool operator ==( const Value &other ) const;
  bool operator !=( const Value &other ) const;

  int compare( const Value &other ) const;

  const char *asCString() const;
  std::string asString() const;
# ifdef JSON_USE_CPPTL
  CppTL::ConstString asConstString() const;
# endif
  Int asInt() const;
  UInt asUInt() const;
  Int64 asInt64() const;
  UInt64 asUInt64() const;
  LargestInt asLargestInt() const;
  LargestUInt asLargestUInt() const;
  float asFloat() const;
  double asDouble() const;
  bool asBool() const;

  bool isNull() const;
  bool isBool() const;
  bool isInt() const;
  bool isUInt() const;
  bool isIntegral() const;
  bool isDouble() const;
  bool isNumeric() const;
  bool isString() const;
  bool isArray() const;
  bool isObject() const;

  bool isConvertibleTo( ValueType other ) const;

  /// Number of values in array or object
  ArrayIndex size() const;

  /// \brief Return true if empty array, empty object, or null;
  /// otherwise, false.
  bool empty() const;

  /// Return isNull()
  bool operator!() const;

  /// Remove all object members and array elements.
  /// \pre type() is arrayValue, objectValue, or nullValue
  /// \post type() is unchanged
  void clear();

  /// Resize the array to size elements.
  /// New elements are initialized to null.
  /// May only be called on nullValue or arrayValue.
  /// \pre type() is arrayValue or nullValue
  /// \post type() is arrayValue
  void resize( ArrayIndex size );

  /// Access an array element (zero based index ).
  /// If the array contains less than index element, then null value are inserted
  /// in the array so that its size is index+1.
  /// (You may need to say 'value[0u]' to get your compiler to distinguish
  ///  this from the operator[] which takes a string.)
  Value &operator[]( ArrayIndex index );

  /// Access an array element (zero based index ).
  /// If the array contains less than index element, then null value are inserted
  /// in the array so that its size is index+1.
  /// (You may need to say 'value[0u]' to get your compiler to distinguish
  ///  this from the operator[] which takes a string.)
  Value &operator[]( int index );

  /// Access an array element (zero based index )
  /// (You may need to say 'value[0u]' to get your compiler to distinguish
  ///  this from the operator[] which takes a string.)
  const Value &operator[]( ArrayIndex index ) const;

  /// Access an array element (zero based index )
  /// (You may need to say 'value[0u]' to get your compiler to distinguish
  ///  this from the operator[] which takes a string.)
  const Value &operator[]( int index ) const;

  /// If the array contains at least index+1 elements, returns the element value,
  /// otherwise returns defaultValue.
  Value get( ArrayIndex index,
             const Value &defaultValue ) const;
  /// Return true if index < size().
  bool isValidIndex( ArrayIndex index ) const;
  /// \brief Append value to array at the end.
  ///
  /// Equivalent to jsonvalue[jsonvalue.size()] = value;
  Value &append( const Value &value );

  /// Access an object value by name, create a null member if it does not exist.
  Value &operator[]( const char *key );
  /// Access an object value by name, returns null if there is no member with that name.
  const Value &operator[]( const char *key ) const;
  /// Access an object value by name, create a null member if it does not exist.
  Value &operator[]( const std::string &key );
  /// Access an object value by name, returns null if there is no member with that name.
  const Value &operator[]( const std::string &key ) const;
  /** \brief Access an object value by name, create a null member if it does not exist.

   * If the object as no entry for that name, then the member name used to store
   * the new entry is not duplicated.
   * Example of use:
   * \code
   * Json::Value object;
   * static const StaticString code("code");
   * object[code] = 1234;
   * \endcode
   */
  Value &operator[]( const StaticString &key );
# ifdef JSON_USE_CPPTL
  /// Access an object value by name, create a null member if it does not exist.
  Value &operator[]( const CppTL::ConstString &key );
  /// Access an object value by name, returns null if there is no member with that name.
  const Value &operator[]( const CppTL::ConstString &key ) const;
# endif
  /// Return the member named key if it exist, defaultValue otherwise.
  Value get( const char *key,
             const Value &defaultValue ) const;
  /// Return the member named key if it exist, defaultValue otherwise.
  Value get( const std::string &key,
             const Value &defaultValue ) const;
# ifdef JSON_USE_CPPTL
  /// Return the member named key if it exist, defaultValue otherwise.
  Value get( const CppTL::ConstString &key,
             const Value &defaultValue ) const;
# endif
  /// \brief Remove and return the named member.
  ///
  /// Do nothing if it did not exist.
  /// \return the removed Value, or null.
  /// \pre type() is objectValue or nullValue
  /// \post type() is unchanged
  Value removeMember( const char* key );
  /// Same as removeMember(const char*)
  Value removeMember( const std::string &key );

  /// Return true if the object has a member named key.
  bool isMember( const char *key ) const;
  /// Return true if the object has a member named key.
  bool isMember( const std::string &key ) const;
# ifdef JSON_USE_CPPTL
  /// Return true if the object has a member named key.
  bool isMember( const CppTL::ConstString &key ) const;
# endif

  /// \brief Return a list of the member names.
  ///
  /// If null, return an empty list.
  /// \pre type() is objectValue or nullValue
  /// \post if type() was nullValue, it remains nullValue
  Members getMemberNames() const;

//# ifdef JSON_USE_CPPTL
//      EnumMemberNames enumMemberNames() const;
//      EnumValues enumValues() const;
//# endif

  /// Comments must be //... or /* ... */
  void setComment( const char *comment,
                   CommentPlacement placement );
  /// Comments must be //... or /* ... */
  void setComment( const std::string &comment,
                   CommentPlacement placement );
  bool hasComment( CommentPlacement placement ) const;
  /// Include delimiters and embedded newlines.
  std::string getComment( CommentPlacement placement ) const;

  std::string toStyledString() const;

  const_iterator begin() const;
  const_iterator end() const;

  iterator begin();
  iterator end();

 private:
  Value &resolveReference( const char *key,
                           bool isStatic );

# ifdef JSON_VALUE_USE_INTERNAL_MAP
  inline bool isItemAvailable() const {
    return itemIsUsed_ == 0;
  }

  inline void setItemUsed( bool isUsed = true ) {
    itemIsUsed_ = isUsed ? 1 : 0;
  }

  inline bool isMemberNameStatic() const {
    return memberNameIsStatic_ == 0;
  }

  inline void setMemberNameIsStatic( bool isStatic ) {
    memberNameIsStatic_ = isStatic ? 1 : 0;
  }
# endif // # ifdef JSON_VALUE_USE_INTERNAL_MAP

 private:
  struct CommentInfo {
    CommentInfo();
    ~CommentInfo();

    void setComment( const char *text );

    char *comment_;
  };

  //struct MemberNamesTransform
  //{
  //   typedef const char *result_type;
  //   const char *operator()( const CZString &name ) const
  //   {
  //      return name.c_str();
  //   }
  //};

  union ValueHolder {
    LargestInt int_;
    LargestUInt uint_;
    double real_;
    bool bool_;
    char *string_;
# ifdef JSON_VALUE_USE_INTERNAL_MAP
    ValueInternalArray *array_;
    ValueInternalMap *map_;
#else
    ObjectValues *map_;
# endif
  } value_;
  ValueType type_ : 8;
  int allocated_ : 1;     // Notes: if declared as bool, bitfield is useless.
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  unsigned int itemIsUsed_ : 1;      // used by the ValueInternalMap container.
  int memberNameIsStatic_ : 1;       // used by the ValueInternalMap container.
# endif
  CommentInfo *comments_;
};


/** \brief Experimental and untested: represents an element of the "path" to access a node.
 */
class PathArgument
{
 public:
  friend class Path;

  PathArgument();
  PathArgument( ArrayIndex index );
  PathArgument( const char *key );
  PathArgument( const std::string &key );

 private:
  enum Kind {
      kindNone = 0,
      kindIndex,
      kindKey
  };
  std::string key_;
  ArrayIndex index_;
  Kind kind_;
};

/** \brief Experimental and untested: represents a "path" to access a node.
 *
 * Syntax:
 * - "." => root node
 * - ".[n]" => elements at index 'n' of root node (an array value)
 * - ".name" => member named 'name' of root node (an object value)
 * - ".name1.name2.name3"
 * - ".[0][1][2].name1[3]"
 * - ".%" => member name is provided as parameter
 * - ".[%]" => index is provied as parameter
 */
class Path
{
 public:
  Path( const std::string &path,
        const PathArgument &a1 = PathArgument(),
        const PathArgument &a2 = PathArgument(),
        const PathArgument &a3 = PathArgument(),
        const PathArgument &a4 = PathArgument(),
        const PathArgument &a5 = PathArgument() );

  const Value &resolve( const Value &root ) const;
  Value resolve( const Value &root,
                 const Value &defaultValue ) const;
  /// Creates the "path" to access the specified node and returns a reference on the node.
  Value &make( Value &root ) const;

 private:
  typedef std::vector<const PathArgument *> InArgs;
  typedef std::vector<PathArgument> Args;

  void makePath( const std::string &path,
                 const InArgs &in );
  void addPathInArg( const std::string &path,
                     const InArgs &in,
                     InArgs::const_iterator &itInArg,
                     PathArgument::Kind kind );
  void invalidPath( const std::string &path,
                    int location );

  Args args_;
};



#ifdef JSON_VALUE_USE_INTERNAL_MAP
/** \brief Allocator to customize Value internal map.
 * Below is an example of a simple implementation (default implementation actually
 * use memory pool for speed).
 * \code
   class DefaultValueMapAllocator : public ValueMapAllocator
   {
   public: // overridden from ValueMapAllocator
      virtual ValueInternalMap *newMap()
      {
         return new ValueInternalMap();
      }

      virtual ValueInternalMap *newMapCopy( const ValueInternalMap &other )
      {
         return new ValueInternalMap( other );
      }

      virtual void destructMap( ValueInternalMap *map )
      {
         delete map;
      }

      virtual ValueInternalLink *allocateMapBuckets( unsigned int size )
      {
         return new ValueInternalLink[size];
      }

      virtual void releaseMapBuckets( ValueInternalLink *links )
      {
         delete [] links;
      }

      virtual ValueInternalLink *allocateMapLink()
      {
         return new ValueInternalLink();
      }

      virtual void releaseMapLink( ValueInternalLink *link )
      {
         delete link;
      }
   };
 * \endcode
 */
class JSON_API ValueMapAllocator
{
 public:
  virtual ~ValueMapAllocator();
  virtual ValueInternalMap *newMap() = 0;
  virtual ValueInternalMap *newMapCopy( const ValueInternalMap &other ) = 0;
  virtual void destructMap( ValueInternalMap *map ) = 0;
  virtual ValueInternalLink *allocateMapBuckets( unsigned int size ) = 0;
  virtual void releaseMapBuckets( ValueInternalLink *links ) = 0;
  virtual ValueInternalLink *allocateMapLink() = 0;
  virtual void releaseMapLink( ValueInternalLink *link ) = 0;
};

/** \brief ValueInternalMap hash-map bucket chain link (for internal use only).
 * \internal previous_ & next_ allows for bidirectional traversal.
 */
class JSON_API ValueInternalLink
{
 public:
  enum { itemPerLink = 6 };  // sizeof(ValueInternalLink) = 128 on 32 bits architecture.
  enum InternalFlags {
      flagAvailable = 0,
      flagUsed = 1
  };

  ValueInternalLink();

  ~ValueInternalLink();

  Value items_[itemPerLink];
  char *keys_[itemPerLink];
  ValueInternalLink *previous_;
  ValueInternalLink *next_;
};


/** \brief A linked page based hash-table implementation used internally by Value.
 * \internal ValueInternalMap is a tradional bucket based hash-table, with a linked
 * list in each bucket to handle collision. There is an addional twist in that
 * each node of the collision linked list is a page containing a fixed amount of
 * value. This provides a better compromise between memory usage and speed.
 *
 * Each bucket is made up of a chained list of ValueInternalLink. The last
 * link of a given bucket can be found in the 'previous_' field of the following bucket.
 * The last link of the last bucket is stored in tailLink_ as it has no following bucket.
 * Only the last link of a bucket may contains 'available' item. The last link always
 * contains at least one element unless is it the bucket one very first link.
 */
class JSON_API ValueInternalMap
{
  friend class ValueIteratorBase;
  friend class Value;
 public:
  typedef unsigned int HashKey;
  typedef unsigned int BucketIndex;

# ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION
  struct IteratorState {
    IteratorState()
      : map_(0)
    , link_(0)
    , itemIndex_(0)
    , bucketIndex_(0) {
    }
    ValueInternalMap *map_;
    ValueInternalLink *link_;
    BucketIndex itemIndex_;
    BucketIndex bucketIndex_;
  };
# endif // ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION

  ValueInternalMap();
  ValueInternalMap( const ValueInternalMap &other );
  ValueInternalMap &operator =( const ValueInternalMap &other );
  ~ValueInternalMap();

  void swap( ValueInternalMap &other );

  BucketIndex size() const;

  void clear();

  bool reserveDelta( BucketIndex growth );

  bool reserve( BucketIndex newItemCount );

  const Value *find( const char *key ) const;

  Value *find( const char *key );

  Value &resolveReference( const char *key,
                           bool isStatic );

  void remove( const char *key );

  void doActualRemove( ValueInternalLink *link,
                       BucketIndex index,
                       BucketIndex bucketIndex );

  ValueInternalLink *&getLastLinkInBucket( BucketIndex bucketIndex );

  Value &setNewItem( const char *key,
                     bool isStatic,
                     ValueInternalLink *link,
                     BucketIndex index );

  Value &unsafeAdd( const char *key,
                    bool isStatic,
                    HashKey hashedKey );

  HashKey hash( const char *key ) const;

  int compare( const ValueInternalMap &other ) const;

 private:
  void makeBeginIterator( IteratorState &it ) const;
  void makeEndIterator( IteratorState &it ) const;
  static bool equals( const IteratorState &x, const IteratorState &other );
  static void increment( IteratorState &iterator );
  static void incrementBucket( IteratorState &iterator );
  static void decrement( IteratorState &iterator );
  static const char *key( const IteratorState &iterator );
  static const char *key( const IteratorState &iterator, bool &isStatic );
  static Value &value( const IteratorState &iterator );
  static int distance( const IteratorState &x, const IteratorState &y );

 private:
  ValueInternalLink *buckets_;
  ValueInternalLink *tailLink_;
  BucketIndex bucketsSize_;
  BucketIndex itemCount_;
};

/** \brief A simplified deque implementation used internally by Value.
* \internal
* It is based on a list of fixed "page", each page contains a fixed number of items.
* Instead of using a linked-list, a array of pointer is used for fast item look-up.
* Look-up for an element is as follow:
* - compute page index: pageIndex = itemIndex / itemsPerPage
* - look-up item in page: pages_[pageIndex][itemIndex % itemsPerPage]
*
* Insertion is amortized constant time (only the array containing the index of pointers
* need to be reallocated when items are appended).
*/
class JSON_API ValueInternalArray
{
  friend class Value;
  friend class ValueIteratorBase;
 public:
  enum { itemsPerPage = 8 };    // should be a power of 2 for fast divide and modulo.
  typedef Value::ArrayIndex ArrayIndex;
  typedef unsigned int PageIndex;

# ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION
  struct IteratorState { // Must be a POD
    IteratorState()
      : array_(0)
    , currentPageIndex_(0)
    , currentItemIndex_(0) {
    }
    ValueInternalArray *array_;
    Value **currentPageIndex_;
    unsigned int currentItemIndex_;
  };
# endif // ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION

  ValueInternalArray();
  ValueInternalArray( const ValueInternalArray &other );
  ValueInternalArray &operator =( const ValueInternalArray &other );
  ~ValueInternalArray();
  void swap( ValueInternalArray &other );

  void clear();
  void resize( ArrayIndex newSize );

  Value &resolveReference( ArrayIndex index );

  Value *find( ArrayIndex index ) const;

  ArrayIndex size() const;

  int compare( const ValueInternalArray &other ) const;

 private:
  static bool equals( const IteratorState &x, const IteratorState &other );
  static void increment( IteratorState &iterator );
  static void decrement( IteratorState &iterator );
  static Value &dereference( const IteratorState &iterator );
  static Value &unsafeDereference( const IteratorState &iterator );
  static int distance( const IteratorState &x, const IteratorState &y );
  static ArrayIndex indexOf( const IteratorState &iterator );
  void makeBeginIterator( IteratorState &it ) const;
  void makeEndIterator( IteratorState &it ) const;
  void makeIterator( IteratorState &it, ArrayIndex index ) const;

  void makeIndexValid( ArrayIndex index );

  Value **pages_;
  ArrayIndex size_;
  PageIndex pageCount_;
};

/** \brief Experimental: do not use. Allocator to customize Value internal array.
 * Below is an example of a simple implementation (actual implementation use
 * memory pool).
   \code
class DefaultValueArrayAllocator : public ValueArrayAllocator
{
public: // overridden from ValueArrayAllocator
virtual ~DefaultValueArrayAllocator()
{
}

virtual ValueInternalArray *newArray()
{
   return new ValueInternalArray();
}

virtual ValueInternalArray *newArrayCopy( const ValueInternalArray &other )
{
   return new ValueInternalArray( other );
}

virtual void destruct( ValueInternalArray *array )
{
   delete array;
}

virtual void reallocateArrayPageIndex( Value **&indexes,
                                       ValueInternalArray::PageIndex &indexCount,
                                       ValueInternalArray::PageIndex minNewIndexCount )
{
   ValueInternalArray::PageIndex newIndexCount = (indexCount*3)/2 + 1;
   if ( minNewIndexCount > newIndexCount )
      newIndexCount = minNewIndexCount;
   void *newIndexes = realloc( indexes, sizeof(Value*) * newIndexCount );
   if ( !newIndexes )
      throw std::bad_alloc();
   indexCount = newIndexCount;
   indexes = static_cast<Value **>( newIndexes );
}
virtual void releaseArrayPageIndex( Value **indexes,
                                    ValueInternalArray::PageIndex indexCount )
{
   if ( indexes )
      free( indexes );
}

virtual Value *allocateArrayPage()
{
   return static_cast<Value *>( malloc( sizeof(Value) * ValueInternalArray::itemsPerPage ) );
}

virtual void releaseArrayPage( Value *value )
{
   if ( value )
      free( value );
}
};
   \endcode
 */
class JSON_API ValueArrayAllocator
{
 public:
  virtual ~ValueArrayAllocator();
  virtual ValueInternalArray *newArray() = 0;
  virtual ValueInternalArray *newArrayCopy( const ValueInternalArray &other ) = 0;
  virtual void destructArray( ValueInternalArray *array ) = 0;
  /** \brief Reallocate array page index.
   * Reallocates an array of pointer on each page.
   * \param indexes [input] pointer on the current index. May be \c NULL.
   *                [output] pointer on the new index of at least
   *                         \a minNewIndexCount pages.
   * \param indexCount [input] current number of pages in the index.
   *                   [output] number of page the reallocated index can handle.
   *                            \b MUST be >= \a minNewIndexCount.
   * \param minNewIndexCount Minimum number of page the new index must be able to
   *                         handle.
   */
  virtual void reallocateArrayPageIndex( Value **&indexes,
                                         ValueInternalArray::PageIndex &indexCount,
                                         ValueInternalArray::PageIndex minNewIndexCount ) = 0;
  virtual void releaseArrayPageIndex( Value **indexes,
                                      ValueInternalArray::PageIndex indexCount ) = 0;
  virtual Value *allocateArrayPage() = 0;
  virtual void releaseArrayPage( Value *value ) = 0;
};
#endif // #ifdef JSON_VALUE_USE_INTERNAL_MAP


/** \brief base class for Value iterators.
 *
 */
class ValueIteratorBase
{
 public:
  typedef unsigned int size_t;
  typedef int difference_type;
  typedef ValueIteratorBase SelfType;

  ValueIteratorBase();
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  explicit ValueIteratorBase( const Value::ObjectValues::iterator &current );
#else
  ValueIteratorBase( const ValueInternalArray::IteratorState &state );
  ValueIteratorBase( const ValueInternalMap::IteratorState &state );
#endif

  bool operator ==( const SelfType &other ) const {
    return isEqual( other );
  }

  bool operator !=( const SelfType &other ) const {
    return !isEqual( other );
  }

  difference_type operator -( const SelfType &other ) const {
    return computeDistance( other );
  }

  /// Return either the index or the member name of the referenced value as a Value.
  Value key() const;

  /// Return the index of the referenced Value. -1 if it is not an arrayValue.
  UInt index() const;

  /// Return the member name of the referenced Value. "" if it is not an objectValue.
  const char *memberName() const;

 protected:
  Value &deref() const;

  void increment();

  void decrement();

  difference_type computeDistance( const SelfType &other ) const;

  bool isEqual( const SelfType &other ) const;

  void copy( const SelfType &other );

 private:
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  Value::ObjectValues::iterator current_;
  // Indicates that iterator is for a null value.
  bool isNull_;
#else
  union {
    ValueInternalArray::IteratorState array_;
    ValueInternalMap::IteratorState map_;
  } iterator_;
  bool isArray_;
#endif
};

/** \brief const iterator for object and array value.
 *
 */
class ValueConstIterator : public ValueIteratorBase
{
  friend class Value;
 public:
  typedef unsigned int size_t;
  typedef int difference_type;
  typedef const Value &reference;
  typedef const Value *pointer;
  typedef ValueConstIterator SelfType;

  ValueConstIterator();
 private:
  /*! \internal Use by Value to create an iterator.
   */
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  explicit ValueConstIterator( const Value::ObjectValues::iterator &current );
#else
  ValueConstIterator( const ValueInternalArray::IteratorState &state );
  ValueConstIterator( const ValueInternalMap::IteratorState &state );
#endif
 public:
  SelfType &operator =( const ValueIteratorBase &other );

  SelfType operator++( int ) {
    SelfType temp( *this );
    ++*this;
    return temp;
  }

  SelfType operator--( int ) {
    SelfType temp( *this );
    --*this;
    return temp;
  }

  SelfType &operator--() {
    decrement();
    return *this;
  }

  SelfType &operator++() {
    increment();
    return *this;
  }

  reference operator *() const {
    return deref();
  }
};


/** \brief Iterator for object and array value.
 */
class ValueIterator : public ValueIteratorBase
{
  friend class Value;
 public:
  typedef unsigned int size_t;
  typedef int difference_type;
  typedef Value &reference;
  typedef Value *pointer;
  typedef ValueIterator SelfType;

  ValueIterator();
  ValueIterator( const ValueConstIterator &other );
  ValueIterator( const ValueIterator &other );
 private:
  /*! \internal Use by Value to create an iterator.
   */
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  explicit ValueIterator( const Value::ObjectValues::iterator &current );
#else
  ValueIterator( const ValueInternalArray::IteratorState &state );
  ValueIterator( const ValueInternalMap::IteratorState &state );
#endif
 public:

  SelfType &operator =( const SelfType &other );

  SelfType operator++( int ) {
    SelfType temp( *this );
    ++*this;
    return temp;
  }

  SelfType operator--( int ) {
    SelfType temp( *this );
    --*this;
    return temp;
  }

  SelfType &operator--() {
    decrement();
    return *this;
  }

  SelfType &operator++() {
    increment();
    return *this;
  }

  reference operator *() const {
    return deref();
  }
};


} // namespace Json


#endif // CPPTL_JSON_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/value.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/reader.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef CPPTL_JSON_READER_H_INCLUDED
# define CPPTL_JSON_READER_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "features.h"
# include "value.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
# include <deque>
# include <stack>
# include <string>
# include <iostream>

namespace Json
{

/** \brief Unserialize a <a HREF="http://www.json.org">JSON</a> document into a Value.
 *
 */
class JSON_API Reader
{
 public:
  typedef char Char;
  typedef const Char *Location;

  /** \brief Constructs a Reader allowing all features
   * for parsing.
   */
  Reader();

  /** \brief Constructs a Reader allowing the specified feature set
   * for parsing.
   */
  Reader( const Features &features );

  /** \brief Read a Value from a <a HREF="http://www.json.org">JSON</a> document.
   * \param document UTF-8 encoded string containing the document to read.
   * \param root [out] Contains the root value of the document if it was
   *             successfully parsed.
   * \param collectComments \c true to collect comment and allow writing them back during
   *                        serialization, \c false to discard comments.
   *                        This parameter is ignored if Features::allowComments_
   *                        is \c false.
   * \return \c true if the document was successfully parsed, \c false if an error occurred.
   */
  bool parse( const std::string &document,
              Value &root,
              bool collectComments = true );

  /** \brief Read a Value from a <a HREF="http://www.json.org">JSON</a> document.
   * \param beginDoc Pointer on the beginning of the UTF-8 encoded string of the document to read.
   * \param endDoc Pointer on the end of the UTF-8 encoded string of the document to read.
   \               Must be >= beginDoc.
   * \param root [out] Contains the root value of the document if it was
   *             successfully parsed.
   * \param collectComments \c true to collect comment and allow writing them back during
   *                        serialization, \c false to discard comments.
   *                        This parameter is ignored if Features::allowComments_
   *                        is \c false.
   * \return \c true if the document was successfully parsed, \c false if an error occurred.
   */
  bool parse( const char *beginDoc, const char *endDoc,
              Value &root,
              bool collectComments = true );

  /// \brief Parse from input stream.
  /// \see Json::operator>>(std::istream&, Json::Value&).
  bool parse( std::istream &is,
              Value &root,
              bool collectComments = true );

  /** \brief Returns a user friendly string that list errors in the parsed document.
   * \return Formatted error message with the list of errors with their location in
   *         the parsed document. An empty string is returned if no error occurred
   *         during parsing.
   * \deprecated Use getFormattedErrorMessages() instead (typo fix).
   */
  JSONCPP_DEPRECATED("Use getFormattedErrorMessages instead")
  std::string getFormatedErrorMessages() const;

  /** \brief Returns a user friendly string that list errors in the parsed document.
   * \return Formatted error message with the list of errors with their location in
   *         the parsed document. An empty string is returned if no error occurred
   *         during parsing.
   */
  std::string getFormattedErrorMessages() const;

 private:
  enum TokenType {
      tokenEndOfStream = 0,
      tokenObjectBegin,
      tokenObjectEnd,
      tokenArrayBegin,
      tokenArrayEnd,
      tokenString,
      tokenNumber,
      tokenTrue,
      tokenFalse,
      tokenNull,
      tokenArraySeparator,
      tokenMemberSeparator,
      tokenComment,
      tokenError
  };

  class Token
{
   public:
    TokenType type_;
    Location start_;
    Location end_;
  };

  class ErrorInfo
{
   public:
    Token token_;
    std::string message_;
    Location extra_;
  };

  typedef std::deque<ErrorInfo> Errors;

  bool expectToken( TokenType type, Token &token, const char *message );
  bool readToken( Token &token );
  void skipSpaces();
  bool match( Location pattern,
              int patternLength );
  bool readComment();
  bool readCStyleComment();
  bool readCppStyleComment();
  bool readString();
  void readNumber();
  bool readValue();
  bool readObject( Token &token );
  bool readArray( Token &token );
  bool decodeNumber( Token &token );
  bool decodeString( Token &token );
  bool decodeString( Token &token, std::string &decoded );
  bool decodeDouble( Token &token );
  bool decodeUnicodeCodePoint( Token &token,
                               Location &current,
                               Location end,
                               unsigned int &unicode );
  bool decodeUnicodeEscapeSequence( Token &token,
                                    Location &current,
                                    Location end,
                                    unsigned int &unicode );
  bool addError( const std::string &message,
                 Token &token,
                 Location extra = 0 );
  bool recoverFromError( TokenType skipUntilToken );
  bool addErrorAndRecover( const std::string &message,
                           Token &token,
                           TokenType skipUntilToken );
  void skipUntilSpace();
  Value &currentValue();
  Char getNextChar();
  void getLocationLineAndColumn( Location location,
                                 int &line,
                                 int &column ) const;
  std::string getLocationLineAndColumn( Location location ) const;
  void addComment( Location begin,
                   Location end,
                   CommentPlacement placement );
  void skipCommentTokens( Token &token );

  typedef std::stack<Value *> Nodes;
  Nodes nodes_;
  Errors errors_;
  std::string document_;
  Location begin_;
  Location end_;
  Location current_;
  Location lastValueEnd_;
  Value *lastValue_;
  std::string commentsBefore_;
  Features features_;
  bool collectComments_;
};

/** \brief Read from 'sin' into 'root'.

 Always keep comments from the input JSON.

 This can be used to read a file into a particular sub-object.
 For example:
 \code
 Json::Value root;
 cin >> root["dir"]["file"];
 cout << root;
 \endcode
 Result:
 \verbatim
 {
 "dir": {
     "file": {
     // The input stream JSON would be nested here.
     }
 }
 }
 \endverbatim
 \throw std::exception on parse error.
 \see Json::operator<<()
*/
std::istream& operator>>( std::istream&, Value& );

} // namespace Json

#endif // CPPTL_JSON_READER_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/reader.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/writer.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_WRITER_H_INCLUDED
# define JSON_WRITER_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
# include "value.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
# include <vector>
# include <string>
# include <iostream>

namespace Json
{

class Value;

/** \brief Abstract class for writers.
 */
class JSON_API Writer
{
 public:
  virtual ~Writer();

  virtual std::string write( const Value &root ) = 0;
};

/** \brief Outputs a Value in <a HREF="http://www.json.org">JSON</a> format without formatting (not human friendly).
 *
 * The JSON document is written in a single line. It is not intended for 'human' consumption,
 * but may be usefull to support feature such as RPC where bandwith is limited.
 * \sa Reader, Value
 */
class JSON_API FastWriter : public Writer
{
 public:
  FastWriter();
  virtual ~FastWriter() {}

  void enableYAMLCompatibility();

 public: // overridden from Writer
  virtual std::string write( const Value &root );

 private:
  void writeValue( const Value &value );

  std::string document_;
  bool yamlCompatiblityEnabled_;
};

/** \brief Writes a Value in <a HREF="http://www.json.org">JSON</a> format in a human friendly way.
 *
 * The rules for line break and indent are as follow:
 * - Object value:
 *     - if empty then print {} without indent and line break
 *     - if not empty the print '{', line break & indent, print one value per line
 *       and then unindent and line break and print '}'.
 * - Array value:
 *     - if empty then print [] without indent and line break
 *     - if the array contains no object value, empty array or some other value types,
 *       and all the values fit on one lines, then print the array on a single line.
 *     - otherwise, it the values do not fit on one line, or the array contains
 *       object or non empty array, then print one value per line.
 *
 * If the Value have comments then they are outputed according to their #CommentPlacement.
 *
 * \sa Reader, Value, Value::setComment()
 */
class JSON_API StyledWriter: public Writer
{
 public:
  StyledWriter();
  virtual ~StyledWriter() {}

 public: // overridden from Writer
  /** \brief Serialize a Value in <a HREF="http://www.json.org">JSON</a> format.
   * \param root Value to serialize.
   * \return String containing the JSON document that represents the root value.
   */
  virtual std::string write( const Value &root );

 private:
  void writeValue( const Value &value );
  void writeArrayValue( const Value &value );
  bool isMultineArray( const Value &value );
  void pushValue( const std::string &value );
  void writeIndent();
  void writeWithIndent( const std::string &value );
  void indent();
  void unindent();
  void writeCommentBeforeValue( const Value &root );
  void writeCommentAfterValueOnSameLine( const Value &root );
  bool hasCommentForValue( const Value &value );
  static std::string normalizeEOL( const std::string &text );

  typedef std::vector<std::string> ChildValues;

  ChildValues childValues_;
  std::string document_;
  std::string indentString_;
  int rightMargin_;
  int indentSize_;
  bool addChildValues_;
};

/** \brief Writes a Value in <a HREF="http://www.json.org">JSON</a> format in a human friendly way,
     to a stream rather than to a string.
 *
 * The rules for line break and indent are as follow:
 * - Object value:
 *     - if empty then print {} without indent and line break
 *     - if not empty the print '{', line break & indent, print one value per line
 *       and then unindent and line break and print '}'.
 * - Array value:
 *     - if empty then print [] without indent and line break
 *     - if the array contains no object value, empty array or some other value types,
 *       and all the values fit on one lines, then print the array on a single line.
 *     - otherwise, it the values do not fit on one line, or the array contains
 *       object or non empty array, then print one value per line.
 *
 * If the Value have comments then they are outputed according to their #CommentPlacement.
 *
 * \param indentation Each level will be indented by this amount extra.
 * \sa Reader, Value, Value::setComment()
 */
class JSON_API StyledStreamWriter
{
 public:
  StyledStreamWriter( std::string indentation="\t" );
  ~StyledStreamWriter() {}

 public:
  /** \brief Serialize a Value in <a HREF="http://www.json.org">JSON</a> format.
   * \param out Stream to write to. (Can be ostringstream, e.g.)
   * \param root Value to serialize.
   * \note There is no point in deriving from Writer, since write() should not return a value.
   */
  void write( std::ostream &out, const Value &root );

 private:
  void writeValue( const Value &value );
  void writeArrayValue( const Value &value );
  bool isMultineArray( const Value &value );
  void pushValue( const std::string &value );
  void writeIndent();
  void writeWithIndent( const std::string &value );
  void indent();
  void unindent();
  void writeCommentBeforeValue( const Value &root );
  void writeCommentAfterValueOnSameLine( const Value &root );
  bool hasCommentForValue( const Value &value );
  static std::string normalizeEOL( const std::string &text );

  typedef std::vector<std::string> ChildValues;

  ChildValues childValues_;
  std::ostream* document_;
  std::string indentString_;
  int rightMargin_;
  std::string indentation_;
  bool addChildValues_;
};

# if defined(JSON_HAS_INT64)
std::string JSON_API valueToString( Int value );
std::string JSON_API valueToString( UInt value );
# endif // if defined(JSON_HAS_INT64)
std::string JSON_API valueToString( LargestInt value );
std::string JSON_API valueToString( LargestUInt value );
std::string JSON_API valueToString( double value );
std::string JSON_API valueToString( bool value );
std::string JSON_API valueToQuotedString( const char *value );

/// \brief Output using the StyledStreamWriter.
/// \see Json::operator>>()
std::ostream& operator<<( std::ostream&, const Value &root );

} // namespace Json



#endif // JSON_WRITER_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/writer.h
// //////////////////////////////////////////////////////////////////////





#endif //ifndef JSON_AMALGATED_H_INCLUDED
//
// end of src/json.h
//


//
// start of src/jsoncustomwriter.h
//
/**********************************************************************
Copyright (c) 2013 by Matt Swain <m.swain@me.com>

The MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

***********************************************************************/

#ifndef PYNE_46Z7LQYFI5HZNASIPCWHVX3X5E
#define PYNE_46Z7LQYFI5HZNASIPCWHVX3X5E

#include <string>

namespace Json
{

/** \brief Writes a Value in <a HREF="http://www.json.org">JSON</a> format with custom formatting.
 *
 * The JSON document is written according to the rules specified in the constructor. Objects and
 * arrays are printed on a single line if they are below a certain length, otherwise they are
 * indented. It is possible to output invalid json if the customizable parameters are specified
 * incorrectly. Set maxWidth to 0 to print output on a single line.
 *
 * \sa Reader, Value
 */
class JSON_API CustomWriter : public Writer
{
 public:
  CustomWriter( std::string opencurly = "{",
                std::string closecurly = "}",
                std::string opensquare = "[",
                std::string closesquare = "]",
                std::string colon = ":",
                std::string comma = ",",
                std::string indent = "  ",
                int maxWidth = 74);
  virtual ~CustomWriter() {}

 public: // overridden from Writer
  virtual std::string write( const Value &root );

 private:
  void writeValue( const Value &value, std::string &doc, bool forceSingleLine );
  bool isMultiline( const Value &value );
  void indent();
  void unindent();

  std::string document_;
  std::string indentString_;
  std::string opencurly_;
  std::string closecurly_;
  std::string opensquare_;
  std::string closesquare_;
  std::string colon_;
  std::string comma_;
  std::string indent_;
  int maxWidth_;
};

}

#endif//
// end of src/jsoncustomwriter.h
//


//
// start of src/material.h
//
/// \brief The ever-important material class and related helpers.
///
/// The material class is effectively a normalized nuclide linked list with
/// associated mass, density, atoms per mol, and metadata.  However, this
/// implementation also contains other functions for mixing materials and generating
/// related materials.

#ifndef PYNE_MR34UE5INRGMZK2QYRDWICFHVM
#define PYNE_MR34UE5INRGMZK2QYRDWICFHVM

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>	// std::ostringstream

#if !defined(JSON_IS_AMALGAMATION)
#define JSON_IS_AMALGAMATION
#endif

#ifndef PYNE_IS_AMALGAMATED
#include "json-forwards.h"
#include "json.h"
#include "h5wrap.h"
#include "utils.h"
#include "nucname.h"
#include "data.h"
#include "decay.h"
#endif

namespace pyne
{
// Set Type Definitions
typedef std::map<int, double> comp_map; ///< Nuclide-mass composition map type
typedef comp_map::iterator comp_iter;   ///< Nuclide-mass composition iter type

#ifdef PYNE_IS_AMALGAMATED
namespace decayers
{
extern comp_map decay(comp_map, double);
}  // namespace decayers
#endif


// These 37 strings are predefined FLUKA materials.
// Materials not on this list requires a MATERIAL card.
static std::string fluka_mat_strings[] = {
  "BLCKHOLE", "VACUUM",   "HYDROGEN", "HELIUM",   "BERYLLIU", "CARBON",
  "NITROGEN", "OXYGEN",   "MAGNESIU", "ALUMINUM", "IRON",     "COPPER",
  "SILVER",   "SILICON",  "GOLD",     "MERCURY",  "LEAD",     "TANTALUM",
  "SODIUM",   "ARGON",    "CALCIUM",  "TIN",      "TUNGSTEN", "TITANIUM",
  "NICKEL",   "WATER",    "POLYSTYR", "PLASCINT", "PMMA",     "BONECOMP",
  "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON", "POLYETHY", "AIR"
};

static int FLUKA_MAT_NUM = 37;

/// Material composed of nuclides.
class Material
{
 protected:

  /// Computes the total mass stored in the composition.
  double get_comp_sum ();

 public:

  // Material Constructors
  Material ();  ///< empty constructor
  /// Constructor from composition map
  /// \param cm composition map
  /// \param m mass value, the mass is set to the sum of the values in the
  ///          composition if \a m is negative.
  /// \param d density value
  /// \param apm atoms per mole
  /// \param attributes initial metadata
  Material(comp_map cm, double m=-1.0, double d=-1.0, double apm=-1.0,
           Json::Value attributes=Json::Value(Json::objectValue));
  /// Constructor from file
  /// \param filename path to file on disk, this file may be either in plaintext
  ///                 or HDF5 format.
  /// \param m mass value, the mass is set to the sum of the values in the
  ///          composition if \a m is negative,
  ///          may be overridden by the value from disk.
  /// \param d density value,
  ///          may be overridden by the value from disk.
  /// \param apm atoms per mole,
  ///          may be overridden by the value from disk.
  /// \param attributes initial metadata,
  ///          may be overridden by the value from disk.
  Material(char * filename, double m=-1.0, double d=-1.0, double apm=-1.0,
           Json::Value attributes=Json::Value(Json::objectValue));
  /// Constructor from file
  /// \param filename path to file on disk, this file may be either in plaintext
  ///                 or HDF5 format.
  /// \param m mass value, the mass is set to the sum of the values in the
  ///          composition if \a m is negative,
  ///          may be overridden by the value from disk.
  /// \param d density value,
  ///          may be overridden by the value from disk.
  /// \param apm atoms per mole,
  ///          may be overridden by the value from disk.
  /// \param attributes initial metadata,
  ///          may be overridden by the value from disk.
  Material(std::string filename, double m=-1.0, double d=-1.0, double apm=-1.0,
           Json::Value attributes=Json::Value(Json::objectValue));
  ~Material (); ///< default destructor

  /// Normalizes the mass values in the composition.
  void norm_comp ();

  // Persistence functions.

  /// Loads the matrial composition from an HDF5 file according to the layout
  /// defined by protocol 0.  This protocol is depratacted.
  /// \param db HDF5 id for the open HDF5 file.
  /// \param datapath Path to the base node for the material in \a db.
  /// \param row The index to read out, may be negative.
  void _load_comp_protocol0(hid_t db, std::string datapath, int row);

  /// Loads the matrial composition from an HDF5 file according to the layout
  /// defined by protocol 1.  This protocol should be used in favor of protocol 0.
  /// \param db HDF5 id for the open HDF5 file.
  /// \param datapath Path to the base node for the material in \a db.
  /// \param row The index to read out, may be negative.
  void _load_comp_protocol1(hid_t db, std::string datapath, int row);

  /// Loads a material from an HDF5 file into this object.
  /// \param filename Path on disk to the HDF5 file.
  /// \param datapath Path to the the material in the file.
  /// \param row The index to read out, may be negative.
  /// \param protocol Flag for layout of material on disk.
  void from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1);

  /// Loads a material from an HDF5 file into this object.
  /// \param filename Path on disk to the HDF5 file.
  /// \param datapath Path to the the material in the file.
  /// \param row The index to read out, may be negative.
  /// \param protocol Flag for layout of material on disk.
  void from_hdf5(std::string filename, std::string datapath="/material",
                 int row=-1, int protocol=1);

  /// Writes this material out to an HDF5 file.
  /// This happens according to protocol 1.
  /// \param filename Path on disk to the HDF5 file.
  /// \param datapath Path to the the material in the file.
  /// \param nucpath Path to the nuclides set in the file.
  /// \param row The index to read out, may be negative. Also note that this is a
  ///            float.  A value of -0.0 indicates that the material should be
  ///            appended to the end of the dataset.
  /// \param chunksize The chunksize for all material data on disk.
  void write_hdf5(char * filename, char * datapath, char * nucpath, float row=-0.0,
                  int chunksize=100);
  /// Writes this material out to an HDF5 file.
  /// This happens according to protocol 1.
  /// \param filename Path on disk to the HDF5 file.
  /// \param datapath Path to the the material in the file.
  /// \param nucpath Path to the nuclides set in the file.
  /// \param row The index to read out, may be negative. Also note that this is a
  ///            float.  A value of -0.0 indicates that the material should be
  ///            appended to the end of the dataset.
  /// \param chunksize The chunksize for all material data on disk.
  void write_hdf5(std::string filename, std::string datapath="/material",
                  std::string nucpath="/nucid", float row=-0.0, int chunksize=100);

  /// Return an mcnp input deck record as a string
  std::string mcnp(std::string frac_type = "mass");
  ///
  /// Return a fluka input deck MATERIAL card as a string
  std::string fluka(int id, std::string frac_type = "mass");
  /// Convenience function to tell whether a given name needs a material card
  bool not_fluka_builtin(std::string fluka_name);
  /// High level call to get details and call material_component(..)
  std::string fluka_material_str(int id);
  /// Intermediate level call to prepare final info and call material_line(..)
  std::string fluka_material_component(int fid, int nucid,
                                       std::string fluka_name);
  /// Format information into a FLUKA material card
  std::string fluka_material_line(int znum, double atomic_mass,
                                  int fid, std::string fluka_name);
  /// Convenience function to format a single fluka field
  std::string fluka_format_field(float field);
  /// Return FLUKA compound card and the material card for the named compound
  /// but not the material cards of the components
  std::string fluka_compound_str(int id, std::string frac_type = "mass");

  /// Reads data from a plaintext file at \a filename into this Material instance.
  void from_text(char * filename);
  /// Reads data from a plaintext file at \a filename into this Material instance.
  void from_text(std::string filename);

  /// Writes the Material out to a simple plaintext file readable by from_text().
  void write_text(char * filename);
  /// Writes the Material out to a simple plaintext file readable by from_text().
  void write_text(std::string filename);

  /// Loads a JSON instance tree into this Material.
  void load_json(Json::Value);
  /// Dumps the Material out to a JSON instance tree.
  Json::Value dump_json();
  /// Reads data from a JSON file at \a filename into this Material instance.
  void from_json(char * filename);
  /// Reads data from a JSON file at \a filename into this Material instance.
  void from_json(std::string filname);
  /// Writes the Material out to a JSON file
  void write_json(char * filename);
  /// Writes the Material out to a JSON file
  void write_json(std::string filename);

  // Fundemental mass stream data
  /// composition, maps nuclides in id form to normalized mass weights.
  comp_map comp;
  double mass;  ///< mass (in arbitrary units) of the Material.
  double density; ///< density (in arbitrary units) of the Material.
  double atoms_per_molecule; ///< The number of atoms per molecule.
  /// container for arbitrary metadata, following the JSON rules.
  Json::Value metadata;

  // Material function definitions
  void normalize ();  ///< Normalizes the mass.
  /// Returns a composition map that has been unnormalized by multiplying each
  /// mass weight by the actual mass of the material.
  comp_map mult_by_mass();
  /// Calculates the atomic weight of this material based on the composition
  /// and the number of atoms per mol.  If \a apm is non-negative then it is
  /// used (and stored on the instance) as the atoms_per_molecule for this calculation.
  /// If \a apm and atoms_per_molecule on this instance are both negative, then the best
  /// guess value calculated from the normailized composition is used here.
  double molecular_mass(double apm=-1.0);
  /// Calculates the activity of a material based on the composition and each
  /// nuclide's mass, decay_const, and atmoic_mass.
  comp_map activity();
  /// Calculates the decay heat of a material based on the composition and
  /// each nuclide's mass, q_val, decay_const, and atomic_mass.
  comp_map decay_heat();
  /// Caclulates the dose per gram using the composition of the the
  /// material, the dose type desired, and the source for dose factors
  ///   dose_type is one of:
  ///     ext_air -- returns mrem/h per g per m^3
  ///     ext_soil -- returns mrem/h per g per m^2
  ///     ingest -- returns mrem per g
  ///     inhale -- returns mrem per g
  ///   source is:
  ///     {EPA=0, DOE=1, GENII=2}, default is EPA
  comp_map dose_per_g(std::string dose_type, int source=0);
  /// Returns a copy of the current material where all natural elements in the
  /// composition are expanded to their natural isotopic abundances.
  Material expand_elements();
  // Returns a copy of the current material where all the isotopes of the elements
  // are added up, atomic-fraction-wise, unless they are in the exception set
  Material collapse_elements(std::set<int> exception_znum);
  // Wrapped version to facilitate calling from python
  Material collapse_elements(int **int_ptr_arry);
  // void print_material( pyne::Material test_mat);
  /// Computes, sets, and returns the mass density when \a num_dens is greater
  /// than or equal zero.  If \a num_dens is negative, this simply returns the
  /// current value of the density member variable.  You may also use / set the
  /// atoms per molecule (atoms_per_molecule) in this function using \a apm.
  double mass_density(double num_dens=-1.0, double apm=-1.0);
  /// Computes and returns the number density of the material using the
  /// mass density if \a mass_dens is greater than or equal to zero.  If
  /// \a mass_dens is negative, the denisty member variable is used instead.
  /// You may also use / set the atoms per molecule (atoms_per_molecule) in this
  /// function using \a apm.
  double number_density(double mass_dens=-1.0, double apm=-1.0);

  // Sub-Stream Computation
  /// Creates a sub-Material with only the nuclides present in \a nucset.
  /// Elements of this set may be either in id form or simple Z numbers.
  Material sub_mat(std::set<int> nucset);
  /// Creates a sub-Material with only the nuclides present in \a nucset.
  /// Elements of this set may be in any form.
  Material sub_mat(std::set<std::string> nucset);

  /// Creates a new Material with the mass weights for all nuclides in \a nucset
  /// set to \a value.
  /// Elements of \a nucset may be either in id form or simple Z numbers.
  Material set_mat(std::set<int> nucset, double value);
  /// Creates a new Material with the mass weights for all nuclides in \a nucset
  /// set to \a value.  Elements of \a nucset may be in any form.
  Material set_mat(std::set<std::string> nucset, double value);

  /// Creates a new Material with the all nuclides in \a nucset removed.
  /// Elements of \a nucset may be either in id form or simple Z numbers.
  Material del_mat(std::set<int> nucset);
  /// Creates a new Material with the all nuclides in \a nucset removed.
  /// Elements of \a nucset may be in any form.
  Material del_mat(std::set<std::string> nucset);

  /// Creates a sub-Material based on a range of id-form integers.
  Material sub_range(int lower=0, int upper=10000000);
  /// Creates a new Material with the mass weights for all nuclides in the id
  /// range set to \a value.
  Material set_range(int lower=0, int upper=10000000, double value=0.0);
  /// Creates a new Material with the all nuclides in the id range removed.
  Material del_range(int lower=0, int upper=10000000);

  /// Creates a sub-Material of only the given element. Assumes element is
  /// id form.
  Material sub_elem(int element);
  /// Creates a sub-Material of only lanthanides.
  Material sub_lan();
  /// Creates a sub-Material of only actinides.
  Material sub_act();
  /// Creates a sub-Material of only transuranics.
  Material sub_tru();
  /// Creates a sub-Material of only minor actinides.
  Material sub_ma();
  /// Creates a sub-Material of only fission products.
  Material sub_fp();

  // Atom fraction functions
  /// Returns a mapping of the nuclides in this material to their atom fractions.
  /// This calculation is based off of the material's molecular weight.
  std::map<int, double> to_atom_frac();
  /// Sets the composition, mass, and atoms_per_molecule of this material to those
  /// calculated from \a atom_fracs, a mapping of nuclides to atom fractions values.
  void from_atom_frac(std::map<int, double> atom_fracs);

  /// Returns a mapping of the nuclides in this material to their atom densities.
  /// This calculation is based off of the material's density.
  std::map<int, double> to_atom_dens();

  // Radioactive Material functions
  /// Returns a list of gamma-rays energies in keV and intensities in
  /// decays/s/atom material unnormalized
  std::vector<std::pair<double, double> > gammas();
  /// Returns a list of x-rays average energies in keV and intensities in
  /// decays/s material unnormalized
  std::vector<std::pair<double, double> > xrays();
  /// Returns a list of photon energies in keV and intensities in
  /// decays/s/atom material unnormalized
  std::vector<std::pair<double, double> > photons(bool norm);
  /// Takes a list of photon energies and intensities and normalizes them
  /// so the sum of the intensities is one
  std::vector<std::pair<double, double> > normalize_radioactivity(
      std::vector<std::pair<double, double> > unnormed);

  /// Decays this material for a given amount of time in seconds
  Material decay(double t);

  // Overloaded Operators
  /// Adds mass to a material instance.
  Material operator+ (double);
  /// Adds two materials together.
  Material operator+ (Material);
  /// Multiplies a material's mass.
  Material operator* (double);
  /// Divides a material's mass.
  Material operator/ (double);
};

/// Converts a Material to a string stream representation for canonical writing.
/// This operator is also defined on inheritors of std::ostream
std::ostream& operator<< (std::ostream& os, Material mat);

/// A stuct for reprensenting fundemental data in a material.
/// Useful for HDF5 representations.
typedef struct material_data {
  double mass;  ///< material mass
  double density; ///< material density
  double atoms_per_mol; ///< material atoms per mole
  double comp[1]; ///< array of material composition mass weights.
} material_data;

/// Custom exception for invalid HDF5 protocol numbers
class MaterialProtocolError: public std::exception
{
  /// marginally helpful error message.
  virtual const char* what() const throw() {
    return "Invalid loading protocol number; please use 0 or 1.";
  }
};

// End pyne namespace
}

#endif  // PYNE_MR34UE5INRGMZK2QYRDWICFHVM
//
// end of src/material.h
//


//
// start of src/tally.h
//
/// \brief The tally class and helper functions.
///
/// The tally class is in essesence a structure containing attributes
/// related to tallies.

#ifndef PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
#define PYNE_IQ4M73STINHJDPRV6KWUZZXOYE

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>

#ifndef PYNE_IS_AMALGAMATED
#include "h5wrap.h"
#include "utils.h"
#endif


namespace pyne
{
class Tally
{
 public:
  std::map<std::string, std::string> rx2fluka;
  std::map<std::string, std::string> rx2mcnp5;
  std::map<std::string, std::string> rx2mcnp6;

  /// Tally Constructors
  Tally (); /// empty constructor

  /// Constructor from passed in vars
  /// \param type the type of tally (flux or current)
  /// \param particle_name the name of the particle type
  /// \param entity the entity id of the tally (eg. surface index,
  ///          volume number)
  /// \param entity_type (volume or surface)
  /// \param entity_name string identifying the entity
  /// \param tally_name string identifying the tally
  /// \param entity_size the physical size of the tally volume
  /// \param normalization the number required to normalize your tally
  Tally(std::string type, std::string particle_name, int entity,
        std::string entity_type, std::string entity_name,
        std::string tally_name = "", double entity_size = 0.0,
        double normalization = 1.0);

  ~Tally (); /// default destructor


  // Create hdf5 datatable for tallies
  hid_t create_dataspace(hid_t file, std::string datapath);

  // create filetype of data
  hid_t create_filetype();

  // create memory type for tally struct
  hid_t create_memtype();

  /// Dummy read method wrapper around c style strings
  /// \param filename the filename of the file to read from
  /// \param datapath _name the name of the region where tallies
  ///          are stored
  /// \param row  the array index of data to access
  void from_hdf5(char * filename, char *datapath, int row = -1);

  /// Main read tally method
  /// \param filename the filename of the file to read from
  /// \param datapath _name the name of the region where tallies
  ///          are stored
  /// \param row  the array index of data to access
  void from_hdf5(std::string filename, std::string datapath, int row = -1);

  /// Dummy write method wrapper around c style strings
  /// \param filename the filename of the file to write to
  /// \param datapath _name the name of the region where tallies
  ///          are to be stored
  void write_hdf5( char * filename, char * datapath);

  /// Main write tally method
  /// \param filename the filename of the file to write to
  /// \param datapath _name the name of the region where tallies
  ///          are to be stored
  void write_hdf5(std::string filename, std::string datapath);

  // mcnp tally
  std::string mcnp(int tally_index = 1, std::string mcnp_version = "mcnp5" );

  // fluka tally
  std::string fluka(std::string unit_number = "-21");

  /// fundamental tally variables
  std::string entity_type; ///< the type of entity (volume,surface)
  std::string entity_name; ///< the name of the entity (optional)
  std::string particle_name; ///< particle name string
  std::string tally_type; ///< type of tally flux or current
  std::string tally_name; ///< name of the tally
  int entity_id; ///< id number of the entity being tallied upon
  double entity_size; ///< the physical size of the entity
  double normalization; ///< the tally normalization
};

/// Converts a Tally to a string stream representation.
std::ostream& operator<< (std::ostream& os, Tally tally);


/// A stuct for reprensenting fundemental data in a tally
/// Maybe Useful for HDF5 representations.
/// following scoptaz's lead here
typedef struct tally_struct {
  int entity_id;
  int entity_type;
  int tally_type;
  const char * particle_name;
  const char * entity_name;
  const char * tally_name;
  double entity_size;
  double normalization;
} tally_struct;

// End pyne namespace
}

#endif
//
// end of src/tally.h
//


//
// start of src/atomic_data.cpp~
//
// // Implements basic nuclear data functions.
// #ifndef PYNE_IS_AMALGAMATED
// #include "atomic_data.h"
// #endif
//
// void pyne::_load_atomic_mass_map_memory() {
//   // header version of atomic weight table data
//   //see if the data table is already loaded
//   if(!atomic_mass_map.empty()) {
//     return;
//   } else {
//     _insert_atomic_mass_map();
//   }
//   //see if the data table is already loaded
//   if(!natural_abund_map.empty()) {
//     return;
//   } else {
//     _insert_abund_map();
//   }
//
//   // calculate the atomic_masses of the elements
//   std::map<int,double> :: iterator it;
//
//   for ( int i = 0 ; i < 92 ; i++ ) {
//     // loop through the natural abundance map
//     double element_atomic_weight = 0.0;
//     int atomic_number = i + 1;
//     for ( it = natural_abund_map.begin() ; it != natural_abund_map.end() ; ++it ){
//       // if the atomic number of the abudance matches the
//       // that of index
//       if(pyne::nucname::znum(it->first) == atomic_number ) {
// 	// take atomic abundance and multiply by mass
// 	// to get the mass of that nuclide / 100 since abundance is in %
// 	element_atomic_weight += (it->second*atomic_mass_map[it->first]/100.0);
//       }
//     }
//     // insert the abundance of the element into the list
//     atomic_mass_map[i*10000000] = element_atomic_weight;
//   }
// }
//
// void pyne::_insert_atomic_mass_map() {
//   atomic_mass_map[10010000] = 1.00782503223;
//   atomic_mass_map[10020000] = 2.01410177812;
//   atomic_mass_map[10030000] = 3.01604927791;
//   atomic_mass_map[20030000] = 3.01602932008;
//   atomic_mass_map[30030000] = 3.030775;
//   atomic_mass_map[10040000] = 4.026431864;
//   atomic_mass_map[20040000] = 4.00260325413;
//   atomic_mass_map[30040000] = 4.027185559;
//   atomic_mass_map[10050000] = 5.035311489;
//   atomic_mass_map[20050000] = 5.012057224;
//   atomic_mass_map[30050000] = 5.0125378;
//   atomic_mass_map[40050000] = 5.03987;
//   atomic_mass_map[10060000] = 6.044955433;
//   atomic_mass_map[20060000] = 6.018885891;
//   atomic_mass_map[30060000] = 6.01512288742;
//   atomic_mass_map[40060000] = 6.019726411;
//   atomic_mass_map[50060000] = 6.0508;
//   atomic_mass_map[10070000] = 7.052749;
//   atomic_mass_map[20070000] = 7.027990654;
//   atomic_mass_map[30070000] = 7.01600343659;
//   atomic_mass_map[40070000] = 7.016928717;
//   atomic_mass_map[50070000] = 7.029712;
//   atomic_mass_map[20080000] = 8.03393439;
//   atomic_mass_map[30080000] = 8.022486246;
//   atomic_mass_map[40080000] = 8.005305102;
//   atomic_mass_map[50080000] = 8.024607326;
//   atomic_mass_map[60080000] = 8.037643042;
//   atomic_mass_map[20090000] = 9.043946494;
//   atomic_mass_map[30090000] = 9.026790191;
//   atomic_mass_map[40090000] = 9.012183065;
//   atomic_mass_map[50090000] = 9.013329649;
//   atomic_mass_map[60090000] = 9.031037204;
//   atomic_mass_map[20100000] = 10.052788655;
//   atomic_mass_map[30100000] = 10.035483453;
//   atomic_mass_map[40100000] = 10.013534695;
//   atomic_mass_map[50100000] = 10.012936949;
//   atomic_mass_map[60100000] = 10.016853307;
//   atomic_mass_map[70100000] = 10.04165363;
//   atomic_mass_map[30110000] = 11.043723581;
//   atomic_mass_map[40110000] = 11.021661081;
//   atomic_mass_map[50110000] = 11.009305355;
//   atomic_mass_map[60110000] = 11.011433563;
//   atomic_mass_map[70110000] = 11.026091034;
//   atomic_mass_map[30120000] = 12.052517322;
//   atomic_mass_map[40120000] = 12.026922082;
//   atomic_mass_map[50120000] = 12.014352658;
//   atomic_mass_map[60120000] = 12.0;
//   atomic_mass_map[70120000] = 12.018613187;
//   atomic_mass_map[80120000] = 12.034261836;
//   atomic_mass_map[30130000] = 13.062631523;
//   atomic_mass_map[40130000] = 13.036134506;
//   atomic_mass_map[50130000] = 13.017780166;
//   atomic_mass_map[60130000] = 13.0033548351;
//   atomic_mass_map[70130000] = 13.005738609;
//   atomic_mass_map[80130000] = 13.024815446;
//   atomic_mass_map[40140000] = 14.04289292;
//   atomic_mass_map[50140000] = 14.025404012;
//   atomic_mass_map[60140000] = 14.0032419884;
//   atomic_mass_map[70140000] = 14.0030740044;
//   atomic_mass_map[80140000] = 14.008596359;
//   atomic_mass_map[90140000] = 14.034315207;
//   atomic_mass_map[40150000] = 15.05342;
//   atomic_mass_map[50150000] = 15.03108768;
//   atomic_mass_map[60150000] = 15.010599256;
//   atomic_mass_map[70150000] = 15.0001088989;
//   atomic_mass_map[80150000] = 15.003065618;
//   atomic_mass_map[90150000] = 15.018042853;
//   atomic_mass_map[40160000] = 16.061672036;
//   atomic_mass_map[50160000] = 16.039841663;
//   atomic_mass_map[60160000] = 16.014701252;
//   atomic_mass_map[70160000] = 16.006101925;
//   atomic_mass_map[80160000] = 15.9949146196;
//   atomic_mass_map[90160000] = 16.011465725;
//   atomic_mass_map[100160000] = 16.025750197;
//   atomic_mass_map[50170000] = 17.046989906;
//   atomic_mass_map[60170000] = 17.02257747;
//   atomic_mass_map[70170000] = 17.008448873;
//   atomic_mass_map[80170000] = 16.9991317565;
//   atomic_mass_map[90170000] = 17.002095237;
//   atomic_mass_map[100170000] = 17.017713963;
//   atomic_mass_map[50180000] = 18.055660189;
//   atomic_mass_map[60180000] = 18.026750708;
//   atomic_mass_map[70180000] = 18.014077565;
//   atomic_mass_map[80180000] = 17.9991596129;
//   atomic_mass_map[90180000] = 18.000937325;
//   atomic_mass_map[100180000] = 18.005708703;
//   atomic_mass_map[110180000] = 18.026878252;
//   atomic_mass_map[50190000] = 19.0631;
//   atomic_mass_map[60190000] = 19.034796372;
//   atomic_mass_map[70190000] = 19.017021603;
//   atomic_mass_map[80190000] = 19.00357797;
//   atomic_mass_map[90190000] = 18.9984031627;
//   atomic_mass_map[100190000] = 19.001880907;
//   atomic_mass_map[110190000] = 19.013880272;
//   atomic_mass_map[120190000] = 19.034169186;
//   atomic_mass_map[50200000] = 20.07207;
//   atomic_mass_map[60200000] = 20.040319754;
//   atomic_mass_map[70200000] = 20.023365807;
//   atomic_mass_map[80200000] = 20.004075354;
//   atomic_mass_map[90200000] = 19.999981252;
//   atomic_mass_map[100200000] = 19.9924401762;
//   atomic_mass_map[110200000] = 20.007354426;
//   atomic_mass_map[120200000] = 20.018850004;
//   atomic_mass_map[50210000] = 21.08129;
//   atomic_mass_map[60210000] = 21.049;
//   atomic_mass_map[70210000] = 21.02710824;
//   atomic_mass_map[80210000] = 21.008654813;
//   atomic_mass_map[90210000] = 20.999948889;
//   atomic_mass_map[100210000] = 20.993846685;
//   atomic_mass_map[110210000] = 20.997654695;
//   atomic_mass_map[120210000] = 21.01171614;
//   atomic_mass_map[130210000] = 21.028975;
//   atomic_mass_map[60220000] = 22.057531496;
//   atomic_mass_map[70220000] = 22.034394934;
//   atomic_mass_map[80220000] = 22.009966057;
//   atomic_mass_map[90220000] = 22.002998813;
//   atomic_mass_map[100220000] = 21.991385114;
//   atomic_mass_map[110220000] = 21.994437411;
//   atomic_mass_map[120220000] = 21.999570648;
//   atomic_mass_map[130220000] = 22.01954;
//   atomic_mass_map[140220000] = 22.03579;
//   atomic_mass_map[60230000] = 23.06889;
//   atomic_mass_map[70230000] = 23.04114;
//   atomic_mass_map[80230000] = 23.015695922;
//   atomic_mass_map[90230000] = 23.003556696;
//   atomic_mass_map[100230000] = 22.994466905;
//   atomic_mass_map[110230000] = 22.989769282;
//   atomic_mass_map[120230000] = 22.994124208;
//   atomic_mass_map[130230000] = 23.007244351;
//   atomic_mass_map[140230000] = 23.02544;
//   atomic_mass_map[70240000] = 24.05039;
//   atomic_mass_map[80240000] = 24.019861;
//   atomic_mass_map[90240000] = 24.008115485;
//   atomic_mass_map[100240000] = 23.99361065;
//   atomic_mass_map[110240000] = 23.99096295;
//   atomic_mass_map[120240000] = 23.985041697;
//   atomic_mass_map[130240000] = 23.999948883;
//   atomic_mass_map[140240000] = 24.011534538;
//   atomic_mass_map[150240000] = 24.03577;
//   atomic_mass_map[70250000] = 25.0601;
//   atomic_mass_map[80250000] = 25.029358986;
//   atomic_mass_map[90250000] = 25.012199229;
//   atomic_mass_map[100250000] = 24.997788707;
//   atomic_mass_map[110250000] = 24.989953969;
//   atomic_mass_map[120250000] = 24.985836976;
//   atomic_mass_map[130250000] = 24.990428102;
//   atomic_mass_map[140250000] = 25.004108808;
//   atomic_mass_map[150250000] = 25.02119;
//   atomic_mass_map[80260000] = 26.03728745;
//   atomic_mass_map[90260000] = 26.020037768;
//   atomic_mass_map[100260000] = 26.000514705;
//   atomic_mass_map[110260000] = 25.992634649;
//   atomic_mass_map[120260000] = 25.982592968;
//   atomic_mass_map[130260000] = 25.986891904;
//   atomic_mass_map[140260000] = 25.992333845;
//   atomic_mass_map[150260000] = 26.01178;
//   atomic_mass_map[160260000] = 26.02907;
//   atomic_mass_map[80270000] = 27.04772;
//   atomic_mass_map[90270000] = 27.026441;
//   atomic_mass_map[100270000] = 27.007553268;
//   atomic_mass_map[110270000] = 26.994076531;
//   atomic_mass_map[120270000] = 26.984340624;
//   atomic_mass_map[130270000] = 26.981538531;
//   atomic_mass_map[140270000] = 26.986704811;
//   atomic_mass_map[150270000] = 26.999224406;
//   atomic_mass_map[160270000] = 27.01828;
//   atomic_mass_map[80280000] = 28.05591;
//   atomic_mass_map[90280000] = 28.035342095;
//   atomic_mass_map[100280000] = 28.012121998;
//   atomic_mass_map[110280000] = 27.998939;
//   atomic_mass_map[120280000] = 27.983876728;
//   atomic_mass_map[130280000] = 27.98191021;
//   atomic_mass_map[140280000] = 27.9769265347;
//   atomic_mass_map[150280000] = 27.992326585;
//   atomic_mass_map[160280000] = 28.004372766;
//   atomic_mass_map[170280000] = 28.02954;
//   atomic_mass_map[90290000] = 29.04254;
//   atomic_mass_map[100290000] = 29.019753;
//   atomic_mass_map[110290000] = 29.002877073;
//   atomic_mass_map[120290000] = 28.988617393;
//   atomic_mass_map[130290000] = 28.98045649;
//   atomic_mass_map[140290000] = 28.9764946649;
//   atomic_mass_map[150290000] = 28.981800794;
//   atomic_mass_map[160290000] = 28.996611456;
//   atomic_mass_map[170290000] = 29.01478;
//   atomic_mass_map[90300000] = 30.05165;
//   atomic_mass_map[100300000] = 30.024734;
//   atomic_mass_map[110300000] = 30.009097932;
//   atomic_mass_map[120300000] = 29.990462926;
//   atomic_mass_map[130300000] = 29.98296022;
//   atomic_mass_map[140300000] = 29.973770136;
//   atomic_mass_map[150300000] = 29.978313753;
//   atomic_mass_map[160300000] = 29.984907033;
//   atomic_mass_map[170300000] = 30.00477;
//   atomic_mass_map[180300000] = 30.02307;
//   atomic_mass_map[90310000] = 31.059709;
//   atomic_mass_map[100310000] = 31.033087;
//   atomic_mass_map[110310000] = 31.013162656;
//   atomic_mass_map[120310000] = 30.996648032;
//   atomic_mass_map[130310000] = 30.983945171;
//   atomic_mass_map[140310000] = 30.975363194;
//   atomic_mass_map[150310000] = 30.9737619984;
//   atomic_mass_map[160310000] = 30.979557007;
//   atomic_mass_map[170310000] = 30.992414203;
//   atomic_mass_map[180310000] = 31.012124;
//   atomic_mass_map[100320000] = 32.03972;
//   atomic_mass_map[110320000] = 32.020193;
//   atomic_mass_map[120320000] = 31.999110239;
//   atomic_mass_map[130320000] = 31.988085239;
//   atomic_mass_map[140320000] = 31.974151539;
//   atomic_mass_map[150320000] = 31.973907643;
//   atomic_mass_map[160320000] = 31.9720711744;
//   atomic_mass_map[170320000] = 31.985684637;
//   atomic_mass_map[180320000] = 31.997637826;
//   atomic_mass_map[190320000] = 32.02265;
//   atomic_mass_map[100330000] = 33.04938;
//   atomic_mass_map[110330000] = 33.02573;
//   atomic_mass_map[120330000] = 33.005327145;
//   atomic_mass_map[130330000] = 32.990908977;
//   atomic_mass_map[140330000] = 32.977976964;
//   atomic_mass_map[150330000] = 32.971725694;
//   atomic_mass_map[160330000] = 32.9714589098;
//   atomic_mass_map[170330000] = 32.977451989;
//   atomic_mass_map[180330000] = 32.989925546;
//   atomic_mass_map[190330000] = 33.00756;
//   atomic_mass_map[100340000] = 34.056728;
//   atomic_mass_map[110340000] = 34.03359;
//   atomic_mass_map[120340000] = 34.008935481;
//   atomic_mass_map[130340000] = 33.996705398;
//   atomic_mass_map[140340000] = 33.978575569;
//   atomic_mass_map[150340000] = 33.973645886;
//   atomic_mass_map[160340000] = 33.967867004;
//   atomic_mass_map[170340000] = 33.973762485;
//   atomic_mass_map[180340000] = 33.98027009;
//   atomic_mass_map[190340000] = 33.99869;
//   atomic_mass_map[200340000] = 34.01487;
//   atomic_mass_map[110350000] = 35.040623;
//   atomic_mass_map[120350000] = 35.01679;
//   atomic_mass_map[130350000] = 34.999764;
//   atomic_mass_map[140350000] = 34.984583476;
//   atomic_mass_map[150350000] = 34.973314062;
//   atomic_mass_map[160350000] = 34.96903231;
//   atomic_mass_map[170350000] = 34.968852682;
//   atomic_mass_map[180350000] = 34.975257586;
//   atomic_mass_map[190350000] = 34.988005407;
//   atomic_mass_map[200350000] = 35.00514;
//   atomic_mass_map[110360000] = 36.049288;
//   atomic_mass_map[120360000] = 36.021879;
//   atomic_mass_map[130360000] = 36.006388;
//   atomic_mass_map[140360000] = 35.986695219;
//   atomic_mass_map[150360000] = 35.978259625;
//   atomic_mass_map[160360000] = 35.967080706;
//   atomic_mass_map[170360000] = 35.968306809;
//   atomic_mass_map[180360000] = 35.967545105;
//   atomic_mass_map[190360000] = 35.98130201;
//   atomic_mass_map[200360000] = 35.993074404;
//   atomic_mass_map[210360000] = 36.01648;
//   atomic_mass_map[110370000] = 37.057051;
//   atomic_mass_map[120370000] = 37.03037;
//   atomic_mass_map[130370000] = 37.010531;
//   atomic_mass_map[140370000] = 36.99292074;
//   atomic_mass_map[150370000] = 36.979606841;
//   atomic_mass_map[160370000] = 36.971125514;
//   atomic_mass_map[170370000] = 36.965902602;
//   atomic_mass_map[180370000] = 36.966776331;
//   atomic_mass_map[190370000] = 36.973375889;
//   atomic_mass_map[200370000] = 36.985897852;
//   atomic_mass_map[210370000] = 37.00374;
//   atomic_mass_map[120380000] = 38.03658;
//   atomic_mass_map[130380000] = 38.017402;
//   atomic_mass_map[140380000] = 37.995523;
//   atomic_mass_map[150380000] = 37.984251583;
//   atomic_mass_map[160380000] = 37.971163328;
//   atomic_mass_map[170380000] = 37.968010436;
//   atomic_mass_map[180380000] = 37.962732106;
//   atomic_mass_map[190380000] = 37.969081117;
//   atomic_mass_map[200380000] = 37.976319224;
//   atomic_mass_map[210380000] = 37.99512;
//   atomic_mass_map[220380000] = 38.01145;
//   atomic_mass_map[120390000] = 39.045384;
//   atomic_mass_map[130390000] = 39.02254;
//   atomic_mass_map[140390000] = 39.002491;
//   atomic_mass_map[150390000] = 38.986227208;
//   atomic_mass_map[160390000] = 38.975134197;
//   atomic_mass_map[170390000] = 38.968008176;
//   atomic_mass_map[180390000] = 38.964313038;
//   atomic_mass_map[190390000] = 38.9637064864;
//   atomic_mass_map[200390000] = 38.970710813;
//   atomic_mass_map[210390000] = 38.984784968;
//   atomic_mass_map[220390000] = 39.00236;
//   atomic_mass_map[120400000] = 40.05218;
//   atomic_mass_map[130400000] = 40.03003;
//   atomic_mass_map[140400000] = 40.005829;
//   atomic_mass_map[150400000] = 39.991331748;
//   atomic_mass_map[160400000] = 39.975482562;
//   atomic_mass_map[170400000] = 39.970415469;
//   atomic_mass_map[180400000] = 39.9623831237;
//   atomic_mass_map[190400000] = 39.963998166;
//   atomic_mass_map[200400000] = 39.962590863;
//   atomic_mass_map[210400000] = 39.977967291;
//   atomic_mass_map[220400000] = 39.990498719;
//   atomic_mass_map[230400000] = 40.01276;
//   atomic_mass_map[130410000] = 41.03638;
//   atomic_mass_map[140410000] = 41.013011;
//   atomic_mass_map[150410000] = 40.994654;
//   atomic_mass_map[160410000] = 40.979593451;
//   atomic_mass_map[170410000] = 40.970684525;
//   atomic_mass_map[180410000] = 40.96450057;
//   atomic_mass_map[190410000] = 40.9618252579;
//   atomic_mass_map[200410000] = 40.962277924;
//   atomic_mass_map[210410000] = 40.969251105;
//   atomic_mass_map[220410000] = 40.983148;
//   atomic_mass_map[230410000] = 41.00021;
//   atomic_mass_map[130420000] = 42.04384;
//   atomic_mass_map[140420000] = 42.01778;
//   atomic_mass_map[150420000] = 42.001084;
//   atomic_mass_map[160420000] = 41.9810651;
//   atomic_mass_map[170420000] = 41.973254804;
//   atomic_mass_map[180420000] = 41.963045736;
//   atomic_mass_map[190420000] = 41.962402306;
//   atomic_mass_map[200420000] = 41.95861783;
//   atomic_mass_map[210420000] = 41.965516532;
//   atomic_mass_map[220420000] = 41.973049034;
//   atomic_mass_map[230420000] = 41.99182;
//   atomic_mass_map[240420000] = 42.0067;
//   atomic_mass_map[130430000] = 43.05147;
//   atomic_mass_map[140430000] = 43.0248;
//   atomic_mass_map[150430000] = 43.005024;
//   atomic_mass_map[160430000] = 42.986907635;
//   atomic_mass_map[170430000] = 42.973888584;
//   atomic_mass_map[180430000] = 42.965636055;
//   atomic_mass_map[190430000] = 42.960734702;
//   atomic_mass_map[200430000] = 42.958766438;
//   atomic_mass_map[210430000] = 42.961150474;
//   atomic_mass_map[220430000] = 42.96852252;
//   atomic_mass_map[230430000] = 42.980766;
//   atomic_mass_map[240430000] = 42.99753;
//   atomic_mass_map[140440000] = 44.03061;
//   atomic_mass_map[150440000] = 44.01121;
//   atomic_mass_map[160440000] = 43.990118848;
//   atomic_mass_map[170440000] = 43.977874543;
//   atomic_mass_map[180440000] = 43.964923815;
//   atomic_mass_map[190440000] = 43.961586985;
//   atomic_mass_map[200440000] = 43.955481561;
//   atomic_mass_map[210440000] = 43.959402875;
//   atomic_mass_map[220440000] = 43.959689949;
//   atomic_mass_map[230440000] = 43.97411;
//   atomic_mass_map[240440000] = 43.98536;
//   atomic_mass_map[250440000] = 44.00715;
//   atomic_mass_map[140450000] = 45.03995;
//   atomic_mass_map[150450000] = 45.01645;
//   atomic_mass_map[160450000] = 44.995717;
//   atomic_mass_map[170450000] = 44.98029;
//   atomic_mass_map[180450000] = 44.968039733;
//   atomic_mass_map[190450000] = 44.960691493;
//   atomic_mass_map[200450000] = 44.95618635;
//   atomic_mass_map[210450000] = 44.955908275;
//   atomic_mass_map[220450000] = 44.958121983;
//   atomic_mass_map[230450000] = 44.96577482;
//   atomic_mass_map[240450000] = 44.97905;
//   atomic_mass_map[250450000] = 44.99449;
//   atomic_mass_map[260450000] = 45.014419;
//   atomic_mass_map[150460000] = 46.02446;
//   atomic_mass_map[160460000] = 46.00004;
//   atomic_mass_map[170460000] = 45.985174;
//   atomic_mass_map[180460000] = 45.968082712;
//   atomic_mass_map[190460000] = 45.961981586;
//   atomic_mass_map[200460000] = 45.953689023;
//   atomic_mass_map[210460000] = 45.955168257;
//   atomic_mass_map[220460000] = 45.952627718;
//   atomic_mass_map[230460000] = 45.960198775;
//   atomic_mass_map[240460000] = 45.968358861;
//   atomic_mass_map[250460000] = 45.98609;
//   atomic_mass_map[260460000] = 46.00063;
//   atomic_mass_map[150470000] = 47.03139;
//   atomic_mass_map[160470000] = 47.00795;
//   atomic_mass_map[170470000] = 46.98916;
//   atomic_mass_map[180470000] = 46.972934865;
//   atomic_mass_map[190470000] = 46.961661614;
//   atomic_mass_map[200470000] = 46.95454243;
//   atomic_mass_map[210470000] = 46.95240374;
//   atomic_mass_map[220470000] = 46.951758787;
//   atomic_mass_map[230470000] = 46.954904914;
//   atomic_mass_map[240470000] = 46.962897359;
//   atomic_mass_map[250470000] = 46.975775;
//   atomic_mass_map[260470000] = 46.99185;
//   atomic_mass_map[270470000] = 47.01057;
//   atomic_mass_map[160480000] = 48.0137;
//   atomic_mass_map[170480000] = 47.99564;
//   atomic_mass_map[180480000] = 47.97591;
//   atomic_mass_map[190480000] = 47.965341186;
//   atomic_mass_map[200480000] = 47.952522765;
//   atomic_mass_map[210480000] = 47.952223611;
//   atomic_mass_map[220480000] = 47.947941979;
//   atomic_mass_map[230480000] = 47.952252223;
//   atomic_mass_map[240480000] = 47.954029061;
//   atomic_mass_map[250480000] = 47.96852;
//   atomic_mass_map[260480000] = 47.98023;
//   atomic_mass_map[270480000] = 48.00093;
//   atomic_mass_map[280480000] = 48.017688;
//   atomic_mass_map[160490000] = 49.022762;
//   atomic_mass_map[170490000] = 49.00123;
//   atomic_mass_map[180490000] = 48.9819;
//   atomic_mass_map[190490000] = 48.968210755;
//   atomic_mass_map[200490000] = 48.955662736;
//   atomic_mass_map[210490000] = 48.950014629;
//   atomic_mass_map[220490000] = 48.947865676;
//   atomic_mass_map[230490000] = 48.948511795;
//   atomic_mass_map[240490000] = 48.951333349;
//   atomic_mass_map[250490000] = 48.959595297;
//   atomic_mass_map[260490000] = 48.973429;
//   atomic_mass_map[270490000] = 48.98891;
//   atomic_mass_map[280490000] = 49.0077;
//   atomic_mass_map[170500000] = 50.00905;
//   atomic_mass_map[180500000] = 49.98613;
//   atomic_mass_map[190500000] = 49.972380017;
//   atomic_mass_map[200500000] = 49.957499217;
//   atomic_mass_map[210500000] = 49.952176274;
//   atomic_mass_map[220500000] = 49.944786889;
//   atomic_mass_map[230500000] = 49.947156014;
//   atomic_mass_map[240500000] = 49.946041833;
//   atomic_mass_map[250500000] = 49.954237781;
//   atomic_mass_map[260500000] = 49.962974929;
//   atomic_mass_map[270500000] = 49.98091;
//   atomic_mass_map[280500000] = 49.99474;
//   atomic_mass_map[170510000] = 51.01554;
//   atomic_mass_map[180510000] = 50.9937;
//   atomic_mass_map[190510000] = 50.975828036;
//   atomic_mass_map[200510000] = 50.960988981;
//   atomic_mass_map[210510000] = 50.953591956;
//   atomic_mass_map[220510000] = 50.946610651;
//   atomic_mass_map[230510000] = 50.943957036;
//   atomic_mass_map[240510000] = 50.944765018;
//   atomic_mass_map[250510000] = 50.948208475;
//   atomic_mass_map[260510000] = 50.956841021;
//   atomic_mass_map[270510000] = 50.970647;
//   atomic_mass_map[280510000] = 50.98611;
//   atomic_mass_map[180520000] = 51.99896;
//   atomic_mass_map[190520000] = 51.98224;
//   atomic_mass_map[200520000] = 51.963217145;
//   atomic_mass_map[210520000] = 51.956879479;
//   atomic_mass_map[220520000] = 51.946893006;
//   atomic_mass_map[230520000] = 51.944773007;
//   atomic_mass_map[240520000] = 51.940506231;
//   atomic_mass_map[250520000] = 51.945563949;
//   atomic_mass_map[260520000] = 51.948113146;
//   atomic_mass_map[270520000] = 51.96351;
//   atomic_mass_map[280520000] = 51.9748;
//   atomic_mass_map[290520000] = 51.99671;
//   atomic_mass_map[180530000] = 53.00729;
//   atomic_mass_map[190530000] = 52.98746;
//   atomic_mass_map[200530000] = 52.96945;
//   atomic_mass_map[210530000] = 52.95909;
//   atomic_mass_map[220530000] = 52.949725933;
//   atomic_mass_map[230530000] = 52.944336741;
//   atomic_mass_map[240530000] = 52.940648147;
//   atomic_mass_map[250530000] = 52.941288891;
//   atomic_mass_map[260530000] = 52.945306428;
//   atomic_mass_map[270530000] = 52.954204077;
//   atomic_mass_map[280530000] = 52.96819;
//   atomic_mass_map[290530000] = 52.98459;
//   atomic_mass_map[190540000] = 53.99463;
//   atomic_mass_map[200540000] = 53.9734;
//   atomic_mass_map[210540000] = 53.96393;
//   atomic_mass_map[220540000] = 53.951049547;
//   atomic_mass_map[230540000] = 53.946438617;
//   atomic_mass_map[240540000] = 53.938879158;
//   atomic_mass_map[250540000] = 53.940357615;
//   atomic_mass_map[260540000] = 53.939608986;
//   atomic_mass_map[270540000] = 53.948459872;
//   atomic_mass_map[280540000] = 53.957892463;
//   atomic_mass_map[290540000] = 53.97666;
//   atomic_mass_map[300540000] = 53.992039;
//   atomic_mass_map[190550000] = 55.00076;
//   atomic_mass_map[200550000] = 54.9803;
//   atomic_mass_map[210550000] = 54.967818892;
//   atomic_mass_map[220550000] = 54.955268064;
//   atomic_mass_map[230550000] = 54.94724215;
//   atomic_mass_map[240550000] = 54.940838434;
//   atomic_mass_map[250550000] = 54.93804391;
//   atomic_mass_map[260550000] = 54.938291994;
//   atomic_mass_map[270550000] = 54.941997202;
//   atomic_mass_map[280550000] = 54.951330632;
//   atomic_mass_map[290550000] = 54.966038;
//   atomic_mass_map[300550000] = 54.98398;
//   atomic_mass_map[190560000] = 56.00851;
//   atomic_mass_map[200560000] = 55.98508;
//   atomic_mass_map[210560000] = 55.97345;
//   atomic_mass_map[220560000] = 55.957911353;
//   atomic_mass_map[230560000] = 55.95048481;
//   atomic_mass_map[240560000] = 55.940653139;
//   atomic_mass_map[250560000] = 55.938903691;
//   atomic_mass_map[260560000] = 55.934936326;
//   atomic_mass_map[270560000] = 55.939838798;
//   atomic_mass_map[280560000] = 55.942128549;
//   atomic_mass_map[290560000] = 55.95895;
//   atomic_mass_map[300560000] = 55.97254;
//   atomic_mass_map[310560000] = 55.99536;
//   atomic_mass_map[200570000] = 56.99262;
//   atomic_mass_map[210570000] = 56.97777;
//   atomic_mass_map[220570000] = 56.963641626;
//   atomic_mass_map[230570000] = 56.952518869;
//   atomic_mass_map[240570000] = 56.943613013;
//   atomic_mass_map[250570000] = 56.938286096;
//   atomic_mass_map[260570000] = 56.935392841;
//   atomic_mass_map[270570000] = 56.936290574;
//   atomic_mass_map[280570000] = 56.939792184;
//   atomic_mass_map[290570000] = 56.949212498;
//   atomic_mass_map[300570000] = 56.96506;
//   atomic_mass_map[310570000] = 56.9832;
//   atomic_mass_map[200580000] = 57.99794;
//   atomic_mass_map[210580000] = 57.98403;
//   atomic_mass_map[220580000] = 57.9666;
//   atomic_mass_map[230580000] = 57.956715235;
//   atomic_mass_map[240580000] = 57.944353129;
//   atomic_mass_map[250580000] = 57.940066645;
//   atomic_mass_map[260580000] = 57.933274431;
//   atomic_mass_map[270580000] = 57.935752073;
//   atomic_mass_map[280580000] = 57.935342414;
//   atomic_mass_map[290580000] = 57.944533046;
//   atomic_mass_map[300580000] = 57.954591062;
//   atomic_mass_map[310580000] = 57.97478;
//   atomic_mass_map[320580000] = 57.99172;
//   atomic_mass_map[210590000] = 58.98894;
//   atomic_mass_map[220590000] = 58.97247;
//   atomic_mass_map[230590000] = 58.959385659;
//   atomic_mass_map[240590000] = 58.948586367;
//   atomic_mass_map[250590000] = 58.940391113;
//   atomic_mass_map[260590000] = 58.934874338;
//   atomic_mass_map[270590000] = 58.933194288;
//   atomic_mass_map[280590000] = 58.934346202;
//   atomic_mass_map[290590000] = 58.939497482;
//   atomic_mass_map[300590000] = 58.949312657;
//   atomic_mass_map[310590000] = 58.96353;
//   atomic_mass_map[320590000] = 58.98249;
//   atomic_mass_map[210600000] = 59.99565;
//   atomic_mass_map[220600000] = 59.97603;
//   atomic_mass_map[230600000] = 59.96431329;
//   atomic_mass_map[240600000] = 59.950076033;
//   atomic_mass_map[250600000] = 59.943136576;
//   atomic_mass_map[260600000] = 59.9340711;
//   atomic_mass_map[270600000] = 59.933816299;
//   atomic_mass_map[280600000] = 59.930785885;
//   atomic_mass_map[290600000] = 59.937364544;
//   atomic_mass_map[300600000] = 59.941842103;
//   atomic_mass_map[310600000] = 59.95729;
//   atomic_mass_map[320600000] = 59.97036;
//   atomic_mass_map[330600000] = 59.99388;
//   atomic_mass_map[210610000] = 61.001;
//   atomic_mass_map[220610000] = 60.98245;
//   atomic_mass_map[230610000] = 60.96725;
//   atomic_mass_map[240610000] = 60.954422585;
//   atomic_mass_map[250610000] = 60.944452544;
//   atomic_mass_map[260610000] = 60.936746244;
//   atomic_mass_map[270610000] = 60.93247662;
//   atomic_mass_map[280610000] = 60.93105557;
//   atomic_mass_map[290610000] = 60.933457612;
//   atomic_mass_map[300610000] = 60.939507188;
//   atomic_mass_map[310610000] = 60.949398984;
//   atomic_mass_map[320610000] = 60.96379;
//   atomic_mass_map[330610000] = 60.98112;
//   atomic_mass_map[220620000] = 61.98651;
//   atomic_mass_map[230620000] = 61.97265;
//   atomic_mass_map[240620000] = 61.956097451;
//   atomic_mass_map[250620000] = 61.947954;
//   atomic_mass_map[260620000] = 61.936791812;
//   atomic_mass_map[270620000] = 61.934058809;
//   atomic_mass_map[280620000] = 61.928345365;
//   atomic_mass_map[290620000] = 61.932595415;
//   atomic_mass_map[300620000] = 61.934333971;
//   atomic_mass_map[310620000] = 61.944190251;
//   atomic_mass_map[320620000] = 61.95502;
//   atomic_mass_map[330620000] = 61.97361;
//   atomic_mass_map[220630000] = 62.99375;
//   atomic_mass_map[230630000] = 62.97639;
//   atomic_mass_map[240630000] = 62.961650731;
//   atomic_mass_map[250630000] = 62.949664675;
//   atomic_mass_map[260630000] = 62.9402727;
//   atomic_mass_map[270630000] = 62.933600217;
//   atomic_mass_map[280630000] = 62.929669626;
//   atomic_mass_map[290630000] = 62.929597723;
//   atomic_mass_map[300630000] = 62.933211474;
//   atomic_mass_map[310630000] = 62.939294196;
//   atomic_mass_map[320630000] = 62.949628;
//   atomic_mass_map[330630000] = 62.9639;
//   atomic_mass_map[230640000] = 63.98264;
//   atomic_mass_map[240640000] = 63.96408;
//   atomic_mass_map[250640000] = 63.95384937;
//   atomic_mass_map[260640000] = 63.940987763;
//   atomic_mass_map[270640000] = 63.935810764;
//   atomic_mass_map[280640000] = 63.927966816;
//   atomic_mass_map[290640000] = 63.929764342;
//   atomic_mass_map[300640000] = 63.929142013;
//   atomic_mass_map[310640000] = 63.93684044;
//   atomic_mass_map[320640000] = 63.941689913;
//   atomic_mass_map[330640000] = 63.957432;
//   atomic_mass_map[340640000] = 63.97109;
//   atomic_mass_map[230650000] = 64.9875;
//   atomic_mass_map[240650000] = 64.96996;
//   atomic_mass_map[250650000] = 64.95601975;
//   atomic_mass_map[260650000] = 64.945011462;
//   atomic_mass_map[270650000] = 64.936462072;
//   atomic_mass_map[280650000] = 64.930085173;
//   atomic_mass_map[290650000] = 64.927789704;
//   atomic_mass_map[300650000] = 64.92924077;
//   atomic_mass_map[310650000] = 64.932734588;
//   atomic_mass_map[320650000] = 64.939368137;
//   atomic_mass_map[330650000] = 64.949611;
//   atomic_mass_map[340650000] = 64.9644;
//   atomic_mass_map[230660000] = 65.99398;
//   atomic_mass_map[240660000] = 65.97366;
//   atomic_mass_map[250660000] = 65.960546834;
//   atomic_mass_map[260660000] = 65.94624996;
//   atomic_mass_map[270660000] = 65.939442945;
//   atomic_mass_map[280660000] = 65.929139334;
//   atomic_mass_map[290660000] = 65.928869032;
//   atomic_mass_map[300660000] = 65.926033809;
//   atomic_mass_map[310660000] = 65.931589401;
//   atomic_mass_map[320660000] = 65.933862126;
//   atomic_mass_map[330660000] = 65.944148779;
//   atomic_mass_map[340660000] = 65.95559;
//   atomic_mass_map[240670000] = 66.98016;
//   atomic_mass_map[250670000] = 66.96424;
//   atomic_mass_map[260670000] = 66.950543395;
//   atomic_mass_map[270670000] = 66.940609627;
//   atomic_mass_map[280670000] = 66.931569414;
//   atomic_mass_map[290670000] = 66.927730314;
//   atomic_mass_map[300670000] = 66.927127746;
//   atomic_mass_map[310670000] = 66.928202547;
//   atomic_mass_map[320670000] = 66.932733861;
//   atomic_mass_map[330670000] = 66.939251112;
//   atomic_mass_map[340670000] = 66.949994;
//   atomic_mass_map[350670000] = 66.96465;
//   atomic_mass_map[240680000] = 67.98403;
//   atomic_mass_map[250680000] = 67.96962;
//   atomic_mass_map[260680000] = 67.95295155;
//   atomic_mass_map[270680000] = 67.944257589;
//   atomic_mass_map[280680000] = 67.931868789;
//   atomic_mass_map[290680000] = 67.929610889;
//   atomic_mass_map[300680000] = 67.924844554;
//   atomic_mass_map[310680000] = 67.927980485;
//   atomic_mass_map[320680000] = 67.928095307;
//   atomic_mass_map[330680000] = 67.93677413;
//   atomic_mass_map[340680000] = 67.941825238;
//   atomic_mass_map[350680000] = 67.958732;
//   atomic_mass_map[250690000] = 68.97366;
//   atomic_mass_map[260690000] = 68.95807;
//   atomic_mass_map[270690000] = 68.946141268;
//   atomic_mass_map[280690000] = 68.935610269;
//   atomic_mass_map[290690000] = 68.929429269;
//   atomic_mass_map[300690000] = 68.926550682;
//   atomic_mass_map[310690000] = 68.925573541;
//   atomic_mass_map[320690000] = 68.927964481;
//   atomic_mass_map[330690000] = 68.932246302;
//   atomic_mass_map[340690000] = 68.939414847;
//   atomic_mass_map[350690000] = 68.950497297;
//   atomic_mass_map[360690000] = 68.96518;
//   atomic_mass_map[250700000] = 69.97937;
//   atomic_mass_map[260700000] = 69.96102;
//   atomic_mass_map[270700000] = 69.94963;
//   atomic_mass_map[280700000] = 69.936431302;
//   atomic_mass_map[290700000] = 69.932392079;
//   atomic_mass_map[300700000] = 69.925319208;
//   atomic_mass_map[310700000] = 69.926021927;
//   atomic_mass_map[320700000] = 69.92424875;
//   atomic_mass_map[330700000] = 69.930926195;
//   atomic_mass_map[340700000] = 69.933515523;
//   atomic_mass_map[350700000] = 69.944792323;
//   atomic_mass_map[360700000] = 69.95604;
//   atomic_mass_map[250710000] = 70.98368;
//   atomic_mass_map[260710000] = 70.96672;
//   atomic_mass_map[270710000] = 70.952366923;
//   atomic_mass_map[280710000] = 70.940518964;
//   atomic_mass_map[290710000] = 70.932676833;
//   atomic_mass_map[300710000] = 70.927719582;
//   atomic_mass_map[310710000] = 70.924702577;
//   atomic_mass_map[320710000] = 70.924952327;
//   atomic_mass_map[330710000] = 70.927113801;
//   atomic_mass_map[340710000] = 70.932209433;
//   atomic_mass_map[350710000] = 70.939342155;
//   atomic_mass_map[360710000] = 70.950265696;
//   atomic_mass_map[370710000] = 70.96532;
//   atomic_mass_map[260720000] = 71.96983;
//   atomic_mass_map[270720000] = 71.95729;
//   atomic_mass_map[280720000] = 71.941785926;
//   atomic_mass_map[290720000] = 71.935820307;
//   atomic_mass_map[300720000] = 71.926842807;
//   atomic_mass_map[310720000] = 71.926367474;
//   atomic_mass_map[320720000] = 71.922075826;
//   atomic_mass_map[330720000] = 71.926752294;
//   atomic_mass_map[340720000] = 71.927140507;
//   atomic_mass_map[350720000] = 71.936588619;
//   atomic_mass_map[360720000] = 71.942092407;
//   atomic_mass_map[370720000] = 71.95908;
//   atomic_mass_map[260730000] = 72.97572;
//   atomic_mass_map[270730000] = 72.96039;
//   atomic_mass_map[280730000] = 72.946206683;
//   atomic_mass_map[290730000] = 72.936674378;
//   atomic_mass_map[300730000] = 72.929582582;
//   atomic_mass_map[310730000] = 72.925174682;
//   atomic_mass_map[320730000] = 72.923458956;
//   atomic_mass_map[330730000] = 72.923829086;
//   atomic_mass_map[340730000] = 72.926754882;
//   atomic_mass_map[350730000] = 72.931671503;
//   atomic_mass_map[360730000] = 72.939289195;
//   atomic_mass_map[370730000] = 72.950529;
//   atomic_mass_map[380730000] = 72.9657;
//   atomic_mass_map[260740000] = 73.97935;
//   atomic_mass_map[270740000] = 73.96515;
//   atomic_mass_map[280740000] = 73.94798;
//   atomic_mass_map[290740000] = 73.939874862;
//   atomic_mass_map[300740000] = 73.929407262;
//   atomic_mass_map[310740000] = 73.926945727;
//   atomic_mass_map[320740000] = 73.921177761;
//   atomic_mass_map[330740000] = 73.923928598;
//   atomic_mass_map[340740000] = 73.922475934;
//   atomic_mass_map[350740000] = 73.929910177;
//   atomic_mass_map[360740000] = 73.93308402;
//   atomic_mass_map[370740000] = 73.944265894;
//   atomic_mass_map[380740000] = 73.95617;
//   atomic_mass_map[270750000] = 74.96876;
//   atomic_mass_map[280750000] = 74.9525;
//   atomic_mass_map[290750000] = 74.941522606;
//   atomic_mass_map[300750000] = 74.932840246;
//   atomic_mass_map[310750000] = 74.926500246;
//   atomic_mass_map[320750000] = 74.92285837;
//   atomic_mass_map[330750000] = 74.921594567;
//   atomic_mass_map[340750000] = 74.92252287;
//   atomic_mass_map[350750000] = 74.925810452;
//   atomic_mass_map[360750000] = 74.930945746;
//   atomic_mass_map[370750000] = 74.938573201;
//   atomic_mass_map[380750000] = 74.94995277;
//   atomic_mass_map[270760000] = 75.97413;
//   atomic_mass_map[280760000] = 75.95533;
//   atomic_mass_map[290760000] = 75.945275026;
//   atomic_mass_map[300760000] = 75.933114957;
//   atomic_mass_map[310760000] = 75.928827626;
//   atomic_mass_map[320760000] = 75.921402726;
//   atomic_mass_map[330760000] = 75.922392015;
//   atomic_mass_map[340760000] = 75.919213704;
//   atomic_mass_map[350760000] = 75.924541576;
//   atomic_mass_map[360760000] = 75.925910271;
//   atomic_mass_map[370760000] = 75.935073032;
//   atomic_mass_map[380760000] = 75.941762761;
//   atomic_mass_map[390760000] = 75.95856;
//   atomic_mass_map[280770000] = 76.96055;
//   atomic_mass_map[290770000] = 76.94792;
//   atomic_mass_map[300770000] = 76.936887199;
//   atomic_mass_map[310770000] = 76.9291543;
//   atomic_mass_map[320770000] = 76.923549843;
//   atomic_mass_map[330770000] = 76.920647563;
//   atomic_mass_map[340770000] = 76.919914154;
//   atomic_mass_map[350770000] = 76.921379198;
//   atomic_mass_map[360770000] = 76.92467;
//   atomic_mass_map[370770000] = 76.9304016;
//   atomic_mass_map[380770000] = 76.937945455;
//   atomic_mass_map[390770000] = 76.949781;
//   atomic_mass_map[280780000] = 77.96336;
//   atomic_mass_map[290780000] = 77.95223;
//   atomic_mass_map[300780000] = 77.938289206;
//   atomic_mass_map[310780000] = 77.931608845;
//   atomic_mass_map[320780000] = 77.922852908;
//   atomic_mass_map[330780000] = 77.921827773;
//   atomic_mass_map[340780000] = 77.91730928;
//   atomic_mass_map[350780000] = 77.921145895;
//   atomic_mass_map[360780000] = 77.920364944;
//   atomic_mass_map[370780000] = 77.928141868;
//   atomic_mass_map[380780000] = 77.93217998;
//   atomic_mass_map[390780000] = 77.94361;
//   atomic_mass_map[400780000] = 77.95566;
//   atomic_mass_map[280790000] = 78.970252;
//   atomic_mass_map[290790000] = 78.95502;
//   atomic_mass_map[300790000] = 78.942638069;
//   atomic_mass_map[310790000] = 78.932852301;
//   atomic_mass_map[320790000] = 78.925360127;
//   atomic_mass_map[330790000] = 78.920948423;
//   atomic_mass_map[340790000] = 78.918499287;
//   atomic_mass_map[350790000] = 78.918337579;
//   atomic_mass_map[360790000] = 78.920082923;
//   atomic_mass_map[370790000] = 78.923989864;
//   atomic_mass_map[380790000] = 78.929707664;
//   atomic_mass_map[390790000] = 78.937351299;
//   atomic_mass_map[400790000] = 78.94948;
//   atomic_mass_map[290800000] = 79.96089;
//   atomic_mass_map[300800000] = 79.944552931;
//   atomic_mass_map[310800000] = 79.936420775;
//   atomic_mass_map[320800000] = 79.925350775;
//   atomic_mass_map[330800000] = 79.922474584;
//   atomic_mass_map[340800000] = 79.916521762;
//   atomic_mass_map[350800000] = 79.918529788;
//   atomic_mass_map[360800000] = 79.916378084;
//   atomic_mass_map[370800000] = 79.922516444;
//   atomic_mass_map[380800000] = 79.924517516;
//   atomic_mass_map[390800000] = 79.934356096;
//   atomic_mass_map[400800000] = 79.9404;
//   atomic_mass_map[290810000] = 80.965872;
//   atomic_mass_map[300810000] = 80.950402619;
//   atomic_mass_map[310810000] = 80.938133843;
//   atomic_mass_map[320810000] = 80.928832943;
//   atomic_mass_map[330810000] = 80.92213228;
//   atomic_mass_map[340810000] = 80.917993023;
//   atomic_mass_map[350810000] = 80.91628969;
//   atomic_mass_map[360810000] = 80.916591181;
//   atomic_mass_map[370810000] = 80.918993936;
//   atomic_mass_map[380810000] = 80.923211395;
//   atomic_mass_map[390810000] = 80.929455556;
//   atomic_mass_map[400810000] = 80.937308892;
//   atomic_mass_map[410810000] = 80.9496;
//   atomic_mass_map[290820000] = 81.972442;
//   atomic_mass_map[300820000] = 81.95426;
//   atomic_mass_map[310820000] = 81.943176533;
//   atomic_mass_map[320820000] = 81.929774033;
//   atomic_mass_map[330820000] = 81.924741233;
//   atomic_mass_map[340820000] = 81.916699497;
//   atomic_mass_map[350820000] = 81.916803246;
//   atomic_mass_map[360820000] = 81.91348273;
//   atomic_mass_map[370820000] = 81.918209024;
//   atomic_mass_map[380820000] = 81.918399855;
//   atomic_mass_map[390820000] = 81.926931401;
//   atomic_mass_map[400820000] = 81.931354;
//   atomic_mass_map[410820000] = 81.94396;
//   atomic_mass_map[300830000] = 82.96056;
//   atomic_mass_map[310830000] = 82.947120301;
//   atomic_mass_map[320830000] = 82.934539101;
//   atomic_mass_map[330830000] = 82.925206901;
//   atomic_mass_map[340830000] = 82.919118569;
//   atomic_mass_map[350830000] = 82.91517562;
//   atomic_mass_map[360830000] = 82.914127164;
//   atomic_mass_map[370830000] = 82.915114183;
//   atomic_mass_map[380830000] = 82.917554374;
//   atomic_mass_map[390830000] = 82.92248525;
//   atomic_mass_map[400830000] = 82.92924215;
//   atomic_mass_map[410830000] = 82.937293731;
//   atomic_mass_map[420830000] = 82.94988;
//   atomic_mass_map[300840000] = 83.965206;
//   atomic_mass_map[310840000] = 83.95246;
//   atomic_mass_map[320840000] = 83.937575092;
//   atomic_mass_map[330840000] = 83.929303292;
//   atomic_mass_map[340840000] = 83.918466763;
//   atomic_mass_map[350840000] = 83.916496419;
//   atomic_mass_map[360840000] = 83.9114977282;
//   atomic_mass_map[370840000] = 83.914375229;
//   atomic_mass_map[380840000] = 83.913419136;
//   atomic_mass_map[390840000] = 83.920672086;
//   atomic_mass_map[400840000] = 83.923326899;
//   atomic_mass_map[410840000] = 83.934492;
//   atomic_mass_map[420840000] = 83.94149;
//   atomic_mass_map[300850000] = 84.97226;
//   atomic_mass_map[310850000] = 84.95699;
//   atomic_mass_map[320850000] = 84.94296966;
//   atomic_mass_map[330850000] = 84.93216366;
//   atomic_mass_map[340850000] = 84.92226076;
//   atomic_mass_map[350850000] = 84.91564576;
//   atomic_mass_map[360850000] = 84.912527262;
//   atomic_mass_map[370850000] = 84.9117897379;
//   atomic_mass_map[380850000] = 84.912932046;
//   atomic_mass_map[390850000] = 84.916433039;
//   atomic_mass_map[400850000] = 84.921444448;
//   atomic_mass_map[410850000] = 84.928845837;
//   atomic_mass_map[420850000] = 84.938260737;
//   atomic_mass_map[430850000] = 84.95058;
//   atomic_mass_map[310860000] = 85.96301;
//   atomic_mass_map[320860000] = 85.94658;
//   atomic_mass_map[330860000] = 85.936701533;
//   atomic_mass_map[340860000] = 85.924311733;
//   atomic_mass_map[350860000] = 85.918805433;
//   atomic_mass_map[360860000] = 85.9106106269;
//   atomic_mass_map[370860000] = 85.911167425;
//   atomic_mass_map[380860000] = 85.909260608;
//   atomic_mass_map[390860000] = 85.91488598;
//   atomic_mass_map[400860000] = 85.916297204;
//   atomic_mass_map[410860000] = 85.925782798;
//   atomic_mass_map[420860000] = 85.931174817;
//   atomic_mass_map[430860000] = 85.94493;
//   atomic_mass_map[310870000] = 86.968245;
//   atomic_mass_map[320870000] = 86.95268;
//   atomic_mass_map[330870000] = 86.940291718;
//   atomic_mass_map[340870000] = 86.928688618;
//   atomic_mass_map[350870000] = 86.920674018;
//   atomic_mass_map[360870000] = 86.91335476;
//   atomic_mass_map[370870000] = 86.909180531;
//   atomic_mass_map[380870000] = 86.908877531;
//   atomic_mass_map[390870000] = 86.910876138;
//   atomic_mass_map[400870000] = 86.914817988;
//   atomic_mass_map[410870000] = 86.920693747;
//   atomic_mass_map[420870000] = 86.928196201;
//   atomic_mass_map[430870000] = 86.938067187;
//   atomic_mass_map[440870000] = 86.95069;
//   atomic_mass_map[320880000] = 87.95691;
//   atomic_mass_map[330880000] = 87.94555;
//   atomic_mass_map[340880000] = 87.931417492;
//   atomic_mass_map[350880000] = 87.924083292;
//   atomic_mass_map[360880000] = 87.914447881;
//   atomic_mass_map[370880000] = 87.911315592;
//   atomic_mass_map[380880000] = 87.905612542;
//   atomic_mass_map[390880000] = 87.909501563;
//   atomic_mass_map[400880000] = 87.91022129;
//   atomic_mass_map[410880000] = 87.91822171;
//   atomic_mass_map[420880000] = 87.921967781;
//   atomic_mass_map[430880000] = 87.933782381;
//   atomic_mass_map[440880000] = 87.9416;
//   atomic_mass_map[320890000] = 88.96379;
//   atomic_mass_map[330890000] = 88.94976;
//   atomic_mass_map[340890000] = 88.93666906;
//   atomic_mass_map[350890000] = 88.92670456;
//   atomic_mass_map[360890000] = 88.917835451;
//   atomic_mass_map[370890000] = 88.912278298;
//   atomic_mass_map[380890000] = 88.907451095;
//   atomic_mass_map[390890000] = 88.905840348;
//   atomic_mass_map[400890000] = 88.908881441;
//   atomic_mass_map[410890000] = 88.913445073;
//   atomic_mass_map[420890000] = 88.919468151;
//   atomic_mass_map[430890000] = 88.927648651;
//   atomic_mass_map[440890000] = 88.93762;
//   atomic_mass_map[450890000] = 88.950584;
//   atomic_mass_map[320900000] = 89.96863;
//   atomic_mass_map[330900000] = 89.95563;
//   atomic_mass_map[340900000] = 89.940096;
//   atomic_mass_map[350900000] = 89.93129285;
//   atomic_mass_map[360900000] = 89.919527931;
//   atomic_mass_map[370900000] = 89.914798453;
//   atomic_mass_map[380900000] = 89.907730037;
//   atomic_mass_map[390900000] = 89.907143942;
//   atomic_mass_map[400900000] = 89.904697659;
//   atomic_mass_map[410900000] = 89.911258449;
//   atomic_mass_map[420900000] = 89.913930861;
//   atomic_mass_map[430900000] = 89.924073921;
//   atomic_mass_map[440900000] = 89.93034438;
//   atomic_mass_map[450900000] = 89.94422;
//   atomic_mass_map[330910000] = 90.96039;
//   atomic_mass_map[340910000] = 90.94596;
//   atomic_mass_map[350910000] = 90.934398619;
//   atomic_mass_map[360910000] = 90.923806311;
//   atomic_mass_map[370910000] = 90.916537169;
//   atomic_mass_map[380910000] = 90.910195442;
//   atomic_mass_map[390910000] = 90.907297442;
//   atomic_mass_map[400910000] = 90.905639587;
//   atomic_mass_map[410910000] = 90.906989658;
//   atomic_mass_map[420910000] = 90.911745312;
//   atomic_mass_map[430910000] = 90.918425397;
//   atomic_mass_map[440910000] = 90.926741859;
//   atomic_mass_map[450910000] = 90.93688;
//   atomic_mass_map[460910000] = 90.95032;
//   atomic_mass_map[330920000] = 91.96674;
//   atomic_mass_map[340920000] = 91.94984;
//   atomic_mass_map[350920000] = 91.939631597;
//   atomic_mass_map[360920000] = 91.926173094;
//   atomic_mass_map[370920000] = 91.919728389;
//   atomic_mass_map[380920000] = 91.91103819;
//   atomic_mass_map[390920000] = 91.908945142;
//   atomic_mass_map[400920000] = 91.905034675;
//   atomic_mass_map[410920000] = 91.907188081;
//   atomic_mass_map[420920000] = 91.906807959;
//   atomic_mass_map[430920000] = 91.915269779;
//   atomic_mass_map[440920000] = 91.920234375;
//   atomic_mass_map[450920000] = 91.932367694;
//   atomic_mass_map[460920000] = 91.94088;
//   atomic_mass_map[340930000] = 92.95629;
//   atomic_mass_map[350930000] = 92.943134;
//   atomic_mass_map[360930000] = 92.931147174;
//   atomic_mass_map[370930000] = 92.922039269;
//   atomic_mass_map[380930000] = 92.914024228;
//   atomic_mass_map[390930000] = 92.909577886;
//   atomic_mass_map[400930000] = 92.906469947;
//   atomic_mass_map[410930000] = 92.906373004;
//   atomic_mass_map[420930000] = 92.906809577;
//   atomic_mass_map[430930000] = 92.910245952;
//   atomic_mass_map[440930000] = 92.917104444;
//   atomic_mass_map[450930000] = 92.925912781;
//   atomic_mass_map[460930000] = 92.936511;
//   atomic_mass_map[470930000] = 92.95033;
//   atomic_mass_map[340940000] = 93.96049;
//   atomic_mass_map[350940000] = 93.9489;
//   atomic_mass_map[360940000] = 93.934140454;
//   atomic_mass_map[370940000] = 93.926394814;
//   atomic_mass_map[380940000] = 93.915355602;
//   atomic_mass_map[390940000] = 93.911590582;
//   atomic_mass_map[400940000] = 93.906310828;
//   atomic_mass_map[410940000] = 93.907278836;
//   atomic_mass_map[420940000] = 93.905084903;
//   atomic_mass_map[430940000] = 93.909653637;
//   atomic_mass_map[440940000] = 93.911342863;
//   atomic_mass_map[450940000] = 93.921730453;
//   atomic_mass_map[460940000] = 93.929037603;
//   atomic_mass_map[470940000] = 93.943734;
//   atomic_mass_map[340950000] = 94.9673;
//   atomic_mass_map[350950000] = 94.95301;
//   atomic_mass_map[360950000] = 94.939710924;
//   atomic_mass_map[370950000] = 94.929260004;
//   atomic_mass_map[380950000] = 94.919352941;
//   atomic_mass_map[390950000] = 94.912816058;
//   atomic_mass_map[400950000] = 94.90803853;
//   atomic_mass_map[410950000] = 94.906832404;
//   atomic_mass_map[420950000] = 94.905838766;
//   atomic_mass_map[430950000] = 94.907653612;
//   atomic_mass_map[440950000] = 94.910405759;
//   atomic_mass_map[450950000] = 94.915897895;
//   atomic_mass_map[460950000] = 94.924889838;
//   atomic_mass_map[470950000] = 94.93602;
//   atomic_mass_map[480950000] = 94.94994;
//   atomic_mass_map[350960000] = 95.95903;
//   atomic_mass_map[360960000] = 95.943016618;
//   atomic_mass_map[370960000] = 95.93413338;
//   atomic_mass_map[380960000] = 95.921706575;
//   atomic_mass_map[390960000] = 95.915896835;
//   atomic_mass_map[400960000] = 95.908271433;
//   atomic_mass_map[410960000] = 95.908097286;
//   atomic_mass_map[420960000] = 95.904676115;
//   atomic_mass_map[430960000] = 95.907868022;
//   atomic_mass_map[440960000] = 95.907590255;
//   atomic_mass_map[450960000] = 95.914453051;
//   atomic_mass_map[460960000] = 95.918215084;
//   atomic_mass_map[470960000] = 95.930743906;
//   atomic_mass_map[480960000] = 95.94034;
//   atomic_mass_map[350970000] = 96.96344;
//   atomic_mass_map[360970000] = 96.949088785;
//   atomic_mass_map[370970000] = 96.937177136;
//   atomic_mass_map[380970000] = 96.92637396;
//   atomic_mass_map[390970000] = 96.918274106;
//   atomic_mass_map[400970000] = 96.910951206;
//   atomic_mass_map[410970000] = 96.908095932;
//   atomic_mass_map[420970000] = 96.906018118;
//   atomic_mass_map[430970000] = 96.906366706;
//   atomic_mass_map[440970000] = 96.90754712;
//   atomic_mass_map[450970000] = 96.911329216;
//   atomic_mass_map[460970000] = 96.916471988;
//   atomic_mass_map[470970000] = 96.923965326;
//   atomic_mass_map[480970000] = 96.9351;
//   atomic_mass_map[490970000] = 96.94934;
//   atomic_mass_map[350980000] = 97.96946;
//   atomic_mass_map[360980000] = 97.95243;
//   atomic_mass_map[370980000] = 97.941686868;
//   atomic_mass_map[380980000] = 97.928688778;
//   atomic_mass_map[390980000] = 97.922382119;
//   atomic_mass_map[400980000] = 97.912728892;
//   atomic_mass_map[410980000] = 97.910326459;
//   atomic_mass_map[420980000] = 97.90540482;
//   atomic_mass_map[430980000] = 97.907212365;
//   atomic_mass_map[440980000] = 97.905286813;
//   atomic_mass_map[450980000] = 97.910707839;
//   atomic_mass_map[460980000] = 97.912698337;
//   atomic_mass_map[470980000] = 97.921559973;
//   atomic_mass_map[480980000] = 97.927389318;
//   atomic_mass_map[490980000] = 97.94214;
//   atomic_mass_map[360990000] = 98.95839;
//   atomic_mass_map[370990000] = 98.945028735;
//   atomic_mass_map[380990000] = 98.932890666;
//   atomic_mass_map[390990000] = 98.924147979;
//   atomic_mass_map[400990000] = 98.916666746;
//   atomic_mass_map[410990000] = 98.911613177;
//   atomic_mass_map[420990000] = 98.907708509;
//   atomic_mass_map[430990000] = 98.906250844;
//   atomic_mass_map[440990000] = 98.905934082;
//   atomic_mass_map[450990000] = 98.908128239;
//   atomic_mass_map[460990000] = 98.911774806;
//   atomic_mass_map[470990000] = 98.917645768;
//   atomic_mass_map[480990000] = 98.924925848;
//   atomic_mass_map[490990000] = 98.93411;
//   atomic_mass_map[500990000] = 98.94853;
//   atomic_mass_map[361000000] = 99.96237;
//   atomic_mass_map[371000000] = 99.95003;
//   atomic_mass_map[381000000] = 99.935769692;
//   atomic_mass_map[391000000] = 99.927714692;
//   atomic_mass_map[401000000] = 99.918000576;
//   atomic_mass_map[411000000] = 99.914327592;
//   atomic_mass_map[421000000] = 99.907471782;
//   atomic_mass_map[431000000] = 99.907653877;
//   atomic_mass_map[441000000] = 99.904214256;
//   atomic_mass_map[451000000] = 99.908117264;
//   atomic_mass_map[461000000] = 99.908504805;
//   atomic_mass_map[471000000] = 99.916115445;
//   atomic_mass_map[481000000] = 99.92034882;
//   atomic_mass_map[491000000] = 99.93095718;
//   atomic_mass_map[501000000] = 99.938504196;
//   atomic_mass_map[361010000] = 100.96873;
//   atomic_mass_map[371010000] = 100.954039;
//   atomic_mass_map[381010000] = 100.940351743;
//   atomic_mass_map[391010000] = 100.930147705;
//   atomic_mass_map[401010000] = 100.921447964;
//   atomic_mass_map[411010000] = 100.915310254;
//   atomic_mass_map[421010000] = 100.910341447;
//   atomic_mass_map[431010000] = 100.907309057;
//   atomic_mass_map[441010000] = 100.905576872;
//   atomic_mass_map[451010000] = 100.906160613;
//   atomic_mass_map[461010000] = 100.908286412;
//   atomic_mass_map[471010000] = 100.912683953;
//   atomic_mass_map[481010000] = 100.918586211;
//   atomic_mass_map[491010000] = 100.92634;
//   atomic_mass_map[501010000] = 100.935259244;
//   atomic_mass_map[371020000] = 101.95952;
//   atomic_mass_map[381020000] = 101.943790979;
//   atomic_mass_map[391020000] = 101.934327687;
//   atomic_mass_map[401020000] = 101.92314093;
//   atomic_mass_map[411020000] = 101.918077197;
//   atomic_mass_map[421020000] = 101.910283414;
//   atomic_mass_map[431020000] = 101.909209733;
//   atomic_mass_map[441020000] = 101.904344096;
//   atomic_mass_map[451020000] = 101.906837373;
//   atomic_mass_map[461020000] = 101.905602187;
//   atomic_mass_map[471020000] = 101.911704708;
//   atomic_mass_map[481020000] = 101.914481967;
//   atomic_mass_map[491020000] = 101.924107138;
//   atomic_mass_map[501020000] = 101.930290753;
//   atomic_mass_map[371030000] = 102.96392;
//   atomic_mass_map[381030000] = 102.94909;
//   atomic_mass_map[391030000] = 102.937243004;
//   atomic_mass_map[401030000] = 102.927190678;
//   atomic_mass_map[411030000] = 102.919457238;
//   atomic_mass_map[421030000] = 102.913078578;
//   atomic_mass_map[431030000] = 102.909176131;
//   atomic_mass_map[441030000] = 102.906318627;
//   atomic_mass_map[451030000] = 102.905497993;
//   atomic_mass_map[461030000] = 102.906080949;
//   atomic_mass_map[471030000] = 102.908963138;
//   atomic_mass_map[481030000] = 102.913416517;
//   atomic_mass_map[491030000] = 102.919881915;
//   atomic_mass_map[501030000] = 102.928105264;
//   atomic_mass_map[511030000] = 102.93969;
//   atomic_mass_map[381040000] = 103.95265;
//   atomic_mass_map[391040000] = 103.94196;
//   atomic_mass_map[401040000] = 103.929435691;
//   atomic_mass_map[411040000] = 103.922892491;
//   atomic_mass_map[421040000] = 103.913734373;
//   atomic_mass_map[431040000] = 103.911424888;
//   atomic_mass_map[441040000] = 103.905427481;
//   atomic_mass_map[451040000] = 103.90664922;
//   atomic_mass_map[461040000] = 103.90403054;
//   atomic_mass_map[471040000] = 103.908623864;
//   atomic_mass_map[481040000] = 103.909856386;
//   atomic_mass_map[491040000] = 103.918214541;
//   atomic_mass_map[501040000] = 103.923105204;
//   atomic_mass_map[511040000] = 103.936477804;
//   atomic_mass_map[381050000] = 104.95855;
//   atomic_mass_map[391050000] = 104.94544;
//   atomic_mass_map[401050000] = 104.934008204;
//   atomic_mass_map[411050000] = 104.924946471;
//   atomic_mass_map[421050000] = 104.916968617;
//   atomic_mass_map[431050000] = 104.911654883;
//   atomic_mass_map[441050000] = 104.907747645;
//   atomic_mass_map[451050000] = 104.905688549;
//   atomic_mass_map[461050000] = 104.905079626;
//   atomic_mass_map[471050000] = 104.906525615;
//   atomic_mass_map[481050000] = 104.909463896;
//   atomic_mass_map[491050000] = 104.914502325;
//   atomic_mass_map[501050000] = 104.921268429;
//   atomic_mass_map[511050000] = 104.931275897;
//   atomic_mass_map[521050000] = 104.943304508;
//   atomic_mass_map[381060000] = 105.962651;
//   atomic_mass_map[391060000] = 105.95056;
//   atomic_mass_map[401060000] = 105.93676;
//   atomic_mass_map[411060000] = 105.928931712;
//   atomic_mass_map[421060000] = 105.918259464;
//   atomic_mass_map[431060000] = 105.914357598;
//   atomic_mass_map[441060000] = 105.907329104;
//   atomic_mass_map[451060000] = 105.907286801;
//   atomic_mass_map[461060000] = 105.903480426;
//   atomic_mass_map[471060000] = 105.906663637;
//   atomic_mass_map[481060000] = 105.906459928;
//   atomic_mass_map[491060000] = 105.913463735;
//   atomic_mass_map[501060000] = 105.916957404;
//   atomic_mass_map[511060000] = 105.928637982;
//   atomic_mass_map[521060000] = 105.937499664;
//   atomic_mass_map[381070000] = 106.968975;
//   atomic_mass_map[391070000] = 106.95452;
//   atomic_mass_map[401070000] = 106.94174;
//   atomic_mass_map[411070000] = 106.931593654;
//   atomic_mass_map[421070000] = 106.922105877;
//   atomic_mass_map[431070000] = 106.915460645;
//   atomic_mass_map[441070000] = 106.909972045;
//   atomic_mass_map[451070000] = 106.906747811;
//   atomic_mass_map[461070000] = 106.905128195;
//   atomic_mass_map[471070000] = 106.905091611;
//   atomic_mass_map[481070000] = 106.906612122;
//   atomic_mass_map[491070000] = 106.910290084;
//   atomic_mass_map[501070000] = 106.915713652;
//   atomic_mass_map[511070000] = 106.924150641;
//   atomic_mass_map[521070000] = 106.935011573;
//   atomic_mass_map[531070000] = 106.94678;
//   atomic_mass_map[391080000] = 107.95996;
//   atomic_mass_map[401080000] = 107.94487;
//   atomic_mass_map[411080000] = 107.936074773;
//   atomic_mass_map[421080000] = 107.92403349;
//   atomic_mass_map[431080000] = 107.918495722;
//   atomic_mass_map[441080000] = 107.910188022;
//   atomic_mass_map[451080000] = 107.908714473;
//   atomic_mass_map[461080000] = 107.90389164;
//   atomic_mass_map[471080000] = 107.905950346;
//   atomic_mass_map[481080000] = 107.90418344;
//   atomic_mass_map[491080000] = 107.909693524;
//   atomic_mass_map[501080000] = 107.911894287;
//   atomic_mass_map[511080000] = 107.922226735;
//   atomic_mass_map[521080000] = 107.929380467;
//   atomic_mass_map[531080000] = 107.943481623;
//   atomic_mass_map[391090000] = 108.964358;
//   atomic_mass_map[401090000] = 108.95041;
//   atomic_mass_map[411090000] = 108.939216;
//   atomic_mass_map[421090000] = 108.92842416;
//   atomic_mass_map[431090000] = 108.920256356;
//   atomic_mass_map[441090000] = 108.913325956;
//   atomic_mass_map[451090000] = 108.908748821;
//   atomic_mass_map[461090000] = 108.905950406;
//   atomic_mass_map[471090000] = 108.904755282;
//   atomic_mass_map[481090000] = 108.904986653;
//   atomic_mass_map[491090000] = 108.907151381;
//   atomic_mass_map[501090000] = 108.91129206;
//   atomic_mass_map[511090000] = 108.918141122;
//   atomic_mass_map[521090000] = 108.927304534;
//   atomic_mass_map[531090000] = 108.938085287;
//   atomic_mass_map[541090000] = 108.950434864;
//   atomic_mass_map[401100000] = 109.95396;
//   atomic_mass_map[411100000] = 109.94403;
//   atomic_mass_map[421100000] = 109.930703673;
//   atomic_mass_map[431100000] = 109.923743534;
//   atomic_mass_map[441100000] = 109.914040696;
//   atomic_mass_map[451100000] = 109.911079429;
//   atomic_mass_map[461100000] = 109.905172199;
//   atomic_mass_map[471100000] = 109.906110226;
//   atomic_mass_map[481100000] = 109.903006606;
//   atomic_mass_map[491100000] = 109.90716981;
//   atomic_mass_map[501100000] = 109.907844835;
//   atomic_mass_map[511100000] = 109.916854287;
//   atomic_mass_map[521100000] = 109.922458091;
//   atomic_mass_map[531100000] = 109.935089034;
//   atomic_mass_map[541100000] = 109.944263102;
//   atomic_mass_map[401110000] = 110.959678;
//   atomic_mass_map[411110000] = 110.94753;
//   atomic_mass_map[421110000] = 110.935654257;
//   atomic_mass_map[431110000] = 110.925901257;
//   atomic_mass_map[441110000] = 110.917569857;
//   atomic_mass_map[451110000] = 110.91164231;
//   atomic_mass_map[461110000] = 110.907689679;
//   atomic_mass_map[471110000] = 110.905295923;
//   atomic_mass_map[481110000] = 110.904182872;
//   atomic_mass_map[491110000] = 110.905108458;
//   atomic_mass_map[501110000] = 110.907740127;
//   atomic_mass_map[511110000] = 110.913218189;
//   atomic_mass_map[521110000] = 110.921000589;
//   atomic_mass_map[531110000] = 110.930269214;
//   atomic_mass_map[541110000] = 110.941607206;
//   atomic_mass_map[401120000] = 111.963703;
//   atomic_mass_map[411120000] = 111.95247;
//   atomic_mass_map[421120000] = 111.93831;
//   atomic_mass_map[431120000] = 111.929945813;
//   atomic_mass_map[441120000] = 111.918809234;
//   atomic_mass_map[451120000] = 111.914403222;
//   atomic_mass_map[461120000] = 111.907329698;
//   atomic_mass_map[471120000] = 111.90704855;
//   atomic_mass_map[481120000] = 111.902762868;
//   atomic_mass_map[491120000] = 111.905537694;
//   atomic_mass_map[501120000] = 111.904823874;
//   atomic_mass_map[511120000] = 111.912399903;
//   atomic_mass_map[521120000] = 111.91672785;
//   atomic_mass_map[531120000] = 111.92800455;
//   atomic_mass_map[541120000] = 111.935558982;
//   atomic_mass_map[551120000] = 111.950308558;
//   atomic_mass_map[411130000] = 112.95651;
//   atomic_mass_map[421130000] = 112.94335;
//   atomic_mass_map[431130000] = 112.932569035;
//   atomic_mass_map[441130000] = 112.922843999;
//   atomic_mass_map[451130000] = 112.915439342;
//   atomic_mass_map[461130000] = 112.910261042;
//   atomic_mass_map[471130000] = 112.906572893;
//   atomic_mass_map[481130000] = 112.904408133;
//   atomic_mass_map[491130000] = 112.904061839;
//   atomic_mass_map[501130000] = 112.905175728;
//   atomic_mass_map[511130000] = 112.909374536;
//   atomic_mass_map[521130000] = 112.915891;
//   atomic_mass_map[531130000] = 112.923650064;
//   atomic_mass_map[541130000] = 112.93322165;
//   atomic_mass_map[551130000] = 112.944429144;
//   atomic_mass_map[411140000] = 113.96201;
//   atomic_mass_map[421140000] = 113.94653;
//   atomic_mass_map[431140000] = 113.93691;
//   atomic_mass_map[441140000] = 113.924613554;
//   atomic_mass_map[451140000] = 113.918718294;
//   atomic_mass_map[461140000] = 113.910368554;
//   atomic_mass_map[471140000] = 113.908823031;
//   atomic_mass_map[481140000] = 113.903365086;
//   atomic_mass_map[491140000] = 113.904917909;
//   atomic_mass_map[501140000] = 113.902782695;
//   atomic_mass_map[511140000] = 113.909290189;
//   atomic_mass_map[521140000] = 113.912089;
//   atomic_mass_map[531140000] = 113.92185;
//   atomic_mass_map[541140000] = 113.927980331;
//   atomic_mass_map[551140000] = 113.941296176;
//   atomic_mass_map[561140000] = 113.950660438;
//   atomic_mass_map[411150000] = 114.96634;
//   atomic_mass_map[421150000] = 114.95196;
//   atomic_mass_map[431150000] = 114.93998;
//   atomic_mass_map[441150000] = 114.928819898;
//   atomic_mass_map[451150000] = 114.920311589;
//   atomic_mass_map[461150000] = 114.913658506;
//   atomic_mass_map[471150000] = 114.908767355;
//   atomic_mass_map[481150000] = 114.905437513;
//   atomic_mass_map[491150000] = 114.903878776;
//   atomic_mass_map[501150000] = 114.903344699;
//   atomic_mass_map[511150000] = 114.906598;
//   atomic_mass_map[521150000] = 114.911902;
//   atomic_mass_map[531150000] = 114.918048;
//   atomic_mass_map[541150000] = 114.926293945;
//   atomic_mass_map[551150000] = 114.93591;
//   atomic_mass_map[561150000] = 114.94737;
//   atomic_mass_map[421160000] = 115.955448;
//   atomic_mass_map[431160000] = 115.94476;
//   atomic_mass_map[441160000] = 115.931219195;
//   atomic_mass_map[451160000] = 115.924058528;
//   atomic_mass_map[461160000] = 115.914296979;
//   atomic_mass_map[471160000] = 115.911386812;
//   atomic_mass_map[481160000] = 115.904763148;
//   atomic_mass_map[491160000] = 115.905259995;
//   atomic_mass_map[501160000] = 115.901742797;
//   atomic_mass_map[511160000] = 115.906793115;
//   atomic_mass_map[521160000] = 115.90846;
//   atomic_mass_map[531160000] = 115.916808658;
//   atomic_mass_map[541160000] = 115.921581112;
//   atomic_mass_map[551160000] = 115.933373;
//   atomic_mass_map[561160000] = 115.94128;
//   atomic_mass_map[571160000] = 115.956304;
//   atomic_mass_map[421170000] = 116.96117;
//   atomic_mass_map[431170000] = 116.94806;
//   atomic_mass_map[441170000] = 116.936103;
//   atomic_mass_map[451170000] = 116.926035391;
//   atomic_mass_map[461170000] = 116.917954721;
//   atomic_mass_map[471170000] = 116.911773934;
//   atomic_mass_map[481170000] = 116.907225956;
//   atomic_mass_map[491170000] = 116.904515678;
//   atomic_mass_map[501170000] = 116.902953983;
//   atomic_mass_map[511170000] = 116.904841508;
//   atomic_mass_map[521170000] = 116.908646298;
//   atomic_mass_map[531170000] = 116.913648312;
//   atomic_mass_map[541170000] = 116.920358761;
//   atomic_mass_map[551170000] = 116.928616726;
//   atomic_mass_map[561170000] = 116.93814057;
//   atomic_mass_map[571170000] = 116.949985;
//   atomic_mass_map[431180000] = 117.95299;
//   atomic_mass_map[441180000] = 117.93853;
//   atomic_mass_map[451180000] = 117.930340208;
//   atomic_mass_map[461180000] = 117.9190667;
//   atomic_mass_map[471180000] = 117.914595487;
//   atomic_mass_map[481180000] = 117.906921869;
//   atomic_mass_map[491180000] = 117.906356616;
//   atomic_mass_map[501180000] = 117.901606574;
//   atomic_mass_map[511180000] = 117.905532139;
//   atomic_mass_map[521180000] = 117.905853629;
//   atomic_mass_map[531180000] = 117.913074;
//   atomic_mass_map[541180000] = 117.91617868;
//   atomic_mass_map[551180000] = 117.926559519;
//   atomic_mass_map[561180000] = 117.93306;
//   atomic_mass_map[571180000] = 117.94673;
//   atomic_mass_map[431190000] = 118.95666;
//   atomic_mass_map[441190000] = 118.94357;
//   atomic_mass_map[451190000] = 118.932556954;
//   atomic_mass_map[461190000] = 118.923340223;
//   atomic_mass_map[471190000] = 118.915570287;
//   atomic_mass_map[481190000] = 118.909846851;
//   atomic_mass_map[491190000] = 118.905850708;
//   atomic_mass_map[501190000] = 118.903311172;
//   atomic_mass_map[511190000] = 118.903945471;
//   atomic_mass_map[521190000] = 118.906407108;
//   atomic_mass_map[531190000] = 118.910074;
//   atomic_mass_map[541190000] = 118.915410714;
//   atomic_mass_map[551190000] = 118.92237733;
//   atomic_mass_map[561190000] = 118.930659686;
//   atomic_mass_map[571190000] = 118.94099;
//   atomic_mass_map[581190000] = 118.95271;
//   atomic_mass_map[431200000] = 119.96187;
//   atomic_mass_map[441200000] = 119.94631;
//   atomic_mass_map[451200000] = 119.93686;
//   atomic_mass_map[461200000] = 119.924551089;
//   atomic_mass_map[471200000] = 119.918784768;
//   atomic_mass_map[481200000] = 119.909868068;
//   atomic_mass_map[491200000] = 119.907966567;
//   atomic_mass_map[501200000] = 119.902201634;
//   atomic_mass_map[511200000] = 119.905079385;
//   atomic_mass_map[521200000] = 119.9040593;
//   atomic_mass_map[531200000] = 119.910087251;
//   atomic_mass_map[541200000] = 119.91178427;
//   atomic_mass_map[551200000] = 119.920677279;
//   atomic_mass_map[561200000] = 119.926045;
//   atomic_mass_map[571200000] = 119.93807;
//   atomic_mass_map[581200000] = 119.94654;
//   atomic_mass_map[441210000] = 120.95164;
//   atomic_mass_map[451210000] = 120.93942;
//   atomic_mass_map[461210000] = 120.928950345;
//   atomic_mass_map[471210000] = 120.920125282;
//   atomic_mass_map[481210000] = 120.912963663;
//   atomic_mass_map[491210000] = 120.907851062;
//   atomic_mass_map[501210000] = 120.904242554;
//   atomic_mass_map[511210000] = 120.903811967;
//   atomic_mass_map[521210000] = 120.904943812;
//   atomic_mass_map[531210000] = 120.90740505;
//   atomic_mass_map[541210000] = 120.911453015;
//   atomic_mass_map[551210000] = 120.917227238;
//   atomic_mass_map[561210000] = 120.92405229;
//   atomic_mass_map[571210000] = 120.93315;
//   atomic_mass_map[581210000] = 120.94335;
//   atomic_mass_map[591210000] = 120.95532;
//   atomic_mass_map[441220000] = 121.95447;
//   atomic_mass_map[451220000] = 121.94399;
//   atomic_mass_map[461220000] = 121.930631696;
//   atomic_mass_map[471220000] = 121.923664449;
//   atomic_mass_map[481220000] = 121.913459053;
//   atomic_mass_map[491220000] = 121.910280738;
//   atomic_mass_map[501220000] = 121.903443774;
//   atomic_mass_map[511220000] = 121.905169948;
//   atomic_mass_map[521220000] = 121.903043455;
//   atomic_mass_map[531220000] = 121.907588841;
//   atomic_mass_map[541220000] = 121.908367658;
//   atomic_mass_map[551220000] = 121.916108146;
//   atomic_mass_map[561220000] = 121.919904;
//   atomic_mass_map[571220000] = 121.93071;
//   atomic_mass_map[581220000] = 121.93787;
//   atomic_mass_map[591220000] = 121.95175;
//   atomic_mass_map[441230000] = 122.95989;
//   atomic_mass_map[451230000] = 122.94685;
//   atomic_mass_map[461230000] = 122.93514;
//   atomic_mass_map[471230000] = 122.925337063;
//   atomic_mass_map[481230000] = 122.916892453;
//   atomic_mass_map[491230000] = 122.910433718;
//   atomic_mass_map[501230000] = 122.905725221;
//   atomic_mass_map[511230000] = 122.904213204;
//   atomic_mass_map[521230000] = 122.904269769;
//   atomic_mass_map[531230000] = 122.905588537;
//   atomic_mass_map[541230000] = 122.908481757;
//   atomic_mass_map[551230000] = 122.912996063;
//   atomic_mass_map[561230000] = 122.918781063;
//   atomic_mass_map[571230000] = 122.9263;
//   atomic_mass_map[581230000] = 122.93528;
//   atomic_mass_map[591230000] = 122.94596;
//   atomic_mass_map[441240000] = 123.96305;
//   atomic_mass_map[451240000] = 123.95151;
//   atomic_mass_map[461240000] = 123.93714;
//   atomic_mass_map[471240000] = 123.928931229;
//   atomic_mass_map[481240000] = 123.917657364;
//   atomic_mass_map[491240000] = 123.913182245;
//   atomic_mass_map[501240000] = 123.905276645;
//   atomic_mass_map[511240000] = 123.905934978;
//   atomic_mass_map[521240000] = 123.902817085;
//   atomic_mass_map[531240000] = 123.906209041;
//   atomic_mass_map[541240000] = 123.905891984;
//   atomic_mass_map[551240000] = 123.912257798;
//   atomic_mass_map[561240000] = 123.915093629;
//   atomic_mass_map[571240000] = 123.924574275;
//   atomic_mass_map[581240000] = 123.93031;
//   atomic_mass_map[591240000] = 123.94294;
//   atomic_mass_map[601240000] = 123.9522;
//   atomic_mass_map[451250000] = 124.95469;
//   atomic_mass_map[461250000] = 124.94179;
//   atomic_mass_map[471250000] = 124.931046;
//   atomic_mass_map[481250000] = 124.921257577;
//   atomic_mass_map[491250000] = 124.913604553;
//   atomic_mass_map[501250000] = 124.907786395;
//   atomic_mass_map[511250000] = 124.905253007;
//   atomic_mass_map[521250000] = 124.90442992;
//   atomic_mass_map[531250000] = 124.904629353;
//   atomic_mass_map[541250000] = 124.906394445;
//   atomic_mass_map[551250000] = 124.909727985;
//   atomic_mass_map[561250000] = 124.914471843;
//   atomic_mass_map[571250000] = 124.920815932;
//   atomic_mass_map[581250000] = 124.92844;
//   atomic_mass_map[591250000] = 124.9377;
//   atomic_mass_map[601250000] = 124.9489;
//   atomic_mass_map[451260000] = 125.95946;
//   atomic_mass_map[461260000] = 125.94416;
//   atomic_mass_map[471260000] = 125.93475;
//   atomic_mass_map[481260000] = 125.922429128;
//   atomic_mass_map[491260000] = 125.916507322;
//   atomic_mass_map[501260000] = 125.907658786;
//   atomic_mass_map[511260000] = 125.907252987;
//   atomic_mass_map[521260000] = 125.903310886;
//   atomic_mass_map[531260000] = 125.905623329;
//   atomic_mass_map[541260000] = 125.904298292;
//   atomic_mass_map[551260000] = 125.909446057;
//   atomic_mass_map[561260000] = 125.911250204;
//   atomic_mass_map[571260000] = 125.919512667;
//   atomic_mass_map[581260000] = 125.923971;
//   atomic_mass_map[591260000] = 125.93524;
//   atomic_mass_map[601260000] = 125.94311;
//   atomic_mass_map[611260000] = 125.95792;
//   atomic_mass_map[461270000] = 126.94907;
//   atomic_mass_map[471270000] = 126.93711;
//   atomic_mass_map[481270000] = 126.926472404;
//   atomic_mass_map[491270000] = 126.917446276;
//   atomic_mass_map[501270000] = 126.910389943;
//   atomic_mass_map[511270000] = 126.906924278;
//   atomic_mass_map[521270000] = 126.905225735;
//   atomic_mass_map[531270000] = 126.904471853;
//   atomic_mass_map[541270000] = 126.905182914;
//   atomic_mass_map[551270000] = 126.907417384;
//   atomic_mass_map[561270000] = 126.911091275;
//   atomic_mass_map[571270000] = 126.916375084;
//   atomic_mass_map[581270000] = 126.922727;
//   atomic_mass_map[591270000] = 126.93071;
//   atomic_mass_map[601270000] = 126.94038;
//   atomic_mass_map[611270000] = 126.95192;
//   atomic_mass_map[461280000] = 127.95183;
//   atomic_mass_map[471280000] = 127.94106;
//   atomic_mass_map[481280000] = 127.927812858;
//   atomic_mass_map[491280000] = 127.920401044;
//   atomic_mass_map[501280000] = 127.910507184;
//   atomic_mass_map[511280000] = 127.909145628;
//   atomic_mass_map[521280000] = 127.904461279;
//   atomic_mass_map[531280000] = 127.905808615;
//   atomic_mass_map[541280000] = 127.903531018;
//   atomic_mass_map[551280000] = 127.907748666;
//   atomic_mass_map[561280000] = 127.908341967;
//   atomic_mass_map[571280000] = 127.915592123;
//   atomic_mass_map[581280000] = 127.918911;
//   atomic_mass_map[591280000] = 127.928791;
//   atomic_mass_map[601280000] = 127.93525;
//   atomic_mass_map[611280000] = 127.9487;
//   atomic_mass_map[621280000] = 127.95842;
//   atomic_mass_map[471290000] = 128.94395;
//   atomic_mass_map[481290000] = 128.93182;
//   atomic_mass_map[491290000] = 128.921805301;
//   atomic_mass_map[501290000] = 128.913464711;
//   atomic_mass_map[511290000] = 128.909146665;
//   atomic_mass_map[521290000] = 128.90659646;
//   atomic_mass_map[531290000] = 128.904983669;
//   atomic_mass_map[541290000] = 128.904780861;
//   atomic_mass_map[551290000] = 128.906065683;
//   atomic_mass_map[561290000] = 128.908680798;
//   atomic_mass_map[571290000] = 128.912694431;
//   atomic_mass_map[581290000] = 128.918102;
//   atomic_mass_map[591290000] = 128.925095;
//   atomic_mass_map[601290000] = 128.933102;
//   atomic_mass_map[611290000] = 128.94323;
//   atomic_mass_map[621290000] = 128.95476;
//   atomic_mass_map[471300000] = 129.950703;
//   atomic_mass_map[481300000] = 129.933940679;
//   atomic_mass_map[491300000] = 129.924976585;
//   atomic_mass_map[501300000] = 129.91397383;
//   atomic_mass_map[511300000] = 129.911662054;
//   atomic_mass_map[521300000] = 129.906222748;
//   atomic_mass_map[531300000] = 129.906670193;
//   atomic_mass_map[541300000] = 129.903509349;
//   atomic_mass_map[551300000] = 129.906709283;
//   atomic_mass_map[561300000] = 129.906320669;
//   atomic_mass_map[571300000] = 129.912369413;
//   atomic_mass_map[581300000] = 129.914736;
//   atomic_mass_map[591300000] = 129.92359;
//   atomic_mass_map[601300000] = 129.928506;
//   atomic_mass_map[611300000] = 129.94053;
//   atomic_mass_map[621300000] = 129.949;
//   atomic_mass_map[631300000] = 129.963689;
//   atomic_mass_map[481310000] = 130.9406;
//   atomic_mass_map[491310000] = 130.926971529;
//   atomic_mass_map[501310000] = 130.917044954;
//   atomic_mass_map[511310000] = 130.911988799;
//   atomic_mass_map[521310000] = 130.908522213;
//   atomic_mass_map[531310000] = 130.906126305;
//   atomic_mass_map[541310000] = 130.905084057;
//   atomic_mass_map[551310000] = 130.905464899;
//   atomic_mass_map[561310000] = 130.906940977;
//   atomic_mass_map[571310000] = 130.91007;
//   atomic_mass_map[581310000] = 130.914429465;
//   atomic_mass_map[591310000] = 130.92023496;
//   atomic_mass_map[601310000] = 130.92724802;
//   atomic_mass_map[611310000] = 130.93567;
//   atomic_mass_map[621310000] = 130.94618;
//   atomic_mass_map[631310000] = 130.957842;
//   atomic_mass_map[481320000] = 131.94604;
//   atomic_mass_map[491320000] = 131.933001273;
//   atomic_mass_map[501320000] = 131.917826725;
//   atomic_mass_map[511320000] = 131.914507691;
//   atomic_mass_map[521320000] = 131.908546716;
//   atomic_mass_map[531320000] = 131.907993514;
//   atomic_mass_map[541320000] = 131.904155086;
//   atomic_mass_map[551320000] = 131.906433914;
//   atomic_mass_map[561320000] = 131.905061128;
//   atomic_mass_map[571320000] = 131.910118979;
//   atomic_mass_map[581320000] = 131.911463775;
//   atomic_mass_map[591320000] = 131.919255;
//   atomic_mass_map[601320000] = 131.923321237;
//   atomic_mass_map[611320000] = 131.93384;
//   atomic_mass_map[621320000] = 131.94087;
//   atomic_mass_map[631320000] = 131.95467;
//   atomic_mass_map[481330000] = 132.95285;
//   atomic_mass_map[491330000] = 132.93831;
//   atomic_mass_map[501330000] = 132.923913404;
//   atomic_mass_map[511330000] = 132.915273198;
//   atomic_mass_map[521330000] = 132.910968766;
//   atomic_mass_map[531330000] = 132.907796968;
//   atomic_mass_map[541330000] = 132.905910751;
//   atomic_mass_map[551330000] = 132.905451961;
//   atomic_mass_map[561330000] = 132.906007351;
//   atomic_mass_map[571330000] = 132.908218;
//   atomic_mass_map[581330000] = 132.911520402;
//   atomic_mass_map[591330000] = 132.916330561;
//   atomic_mass_map[601330000] = 132.922348;
//   atomic_mass_map[611330000] = 132.929782;
//   atomic_mass_map[621330000] = 132.93856;
//   atomic_mass_map[631330000] = 132.94929;
//   atomic_mass_map[641330000] = 132.96133;
//   atomic_mass_map[491340000] = 133.94454;
//   atomic_mass_map[501340000] = 133.928682064;
//   atomic_mass_map[511340000] = 133.920535676;
//   atomic_mass_map[521340000] = 133.911393959;
//   atomic_mass_map[531340000] = 133.909758836;
//   atomic_mass_map[541340000] = 133.905394664;
//   atomic_mass_map[551340000] = 133.906718503;
//   atomic_mass_map[561340000] = 133.904508182;
//   atomic_mass_map[571340000] = 133.908514011;
//   atomic_mass_map[581340000] = 133.908928142;
//   atomic_mass_map[591340000] = 133.915696729;
//   atomic_mass_map[601340000] = 133.91879021;
//   atomic_mass_map[611340000] = 133.928353;
//   atomic_mass_map[621340000] = 133.93411;
//   atomic_mass_map[631340000] = 133.9464;
//   atomic_mass_map[641340000] = 133.95566;
//   atomic_mass_map[491350000] = 134.95005;
//   atomic_mass_map[501350000] = 134.934908606;
//   atomic_mass_map[511350000] = 134.925185106;
//   atomic_mass_map[521350000] = 134.916555706;
//   atomic_mass_map[531350000] = 134.910048847;
//   atomic_mass_map[541350000] = 134.90722778;
//   atomic_mass_map[551350000] = 134.905977049;
//   atomic_mass_map[561350000] = 134.905688375;
//   atomic_mass_map[571350000] = 134.906984363;
//   atomic_mass_map[581350000] = 134.909160599;
//   atomic_mass_map[591350000] = 134.913111775;
//   atomic_mass_map[601350000] = 134.918181321;
//   atomic_mass_map[611350000] = 134.924823;
//   atomic_mass_map[621350000] = 134.93252;
//   atomic_mass_map[631350000] = 134.94187;
//   atomic_mass_map[641350000] = 134.95245;
//   atomic_mass_map[651350000] = 134.96476;
//   atomic_mass_map[501360000] = 135.93999;
//   atomic_mass_map[511360000] = 135.930745879;
//   atomic_mass_map[521360000] = 135.920100608;
//   atomic_mass_map[531360000] = 135.914604412;
//   atomic_mass_map[541360000] = 135.907214484;
//   atomic_mass_map[551360000] = 135.907311358;
//   atomic_mass_map[561360000] = 135.904575727;
//   atomic_mass_map[571360000] = 135.907634962;
//   atomic_mass_map[581360000] = 135.907129205;
//   atomic_mass_map[591360000] = 135.912677456;
//   atomic_mass_map[601360000] = 135.914976064;
//   atomic_mass_map[611360000] = 135.923584586;
//   atomic_mass_map[621360000] = 135.928275556;
//   atomic_mass_map[631360000] = 135.93962;
//   atomic_mass_map[641360000] = 135.9473;
//   atomic_mass_map[651360000] = 135.96129;
//   atomic_mass_map[501370000] = 136.94655;
//   atomic_mass_map[511370000] = 136.935555;
//   atomic_mass_map[521370000] = 136.925598852;
//   atomic_mass_map[531370000] = 136.918028188;
//   atomic_mass_map[541370000] = 136.911557781;
//   atomic_mass_map[551370000] = 136.907089231;
//   atomic_mass_map[561370000] = 136.905827141;
//   atomic_mass_map[571370000] = 136.906450385;
//   atomic_mass_map[581370000] = 136.907762364;
//   atomic_mass_map[591370000] = 136.91067915;
//   atomic_mass_map[601370000] = 136.914562448;
//   atomic_mass_map[611370000] = 136.920479523;
//   atomic_mass_map[621370000] = 136.926970517;
//   atomic_mass_map[631370000] = 136.93546;
//   atomic_mass_map[641370000] = 136.94502;
//   atomic_mass_map[651370000] = 136.95602;
//   atomic_mass_map[501380000] = 137.95184;
//   atomic_mass_map[511380000] = 137.94145;
//   atomic_mass_map[521380000] = 137.929472226;
//   atomic_mass_map[531380000] = 137.922726402;
//   atomic_mass_map[541380000] = 137.914146273;
//   atomic_mass_map[551380000] = 137.911017089;
//   atomic_mass_map[561380000] = 137.905246995;
//   atomic_mass_map[571380000] = 137.907114919;
//   atomic_mass_map[581380000] = 137.905991089;
//   atomic_mass_map[591380000] = 137.910754405;
//   atomic_mass_map[601380000] = 137.911949718;
//   atomic_mass_map[611380000] = 137.919548077;
//   atomic_mass_map[621380000] = 137.923243991;
//   atomic_mass_map[631380000] = 137.933709;
//   atomic_mass_map[641380000] = 137.94025;
//   atomic_mass_map[651380000] = 137.95312;
//   atomic_mass_map[661380000] = 137.9625;
//   atomic_mass_map[511390000] = 138.94655;
//   atomic_mass_map[521390000] = 138.935367194;
//   atomic_mass_map[531390000] = 138.926506206;
//   atomic_mass_map[541390000] = 138.918792203;
//   atomic_mass_map[551390000] = 138.913363758;
//   atomic_mass_map[561390000] = 138.908841099;
//   atomic_mass_map[571390000] = 138.906356256;
//   atomic_mass_map[581390000] = 138.906655111;
//   atomic_mass_map[591390000] = 138.908940762;
//   atomic_mass_map[601390000] = 138.911953649;
//   atomic_mass_map[611390000] = 138.91679967;
//   atomic_mass_map[621390000] = 138.922296635;
//   atomic_mass_map[631390000] = 138.92979231;
//   atomic_mass_map[641390000] = 138.93813;
//   atomic_mass_map[651390000] = 138.94833;
//   atomic_mass_map[661390000] = 138.95959;
//   atomic_mass_map[511400000] = 139.95283;
//   atomic_mass_map[521400000] = 139.939498567;
//   atomic_mass_map[531400000] = 139.931727;
//   atomic_mass_map[541400000] = 139.921645817;
//   atomic_mass_map[551400000] = 139.917283063;
//   atomic_mass_map[561400000] = 139.91060573;
//   atomic_mass_map[571400000] = 139.909480635;
//   atomic_mass_map[581400000] = 139.905443107;
//   atomic_mass_map[591400000] = 139.909080275;
//   atomic_mass_map[601400000] = 139.909549849;
//   atomic_mass_map[611400000] = 139.916039639;
//   atomic_mass_map[621400000] = 139.918994717;
//   atomic_mass_map[631400000] = 139.928087637;
//   atomic_mass_map[641400000] = 139.933674;
//   atomic_mass_map[651400000] = 139.945805049;
//   atomic_mass_map[661400000] = 139.95402;
//   atomic_mass_map[671400000] = 139.968589;
//   atomic_mass_map[521410000] = 140.9458;
//   atomic_mass_map[531410000] = 140.93569;
//   atomic_mass_map[541410000] = 140.926787184;
//   atomic_mass_map[551410000] = 140.920045491;
//   atomic_mass_map[561410000] = 140.914403333;
//   atomic_mass_map[571410000] = 140.91096603;
//   atomic_mass_map[581410000] = 140.908280674;
//   atomic_mass_map[591410000] = 140.907657568;
//   atomic_mass_map[601410000] = 140.909614653;
//   atomic_mass_map[611410000] = 140.913555084;
//   atomic_mass_map[621410000] = 140.918481636;
//   atomic_mass_map[631410000] = 140.924931754;
//   atomic_mass_map[641410000] = 140.932126;
//   atomic_mass_map[651410000] = 140.941448;
//   atomic_mass_map[661410000] = 140.95128;
//   atomic_mass_map[671410000] = 140.963108;
//   atomic_mass_map[521420000] = 141.95022;
//   atomic_mass_map[531420000] = 141.941202;
//   atomic_mass_map[541420000] = 141.929973098;
//   atomic_mass_map[551420000] = 141.924295991;
//   atomic_mass_map[561420000] = 141.916432382;
//   atomic_mass_map[571420000] = 141.91409089;
//   atomic_mass_map[581420000] = 141.909250375;
//   atomic_mass_map[591420000] = 141.910049607;
//   atomic_mass_map[601420000] = 141.907728996;
//   atomic_mass_map[611420000] = 141.912890411;
//   atomic_mass_map[621420000] = 141.915204385;
//   atomic_mass_map[631420000] = 141.92344169;
//   atomic_mass_map[641420000] = 141.928116;
//   atomic_mass_map[651420000] = 141.939280859;
//   atomic_mass_map[661420000] = 141.946194;
//   atomic_mass_map[671420000] = 141.96001;
//   atomic_mass_map[681420000] = 141.9701;
//   atomic_mass_map[521430000] = 142.95676;
//   atomic_mass_map[531430000] = 142.94565;
//   atomic_mass_map[541430000] = 142.935369554;
//   atomic_mass_map[551430000] = 142.927349327;
//   atomic_mass_map[561430000] = 142.920625293;
//   atomic_mass_map[571430000] = 142.916079512;
//   atomic_mass_map[581430000] = 142.91239212;
//   atomic_mass_map[591430000] = 142.910822796;
//   atomic_mass_map[601430000] = 142.909819989;
//   atomic_mass_map[611430000] = 142.910938262;
//   atomic_mass_map[621430000] = 142.914635303;
//   atomic_mass_map[631430000] = 142.920298681;
//   atomic_mass_map[641430000] = 142.926750682;
//   atomic_mass_map[651430000] = 142.935137335;
//   atomic_mass_map[661430000] = 142.943994335;
//   atomic_mass_map[671430000] = 142.95486;
//   atomic_mass_map[681430000] = 142.96662;
//   atomic_mass_map[531440000] = 143.95139;
//   atomic_mass_map[541440000] = 143.938945079;
//   atomic_mass_map[551440000] = 143.932076313;
//   atomic_mass_map[561440000] = 143.922954866;
//   atomic_mass_map[571440000] = 143.919645589;
//   atomic_mass_map[581440000] = 143.913652939;
//   atomic_mass_map[591440000] = 143.913310859;
//   atomic_mass_map[601440000] = 143.910092974;
//   atomic_mass_map[611440000] = 143.912596396;
//   atomic_mass_map[621440000] = 143.912006466;
//   atomic_mass_map[631440000] = 143.918819553;
//   atomic_mass_map[641440000] = 143.922963;
//   atomic_mass_map[651440000] = 143.933045;
//   atomic_mass_map[661440000] = 143.939269515;
//   atomic_mass_map[671440000] = 143.952109715;
//   atomic_mass_map[681440000] = 143.9607;
//   atomic_mass_map[691440000] = 143.976283;
//   atomic_mass_map[531450000] = 144.95605;
//   atomic_mass_map[541450000] = 144.944719634;
//   atomic_mass_map[551450000] = 144.935527435;
//   atomic_mass_map[561450000] = 144.9275184;
//   atomic_mass_map[571450000] = 144.921808068;
//   atomic_mass_map[581450000] = 144.917265228;
//   atomic_mass_map[591450000] = 144.914518156;
//   atomic_mass_map[601450000] = 144.912579322;
//   atomic_mass_map[611450000] = 144.912755935;
//   atomic_mass_map[621450000] = 144.913417339;
//   atomic_mass_map[631450000] = 144.916272629;
//   atomic_mass_map[641450000] = 144.921712821;
//   atomic_mass_map[651450000] = 144.928821947;
//   atomic_mass_map[661450000] = 144.937473995;
//   atomic_mass_map[671450000] = 144.947267395;
//   atomic_mass_map[681450000] = 144.95805;
//   atomic_mass_map[691450000] = 144.970389;
//   atomic_mass_map[541460000] = 145.948518249;
//   atomic_mass_map[551460000] = 145.940344271;
//   atomic_mass_map[561460000] = 145.930283712;
//   atomic_mass_map[571460000] = 145.925875174;
//   atomic_mass_map[581460000] = 145.918802412;
//   atomic_mass_map[591460000] = 145.91767985;
//   atomic_mass_map[601460000] = 145.913122628;
//   atomic_mass_map[611460000] = 145.914702396;
//   atomic_mass_map[621460000] = 145.913046991;
//   atomic_mass_map[631460000] = 145.917211039;
//   atomic_mass_map[641460000] = 145.918318817;
//   atomic_mass_map[651460000] = 145.927252984;
//   atomic_mass_map[661460000] = 145.93284453;
//   atomic_mass_map[671460000] = 145.944993506;
//   atomic_mass_map[681460000] = 145.95241836;
//   atomic_mass_map[691460000] = 145.966837;
//   atomic_mass_map[541470000] = 146.95426;
//   atomic_mass_map[551470000] = 146.944156156;
//   atomic_mass_map[561470000] = 146.9353039;
//   atomic_mass_map[571470000] = 146.9284178;
//   atomic_mass_map[581470000] = 146.922689913;
//   atomic_mass_map[591470000] = 146.919007518;
//   atomic_mass_map[601470000] = 146.916106136;
//   atomic_mass_map[611470000] = 146.915144988;
//   atomic_mass_map[621470000] = 146.914904435;
//   atomic_mass_map[631470000] = 146.916752659;
//   atomic_mass_map[641470000] = 146.919101384;
//   atomic_mass_map[651470000] = 146.924054788;
//   atomic_mass_map[661470000] = 146.931082715;
//   atomic_mass_map[671470000] = 146.940142296;
//   atomic_mass_map[681470000] = 146.949964459;
//   atomic_mass_map[691470000] = 146.961379891;
//   atomic_mass_map[541480000] = 147.95813;
//   atomic_mass_map[551480000] = 147.949225137;
//   atomic_mass_map[561480000] = 147.938170578;
//   atomic_mass_map[571480000] = 147.9326794;
//   atomic_mass_map[581480000] = 147.924424225;
//   atomic_mass_map[591480000] = 147.922130083;
//   atomic_mass_map[601480000] = 147.916899294;
//   atomic_mass_map[611480000] = 147.917481945;
//   atomic_mass_map[621480000] = 147.914829226;
//   atomic_mass_map[631480000] = 147.918089243;
//   atomic_mass_map[641480000] = 147.918121511;
//   atomic_mass_map[651480000] = 147.924281552;
//   atomic_mass_map[661480000] = 147.927156571;
//   atomic_mass_map[671480000] = 147.937743928;
//   atomic_mass_map[681480000] = 147.944735029;
//   atomic_mass_map[691480000] = 147.958384029;
//   atomic_mass_map[701480000] = 147.96758;
//   atomic_mass_map[551490000] = 148.95302;
//   atomic_mass_map[561490000] = 148.94308;
//   atomic_mass_map[571490000] = 148.93535126;
//   atomic_mass_map[581490000] = 148.9284269;
//   atomic_mass_map[591490000] = 148.9237361;
//   atomic_mass_map[601490000] = 148.920154849;
//   atomic_mass_map[611490000] = 148.918342277;
//   atomic_mass_map[621490000] = 148.917192062;
//   atomic_mass_map[631490000] = 148.917937763;
//   atomic_mass_map[641490000] = 148.919348117;
//   atomic_mass_map[651490000] = 148.923253549;
//   atomic_mass_map[661490000] = 148.927321692;
//   atomic_mass_map[671490000] = 148.933802646;
//   atomic_mass_map[681490000] = 148.942306;
//   atomic_mass_map[691490000] = 148.95289;
//   atomic_mass_map[701490000] = 148.96436;
//   atomic_mass_map[551500000] = 149.95833;
//   atomic_mass_map[561500000] = 149.94605;
//   atomic_mass_map[571500000] = 149.93947;
//   atomic_mass_map[581500000] = 149.930384042;
//   atomic_mass_map[591500000] = 149.926676502;
//   atomic_mass_map[601500000] = 149.920902249;
//   atomic_mass_map[611500000] = 149.920990941;
//   atomic_mass_map[621500000] = 149.917282919;
//   atomic_mass_map[631500000] = 149.919707671;
//   atomic_mass_map[641500000] = 149.918664422;
//   atomic_mass_map[651500000] = 149.923664941;
//   atomic_mass_map[661500000] = 149.925593264;
//   atomic_mass_map[671500000] = 149.933498408;
//   atomic_mass_map[681500000] = 149.937915567;
//   atomic_mass_map[691500000] = 149.95009;
//   atomic_mass_map[701500000] = 149.95852;
//   atomic_mass_map[711500000] = 149.973548;
//   atomic_mass_map[551510000] = 150.96258;
//   atomic_mass_map[561510000] = 150.95127;
//   atomic_mass_map[571510000] = 150.94232;
//   atomic_mass_map[581510000] = 150.9342722;
//   atomic_mass_map[591510000] = 150.928309285;
//   atomic_mass_map[601510000] = 150.923840289;
//   atomic_mass_map[611510000] = 150.921217539;
//   atomic_mass_map[621510000] = 150.919939796;
//   atomic_mass_map[631510000] = 150.919857803;
//   atomic_mass_map[641510000] = 150.92035595;
//   atomic_mass_map[651510000] = 150.923109599;
//   atomic_mass_map[661510000] = 150.926191564;
//   atomic_mass_map[671510000] = 150.931698345;
//   atomic_mass_map[681510000] = 150.937448567;
//   atomic_mass_map[691510000] = 150.945487875;
//   atomic_mass_map[701510000] = 150.955402497;
//   atomic_mass_map[711510000] = 150.967677;
//   atomic_mass_map[561520000] = 151.95481;
//   atomic_mass_map[571520000] = 151.94682;
//   atomic_mass_map[581520000] = 151.9366;
//   atomic_mass_map[591520000] = 151.9315529;
//   atomic_mass_map[601520000] = 151.924692216;
//   atomic_mass_map[611520000] = 151.923506181;
//   atomic_mass_map[621520000] = 151.919739721;
//   atomic_mass_map[631520000] = 151.921752184;
//   atomic_mass_map[641520000] = 151.919799494;
//   atomic_mass_map[651520000] = 151.924082936;
//   atomic_mass_map[661520000] = 151.924725286;
//   atomic_mass_map[671520000] = 151.931723623;
//   atomic_mass_map[681520000] = 151.935057085;
//   atomic_mass_map[691520000] = 151.944422;
//   atomic_mass_map[701520000] = 151.9502727;
//   atomic_mass_map[711520000] = 151.96412;
//   atomic_mass_map[561530000] = 152.96036;
//   atomic_mass_map[571530000] = 152.95036;
//   atomic_mass_map[581530000] = 152.94093;
//   atomic_mass_map[591530000] = 152.933903539;
//   atomic_mass_map[601530000] = 152.927717978;
//   atomic_mass_map[611530000] = 152.924156686;
//   atomic_mass_map[621530000] = 152.92210465;
//   atomic_mass_map[631530000] = 152.921238003;
//   atomic_mass_map[641530000] = 152.921758027;
//   atomic_mass_map[651530000] = 152.923442403;
//   atomic_mass_map[661530000] = 152.925772378;
//   atomic_mass_map[671530000] = 152.930206429;
//   atomic_mass_map[681530000] = 152.93508044;
//   atomic_mass_map[691530000] = 152.942040101;
//   atomic_mass_map[701530000] = 152.94932;
//   atomic_mass_map[711530000] = 152.958751054;
//   atomic_mass_map[721530000] = 152.97069;
//   atomic_mass_map[571540000] = 153.95517;
//   atomic_mass_map[581540000] = 153.9438;
//   atomic_mass_map[591540000] = 153.937525741;
//   atomic_mass_map[601540000] = 153.929484894;
//   atomic_mass_map[611540000] = 153.926471531;
//   atomic_mass_map[621540000] = 153.922216861;
//   atomic_mass_map[631540000] = 153.922986962;
//   atomic_mass_map[641540000] = 153.92087406;
//   atomic_mass_map[651540000] = 153.924684767;
//   atomic_mass_map[661540000] = 153.924429277;
//   atomic_mass_map[671540000] = 153.930606834;
//   atomic_mass_map[681540000] = 153.932790842;
//   atomic_mass_map[691540000] = 153.941570033;
//   atomic_mass_map[701540000] = 153.946395655;
//   atomic_mass_map[711540000] = 153.957364;
//   atomic_mass_map[721540000] = 153.96486;
//   atomic_mass_map[571550000] = 154.95901;
//   atomic_mass_map[581550000] = 154.94855;
//   atomic_mass_map[591550000] = 154.940509281;
//   atomic_mass_map[601550000] = 154.933135691;
//   atomic_mass_map[611550000] = 154.928137048;
//   atomic_mass_map[621550000] = 154.924647748;
//   atomic_mass_map[631550000] = 154.922901107;
//   atomic_mass_map[641550000] = 154.922630473;
//   atomic_mass_map[651550000] = 154.923510547;
//   atomic_mass_map[661550000] = 154.925759086;
//   atomic_mass_map[671550000] = 154.929104049;
//   atomic_mass_map[681550000] = 154.933215911;
//   atomic_mass_map[691550000] = 154.93920966;
//   atomic_mass_map[701550000] = 154.945783132;
//   atomic_mass_map[711550000] = 154.9543206;
//   atomic_mass_map[721550000] = 154.96311;
//   atomic_mass_map[731550000] = 154.974245;
//   atomic_mass_map[581560000] = 155.95183;
//   atomic_mass_map[591560000] = 155.94464;
//   atomic_mass_map[601560000] = 155.935078894;
//   atomic_mass_map[611560000] = 155.931117516;
//   atomic_mass_map[621560000] = 155.925536067;
//   atomic_mass_map[631560000] = 155.924760494;
//   atomic_mass_map[641560000] = 155.922131241;
//   atomic_mass_map[651560000] = 155.924755181;
//   atomic_mass_map[661560000] = 155.924284713;
//   atomic_mass_map[671560000] = 155.929706112;
//   atomic_mass_map[681560000] = 155.931067313;
//   atomic_mass_map[691560000] = 155.938991573;
//   atomic_mass_map[701560000] = 155.942824698;
//   atomic_mass_map[711560000] = 155.953032522;
//   atomic_mass_map[721560000] = 155.959347805;
//   atomic_mass_map[731560000] = 155.97203;
//   atomic_mass_map[581570000] = 156.95705;
//   atomic_mass_map[591570000] = 156.94789;
//   atomic_mass_map[601570000] = 156.939386061;
//   atomic_mass_map[611570000] = 156.933121393;
//   atomic_mass_map[621570000] = 156.928418698;
//   atomic_mass_map[631570000] = 156.925433446;
//   atomic_mass_map[641570000] = 156.923968569;
//   atomic_mass_map[651570000] = 156.924033028;
//   atomic_mass_map[661570000] = 156.92547066;
//   atomic_mass_map[671570000] = 156.928254427;
//   atomic_mass_map[681570000] = 156.931948658;
//   atomic_mass_map[691570000] = 156.93694412;
//   atomic_mass_map[701570000] = 156.942645349;
//   atomic_mass_map[711570000] = 156.950126667;
//   atomic_mass_map[721570000] = 156.958236;
//   atomic_mass_map[731570000] = 156.968176167;
//   atomic_mass_map[741570000] = 156.97884;
//   atomic_mass_map[591580000] = 157.95241;
//   atomic_mass_map[601580000] = 157.94197;
//   atomic_mass_map[611580000] = 157.936565144;
//   atomic_mass_map[621580000] = 157.929951004;
//   atomic_mass_map[631580000] = 157.927798606;
//   atomic_mass_map[641580000] = 157.924112348;
//   atomic_mass_map[651580000] = 157.925420947;
//   atomic_mass_map[661580000] = 157.924415875;
//   atomic_mass_map[671580000] = 157.928945969;
//   atomic_mass_map[681580000] = 157.929893474;
//   atomic_mass_map[691580000] = 157.936979525;
//   atomic_mass_map[701580000] = 157.939870549;
//   atomic_mass_map[711580000] = 157.949315507;
//   atomic_mass_map[721580000] = 157.954801092;
//   atomic_mass_map[731580000] = 157.966541;
//   atomic_mass_map[741580000] = 157.974562;
//   atomic_mass_map[591590000] = 158.95589;
//   atomic_mass_map[601590000] = 158.94653;
//   atomic_mass_map[611590000] = 158.939286502;
//   atomic_mass_map[621590000] = 158.933217226;
//   atomic_mass_map[631590000] = 158.929100067;
//   atomic_mass_map[641590000] = 158.926396969;
//   atomic_mass_map[651590000] = 158.92535471;
//   atomic_mass_map[661590000] = 158.925746958;
//   atomic_mass_map[671590000] = 158.927719703;
//   atomic_mass_map[681590000] = 158.93069181;
//   atomic_mass_map[691590000] = 158.934975;
//   atomic_mass_map[701590000] = 158.940054623;
//   atomic_mass_map[711590000] = 158.946635615;
//   atomic_mass_map[721590000] = 158.953995669;
//   atomic_mass_map[731590000] = 158.963022556;
//   atomic_mass_map[741590000] = 158.972638;
//   atomic_mass_map[751590000] = 158.984179;
//   atomic_mass_map[601600000] = 159.9494;
//   atomic_mass_map[611600000] = 159.9431;
//   atomic_mass_map[621600000] = 159.935335311;
//   atomic_mass_map[631600000] = 159.931850939;
//   atomic_mass_map[641600000] = 159.927062411;
//   atomic_mass_map[651600000] = 159.927175556;
//   atomic_mass_map[661600000] = 159.925204646;
//   atomic_mass_map[671600000] = 159.928736606;
//   atomic_mass_map[681600000] = 159.92907713;
//   atomic_mass_map[691600000] = 159.935263106;
//   atomic_mass_map[701600000] = 159.937557406;
//   atomic_mass_map[711600000] = 159.946033;
//   atomic_mass_map[721600000] = 159.950690749;
//   atomic_mass_map[731600000] = 159.96148751;
//   atomic_mass_map[741600000] = 159.968462584;
//   atomic_mass_map[751600000] = 159.981823;
//   atomic_mass_map[601610000] = 160.95428;
//   atomic_mass_map[611610000] = 160.94607;
//   atomic_mass_map[621610000] = 160.93916017;
//   atomic_mass_map[631610000] = 160.933664091;
//   atomic_mass_map[641610000] = 160.929677476;
//   atomic_mass_map[651610000] = 160.927577825;
//   atomic_mass_map[661610000] = 160.926940492;
//   atomic_mass_map[671610000] = 160.927861547;
//   atomic_mass_map[681610000] = 160.930004599;
//   atomic_mass_map[691610000] = 160.933549;
//   atomic_mass_map[701610000] = 160.937907138;
//   atomic_mass_map[711610000] = 160.943572;
//   atomic_mass_map[721610000] = 160.95027837;
//   atomic_mass_map[731610000] = 160.958452265;
//   atomic_mass_map[741610000] = 160.967197;
//   atomic_mass_map[751610000] = 160.977572951;
//   atomic_mass_map[761610000] = 160.989029;
//   atomic_mass_map[611620000] = 161.95022;
//   atomic_mass_map[621620000] = 161.94146;
//   atomic_mass_map[631620000] = 161.936988761;
//   atomic_mass_map[641620000] = 161.930993017;
//   atomic_mass_map[651620000] = 161.92949536;
//   atomic_mass_map[661620000] = 161.926805573;
//   atomic_mass_map[671620000] = 161.929102274;
//   atomic_mass_map[681620000] = 161.928788364;
//   atomic_mass_map[691620000] = 161.934002277;
//   atomic_mass_map[701620000] = 161.935773604;
//   atomic_mass_map[711620000] = 161.943282776;
//   atomic_mass_map[721620000] = 161.947214837;
//   atomic_mass_map[731620000] = 161.957294083;
//   atomic_mass_map[741620000] = 161.963499142;
//   atomic_mass_map[751620000] = 161.975844;
//   atomic_mass_map[761620000] = 161.984431;
//   atomic_mass_map[611630000] = 162.95357;
//   atomic_mass_map[621630000] = 162.94555;
//   atomic_mass_map[631630000] = 162.939195675;
//   atomic_mass_map[641630000] = 162.934176855;
//   atomic_mass_map[651630000] = 162.930654659;
//   atomic_mass_map[661630000] = 162.928738284;
//   atomic_mass_map[671630000] = 162.928741027;
//   atomic_mass_map[681630000] = 162.930040797;
//   atomic_mass_map[691630000] = 162.932659172;
//   atomic_mass_map[701630000] = 162.936339632;
//   atomic_mass_map[711630000] = 162.941179;
//   atomic_mass_map[721630000] = 162.947112946;
//   atomic_mass_map[731630000] = 162.95433711;
//   atomic_mass_map[741630000] = 162.962524342;
//   atomic_mass_map[751630000] = 162.97207986;
//   atomic_mass_map[761630000] = 162.98241;
//   atomic_mass_map[621640000] = 163.94836;
//   atomic_mass_map[631640000] = 163.942744;
//   atomic_mass_map[641640000] = 163.93583;
//   atomic_mass_map[651640000] = 163.933357961;
//   atomic_mass_map[661640000] = 163.929181874;
//   atomic_mass_map[671640000] = 163.930240273;
//   atomic_mass_map[681640000] = 163.929208791;
//   atomic_mass_map[691640000] = 163.933543614;
//   atomic_mass_map[701640000] = 163.934494934;
//   atomic_mass_map[711640000] = 163.941339;
//   atomic_mass_map[721640000] = 163.944370845;
//   atomic_mass_map[731640000] = 163.953534;
//   atomic_mass_map[741640000] = 163.958960683;
//   atomic_mass_map[751640000] = 163.97045287;
//   atomic_mass_map[761640000] = 163.978021712;
//   atomic_mass_map[771640000] = 163.991909;
//   atomic_mass_map[621650000] = 164.95297;
//   atomic_mass_map[631650000] = 164.945587;
//   atomic_mass_map[641650000] = 164.93936;
//   atomic_mass_map[651650000] = 164.93498;
//   atomic_mass_map[661650000] = 164.931710456;
//   atomic_mass_map[671650000] = 164.930328835;
//   atomic_mass_map[681650000] = 164.930734496;
//   atomic_mass_map[691650000] = 164.93244314;
//   atomic_mass_map[701650000] = 164.935270241;
//   atomic_mass_map[711650000] = 164.939406758;
//   atomic_mass_map[721650000] = 164.944567;
//   atomic_mass_map[731650000] = 164.950780572;
//   atomic_mass_map[741650000] = 164.958280788;
//   atomic_mass_map[751650000] = 164.96710266;
//   atomic_mass_map[761650000] = 164.976602;
//   atomic_mass_map[771650000] = 164.987501;
//   atomic_mass_map[631660000] = 165.94962;
//   atomic_mass_map[641660000] = 165.94146;
//   atomic_mass_map[651660000] = 165.937859521;
//   atomic_mass_map[661660000] = 165.932813863;
//   atomic_mass_map[671660000] = 165.932290927;
//   atomic_mass_map[681660000] = 165.93029953;
//   atomic_mass_map[691660000] = 165.9335606;
//   atomic_mass_map[701660000] = 165.933874749;
//   atomic_mass_map[711660000] = 165.939859;
//   atomic_mass_map[721660000] = 165.94218;
//   atomic_mass_map[731660000] = 165.950512;
//   atomic_mass_map[741660000] = 165.95503121;
//   atomic_mass_map[751660000] = 165.965760821;
//   atomic_mass_map[761660000] = 165.972692476;
//   atomic_mass_map[771660000] = 165.985664;
//   atomic_mass_map[781660000] = 165.994855;
//   atomic_mass_map[631670000] = 166.95289;
//   atomic_mass_map[641670000] = 166.94545;
//   atomic_mass_map[651670000] = 166.93996;
//   atomic_mass_map[661670000] = 166.935661379;
//   atomic_mass_map[671670000] = 166.93313855;
//   atomic_mass_map[681670000] = 166.932054617;
//   atomic_mass_map[691670000] = 166.932856192;
//   atomic_mass_map[701670000] = 166.934952988;
//   atomic_mass_map[711670000] = 166.93827;
//   atomic_mass_map[721670000] = 166.9426;
//   atomic_mass_map[731670000] = 166.948093;
//   atomic_mass_map[741670000] = 166.95480456;
//   atomic_mass_map[751670000] = 166.962595;
//   atomic_mass_map[761670000] = 166.971548768;
//   atomic_mass_map[771670000] = 166.981666315;
//   atomic_mass_map[781670000] = 166.992695;
//   atomic_mass_map[641680000] = 167.94808;
//   atomic_mass_map[651680000] = 167.9434;
//   atomic_mass_map[661680000] = 167.937133273;
//   atomic_mass_map[671680000] = 167.935522173;
//   atomic_mass_map[681680000] = 167.932376688;
//   atomic_mass_map[691680000] = 167.934177409;
//   atomic_mass_map[701680000] = 167.933889602;
//   atomic_mass_map[711680000] = 167.938735635;
//   atomic_mass_map[721680000] = 167.940568;
//   atomic_mass_map[731680000] = 167.948047;
//   atomic_mass_map[741680000] = 167.951805537;
//   atomic_mass_map[751680000] = 167.961572608;
//   atomic_mass_map[761680000] = 167.967807751;
//   atomic_mass_map[771680000] = 167.979906642;
//   atomic_mass_map[781680000] = 167.988128665;
//   atomic_mass_map[641690000] = 168.9526;
//   atomic_mass_map[651690000] = 168.94597;
//   atomic_mass_map[661690000] = 168.940313531;
//   atomic_mass_map[671690000] = 168.936878189;
//   atomic_mass_map[681690000] = 168.93459685;
//   atomic_mass_map[691690000] = 168.934217889;
//   atomic_mass_map[701690000] = 168.935182512;
//   atomic_mass_map[711690000] = 168.937644149;
//   atomic_mass_map[721690000] = 168.941259;
//   atomic_mass_map[731690000] = 168.946011;
//   atomic_mass_map[741690000] = 168.951778845;
//   atomic_mass_map[751690000] = 168.958766063;
//   atomic_mass_map[761690000] = 168.967017574;
//   atomic_mass_map[771690000] = 168.976298488;
//   atomic_mass_map[781690000] = 168.986567;
//   atomic_mass_map[791690000] = 168.99808;
//   atomic_mass_map[651700000] = 169.94984;
//   atomic_mass_map[661700000] = 169.94239;
//   atomic_mass_map[671700000] = 169.939624846;
//   atomic_mass_map[681700000] = 169.93547023;
//   atomic_mass_map[691700000] = 169.935806032;
//   atomic_mass_map[701700000] = 169.934766376;
//   atomic_mass_map[711700000] = 169.938478365;
//   atomic_mass_map[721700000] = 169.939609;
//   atomic_mass_map[731700000] = 169.946175;
//   atomic_mass_map[741700000] = 169.949231588;
//   atomic_mass_map[751700000] = 169.958220071;
//   atomic_mass_map[761700000] = 169.963578454;
//   atomic_mass_map[771700000] = 169.974922;
//   atomic_mass_map[781700000] = 169.982496345;
//   atomic_mass_map[791700000] = 169.995972;
//   atomic_mass_map[651710000] = 170.95273;
//   atomic_mass_map[661710000] = 170.94612;
//   atomic_mass_map[671710000] = 170.941471022;
//   atomic_mass_map[681710000] = 170.938035681;
//   atomic_mass_map[691710000] = 170.936433871;
//   atomic_mass_map[701710000] = 170.936330208;
//   atomic_mass_map[711710000] = 170.93791696;
//   atomic_mass_map[721710000] = 170.940492;
//   atomic_mass_map[731710000] = 170.944476;
//   atomic_mass_map[741710000] = 170.949451;
//   atomic_mass_map[751710000] = 170.955716;
//   atomic_mass_map[761710000] = 170.963173969;
//   atomic_mass_map[771710000] = 170.971640252;
//   atomic_mass_map[781710000] = 170.981245248;
//   atomic_mass_map[791710000] = 170.991875791;
//   atomic_mass_map[801710000] = 171.003529;
//   atomic_mass_map[661720000] = 171.94846;
//   atomic_mass_map[671720000] = 171.94473;
//   atomic_mass_map[681720000] = 171.939361858;
//   atomic_mass_map[691720000] = 171.938405521;
//   atomic_mass_map[701720000] = 171.936385872;
//   atomic_mass_map[711720000] = 171.939089103;
//   atomic_mass_map[721720000] = 171.939449716;
//   atomic_mass_map[731720000] = 171.944895;
//   atomic_mass_map[741720000] = 171.947292;
//   atomic_mass_map[751720000] = 171.955419665;
//   atomic_mass_map[761720000] = 171.960017317;
//   atomic_mass_map[771720000] = 171.970607036;
//   atomic_mass_map[781720000] = 171.977350921;
//   atomic_mass_map[791720000] = 171.989942284;
//   atomic_mass_map[801720000] = 171.998808967;
//   atomic_mass_map[661730000] = 172.95283;
//   atomic_mass_map[671730000] = 172.94702;
//   atomic_mass_map[681730000] = 172.9424;
//   atomic_mass_map[691730000] = 172.939608371;
//   atomic_mass_map[701730000] = 172.938215136;
//   atomic_mass_map[711730000] = 172.938934029;
//   atomic_mass_map[721730000] = 172.940513;
//   atomic_mass_map[731730000] = 172.94375;
//   atomic_mass_map[741730000] = 172.947689;
//   atomic_mass_map[751730000] = 172.953243;
//   atomic_mass_map[761730000] = 172.959808465;
//   atomic_mass_map[771730000] = 172.967505636;
//   atomic_mass_map[781730000] = 172.976443058;
//   atomic_mass_map[791730000] = 172.986240924;
//   atomic_mass_map[801730000] = 172.997091;
//   atomic_mass_map[671740000] = 173.95095;
//   atomic_mass_map[681740000] = 173.94423;
//   atomic_mass_map[691740000] = 173.942172953;
//   atomic_mass_map[701740000] = 173.938866437;
//   atomic_mass_map[711740000] = 173.940340854;
//   atomic_mass_map[721740000] = 173.940046141;
//   atomic_mass_map[731740000] = 173.944454;
//   atomic_mass_map[741740000] = 173.946079;
//   atomic_mass_map[751740000] = 173.953115;
//   atomic_mass_map[761740000] = 173.957063507;
//   atomic_mass_map[771740000] = 173.966861045;
//   atomic_mass_map[781740000] = 173.972819528;
//   atomic_mass_map[791740000] = 173.984717;
//   atomic_mass_map[801740000] = 173.992864748;
//   atomic_mass_map[671750000] = 174.95362;
//   atomic_mass_map[681750000] = 174.94777;
//   atomic_mass_map[691750000] = 174.9438412;
//   atomic_mass_map[701750000] = 174.941280797;
//   atomic_mass_map[711750000] = 174.940775191;
//   atomic_mass_map[721750000] = 174.941509187;
//   atomic_mass_map[731750000] = 174.943737;
//   atomic_mass_map[741750000] = 174.946717;
//   atomic_mass_map[751750000] = 174.951381;
//   atomic_mass_map[761750000] = 174.956945313;
//   atomic_mass_map[771750000] = 174.964149589;
//   atomic_mass_map[781750000] = 174.972409721;
//   atomic_mass_map[791750000] = 174.981303712;
//   atomic_mass_map[801750000] = 174.991440747;
//   atomic_mass_map[681760000] = 175.94994;
//   atomic_mass_map[691760000] = 175.946999449;
//   atomic_mass_map[701760000] = 175.942576447;
//   atomic_mass_map[711760000] = 175.94268968;
//   atomic_mass_map[721760000] = 175.941407628;
//   atomic_mass_map[731760000] = 175.944857;
//   atomic_mass_map[741760000] = 175.945634;
//   atomic_mass_map[751760000] = 175.951623;
//   atomic_mass_map[761760000] = 175.954806;
//   atomic_mass_map[771760000] = 175.963650389;
//   atomic_mass_map[781760000] = 175.968938362;
//   atomic_mass_map[791760000] = 175.980250432;
//   atomic_mass_map[801760000] = 175.987360863;
//   atomic_mass_map[811760000] = 176.000624028;
//   atomic_mass_map[681770000] = 176.95399;
//   atomic_mass_map[691770000] = 176.94904;
//   atomic_mass_map[701770000] = 176.945265586;
//   atomic_mass_map[711770000] = 176.943761525;
//   atomic_mass_map[721770000] = 176.943227717;
//   atomic_mass_map[731770000] = 176.944479469;
//   atomic_mass_map[741770000] = 176.946643;
//   atomic_mass_map[751770000] = 176.950328;
//   atomic_mass_map[761770000] = 176.954965628;
//   atomic_mass_map[771770000] = 176.9613015;
//   atomic_mass_map[781770000] = 176.968469537;
//   atomic_mass_map[791770000] = 176.976870439;
//   atomic_mass_map[801770000] = 176.986277319;
//   atomic_mass_map[811770000] = 176.996430829;
//   atomic_mass_map[691780000] = 177.95264;
//   atomic_mass_map[701780000] = 177.946651444;
//   atomic_mass_map[711780000] = 177.945958018;
//   atomic_mass_map[721780000] = 177.943705833;
//   atomic_mass_map[731780000] = 177.945678;
//   atomic_mass_map[741780000] = 177.945883303;
//   atomic_mass_map[751780000] = 177.950989;
//   atomic_mass_map[761780000] = 177.953253627;
//   atomic_mass_map[771780000] = 177.961082;
//   atomic_mass_map[781780000] = 177.965649523;
//   atomic_mass_map[791780000] = 177.97603192;
//   atomic_mass_map[801780000] = 177.982483769;
//   atomic_mass_map[811780000] = 177.994854;
//   atomic_mass_map[821780000] = 178.003831243;
//   atomic_mass_map[691790000] = 178.95534;
//   atomic_mass_map[701790000] = 178.95004;
//   atomic_mass_map[711790000] = 178.947330935;
//   atomic_mass_map[721790000] = 178.945823212;
//   atomic_mass_map[731790000] = 178.945936555;
//   atomic_mass_map[741790000] = 178.947077043;
//   atomic_mass_map[751790000] = 178.94998913;
//   atomic_mass_map[761790000] = 178.953816881;
//   atomic_mass_map[771790000] = 178.959120178;
//   atomic_mass_map[781790000] = 178.965358953;
//   atomic_mass_map[791790000] = 178.973173654;
//   atomic_mass_map[801790000] = 178.981831329;
//   atomic_mass_map[811790000] = 178.991110696;
//   atomic_mass_map[821790000] = 179.002201112;
//   atomic_mass_map[701800000] = 179.95212;
//   atomic_mass_map[711800000] = 179.94988825;
//   atomic_mass_map[721800000] = 179.946557042;
//   atomic_mass_map[731800000] = 179.947464832;
//   atomic_mass_map[741800000] = 179.946710805;
//   atomic_mass_map[751800000] = 179.950791568;
//   atomic_mass_map[761800000] = 179.952375485;
//   atomic_mass_map[771800000] = 179.959229446;
//   atomic_mass_map[781800000] = 179.963031955;
//   atomic_mass_map[791800000] = 179.972523397;
//   atomic_mass_map[801800000] = 179.978260335;
//   atomic_mass_map[811800000] = 179.990056524;
//   atomic_mass_map[821800000] = 179.997928286;
//   atomic_mass_map[701810000] = 180.95589;
//   atomic_mass_map[711810000] = 180.951908;
//   atomic_mass_map[721810000] = 180.949108338;
//   atomic_mass_map[731810000] = 180.947995769;
//   atomic_mass_map[741810000] = 180.948197783;
//   atomic_mass_map[751810000] = 180.950057673;
//   atomic_mass_map[761810000] = 180.953247188;
//   atomic_mass_map[771810000] = 180.957625456;
//   atomic_mass_map[781810000] = 180.963097588;
//   atomic_mass_map[791810000] = 180.970079047;
//   atomic_mass_map[801810000] = 180.977819371;
//   atomic_mass_map[811810000] = 180.986259979;
//   atomic_mass_map[821810000] = 180.99665331;
//   atomic_mass_map[711820000] = 181.95504;
//   atomic_mass_map[721820000] = 181.950561185;
//   atomic_mass_map[731820000] = 181.950151853;
//   atomic_mass_map[741820000] = 181.948203945;
//   atomic_mass_map[751820000] = 181.951209869;
//   atomic_mass_map[761820000] = 181.952110187;
//   atomic_mass_map[771820000] = 181.958076296;
//   atomic_mass_map[781820000] = 181.961171823;
//   atomic_mass_map[791820000] = 181.969617874;
//   atomic_mass_map[801820000] = 181.974689351;
//   atomic_mass_map[811820000] = 181.985713159;
//   atomic_mass_map[821820000] = 181.992672466;
//   atomic_mass_map[711830000] = 182.957363;
//   atomic_mass_map[721830000] = 182.953530444;
//   atomic_mass_map[731830000] = 182.95137262;
//   atomic_mass_map[741830000] = 182.950222748;
//   atomic_mass_map[751830000] = 182.950819638;
//   atomic_mass_map[761830000] = 182.953124719;
//   atomic_mass_map[771830000] = 182.956839968;
//   atomic_mass_map[781830000] = 182.961596804;
//   atomic_mass_map[791830000] = 182.967590635;
//   atomic_mass_map[801830000] = 182.974444783;
//   atomic_mass_map[811830000] = 182.982192839;
//   atomic_mass_map[821830000] = 182.991872098;
//   atomic_mass_map[711840000] = 183.96091;
//   atomic_mass_map[721840000] = 183.955446243;
//   atomic_mass_map[731840000] = 183.954007694;
//   atomic_mass_map[741840000] = 183.950930916;
//   atomic_mass_map[751840000] = 183.952522809;
//   atomic_mass_map[761840000] = 183.952488536;
//   atomic_mass_map[771840000] = 183.957476;
//   atomic_mass_map[781840000] = 183.959915113;
//   atomic_mass_map[791840000] = 183.967451524;
//   atomic_mass_map[801840000] = 183.971713528;
//   atomic_mass_map[811840000] = 183.981885851;
//   atomic_mass_map[821840000] = 183.988135701;
//   atomic_mass_map[831840000] = 184.001274756;
//   atomic_mass_map[711850000] = 184.96362;
//   atomic_mass_map[721850000] = 184.958862;
//   atomic_mass_map[731850000] = 184.955559084;
//   atomic_mass_map[741850000] = 184.953418974;
//   atomic_mass_map[751850000] = 184.952954486;
//   atomic_mass_map[761850000] = 184.954041741;
//   atomic_mass_map[771850000] = 184.956698;
//   atomic_mass_map[781850000] = 184.960613659;
//   atomic_mass_map[791850000] = 184.965789569;
//   atomic_mass_map[801850000] = 184.971899388;
//   atomic_mass_map[811850000] = 184.97878905;
//   atomic_mass_map[821850000] = 184.987610004;
//   atomic_mass_map[831850000] = 184.9976;
//   atomic_mass_map[721860000] = 185.960897;
//   atomic_mass_map[731860000] = 185.958550667;
//   atomic_mass_map[741860000] = 185.954362771;
//   atomic_mass_map[751860000] = 185.954985595;
//   atomic_mass_map[761860000] = 185.953835044;
//   atomic_mass_map[771860000] = 185.957944153;
//   atomic_mass_map[781860000] = 185.959350812;
//   atomic_mass_map[791860000] = 185.965952703;
//   atomic_mass_map[801860000] = 185.969362272;
//   atomic_mass_map[811860000] = 185.978650841;
//   atomic_mass_map[821860000] = 185.984238331;
//   atomic_mass_map[831860000] = 185.996643681;
//   atomic_mass_map[841860000] = 186.004393241;
//   atomic_mass_map[721870000] = 186.96477;
//   atomic_mass_map[731870000] = 186.960386;
//   atomic_mass_map[741870000] = 186.957158841;
//   atomic_mass_map[751870000] = 186.955750071;
//   atomic_mass_map[761870000] = 186.955747422;
//   atomic_mass_map[771870000] = 186.957542;
//   atomic_mass_map[781870000] = 186.960616976;
//   atomic_mass_map[791870000] = 186.964543155;
//   atomic_mass_map[801870000] = 186.969814428;
//   atomic_mass_map[811870000] = 186.975906348;
//   atomic_mass_map[821870000] = 186.983910895;
//   atomic_mass_map[831870000] = 186.993147184;
//   atomic_mass_map[841870000] = 187.003041054;
//   atomic_mass_map[721880000] = 187.96685;
//   atomic_mass_map[731880000] = 187.963916;
//   atomic_mass_map[741880000] = 187.958486177;
//   atomic_mass_map[751880000] = 187.95811151;
//   atomic_mass_map[761880000] = 187.955835174;
//   atomic_mass_map[771880000] = 187.958828095;
//   atomic_mass_map[781880000] = 187.959388889;
//   atomic_mass_map[791880000] = 187.965349392;
//   atomic_mass_map[801880000] = 187.967566887;
//   atomic_mass_map[811880000] = 187.976020886;
//   atomic_mass_map[821880000] = 187.980874815;
//   atomic_mass_map[831880000] = 187.992286943;
//   atomic_mass_map[841880000] = 187.999415569;
//   atomic_mass_map[721890000] = 188.97084;
//   atomic_mass_map[731890000] = 188.96583;
//   atomic_mass_map[741890000] = 188.961763;
//   atomic_mass_map[751890000] = 188.95922602;
//   atomic_mass_map[761890000] = 188.958144162;
//   atomic_mass_map[771890000] = 188.958715028;
//   atomic_mass_map[781890000] = 188.960831282;
//   atomic_mass_map[791890000] = 188.963948286;
//   atomic_mass_map[801890000] = 188.968194927;
//   atomic_mass_map[811890000] = 188.973588471;
//   atomic_mass_map[821890000] = 188.980807;
//   atomic_mass_map[831890000] = 188.989194984;
//   atomic_mass_map[841890000] = 188.998473429;
//   atomic_mass_map[731900000] = 189.96939;
//   atomic_mass_map[741900000] = 189.963090589;
//   atomic_mass_map[751900000] = 189.96174426;
//   atomic_mass_map[761900000] = 189.958443702;
//   atomic_mass_map[771900000] = 189.960541215;
//   atomic_mass_map[781900000] = 189.959929707;
//   atomic_mass_map[791900000] = 189.96469839;
//   atomic_mass_map[801900000] = 189.966322735;
//   atomic_mass_map[811900000] = 189.973828;
//   atomic_mass_map[821900000] = 189.978081999;
//   atomic_mass_map[831900000] = 189.988621561;
//   atomic_mass_map[841900000] = 189.995100569;
//   atomic_mass_map[731910000] = 190.97156;
//   atomic_mass_map[741910000] = 190.966531;
//   atomic_mass_map[751910000] = 190.963121551;
//   atomic_mass_map[761910000] = 190.960926361;
//   atomic_mass_map[771910000] = 190.960589293;
//   atomic_mass_map[781910000] = 190.961672912;
//   atomic_mass_map[791910000] = 190.963702248;
//   atomic_mass_map[801910000] = 190.96715716;
//   atomic_mass_map[811910000] = 190.971784199;
//   atomic_mass_map[821910000] = 190.978276;
//   atomic_mass_map[831910000] = 190.985786556;
//   atomic_mass_map[841910000] = 190.994558462;
//   atomic_mass_map[851910000] = 191.004147909;
//   atomic_mass_map[731920000] = 191.97514;
//   atomic_mass_map[741920000] = 191.96817;
//   atomic_mass_map[751920000] = 191.966088;
//   atomic_mass_map[761920000] = 191.961476998;
//   atomic_mass_map[771920000] = 191.962600247;
//   atomic_mass_map[781920000] = 191.961038746;
//   atomic_mass_map[791920000] = 191.964813694;
//   atomic_mass_map[801920000] = 191.965634702;
//   atomic_mass_map[811920000] = 191.972225;
//   atomic_mass_map[821920000] = 191.975775008;
//   atomic_mass_map[831920000] = 191.985469058;
//   atomic_mass_map[841920000] = 191.991335925;
//   atomic_mass_map[851920000] = 192.003151707;
//   atomic_mass_map[741930000] = 192.97178;
//   atomic_mass_map[751930000] = 192.967541;
//   atomic_mass_map[761930000] = 192.96414787;
//   atomic_mass_map[771930000] = 192.962921587;
//   atomic_mass_map[781930000] = 192.96298238;
//   atomic_mass_map[791930000] = 192.964137257;
//   atomic_mass_map[801930000] = 192.966652658;
//   atomic_mass_map[811930000] = 192.970501998;
//   atomic_mass_map[821930000] = 192.976173234;
//   atomic_mass_map[831930000] = 192.982959813;
//   atomic_mass_map[841930000] = 192.991025678;
//   atomic_mass_map[851930000] = 192.999927418;
//   atomic_mass_map[861930000] = 193.009707979;
//   atomic_mass_map[741940000] = 193.97367;
//   atomic_mass_map[751940000] = 193.97076;
//   atomic_mass_map[761940000] = 193.96517724;
//   atomic_mass_map[771940000] = 193.965073536;
//   atomic_mass_map[781940000] = 193.96268085;
//   atomic_mass_map[791940000] = 193.965417754;
//   atomic_mass_map[801940000] = 193.965449112;
//   atomic_mass_map[811940000] = 193.971081412;
//   atomic_mass_map[821940000] = 193.97401225;
//   atomic_mass_map[831940000] = 193.982785;
//   atomic_mass_map[841940000] = 193.988186101;
//   atomic_mass_map[851940000] = 193.999236139;
//   atomic_mass_map[861940000] = 194.006144389;
//   atomic_mass_map[751950000] = 194.97254;
//   atomic_mass_map[761950000] = 194.968318;
//   atomic_mass_map[771950000] = 194.96597473;
//   atomic_mass_map[781950000] = 194.964791719;
//   atomic_mass_map[791950000] = 194.965035225;
//   atomic_mass_map[801950000] = 194.96672054;
//   atomic_mass_map[811950000] = 194.969774268;
//   atomic_mass_map[821950000] = 194.974542922;
//   atomic_mass_map[831950000] = 194.980648781;
//   atomic_mass_map[841950000] = 194.988125532;
//   atomic_mass_map[851950000] = 194.996268546;
//   atomic_mass_map[861950000] = 195.005421673;
//   atomic_mass_map[751960000] = 195.9758;
//   atomic_mass_map[761960000] = 195.969641463;
//   atomic_mass_map[771960000] = 195.968397112;
//   atomic_mass_map[781960000] = 195.964952091;
//   atomic_mass_map[791960000] = 195.966569908;
//   atomic_mass_map[801960000] = 195.96583256;
//   atomic_mass_map[811960000] = 195.970481193;
//   atomic_mass_map[821960000] = 195.972774233;
//   atomic_mass_map[831960000] = 195.980666509;
//   atomic_mass_map[841960000] = 195.985525902;
//   atomic_mass_map[851960000] = 195.995799619;
//   atomic_mass_map[861960000] = 196.002115998;
//   atomic_mass_map[751970000] = 196.97799;
//   atomic_mass_map[761970000] = 196.97283;
//   atomic_mass_map[771970000] = 196.969655415;
//   atomic_mass_map[781970000] = 196.967340687;
//   atomic_mass_map[791970000] = 196.966568786;
//   atomic_mass_map[801970000] = 196.967212847;
//   atomic_mass_map[811970000] = 196.96957589;
//   atomic_mass_map[821970000] = 196.973431166;
//   atomic_mass_map[831970000] = 196.978865099;
//   atomic_mass_map[841970000] = 196.985659522;
//   atomic_mass_map[851970000] = 196.993189187;
//   atomic_mass_map[861970000] = 197.00158462;
//   atomic_mass_map[751980000] = 197.9816;
//   atomic_mass_map[761980000] = 197.97441;
//   atomic_mass_map[771980000] = 197.97228;
//   atomic_mass_map[781980000] = 197.96789492;
//   atomic_mass_map[791980000] = 197.96824242;
//   atomic_mass_map[801980000] = 197.966768602;
//   atomic_mass_map[811980000] = 197.970483065;
//   atomic_mass_map[821980000] = 197.972034077;
//   atomic_mass_map[831980000] = 197.979206;
//   atomic_mass_map[841980000] = 197.983389132;
//   atomic_mass_map[851980000] = 197.992784;
//   atomic_mass_map[861980000] = 197.998679156;
//   atomic_mass_map[761990000] = 198.97801;
//   atomic_mass_map[771990000] = 198.973805301;
//   atomic_mass_map[781990000] = 198.970595224;
//   atomic_mass_map[791990000] = 198.968765282;
//   atomic_mass_map[801990000] = 198.968280643;
//   atomic_mass_map[811990000] = 198.969877;
//   atomic_mass_map[821990000] = 198.97291268;
//   atomic_mass_map[831990000] = 198.97767301;
//   atomic_mass_map[841990000] = 198.983667115;
//   atomic_mass_map[851990000] = 198.990527654;
//   atomic_mass_map[861990000] = 198.998389856;
//   atomic_mass_map[871990000] = 199.007258594;
//   atomic_mass_map[762000000] = 199.97984;
//   atomic_mass_map[772000000] = 199.9768;
//   atomic_mass_map[782000000] = 199.971442807;
//   atomic_mass_map[792000000] = 199.970756456;
//   atomic_mass_map[802000000] = 199.96832659;
//   atomic_mass_map[812000000] = 199.970963258;
//   atomic_mass_map[822000000] = 199.971818893;
//   atomic_mass_map[832000000] = 199.978131179;
//   atomic_mass_map[842000000] = 199.98179879;
//   atomic_mass_map[852000000] = 199.990351015;
//   atomic_mass_map[862000000] = 199.995690431;
//   atomic_mass_map[872000000] = 200.006586003;
//   atomic_mass_map[762010000] = 200.98364;
//   atomic_mass_map[772010000] = 200.97864;
//   atomic_mass_map[782010000] = 200.974513112;
//   atomic_mass_map[792010000] = 200.971657484;
//   atomic_mass_map[802010000] = 200.970302839;
//   atomic_mass_map[812010000] = 200.970822212;
//   atomic_mass_map[822010000] = 200.972882917;
//   atomic_mass_map[832010000] = 200.977010416;
//   atomic_mass_map[842010000] = 200.982259806;
//   atomic_mass_map[852010000] = 200.988417147;
//   atomic_mass_map[862010000] = 200.99562801;
//   atomic_mass_map[872010000] = 201.003866593;
//   atomic_mass_map[882010000] = 201.012712;
//   atomic_mass_map[762020000] = 201.98595;
//   atomic_mass_map[772020000] = 201.98199;
//   atomic_mass_map[782020000] = 201.975639;
//   atomic_mass_map[792020000] = 201.973856;
//   atomic_mass_map[802020000] = 201.9706434;
//   atomic_mass_map[812020000] = 201.972102441;
//   atomic_mass_map[822020000] = 201.972152026;
//   atomic_mass_map[832020000] = 201.977733861;
//   atomic_mass_map[842020000] = 201.980757657;
//   atomic_mass_map[852020000] = 201.988630302;
//   atomic_mass_map[862020000] = 201.993264277;
//   atomic_mass_map[872020000] = 202.00332;
//   atomic_mass_map[882020000] = 202.009759744;
//   atomic_mass_map[772030000] = 202.98423;
//   atomic_mass_map[782030000] = 202.97893;
//   atomic_mass_map[792030000] = 202.975154446;
//   atomic_mass_map[802030000] = 202.972872837;
//   atomic_mass_map[812030000] = 202.972344593;
//   atomic_mass_map[822030000] = 202.973391126;
//   atomic_mass_map[832030000] = 202.976892779;
//   atomic_mass_map[842030000] = 202.981416133;
//   atomic_mass_map[852030000] = 202.98694299;
//   atomic_mass_map[862030000] = 202.993387741;
//   atomic_mass_map[872030000] = 203.000940729;
//   atomic_mass_map[882030000] = 203.009304274;
//   atomic_mass_map[772040000] = 203.9896;
//   atomic_mass_map[782040000] = 203.98076;
//   atomic_mass_map[792040000] = 203.977831;
//   atomic_mass_map[802040000] = 203.973493981;
//   atomic_mass_map[812040000] = 203.9738639;
//   atomic_mass_map[822040000] = 203.973043981;
//   atomic_mass_map[832040000] = 203.977836062;
//   atomic_mass_map[842040000] = 203.980310339;
//   atomic_mass_map[852040000] = 203.987251198;
//   atomic_mass_map[862040000] = 203.991429922;
//   atomic_mass_map[872040000] = 204.000652059;
//   atomic_mass_map[882040000] = 204.006491866;
//   atomic_mass_map[782050000] = 204.98608;
//   atomic_mass_map[792050000] = 204.97985;
//   atomic_mass_map[802050000] = 204.976073417;
//   atomic_mass_map[812050000] = 204.974427801;
//   atomic_mass_map[822050000] = 204.974482157;
//   atomic_mass_map[832050000] = 204.977386694;
//   atomic_mass_map[842050000] = 204.981203091;
//   atomic_mass_map[852050000] = 204.986075861;
//   atomic_mass_map[862050000] = 204.99171884;
//   atomic_mass_map[872050000] = 204.998593858;
//   atomic_mass_map[882050000] = 205.006268245;
//   atomic_mass_map[782060000] = 205.98966;
//   atomic_mass_map[792060000] = 205.98474;
//   atomic_mass_map[802060000] = 205.977514401;
//   atomic_mass_map[812060000] = 205.976110589;
//   atomic_mass_map[822060000] = 205.974465683;
//   atomic_mass_map[832060000] = 205.978499317;
//   atomic_mass_map[842060000] = 205.980473991;
//   atomic_mass_map[852060000] = 205.986656933;
//   atomic_mass_map[862060000] = 205.990214219;
//   atomic_mass_map[872060000] = 205.998666132;
//   atomic_mass_map[882060000] = 206.003828053;
//   atomic_mass_map[892060000] = 206.014452;
//   atomic_mass_map[792070000] = 206.9884;
//   atomic_mass_map[802070000] = 206.9823;
//   atomic_mass_map[812070000] = 206.977419712;
//   atomic_mass_map[822070000] = 206.975897297;
//   atomic_mass_map[832070000] = 206.978471022;
//   atomic_mass_map[842070000] = 206.981593762;
//   atomic_mass_map[852070000] = 206.985800332;
//   atomic_mass_map[862070000] = 206.990730254;
//   atomic_mass_map[872070000] = 206.996946478;
//   atomic_mass_map[882070000] = 207.003799159;
//   atomic_mass_map[892070000] = 207.01196583;
//   atomic_mass_map[792080000] = 207.99345;
//   atomic_mass_map[802080000] = 207.985759;
//   atomic_mass_map[812080000] = 207.982019001;
//   atomic_mass_map[822080000] = 207.976652481;
//   atomic_mass_map[832080000] = 207.97974253;
//   atomic_mass_map[842080000] = 207.981246092;
//   atomic_mass_map[852080000] = 207.986613302;
//   atomic_mass_map[862080000] = 207.989634687;
//   atomic_mass_map[872080000] = 207.997138132;
//   atomic_mass_map[882080000] = 208.001841122;
//   atomic_mass_map[892080000] = 208.011550406;
//   atomic_mass_map[902080000] = 208.017900275;
//   atomic_mass_map[792090000] = 208.99735;
//   atomic_mass_map[802090000] = 208.99072;
//   atomic_mass_map[812090000] = 208.985359353;
//   atomic_mass_map[822090000] = 208.981090461;
//   atomic_mass_map[832090000] = 208.980399068;
//   atomic_mass_map[842090000] = 208.982430836;
//   atomic_mass_map[852090000] = 208.986170215;
//   atomic_mass_map[862090000] = 208.99041451;
//   atomic_mass_map[872090000] = 208.995954932;
//   atomic_mass_map[882090000] = 209.001990455;
//   atomic_mass_map[892090000] = 209.009494762;
//   atomic_mass_map[902090000] = 209.017752974;
//   atomic_mass_map[792100000] = 210.0025;
//   atomic_mass_map[802100000] = 209.99424;
//   atomic_mass_map[812100000] = 209.990073958;
//   atomic_mass_map[822100000] = 209.984188861;
//   atomic_mass_map[832100000] = 209.984120705;
//   atomic_mass_map[842100000] = 209.982874076;
//   atomic_mass_map[852100000] = 209.987147898;
//   atomic_mass_map[862100000] = 209.989689107;
//   atomic_mass_map[872100000] = 209.996422409;
//   atomic_mass_map[882100000] = 210.000494169;
//   atomic_mass_map[892100000] = 210.009436052;
//   atomic_mass_map[902100000] = 210.015093642;
//   atomic_mass_map[802110000] = 210.99933;
//   atomic_mass_map[812110000] = 210.993475;
//   atomic_mass_map[822110000] = 210.988737125;
//   atomic_mass_map[832110000] = 210.987269741;
//   atomic_mass_map[842110000] = 210.986653562;
//   atomic_mass_map[852110000] = 210.987496614;
//   atomic_mass_map[862110000] = 210.99060111;
//   atomic_mass_map[872110000] = 210.995555723;
//   atomic_mass_map[882110000] = 211.000893214;
//   atomic_mass_map[892110000] = 211.007731899;
//   atomic_mass_map[902110000] = 211.014929468;
//   atomic_mass_map[802120000] = 212.00296;
//   atomic_mass_map[812120000] = 211.998335;
//   atomic_mass_map[822120000] = 211.991897703;
//   atomic_mass_map[832120000] = 211.991286026;
//   atomic_mass_map[842120000] = 211.988868376;
//   atomic_mass_map[852120000] = 211.990737688;
//   atomic_mass_map[862120000] = 211.990703919;
//   atomic_mass_map[872120000] = 211.996225652;
//   atomic_mass_map[882120000] = 211.999786715;
//   atomic_mass_map[892120000] = 212.007813171;
//   atomic_mass_map[902120000] = 212.012987595;
//   atomic_mass_map[912120000] = 212.023202993;
//   atomic_mass_map[802130000] = 213.00823;
//   atomic_mass_map[812130000] = 213.001915;
//   atomic_mass_map[822130000] = 212.996562894;
//   atomic_mass_map[832130000] = 212.994385067;
//   atomic_mass_map[842130000] = 212.992857623;
//   atomic_mass_map[852130000] = 212.992936978;
//   atomic_mass_map[862130000] = 212.993883069;
//   atomic_mass_map[872130000] = 212.996186041;
//   atomic_mass_map[882130000] = 213.000384258;
//   atomic_mass_map[892130000] = 213.00660902;
//   atomic_mass_map[902130000] = 213.013009223;
//   atomic_mass_map[912130000] = 213.021109239;
//   atomic_mass_map[802140000] = 214.012;
//   atomic_mass_map[812140000] = 214.00694;
//   atomic_mass_map[822140000] = 213.999805917;
//   atomic_mass_map[832140000] = 213.998711926;
//   atomic_mass_map[842140000] = 213.995201683;
//   atomic_mass_map[852140000] = 213.996372066;
//   atomic_mass_map[862140000] = 213.995362956;
//   atomic_mass_map[872140000] = 213.998971333;
//   atomic_mass_map[882140000] = 214.000099727;
//   atomic_mass_map[892140000] = 214.006918488;
//   atomic_mass_map[902140000] = 214.011500159;
//   atomic_mass_map[912140000] = 214.020918483;
//   atomic_mass_map[802150000] = 215.0174;
//   atomic_mass_map[812150000] = 215.01064;
//   atomic_mass_map[822150000] = 215.004741;
//   atomic_mass_map[832150000] = 215.001769823;
//   atomic_mass_map[842150000] = 214.999420146;
//   atomic_mass_map[852150000] = 214.998652849;
//   atomic_mass_map[862150000] = 214.99874589;
//   atomic_mass_map[872150000] = 215.000341838;
//   atomic_mass_map[882150000] = 215.00272042;
//   atomic_mass_map[892150000] = 215.006474511;
//   atomic_mass_map[902150000] = 215.011724774;
//   atomic_mass_map[912150000] = 215.019182929;
//   atomic_mass_map[802160000] = 216.02132;
//   atomic_mass_map[812160000] = 216.0158;
//   atomic_mass_map[822160000] = 216.00803;
//   atomic_mass_map[832160000] = 216.006305989;
//   atomic_mass_map[842160000] = 216.001915186;
//   atomic_mass_map[852160000] = 216.002423555;
//   atomic_mass_map[862160000] = 216.000271859;
//   atomic_mass_map[872160000] = 216.003189865;
//   atomic_mass_map[882160000] = 216.003533423;
//   atomic_mass_map[892160000] = 216.008743482;
//   atomic_mass_map[902160000] = 216.011056007;
//   atomic_mass_map[912160000] = 216.019108913;
//   atomic_mass_map[812170000] = 217.01966;
//   atomic_mass_map[822170000] = 217.01314;
//   atomic_mass_map[832170000] = 217.009372;
//   atomic_mass_map[842170000] = 217.006318159;
//   atomic_mass_map[852170000] = 217.004719222;
//   atomic_mass_map[862170000] = 217.003928018;
//   atomic_mass_map[872170000] = 217.004632282;
//   atomic_mass_map[882170000] = 217.006320726;
//   atomic_mass_map[892170000] = 217.009343872;
//   atomic_mass_map[902170000] = 217.013116701;
//   atomic_mass_map[912170000] = 217.018325363;
//   atomic_mass_map[922170000] = 217.024661;
//   atomic_mass_map[812180000] = 218.02479;
//   atomic_mass_map[822180000] = 218.01659;
//   atomic_mass_map[832180000] = 218.014188;
//   atomic_mass_map[842180000] = 218.008973546;
//   atomic_mass_map[852180000] = 218.008694723;
//   atomic_mass_map[862180000] = 218.005601586;
//   atomic_mass_map[872180000] = 218.007578653;
//   atomic_mass_map[882180000] = 218.007140631;
//   atomic_mass_map[892180000] = 218.01164164;
//   atomic_mass_map[902180000] = 218.013276331;
//   atomic_mass_map[912180000] = 218.020058579;
//   atomic_mass_map[922180000] = 218.023523472;
//   atomic_mass_map[822190000] = 219.02177;
//   atomic_mass_map[832190000] = 219.01748;
//   atomic_mass_map[842190000] = 219.013614;
//   atomic_mass_map[852190000] = 219.011161848;
//   atomic_mass_map[862190000] = 219.009480361;
//   atomic_mass_map[872190000] = 219.009252427;
//   atomic_mass_map[882190000] = 219.010085484;
//   atomic_mass_map[892190000] = 219.01242073;
//   atomic_mass_map[902190000] = 219.015537481;
//   atomic_mass_map[912190000] = 219.019904029;
//   atomic_mass_map[922190000] = 219.02499913;
//   atomic_mass_map[932190000] = 219.03143;
//   atomic_mass_map[822200000] = 220.02541;
//   atomic_mass_map[832200000] = 220.02235;
//   atomic_mass_map[842200000] = 220.016386;
//   atomic_mass_map[852200000] = 220.015433;
//   atomic_mass_map[862200000] = 220.01139413;
//   atomic_mass_map[872200000] = 220.012327702;
//   atomic_mass_map[882200000] = 220.011025872;
//   atomic_mass_map[892200000] = 220.014754854;
//   atomic_mass_map[902200000] = 220.015748148;
//   atomic_mass_map[912200000] = 220.021705;
//   atomic_mass_map[922200000] = 220.02462;
//   atomic_mass_map[932200000] = 220.03254;
//   atomic_mass_map[832210000] = 221.02587;
//   atomic_mass_map[842210000] = 221.021228;
//   atomic_mass_map[852210000] = 221.018017;
//   atomic_mass_map[862210000] = 221.015537076;
//   atomic_mass_map[872210000] = 221.014255163;
//   atomic_mass_map[882210000] = 221.01391768;
//   atomic_mass_map[892210000] = 221.015591579;
//   atomic_mass_map[902210000] = 221.018184072;
//   atomic_mass_map[912210000] = 221.021874941;
//   atomic_mass_map[922210000] = 221.026284;
//   atomic_mass_map[932210000] = 221.03204;
//   atomic_mass_map[832220000] = 222.03078;
//   atomic_mass_map[842220000] = 222.02414;
//   atomic_mass_map[852220000] = 222.022494;
//   atomic_mass_map[862220000] = 222.017578246;
//   atomic_mass_map[872220000] = 222.017552192;
//   atomic_mass_map[882220000] = 222.015374846;
//   atomic_mass_map[892220000] = 222.017844181;
//   atomic_mass_map[902220000] = 222.018468521;
//   atomic_mass_map[912220000] = 222.023784;
//   atomic_mass_map[922220000] = 222.026003;
//   atomic_mass_map[932220000] = 222.0333;
//   atomic_mass_map[832230000] = 223.0345;
//   atomic_mass_map[842230000] = 223.02907;
//   atomic_mass_map[852230000] = 223.025151;
//   atomic_mass_map[862230000] = 223.021889286;
//   atomic_mass_map[872230000] = 223.019736013;
//   atomic_mass_map[882230000] = 223.018502327;
//   atomic_mass_map[892230000] = 223.019137747;
//   atomic_mass_map[902230000] = 223.020811853;
//   atomic_mass_map[912230000] = 223.023962614;
//   atomic_mass_map[922230000] = 223.027739187;
//   atomic_mass_map[932230000] = 223.03285;
//   atomic_mass_map[832240000] = 224.03947;
//   atomic_mass_map[842240000] = 224.03211;
//   atomic_mass_map[852240000] = 224.029749;
//   atomic_mass_map[862240000] = 224.024095805;
//   atomic_mass_map[872240000] = 224.023398;
//   atomic_mass_map[882240000] = 224.020211968;
//   atomic_mass_map[892240000] = 224.021723163;
//   atomic_mass_map[902240000] = 224.021464382;
//   atomic_mass_map[912240000] = 224.025617614;
//   atomic_mass_map[922240000] = 224.027605163;
//   atomic_mass_map[932240000] = 224.03422;
//   atomic_mass_map[842250000] = 225.03707;
//   atomic_mass_map[852250000] = 225.03263;
//   atomic_mass_map[862250000] = 225.028485574;
//   atomic_mass_map[872250000] = 225.025572682;
//   atomic_mass_map[882250000] = 225.023611857;
//   atomic_mass_map[892250000] = 225.023229987;
//   atomic_mass_map[902250000] = 225.023951363;
//   atomic_mass_map[912250000] = 225.026131009;
//   atomic_mass_map[922250000] = 225.02939135;
//   atomic_mass_map[932250000] = 225.033910892;
//   atomic_mass_map[842260000] = 226.04031;
//   atomic_mass_map[852260000] = 226.03716;
//   atomic_mass_map[862260000] = 226.030861382;
//   atomic_mass_map[872260000] = 226.029566;
//   atomic_mass_map[882260000] = 226.02541033;
//   atomic_mass_map[892260000] = 226.026098383;
//   atomic_mass_map[902260000] = 226.024903383;
//   atomic_mass_map[912260000] = 226.027948082;
//   atomic_mass_map[922260000] = 226.029339101;
//   atomic_mass_map[932260000] = 226.035188;
//   atomic_mass_map[842270000] = 227.04539;
//   atomic_mass_map[852270000] = 227.04024;
//   atomic_mass_map[862270000] = 227.035304396;
//   atomic_mass_map[872270000] = 227.031869;
//   atomic_mass_map[882270000] = 227.029178349;
//   atomic_mass_map[892270000] = 227.027752283;
//   atomic_mass_map[902270000] = 227.027704227;
//   atomic_mass_map[912270000] = 227.028805351;
//   atomic_mass_map[922270000] = 227.031156773;
//   atomic_mass_map[932270000] = 227.034957129;
//   atomic_mass_map[852280000] = 228.04475;
//   atomic_mass_map[862280000] = 228.037835418;
//   atomic_mass_map[872280000] = 228.035823175;
//   atomic_mass_map[882280000] = 228.031070728;
//   atomic_mass_map[892280000] = 228.031021547;
//   atomic_mass_map[902280000] = 228.028741272;
//   atomic_mass_map[912280000] = 228.031051673;
//   atomic_mass_map[922280000] = 228.031371492;
//   atomic_mass_map[932280000] = 228.036066866;
//   atomic_mass_map[942280000] = 228.038732491;
//   atomic_mass_map[852290000] = 229.04812;
//   atomic_mass_map[862290000] = 229.042257277;
//   atomic_mass_map[872290000] = 229.038298;
//   atomic_mass_map[882290000] = 229.034942262;
//   atomic_mass_map[892290000] = 229.032955884;
//   atomic_mass_map[902290000] = 229.031762713;
//   atomic_mass_map[912290000] = 229.032097236;
//   atomic_mass_map[922290000] = 229.03350628;
//   atomic_mass_map[932290000] = 229.036264139;
//   atomic_mass_map[942290000] = 229.040143614;
//   atomic_mass_map[862300000] = 230.04514;
//   atomic_mass_map[872300000] = 230.042416349;
//   atomic_mass_map[882300000] = 230.037054781;
//   atomic_mass_map[892300000] = 230.036327;
//   atomic_mass_map[902300000] = 230.03313413;
//   atomic_mass_map[912300000] = 230.034541047;
//   atomic_mass_map[922300000] = 230.033940096;
//   atomic_mass_map[932300000] = 230.037827926;
//   atomic_mass_map[942300000] = 230.039650283;
//   atomic_mass_map[952300000] = 230.04609;
//   atomic_mass_map[862310000] = 231.04987;
//   atomic_mass_map[872310000] = 231.045158;
//   atomic_mass_map[882310000] = 231.041027087;
//   atomic_mass_map[892310000] = 231.038393;
//   atomic_mass_map[902310000] = 231.036304628;
//   atomic_mass_map[912310000] = 231.035884243;
//   atomic_mass_map[922310000] = 231.036293861;
//   atomic_mass_map[932310000] = 231.038245364;
//   atomic_mass_map[942310000] = 231.041101511;
//   atomic_mass_map[952310000] = 231.04556;
//   atomic_mass_map[872320000] = 232.049368;
//   atomic_mass_map[882320000] = 232.04347527;
//   atomic_mass_map[892320000] = 232.042034;
//   atomic_mass_map[902320000] = 232.03805576;
//   atomic_mass_map[912320000] = 232.038591737;
//   atomic_mass_map[922320000] = 232.037156297;
//   atomic_mass_map[932320000] = 232.040109;
//   atomic_mass_map[942320000] = 232.041184582;
//   atomic_mass_map[952320000] = 232.04645;
//   atomic_mass_map[962320000] = 232.049817;
//   atomic_mass_map[872330000] = 233.05264;
//   atomic_mass_map[882330000] = 233.047582;
//   atomic_mass_map[892330000] = 233.044346;
//   atomic_mass_map[902330000] = 233.041582278;
//   atomic_mass_map[912330000] = 233.040247222;
//   atomic_mass_map[922330000] = 233.039635525;
//   atomic_mass_map[932330000] = 233.040740989;
//   atomic_mass_map[942330000] = 233.042997716;
//   atomic_mass_map[952330000] = 233.046445;
//   atomic_mass_map[962330000] = 233.050770001;
//   atomic_mass_map[882340000] = 234.050342;
//   atomic_mass_map[892340000] = 234.048139;
//   atomic_mass_map[902340000] = 234.043601407;
//   atomic_mass_map[912340000] = 234.043307161;
//   atomic_mass_map[922340000] = 234.040952306;
//   atomic_mass_map[932340000] = 234.042895256;
//   atomic_mass_map[942340000] = 234.043317387;
//   atomic_mass_map[952340000] = 234.047731;
//   atomic_mass_map[962340000] = 234.050160237;
//   atomic_mass_map[972340000] = 234.057267;
//   atomic_mass_map[882350000] = 235.05497;
//   atomic_mass_map[892350000] = 235.05084;
//   atomic_mass_map[902350000] = 235.047255;
//   atomic_mass_map[912350000] = 235.045399;
//   atomic_mass_map[922350000] = 235.043930131;
//   atomic_mass_map[932350000] = 235.044063487;
//   atomic_mass_map[942350000] = 235.045286206;
//   atomic_mass_map[952350000] = 235.047908245;
//   atomic_mass_map[962350000] = 235.051542;
//   atomic_mass_map[972350000] = 235.05658;
//   atomic_mass_map[892360000] = 236.054988;
//   atomic_mass_map[902360000] = 236.049657;
//   atomic_mass_map[912360000] = 236.048668;
//   atomic_mass_map[922360000] = 236.04556821;
//   atomic_mass_map[932360000] = 236.046569744;
//   atomic_mass_map[942360000] = 236.046058109;
//   atomic_mass_map[952360000] = 236.049428;
//   atomic_mass_map[962360000] = 236.051374477;
//   atomic_mass_map[972360000] = 236.05748;
//   atomic_mass_map[892370000] = 237.05827;
//   atomic_mass_map[902370000] = 237.053629;
//   atomic_mass_map[912370000] = 237.051023;
//   atomic_mass_map[922370000] = 237.048730378;
//   atomic_mass_map[932370000] = 237.048173649;
//   atomic_mass_map[942370000] = 237.04840983;
//   atomic_mass_map[952370000] = 237.049996;
//   atomic_mass_map[962370000] = 237.052869294;
//   atomic_mass_map[972370000] = 237.0571;
//   atomic_mass_map[982370000] = 237.062197789;
//   atomic_mass_map[902380000] = 238.056496;
//   atomic_mass_map[912380000] = 238.054637;
//   atomic_mass_map[922380000] = 238.050788423;
//   atomic_mass_map[932380000] = 238.050946611;
//   atomic_mass_map[942380000] = 238.049560111;
//   atomic_mass_map[952380000] = 238.051984542;
//   atomic_mass_map[962380000] = 238.05308142;
//   atomic_mass_map[972380000] = 238.058204;
//   atomic_mass_map[982380000] = 238.06149;
//   atomic_mass_map[902390000] = 239.06077;
//   atomic_mass_map[912390000] = 239.05726;
//   atomic_mass_map[922390000] = 239.054293475;
//   atomic_mass_map[932390000] = 239.052939241;
//   atomic_mass_map[942390000] = 239.052163591;
//   atomic_mass_map[952390000] = 239.053024699;
//   atomic_mass_map[962390000] = 239.054910117;
//   atomic_mass_map[972390000] = 239.058241;
//   atomic_mass_map[982390000] = 239.062529;
//   atomic_mass_map[992390000] = 239.06823;
//   atomic_mass_map[912400000] = 240.06098;
//   atomic_mass_map[922400000] = 240.056593356;
//   atomic_mass_map[932400000] = 240.056165427;
//   atomic_mass_map[942400000] = 240.05381375;
//   atomic_mass_map[952400000] = 240.055300384;
//   atomic_mass_map[962400000] = 240.055529681;
//   atomic_mass_map[972400000] = 240.059759;
//   atomic_mass_map[982400000] = 240.062255728;
//   atomic_mass_map[992400000] = 240.06892;
//   atomic_mass_map[912410000] = 241.06408;
//   atomic_mass_map[922410000] = 241.06033;
//   atomic_mass_map[932410000] = 241.058252636;
//   atomic_mass_map[942410000] = 241.056851661;
//   atomic_mass_map[952410000] = 241.056829349;
//   atomic_mass_map[962410000] = 241.05765317;
//   atomic_mass_map[972410000] = 241.060155;
//   atomic_mass_map[982410000] = 241.06369;
//   atomic_mass_map[992410000] = 241.06856;
//   atomic_mass_map[1002410000] = 241.07421;
//   atomic_mass_map[922420000] = 242.062932;
//   atomic_mass_map[932420000] = 242.061641379;
//   atomic_mass_map[942420000] = 242.058742809;
//   atomic_mass_map[952420000] = 242.059549364;
//   atomic_mass_map[962420000] = 242.058836039;
//   atomic_mass_map[972420000] = 242.061982;
//   atomic_mass_map[982420000] = 242.063754274;
//   atomic_mass_map[992420000] = 242.069567;
//   atomic_mass_map[1002420000] = 242.07343;
//   atomic_mass_map[922430000] = 243.06699;
//   atomic_mass_map[932430000] = 243.06428;
//   atomic_mass_map[942430000] = 243.062003595;
//   atomic_mass_map[952430000] = 243.061381302;
//   atomic_mass_map[962430000] = 243.061389325;
//   atomic_mass_map[972430000] = 243.063007791;
//   atomic_mass_map[982430000] = 243.065477;
//   atomic_mass_map[992430000] = 243.06951;
//   atomic_mass_map[1002430000] = 243.074465;
//   atomic_mass_map[932440000] = 244.06785;
//   atomic_mass_map[942440000] = 244.06420526;
//   atomic_mass_map[952440000] = 244.064285052;
//   atomic_mass_map[962440000] = 244.062752783;
//   atomic_mass_map[972440000] = 244.065180979;
//   atomic_mass_map[982440000] = 244.066000801;
//   atomic_mass_map[992440000] = 244.070883;
//   atomic_mass_map[1002440000] = 244.074038;
//   atomic_mass_map[932450000] = 245.0708;
//   atomic_mass_map[942450000] = 245.067826058;
//   atomic_mass_map[952450000] = 245.066454834;
//   atomic_mass_map[962450000] = 245.065491454;
//   atomic_mass_map[972450000] = 245.066361821;
//   atomic_mass_map[982450000] = 245.068048677;
//   atomic_mass_map[992450000] = 245.071249;
//   atomic_mass_map[1002450000] = 245.075349;
//   atomic_mass_map[1012450000] = 245.080808;
//   atomic_mass_map[942460000] = 246.070205458;
//   atomic_mass_map[952460000] = 246.069775;
//   atomic_mass_map[962460000] = 246.067223841;
//   atomic_mass_map[972460000] = 246.068673126;
//   atomic_mass_map[982460000] = 246.068805531;
//   atomic_mass_map[992460000] = 246.072896;
//   atomic_mass_map[1002460000] = 246.07535047;
//   atomic_mass_map[1012460000] = 246.081713;
//   atomic_mass_map[942470000] = 247.07419;
//   atomic_mass_map[952470000] = 247.072093;
//   atomic_mass_map[962470000] = 247.070354131;
//   atomic_mass_map[972470000] = 247.070307302;
//   atomic_mass_map[982470000] = 247.070965463;
//   atomic_mass_map[992470000] = 247.073622017;
//   atomic_mass_map[1002470000] = 247.076945;
//   atomic_mass_map[1012470000] = 247.081522;
//   atomic_mass_map[952480000] = 248.075753;
//   atomic_mass_map[962480000] = 248.072349862;
//   atomic_mass_map[972480000] = 248.073088;
//   atomic_mass_map[982480000] = 248.072185066;
//   atomic_mass_map[992480000] = 248.075471;
//   atomic_mass_map[1002480000] = 248.077186463;
//   atomic_mass_map[1012480000] = 248.082823;
//   atomic_mass_map[1022480000] = 248.08655;
//   atomic_mass_map[952490000] = 249.07848;
//   atomic_mass_map[962490000] = 249.075954767;
//   atomic_mass_map[972490000] = 249.074987676;
//   atomic_mass_map[982490000] = 249.074853903;
//   atomic_mass_map[992490000] = 249.076411;
//   atomic_mass_map[1002490000] = 249.07892755;
//   atomic_mass_map[1012490000] = 249.082914;
//   atomic_mass_map[1022490000] = 249.087797;
//   atomic_mass_map[962500000] = 250.078358313;
//   atomic_mass_map[972500000] = 250.078316698;
//   atomic_mass_map[982500000] = 250.076406244;
//   atomic_mass_map[992500000] = 250.078612;
//   atomic_mass_map[1002500000] = 250.079521034;
//   atomic_mass_map[1012500000] = 250.084415;
//   atomic_mass_map[1022500000] = 250.087562;
//   atomic_mass_map[962510000] = 251.082286441;
//   atomic_mass_map[972510000] = 251.080762009;
//   atomic_mass_map[982510000] = 251.079588625;
//   atomic_mass_map[992510000] = 251.079993586;
//   atomic_mass_map[1002510000] = 251.08153989;
//   atomic_mass_map[1012510000] = 251.084774376;
//   atomic_mass_map[1022510000] = 251.088944;
//   atomic_mass_map[1032510000] = 251.09418;
//   atomic_mass_map[962520000] = 252.08487;
//   atomic_mass_map[972520000] = 252.084311;
//   atomic_mass_map[982520000] = 252.081627199;
//   atomic_mass_map[992520000] = 252.082979865;
//   atomic_mass_map[1002520000] = 252.08246706;
//   atomic_mass_map[1012520000] = 252.086432;
//   atomic_mass_map[1022520000] = 252.088966908;
//   atomic_mass_map[1032520000] = 252.095264;
//   atomic_mass_map[972530000] = 253.08688;
//   atomic_mass_map[982530000] = 253.085134499;
//   atomic_mass_map[992530000] = 253.084825715;
//   atomic_mass_map[1002530000] = 253.085184571;
//   atomic_mass_map[1012530000] = 253.087144;
//   atomic_mass_map[1022530000] = 253.090564103;
//   atomic_mass_map[1032530000] = 253.095091;
//   atomic_mass_map[1042530000] = 253.100438;
//   atomic_mass_map[972540000] = 254.0906;
//   atomic_mass_map[982540000] = 254.087324263;
//   atomic_mass_map[992540000] = 254.088022199;
//   atomic_mass_map[1002540000] = 254.086854397;
//   atomic_mass_map[1012540000] = 254.089592;
//   atomic_mass_map[1022540000] = 254.090955661;
//   atomic_mass_map[1032540000] = 254.096483;
//   atomic_mass_map[1042540000] = 254.100053;
//   atomic_mass_map[982550000] = 255.091048;
//   atomic_mass_map[992550000] = 255.090274958;
//   atomic_mass_map[1002550000] = 255.089964038;
//   atomic_mass_map[1012550000] = 255.091084149;
//   atomic_mass_map[1022550000] = 255.093191404;
//   atomic_mass_map[1032550000] = 255.096562404;
//   atomic_mass_map[1042550000] = 255.101268;
//   atomic_mass_map[1052550000] = 255.107068;
//   atomic_mass_map[982560000] = 256.093442;
//   atomic_mass_map[992560000] = 256.093599;
//   atomic_mass_map[1002560000] = 256.091774469;
//   atomic_mass_map[1012560000] = 256.093889;
//   atomic_mass_map[1022560000] = 256.09428287;
//   atomic_mass_map[1032560000] = 256.09849403;
//   atomic_mass_map[1042560000] = 256.101152218;
//   atomic_mass_map[1052560000] = 256.10789;
//   atomic_mass_map[992570000] = 257.095979;
//   atomic_mass_map[1002570000] = 257.095106078;
//   atomic_mass_map[1012570000] = 257.095542387;
//   atomic_mass_map[1022570000] = 257.09688783;
//   atomic_mass_map[1032570000] = 257.099418;
//   atomic_mass_map[1042570000] = 257.10291812;
//   atomic_mass_map[1052570000] = 257.107578;
//   atomic_mass_map[992580000] = 258.09952;
//   atomic_mass_map[1002580000] = 258.097077;
//   atomic_mass_map[1012580000] = 258.098431496;
//   atomic_mass_map[1022580000] = 258.098207;
//   atomic_mass_map[1032580000] = 258.101755;
//   atomic_mass_map[1042580000] = 258.103427679;
//   atomic_mass_map[1052580000] = 258.109285;
//   atomic_mass_map[1062580000] = 258.112983;
//   atomic_mass_map[1002590000] = 259.100597;
//   atomic_mass_map[1012590000] = 259.100511;
//   atomic_mass_map[1022590000] = 259.101032;
//   atomic_mass_map[1032590000] = 259.102902;
//   atomic_mass_map[1042590000] = 259.105596;
//   atomic_mass_map[1052590000] = 259.109491866;
//   atomic_mass_map[1062590000] = 259.114396;
//   atomic_mass_map[1002600000] = 260.102809;
//   atomic_mass_map[1012600000] = 260.103653;
//   atomic_mass_map[1022600000] = 260.102644;
//   atomic_mass_map[1032600000] = 260.105505;
//   atomic_mass_map[1042600000] = 260.106441;
//   atomic_mass_map[1052600000] = 260.111297;
//   atomic_mass_map[1062600000] = 260.114384105;
//   atomic_mass_map[1072600000] = 260.121659;
//   atomic_mass_map[1012610000] = 261.105828;
//   atomic_mass_map[1022610000] = 261.105697;
//   atomic_mass_map[1032610000] = 261.106884;
//   atomic_mass_map[1042610000] = 261.108773401;
//   atomic_mass_map[1052610000] = 261.111917;
//   atomic_mass_map[1062610000] = 261.115949461;
//   atomic_mass_map[1072610000] = 261.121455;
//   atomic_mass_map[1012620000] = 262.109101;
//   atomic_mass_map[1022620000] = 262.107464;
//   atomic_mass_map[1032620000] = 262.109612;
//   atomic_mass_map[1042620000] = 262.109925;
//   atomic_mass_map[1052620000] = 262.114072;
//   atomic_mass_map[1062620000] = 262.116336679;
//   atomic_mass_map[1072620000] = 262.122967;
//   atomic_mass_map[1022630000] = 263.110715;
//   atomic_mass_map[1032630000] = 263.111359;
//   atomic_mass_map[1042630000] = 263.112495;
//   atomic_mass_map[1052630000] = 263.11499;
//   atomic_mass_map[1062630000] = 263.118294;
//   atomic_mass_map[1072630000] = 263.122916;
//   atomic_mass_map[1082630000] = 263.128522;
//   atomic_mass_map[1022640000] = 264.112734;
//   atomic_mass_map[1032640000] = 264.114201;
//   atomic_mass_map[1042640000] = 264.113879;
//   atomic_mass_map[1052640000] = 264.117405;
//   atomic_mass_map[1062640000] = 264.118931;
//   atomic_mass_map[1072640000] = 264.124593;
//   atomic_mass_map[1082640000] = 264.128356917;
//   atomic_mass_map[1032650000] = 265.116193;
//   atomic_mass_map[1042650000] = 265.116684;
//   atomic_mass_map[1052650000] = 265.118613;
//   atomic_mass_map[1062650000] = 265.121093;
//   atomic_mass_map[1072650000] = 265.124914;
//   atomic_mass_map[1082650000] = 265.129792986;
//   atomic_mass_map[1092650000] = 265.135996;
//   atomic_mass_map[1032660000] = 266.119831;
//   atomic_mass_map[1042660000] = 266.118172;
//   atomic_mass_map[1052660000] = 266.121029;
//   atomic_mass_map[1062660000] = 266.121975;
//   atomic_mass_map[1072660000] = 266.126794;
//   atomic_mass_map[1082660000] = 266.1300464;
//   atomic_mass_map[1092660000] = 266.137374;
//   atomic_mass_map[1042670000] = 267.121789;
//   atomic_mass_map[1052670000] = 267.122465;
//   atomic_mass_map[1062670000] = 267.124357;
//   atomic_mass_map[1072670000] = 267.127501;
//   atomic_mass_map[1082670000] = 267.131673;
//   atomic_mass_map[1092670000] = 267.137189;
//   atomic_mass_map[1102670000] = 267.143768;
//   atomic_mass_map[1042680000] = 268.123968;
//   atomic_mass_map[1052680000] = 268.125671;
//   atomic_mass_map[1062680000] = 268.125392;
//   atomic_mass_map[1072680000] = 268.129692;
//   atomic_mass_map[1082680000] = 268.131865;
//   atomic_mass_map[1092680000] = 268.138649;
//   atomic_mass_map[1102680000] = 268.143478;
//   atomic_mass_map[1052690000] = 269.127911;
//   atomic_mass_map[1062690000] = 269.128627;
//   atomic_mass_map[1072690000] = 269.130416;
//   atomic_mass_map[1082690000] = 269.133753;
//   atomic_mass_map[1092690000] = 269.138822;
//   atomic_mass_map[1102690000] = 269.144752124;
//   atomic_mass_map[1052700000] = 270.131356;
//   atomic_mass_map[1062700000] = 270.130427;
//   atomic_mass_map[1072700000] = 270.133363;
//   atomic_mass_map[1082700000] = 270.13429;
//   atomic_mass_map[1092700000] = 270.140327;
//   atomic_mass_map[1102700000] = 270.144584153;
//   atomic_mass_map[1062710000] = 271.133933;
//   atomic_mass_map[1072710000] = 271.135256;
//   atomic_mass_map[1082710000] = 271.13717;
//   atomic_mass_map[1092710000] = 271.140744;
//   atomic_mass_map[1102710000] = 271.145946;
//   atomic_mass_map[1062720000] = 272.13589;
//   atomic_mass_map[1072720000] = 272.138264;
//   atomic_mass_map[1082720000] = 272.138495;
//   atomic_mass_map[1092720000] = 272.143406;
//   atomic_mass_map[1102720000] = 272.14602;
//   atomic_mass_map[1112720000] = 272.153273;
//   atomic_mass_map[1062730000] = 273.13958;
//   atomic_mass_map[1072730000] = 273.14024;
//   atomic_mass_map[1082730000] = 273.141679;
//   atomic_mass_map[1092730000] = 273.144399;
//   atomic_mass_map[1102730000] = 273.14856;
//   atomic_mass_map[1112730000] = 273.153127;
//   atomic_mass_map[1072740000] = 274.143548;
//   atomic_mass_map[1082740000] = 274.143304;
//   atomic_mass_map[1092740000] = 274.147245;
//   atomic_mass_map[1102740000] = 274.149411;
//   atomic_mass_map[1112740000] = 274.155253;
//   atomic_mass_map[1072750000] = 275.14567;
//   atomic_mass_map[1082750000] = 275.146668;
//   atomic_mass_map[1092750000] = 275.14882;
//   atomic_mass_map[1102750000] = 275.152033;
//   atomic_mass_map[1112750000] = 275.155939;
//   atomic_mass_map[1082760000] = 276.148455;
//   atomic_mass_map[1092760000] = 276.151594;
//   atomic_mass_map[1102760000] = 276.153025;
//   atomic_mass_map[1112760000] = 276.158334;
//   atomic_mass_map[1122760000] = 276.16141;
//   atomic_mass_map[1082770000] = 277.151899;
//   atomic_mass_map[1092770000] = 277.153268;
//   atomic_mass_map[1102770000] = 277.155914;
//   atomic_mass_map[1112770000] = 277.159069;
//   atomic_mass_map[1122770000] = 277.163641;
//   atomic_mass_map[1092780000] = 278.156307;
//   atomic_mass_map[1102780000] = 278.15704;
//   atomic_mass_map[1112780000] = 278.161493;
//   atomic_mass_map[1122780000] = 278.164156;
//   atomic_mass_map[1132780000] = 278.170578;
//   atomic_mass_map[1092790000] = 279.158075;
//   atomic_mass_map[1102790000] = 279.160097;
//   atomic_mass_map[1112790000] = 279.162722;
//   atomic_mass_map[1122790000] = 279.166542;
//   atomic_mass_map[1132790000] = 279.17095;
//   atomic_mass_map[1102800000] = 280.161311;
//   atomic_mass_map[1112800000] = 280.165138;
//   atomic_mass_map[1122800000] = 280.167147;
//   atomic_mass_map[1132800000] = 280.17293;
//   atomic_mass_map[1102810000] = 281.164511;
//   atomic_mass_map[1112810000] = 281.16636;
//   atomic_mass_map[1122810000] = 281.169746;
//   atomic_mass_map[1132810000] = 281.17348;
//   atomic_mass_map[1112820000] = 282.169119;
//   atomic_mass_map[1122820000] = 282.170496;
//   atomic_mass_map[1132820000] = 282.175672;
//   atomic_mass_map[1112830000] = 283.170544;
//   atomic_mass_map[1122830000] = 283.17327;
//   atomic_mass_map[1132830000] = 283.176571;
//   atomic_mass_map[1122840000] = 284.174156;
//   atomic_mass_map[1132840000] = 284.178727;
//   atomic_mass_map[1122850000] = 285.177117;
//   atomic_mass_map[1132850000] = 285.179727;
//   atomic_mass_map[1142850000] = 285.183643;
//   atomic_mass_map[1132860000] = 286.182208;
//   atomic_mass_map[1142860000] = 286.184235;
//   atomic_mass_map[1132870000] = 287.183389;
//   atomic_mass_map[1142870000] = 287.186783;
//   atomic_mass_map[1152870000] = 287.190704;
//   atomic_mass_map[1142880000] = 288.187572;
//   atomic_mass_map[1152880000] = 288.19274;
//   atomic_mass_map[1142890000] = 289.190419;
//   atomic_mass_map[1152890000] = 289.193627;
//   atomic_mass_map[1162890000] = 289.198162;
//   atomic_mass_map[1152900000] = 290.195975;
//   atomic_mass_map[1162900000] = 290.198638;
//   atomic_mass_map[1152910000] = 291.197071;
//   atomic_mass_map[1162910000] = 291.201077;
//   atomic_mass_map[1172910000] = 291.205535;
//   atomic_mass_map[1162920000] = 292.201742;
//   atomic_mass_map[1172920000] = 292.207463;
//   atomic_mass_map[1162930000] = 293.204487;
//   atomic_mass_map[1172930000] = 293.208236;
//   atomic_mass_map[1182930000] = 293.213562;
//   atomic_mass_map[1172940000] = 294.210462;
//   atomic_mass_map[1182940000] = 294.213921;
//   atomic_mass_map[1182950000] = 295.21624;
// }
//
// void pyne::_insert_abund_map() {
//   natural_abund_map[10010000] = 99.9885;
//   natural_abund_map[10020000] = 0.0115;
//   natural_abund_map[20030000] = 0.000134;
//   natural_abund_map[20040000] = 99.999866;
//   natural_abund_map[30060000] = 7.59;
//   natural_abund_map[30070000] = 92.41;
//   natural_abund_map[40090000] = 100.0;
//   natural_abund_map[50100000] = 19.9;
//   natural_abund_map[50110000] = 80.1;
//   natural_abund_map[60120000] = 98.93;
//   natural_abund_map[60130000] = 1.07;
//   natural_abund_map[70140000] = 99.636;
//   natural_abund_map[70150000] = 0.364;
//   natural_abund_map[80160000] = 99.757;
//   natural_abund_map[80170000] = 0.038;
//   natural_abund_map[80180000] = 0.205;
//   natural_abund_map[90190000] = 100.0;
//   natural_abund_map[100200000] = 90.48;
//   natural_abund_map[100210000] = 0.27;
//   natural_abund_map[100220000] = 9.25;
//   natural_abund_map[110230000] = 100.0;
//   natural_abund_map[120240000] = 78.99;
//   natural_abund_map[120250000] = 10.00;
//   natural_abund_map[120260000] = 11.01;
//   natural_abund_map[130270000] = 100.0;
//   natural_abund_map[140280000] = 92.223;
//   natural_abund_map[140290000] = 4.685;
//   natural_abund_map[140300000] = 3.092;
//   natural_abund_map[150310000] = 100.0;
//   natural_abund_map[160320000] = 94.99;
//   natural_abund_map[160330000] = 0.75;
//   natural_abund_map[160340000] = 4.25;
//   natural_abund_map[160360000] = 0.01;
//   natural_abund_map[170350000] = 75.76;
//   natural_abund_map[170370000] = 24.24;
//   natural_abund_map[180360000] = 0.3336;
//   natural_abund_map[180380000] = 0.0629;
//   natural_abund_map[180400000] = 99.6035;
//   natural_abund_map[190390000] = 93.2581;
//   natural_abund_map[190400000] = 0.0117;
//   natural_abund_map[190410000] = 6.7302;
//   natural_abund_map[200400000] = 96.941;
//   natural_abund_map[200420000] = 0.647;
//   natural_abund_map[200430000] = 0.135;
//   natural_abund_map[200440000] = 2.086;
//   natural_abund_map[200460000] = 0.004;
//   natural_abund_map[200480000] = 0.187;
//   natural_abund_map[210450000] = 100.0;
//   natural_abund_map[220460000] = 8.25;
//   natural_abund_map[220470000] = 7.44;
//   natural_abund_map[220480000] = 73.72;
//   natural_abund_map[220490000] = 5.41;
//   natural_abund_map[220500000] = 5.18;
//   natural_abund_map[230500000] = 0.250;
//   natural_abund_map[230510000] = 99.750;
//   natural_abund_map[240500000] = 4.345;
//   natural_abund_map[240520000] = 83.789;
//   natural_abund_map[240530000] = 9.501;
//   natural_abund_map[240540000] = 2.365;
//   natural_abund_map[250550000] = 100.0;
//   natural_abund_map[260540000] = 5.845;
//   natural_abund_map[260560000] = 91.754;
//   natural_abund_map[260570000] = 2.119;
//   natural_abund_map[260580000] = 0.282;
//   natural_abund_map[270590000] = 100.0;
//   natural_abund_map[280580000] = 68.077;
//   natural_abund_map[280600000] = 26.223;
//   natural_abund_map[280610000] = 1.1399;
//   natural_abund_map[280620000] = 3.6346;
//   natural_abund_map[280640000] = 0.9255;
//   natural_abund_map[290630000] = 69.15;
//   natural_abund_map[290650000] = 30.85;
//   natural_abund_map[300640000] = 49.17;
//   natural_abund_map[300660000] = 27.73;
//   natural_abund_map[300670000] = 4.04;
//   natural_abund_map[300680000] = 18.45;
//   natural_abund_map[300700000] = 0.61;
//   natural_abund_map[310690000] = 60.108;
//   natural_abund_map[310710000] = 39.892;
//   natural_abund_map[320700000] = 20.57;
//   natural_abund_map[320720000] = 27.45;
//   natural_abund_map[320730000] = 7.75;
//   natural_abund_map[320740000] = 36.50;
//   natural_abund_map[320760000] = 7.73;
//   natural_abund_map[330750000] = 100.0;
//   natural_abund_map[340740000] = 0.89;
//   natural_abund_map[340760000] = 9.37;
//   natural_abund_map[340770000] = 7.63;
//   natural_abund_map[340780000] = 23.77;
//   natural_abund_map[340800000] = 49.61;
//   natural_abund_map[340820000] = 8.73;
//   natural_abund_map[350790000] = 50.69;
//   natural_abund_map[350810000] = 49.31;
//   natural_abund_map[360780000] = 0.355;
//   natural_abund_map[360800000] = 2.286;
//   natural_abund_map[360820000] = 11.593;
//   natural_abund_map[360830000] = 11.500;
//   natural_abund_map[360840000] = 56.987;
//   natural_abund_map[360860000] = 17.279;
//   natural_abund_map[370850000] = 72.17;
//   natural_abund_map[370870000] = 27.83;
//   natural_abund_map[380840000] = 0.56;
//   natural_abund_map[380860000] = 9.86;
//   natural_abund_map[380870000] = 7.00;
//   natural_abund_map[380880000] = 82.58;
//   natural_abund_map[390890000] = 100.0;
//   natural_abund_map[400900000] = 51.45;
//   natural_abund_map[400910000] = 11.22;
//   natural_abund_map[400920000] = 17.15;
//   natural_abund_map[400940000] = 17.38;
//   natural_abund_map[400960000] = 2.80;
//   natural_abund_map[410930000] = 100.0;
//   natural_abund_map[420920000] = 14.53;
//   natural_abund_map[420940000] = 9.15;
//   natural_abund_map[420950000] = 15.84;
//   natural_abund_map[420960000] = 16.67;
//   natural_abund_map[420970000] = 9.60;
//   natural_abund_map[420980000] = 24.39;
//   natural_abund_map[421000000] = 9.82;
//   natural_abund_map[440960000] = 5.54;
//   natural_abund_map[440980000] = 1.87;
//   natural_abund_map[440990000] = 12.76;
//   natural_abund_map[441000000] = 12.60;
//   natural_abund_map[441010000] = 17.06;
//   natural_abund_map[441020000] = 31.55;
//   natural_abund_map[441040000] = 18.62;
//   natural_abund_map[451030000] = 100.0;
//   natural_abund_map[461020000] = 1.02;
//   natural_abund_map[461040000] = 11.14;
//   natural_abund_map[461050000] = 22.33;
//   natural_abund_map[461060000] = 27.33;
//   natural_abund_map[461080000] = 26.46;
//   natural_abund_map[461100000] = 11.72;
//   natural_abund_map[471070000] = 51.839;
//   natural_abund_map[471090000] = 48.161;
//   natural_abund_map[481060000] = 1.25;
//   natural_abund_map[481080000] = 0.89;
//   natural_abund_map[481100000] = 12.49;
//   natural_abund_map[481110000] = 12.80;
//   natural_abund_map[481120000] = 24.13;
//   natural_abund_map[481130000] = 12.22;
//   natural_abund_map[481140000] = 28.73;
//   natural_abund_map[481160000] = 7.49;
//   natural_abund_map[491130000] = 4.29;
//   natural_abund_map[491150000] = 95.71;
//   natural_abund_map[501120000] = 0.97;
//   natural_abund_map[501140000] = 0.66;
//   natural_abund_map[501150000] = 0.34;
//   natural_abund_map[501160000] = 14.54;
//   natural_abund_map[501170000] = 7.68;
//   natural_abund_map[501180000] = 24.22;
//   natural_abund_map[501190000] = 8.59;
//   natural_abund_map[501200000] = 32.58;
//   natural_abund_map[501220000] = 4.63;
//   natural_abund_map[501240000] = 5.79;
//   natural_abund_map[511210000] = 57.21;
//   natural_abund_map[511230000] = 42.79;
//   natural_abund_map[521200000] = 0.09;
//   natural_abund_map[521220000] = 2.55;
//   natural_abund_map[521230000] = 0.89;
//   natural_abund_map[521240000] = 4.74;
//   natural_abund_map[521250000] = 7.07;
//   natural_abund_map[521260000] = 18.84;
//   natural_abund_map[521280000] = 31.74;
//   natural_abund_map[521300000] = 34.08;
//   natural_abund_map[531270000] = 100.0;
//   natural_abund_map[541240000] = 0.0952;
//   natural_abund_map[541260000] = 0.0890;
//   natural_abund_map[541280000] = 1.9102;
//   natural_abund_map[541290000] = 26.4006;
//   natural_abund_map[541300000] = 4.0710;
//   natural_abund_map[541310000] = 21.2324;
//   natural_abund_map[541320000] = 26.9086;
//   natural_abund_map[541340000] = 10.4357;
//   natural_abund_map[541360000] = 8.8573;
//   natural_abund_map[551330000] = 100.0;
//   natural_abund_map[561300000] = 0.106;
//   natural_abund_map[561320000] = 0.101;
//   natural_abund_map[561340000] = 2.417;
//   natural_abund_map[561350000] = 6.592;
//   natural_abund_map[561360000] = 7.854;
//   natural_abund_map[561370000] = 11.232;
//   natural_abund_map[561380000] = 71.698;
//   natural_abund_map[571380000] = 0.08881;
//   natural_abund_map[571390000] = 99.91119;
//   natural_abund_map[581360000] = 0.185;
//   natural_abund_map[581380000] = 0.251;
//   natural_abund_map[581400000] = 88.450;
//   natural_abund_map[581420000] = 11.114;
//   natural_abund_map[591410000] = 100.0;
//   natural_abund_map[601420000] = 27.152;
//   natural_abund_map[601430000] = 12.174;
//   natural_abund_map[601440000] = 23.798;
//   natural_abund_map[601450000] = 8.293;
//   natural_abund_map[601460000] = 17.189;
//   natural_abund_map[601480000] = 5.756;
//   natural_abund_map[601500000] = 5.638;
//   natural_abund_map[621440000] = 3.07;
//   natural_abund_map[621470000] = 14.99;
//   natural_abund_map[621480000] = 11.24;
//   natural_abund_map[621490000] = 13.82;
//   natural_abund_map[621500000] = 7.38;
//   natural_abund_map[621520000] = 26.75;
//   natural_abund_map[621540000] = 22.75;
//   natural_abund_map[631510000] = 47.81;
//   natural_abund_map[631530000] = 52.19;
//   natural_abund_map[641520000] = 0.20;
//   natural_abund_map[641540000] = 2.18;
//   natural_abund_map[641550000] = 14.80;
//   natural_abund_map[641560000] = 20.47;
//   natural_abund_map[641570000] = 15.65;
//   natural_abund_map[641580000] = 24.84;
//   natural_abund_map[641600000] = 21.86;
//   natural_abund_map[651590000] = 100.0;
//   natural_abund_map[661560000] = 0.056;
//   natural_abund_map[661580000] = 0.095;
//   natural_abund_map[661600000] = 2.329;
//   natural_abund_map[661610000] = 18.889;
//   natural_abund_map[661620000] = 25.475;
//   natural_abund_map[661630000] = 24.896;
//   natural_abund_map[661640000] = 28.260;
//   natural_abund_map[671650000] = 100.0;
//   natural_abund_map[681620000] = 0.139;
//   natural_abund_map[681640000] = 1.601;
//   natural_abund_map[681660000] = 33.503;
//   natural_abund_map[681670000] = 22.869;
//   natural_abund_map[681680000] = 26.978;
//   natural_abund_map[681700000] = 14.910;
//   natural_abund_map[691690000] = 100.0;
//   natural_abund_map[701680000] = 0.123;
//   natural_abund_map[701700000] = 2.982;
//   natural_abund_map[701710000] = 14.09;
//   natural_abund_map[701720000] = 21.68;
//   natural_abund_map[701730000] = 16.103;
//   natural_abund_map[701740000] = 32.026;
//   natural_abund_map[701760000] = 12.996;
//   natural_abund_map[711750000] = 97.401;
//   natural_abund_map[711760000] = 2.599;
//   natural_abund_map[721740000] = 0.16;
//   natural_abund_map[721760000] = 5.26;
//   natural_abund_map[721770000] = 18.60;
//   natural_abund_map[721780000] = 27.28;
//   natural_abund_map[721790000] = 13.62;
//   natural_abund_map[721800000] = 35.08;
//   natural_abund_map[731800000] = 0.01201;
//   natural_abund_map[731810000] = 99.98799;
//   natural_abund_map[741800000] = 0.12;
//   natural_abund_map[741820000] = 26.50;
//   natural_abund_map[741830000] = 14.31;
//   natural_abund_map[741840000] = 30.64;
//   natural_abund_map[741860000] = 28.43;
//   natural_abund_map[751850000] = 37.40;
//   natural_abund_map[751870000] = 62.60;
//   natural_abund_map[761840000] = 0.02;
//   natural_abund_map[761860000] = 1.59;
//   natural_abund_map[761870000] = 1.96;
//   natural_abund_map[761880000] = 13.24;
//   natural_abund_map[761890000] = 16.15;
//   natural_abund_map[761900000] = 26.26;
//   natural_abund_map[761920000] = 40.78;
//   natural_abund_map[771910000] = 37.3;
//   natural_abund_map[771930000] = 62.7;
//   natural_abund_map[781900000] = 0.012;
//   natural_abund_map[781920000] = 0.782;
//   natural_abund_map[781940000] = 32.86;
//   natural_abund_map[781950000] = 33.78;
//   natural_abund_map[781960000] = 25.21;
//   natural_abund_map[781980000] = 7.356;
//   natural_abund_map[791970000] = 100.0;
//   natural_abund_map[801960000] = 0.15;
//   natural_abund_map[801980000] = 9.97;
//   natural_abund_map[801990000] = 16.87;
//   natural_abund_map[802000000] = 23.10;
//   natural_abund_map[802010000] = 13.18;
//   natural_abund_map[802020000] = 29.86;
//   natural_abund_map[802040000] = 6.87;
//   natural_abund_map[812030000] = 29.52;
//   natural_abund_map[812050000] = 70.48;
//   natural_abund_map[822040000] = 1.4;
//   natural_abund_map[822060000] = 24.1;
//   natural_abund_map[822070000] = 22.1;
//   natural_abund_map[822080000] = 52.4;
//   natural_abund_map[832090000] = 100.0;
//   natural_abund_map[902320000] = 100.0;
//   natural_abund_map[912310000] = 100.0;
//   natural_abund_map[922340000] = 0.0054;
//   natural_abund_map[922350000] = 0.7204;
//   natural_abund_map[922380000] = 99.2742;
// }
// //
// end of src/atomic_data.cpp~
//


//
// start of src/atomic_data.h
//
/// \/file atomic_nuclear_data.h
/// \/author Andrew Davis (andrew.davis@wisc.edu)
///
/// \/brief Implements all the fundamental atomic & nuclear data data
#include <map>

namespace pyne
{
/// main function to be called when you wish to load the nuclide data
/// into memory
void _load_atomic_mass_map_memory();
/// function to create mapping from nuclides in id form
/// to their atomic masses

void _insert_atomic_mass_map();

/// function to create mapping from nuclides in id form
/// to their natural abundances
void _insert_abund_map();

/// Mapping from nuclides in id form to their natural abundances
extern std::map<int,double> natural_abund_map;

/// Mapping from nuclides in id form to their atomic masses.
extern std::map<int,double> atomic_mass_map;

/// Mapping from nuclides in id form to the associated error in
/// abdundance
extern std::map<int,double> atomic_mass_error_map;
} // namespace pyne
//
// end of src/atomic_data.h
//


//
// start of src/atomic_data.h~
//
// /// \/file atomic_nuclear_data.h
// /// \/author Andrew Davis (andrew.davis@wisc.edu)
// ///
// /// \/brief Impliments all the fundamental atomic & nuclear data data
// #include <map>
//
// namespace pyne
// {
//   /// main function to be called when you whish to load the nuclide data
//   /// into memory
//   void _load_atomic_mass_map_memory();
//   /// function to create mapping from nuclides in id form
//   /// to their atomic masses
//
//   void _insert_atomic_mass_map();
//
//   /// function to create mapping from nuclides in id form
//   /// to their natural abundances
//   void _insert_abund_map();
//
//   /// Mapping from nuclides in id form to their natural abundances
//   extern std::map<int,double> natural_abund_map;
//
//   /// Mapping from nuclides in id form to their atomic masses.
//   extern std::map<int,double> atomic_mass_map;
//
//   /// Mapping from nuclides in id form to the associated error in
//   /// abdundance
//   extern std::map<int,double> atomic_mass_error_map;
// } // namespace pyne
// //
// end of src/atomic_data.h~
//


#endif  // PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY
