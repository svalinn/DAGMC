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
//   src/enrichment_cascade.h
//   src/enrichment.h
//   src/enrichment_symbolic.h
//   src/tally.h
//   src/decay.h
//   src/atomic_data.h

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
// http://sigma.mcs.anl.gov/moab-library/
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
// start of src/enrichment_cascade.h
//
/// \brief A container representing enrichment cascades.

#ifndef PYNE_3QGDWZMLZBHDHI424JL52SQHN4
#define PYNE_3QGDWZMLZBHDHI424JL52SQHN4

#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#include "nucname.h"
#include "data.h"
#include "material.h"
#endif

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

namespace pyne
{
namespace enrichment
{

/// A set of physical parameters used to specify an enrichment cascade.
class Cascade
{

 public:

  ///  default constructors
  Cascade();

  /// default destructor
  ~Cascade();

  // Attributes
  double alpha; ///< stage separation factor
  double Mstar; ///< mass separation factor

  int j; ///< Component to enrich (U-235), id form
  int k; ///< Component to de-enrich, or strip (U-238), id form

  double N; ///< number of enriching stages
  double M; ///< number of stripping stages

  double x_feed_j; ///< enrichment of the #j-th isotope in the feed stream
  double x_prod_j; ///< enrichment of the #j-th isotope in the product stream
  double x_tail_j; ///< enrichment of the #j-th isotope in the tails stream

  pyne::Material mat_feed; ///< feed material
  pyne::Material mat_prod; ///< product material
  pyne::Material mat_tail; ///< tails material

  double l_t_per_feed; ///< Total flow rate per feed rate.
  double swu_per_feed; ///< This is the SWU for 1 kg of Feed material.
  double swu_per_prod; ///< This is the SWU for 1 kg of Product material.

  // member functions
  void _reset_xjs();  ///< Sets #x_feed_j to #j-th value of #mat_feed.
};

// end enrichment
}
// end pyne
}

#endif


//
// end of src/enrichment_cascade.h
//


//
// start of src/enrichment.h
//
/// \file enrichment.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Top-level enrichment functionality.

#ifndef PYNE_B3ANNCKDQ5HEJLI33RPZPDNX6A
#define PYNE_B3ANNCKDQ5HEJLI33RPZPDNX6A

#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_symbolic.h"
#endif

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

namespace pyne
{
//! Enrichment Component Class and Functions
namespace enrichment
{

/// Greates a cascade instance with default values for a uranium enrichment.
Cascade _fill_default_uranium_cascade();
/// a cascade instance with default values for a uranium enrichment.
extern Cascade default_uranium_cascade;

/// Computes the feed per product mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double feed_per_prod(double x_feed, double x_prod, double x_tail);
/// Computes the feed per tails mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double feed_per_tail(double x_feed, double x_prod, double x_tail);
/// Computes the product per feed mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double prod_per_feed(double x_feed, double x_prod, double x_tail);
/// Computes the product per tails mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double prod_per_tail(double x_feed, double x_prod, double x_tail);
/// Computes the tails per feed mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double tail_per_feed(double x_feed, double x_prod, double x_tail);
/// Computes the tails per product mass ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double tail_per_prod(double x_feed, double x_prod, double x_tail);
/// Computes the value or separation potential of an assay \a x.
double value_func(double x);
/// Computes the swu per feed ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double swu_per_feed(double x_feed, double x_prod, double x_tail);
/// Computes the swu per product ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double swu_per_prod(double x_feed, double x_prod, double x_tail);
/// Computes the swu per tails ratio for feed, product, and tails
/// enrichments \a x_feed, \a x_prod, and \a x_tails.
double swu_per_tail(double x_feed, double x_prod, double x_tail);

/// Computes the nuclide-specific stage separation factor from the
/// overall stage separation factor \a alpha, the key mass \a Mstar,
/// and the nulide's atomic mass \a M_i.
double alphastar_i(double alpha, double Mstar, double M_i);

/// \{ Numeric Solver Functions
/// \{
/// These functions enable the numeric cascade solver.

/// Finds the total flow rate (L) over the feed flow rate (F), the number of
/// enriching stages (N), and the number of stripping stages (M).
/// \param orig_casc Original input cascade.
/// \param tolerance Maximum numerical error allowed in L/F, N, and M.
/// \param max_iter Maximum number of iterations for to perform.
/// \return A cascade whose N & M coorespond to the L/F value.
Cascade solve_numeric(Cascade & orig_casc, double tolerance=1.0E-7,
                      int max_iter=100);
/// So,ves for valid stage numbers N &nd M of a casc.
/// \param casc Cascade instance, modified in-place.
/// \param tolerance Maximum numerical error allowed in N and M.
void _recompute_nm(Cascade & casc, double tolerance=1.0E-7);
/// This function takes a given initial guess number of enriching and stripping
/// stages for a given composition of fuel with a given jth key component, knowing
/// the values that are desired in both Product and Tails streams.  Having this it
/// solves for what the actual N and M stage numbers are and also what the product
/// and waste streams compositions are.
/// \param casc Cascade instance, modified in-place.
void _recompute_prod_tail_mats(Cascade & casc);
/// This function solves the whole system of equations.  It uses
/// _recompute_prod_tail_mats() to find the roots for the enriching and stripping
/// stage numbers.  It then checks to see if the product and waste streams meet
/// their target enrichments for the jth component like they should.  If they don't
/// then it trys other values of N and M varied by the Secant method.
/// \param casc Input cascade.
/// \param tolerance Maximum numerical error allowed in L/F, N, and M.
/// \param max_iter Maximum number of iterations for to perform.
/// \return A cascade whose N & M coorespond to the L/F value.
Cascade _norm_comp_secant(Cascade & casc, double tolerance=1.0E-7, int max_iter=100);
/// Solves for a stage separative power relevant to the ith component
/// per unit of flow G.  This is taken from Equation 31 divided by G
/// from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
/// 'On a Criterion Efficiency for Multi-Isotope Mixtures Separation',
/// Separation Science and Technology, 34:3, 343 - 357"
/// To link to this article: DOI: 10.1081/SS-100100654
/// URL: http://dx.doi.org/10.1081/SS-100100654
/// \param casc Input cascade.
/// \param i nuclide in id form.
double _deltaU_i_OverG(Cascade & casc, int i);
/// \}

/// \name Multicomponent Functions
/// \{
/// Finds a value of Mstar by minimzing the seperative power.
/// Note that Mstar on \a orig_casc represents an intial guess at what Mstar might
/// be. This is the final function that actually solves for an optimized M* that
/// makes the cascade!
/// \param orig_casc Original input cascade.
/// \param solver flag for solver to use, may be 'symbolic' or 'numeric'.
/// \param tolerance Maximum numerical error allowed in L/F, N, and M.
/// \param max_iter Maximum number of iterations for to perform.
/// \return A cascade whose N & M coorespond to the L/F value.
Cascade multicomponent(Cascade & orig_casc, char * solver,
                       double tolerance=1.0E-7, int max_iter=100);
Cascade multicomponent(Cascade & orig_casc, std::string solver="symbolic",
                       double tolerance=1.0E-7, int max_iter=100);
/// \}

/// Custom exception for when an enrichment solver has entered an infinite loop.
class EnrichmentInfiniteLoopError: public std::exception
{
  /// Returns a moderately helpful error message.
  virtual const char* what() const throw() {
    return "Inifinite loop found while calculating enrichment cascade.";
  }
};

/// Custom exception for when an enrichment solver has reached its maximum
/// number of iterations.
class EnrichmentIterationLimit: public std::exception
{
  /// Returns a moderately helpful error message.
  virtual const char* what() const throw() {
    return "Iteration limit hit durring enrichment calculation.";
  }
};

/// Custom exception for when an enrichment solver iteration has produced a NaN.
class EnrichmentIterationNaN: public std::exception
{
  /// Returns a moderately helpful error message.
  virtual const char* what() const throw() {
    return "Iteration has hit a point where some values are not-a-number.";
  }
};

// end enrichment
}
// end pyne
}

#endif
//
// end of src/enrichment.h
//


//
// start of src/enrichment_symbolic.h
//

/// \file enrichment_symbolic.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief A multicomponent enrichment cascade solver using
///     a symbolic solution to the mass flow rate equations.

/*********************************************************/
/***            Symbolic Enrichment Functions          ***/
/*** WARNING: This file is auto-generated.             ***/
/***                  DO NOT MODIFY!!!                 ***/
/*********************************************************/

#ifndef PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM
#define PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM

#include <math.h>

#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_cascade.h"
#endif

namespace pyne
{
namespace enrichment
{

/// A multicomponent enrichment cascade solver using
/// a symbolic solution to the mass flow rate equations.
/// \param orig_casc The original state of the cascade.
/// \return A cascade solved for new N, M, and total flow
///         rates.
Cascade solve_symbolic(Cascade & orig_casc);

// end enrichment
}
// end pyne
}

#endif
//
// end of src/enrichment_symbolic.h
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
// start of src/decay.h
//
#ifndef PYNE_DECAY_IS_DUMMY
#ifndef PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
#define PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//              WARNING
// This file has been auto generated
// Do not modify directly. You have
// been warned. This is that warning
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include <map>
//#include <cmath>

#ifndef PYNE_IS_AMALGAMATED
#include "data.h"
#include "nucname.h"
#endif

namespace pyne
{
namespace decayers
{

extern const int all_nucs[3979];

std::map<int, double> decay(std::map<int, double> comp, double t);

}  // namespace decayers
}  // namespace pyne

#endif  // PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
#endif  // PYNE_DECAY_IS_DUMMY//
// end of src/decay.h
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


#endif  // PYNE_52BMSKGZ3FHG3NQI566D4I2ZLY
