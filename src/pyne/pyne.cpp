// This file is composed of the following original files:

//   license.txt
//   src/utils.cpp
//   src/state_map.cpp
//   src/nucname.cpp
//   src/rxname.cpp
//   src/particle.cpp
//   src/data.cpp
//   src/jsoncpp.cpp
//   src/jsoncustomwriter.cpp
//   src/material.cpp
//   src/material_library.cpp
//   src/tally.cpp
//   src/atomic_data.cpp
//   src/measure.cpp
//   src/source_sampling.cpp

// PyNE amalgated source http://pyne.io/
#include "pyne.h"

//
// start of license.txt
//
// Copyright 2011-2019, the PyNE Development Team. All rights reserved.
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
// https://press3.mcs.anl.gov/sigma/moab-library
//
// -------------------------------------------------------------------------------
// The files in fortranformat/ are covered by:
//
// The MIT License. Copyright (c) 2011 Brendan Arnold
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// https://bitbucket.org/brendanarnold/py-fortranformat/src/
// //
// end of license.txt
//


//
// start of src/utils.cpp
//
// General Library
#ifndef PYNE_IS_AMALGAMATED
extern "C" double endftod_(char* str, int len);
#endif
#include <iomanip>

#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#endif


// PyNE Globals

std::string pyne::PYNE_DATA = "";
std::string pyne::NUC_DATA_PATH = "";
std::string pyne::VERSION = "0.7.1";

void pyne::pyne_start() {
#if defined __WIN_MSVC__
  char* tmpPYNE_DATA;
  size_t lenPYNE_DATA;
  errno_t errPYNE_DATA = _dupenv_s(&tmpPYNE_DATA, &lenPYNE_DATA, "PYNE_DATA");
  if (errPYNE_DATA)
    tmpPYNE_DATA = (char*) "<NOT_FOUND>";
  PYNE_DATA = (std::string) tmpPYNE_DATA;

  char* tmpNUC_DATA_PATH;
  size_t lenNUC_DATA_PATH;
  errno_t errNUC_DATA_PATH = _dupenv_s(&tmpNUC_DATA_PATH, &lenNUC_DATA_PATH, "NUC_DATA_PATH");
  if (errPYNE_DATA)
    tmpNUC_DATA_PATH = (char*) "<NOT_FOUND>";
  NUC_DATA_PATH = (std::string) tmpNUC_DATA_PATH;
#else
  char* tmppath;
  tmppath = getenv("PYNE_DATA");
  if (tmppath == NULL)
    tmppath = (char*) "<NOT_FOUND>";
  PYNE_DATA = std::string(tmppath);

  tmppath = getenv("NUC_DATA_PATH");
  if (tmppath == NULL)
    tmppath = (char*) "<NOT_FOUND>";
  NUC_DATA_PATH = std::string(tmppath);
#endif
  return;
}



// String Transformations
std::string pyne::to_str(int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(unsigned int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(double t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(bool t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}


int pyne::to_int(std::string s) {
  return atoi(s.c_str());
}

double pyne::to_dbl(std::string s) {
  return strtod(s.c_str(), NULL);
}

double pyne::endftod_cpp(char* s) {
  // Converts string from ENDF only handles "E-less" format but is 5x faster
  int pos, mant, exp;
  double v, dbl_exp;

  mant = exp = 0;
  if (s[2] == '.') {
    // Convert an ENDF float
    if (s[9] == '+' or s[9] == '-') {
      // All these factors of ten are from place values.
      mant = s[8] + 10 * s[7] + 100 * s[6] + 1000 * s[5] + 10000 * s[4] + \
             100000 * s[3] + 1000000 * s[1] - 1111111 * '0';
      exp = s[10] - '0';
      // Make the right power of 10.
      dbl_exp = exp & 01 ? 10. : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 100. : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e4 : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e8 : 1;
      // Adjust for powers of ten from treating mantissa as an integer.
      dbl_exp = (s[9] == '-' ? 1 / dbl_exp : dbl_exp) * 1.0e-6;
      // Get mantissa sign, apply exponent.
      v = mant * (s[0] == '-' ? -1 : 1) * dbl_exp;
    } else {
      mant = s[7] + 10 * s[6] + 100 * s[5] + 1000 * s[4] + 10000 * s[3] + \
             100000 * s[1] - 111111 * '0';
      exp = s[10] + 10 * s[9] - 11 * '0';
      dbl_exp = exp & 01 ? 10. : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 100. : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e4 : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e8 : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e16 : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e32 : 1;
      dbl_exp *= (exp >>= 1) & 01 ? 1.0e64 : 1;
      dbl_exp = (s[8] == '-' ? 1 / dbl_exp : dbl_exp) * 1.0e-5;
      v = mant * (s[0] == '-' ? -1 : 1) * dbl_exp;
    }
  }

  // Convert an ENDF int to float; we start from the last char in the field and
  // move forward until we hit a non-digit.
  else {
    v = 0;
    mant = 1; // Here we use mant for the place value about to be read in.
    pos = 10;
    while (s[pos] != '-' and s[pos] != '+' and s[pos] != ' ' and pos > 0) {
      v += mant * (s[pos] - '0');
      mant *= 10;
      pos--;
    }
    v *= (s[pos] == '-' ? -1 : 1);
  }
  return v;
}

double pyne::endftod_f(char* s) {
#ifdef PYNE_IS_AMALGAMATED
  return endftod_cpp(s);
#else
  return endftod_(s, 12);
#endif
}

double (*pyne::endftod)(char* s) = &pyne::endftod_f;

void pyne::use_fast_endftod() {
  pyne::endftod = &pyne::endftod_cpp;
}

std::string pyne::to_upper(std::string s) {
  // change each element of the string to upper case.
  for (unsigned int i = 0; i < s.length(); i++)
    s[i] = toupper(s[i]);
  return s;
}

std::string pyne::to_lower(std::string s) {
  // change each element of the string to lower case
  for (unsigned int i = 0; i < s.length(); i++)
    s[i] = tolower(s[i]);
  return s;
}

std::string pyne::comment_line_wrapping(std::string line,
                                        std::string comment_prefix,
                                        int line_length) {
  std::ostringstream oss;

  line_length -= comment_prefix.length();

  // Include as is if short enough
  while (line.length() > line_length) {
    oss << comment_prefix << line.substr(0, line_length) << std::endl;
    line.erase(0, line_length);
  }

  if (line.length() > 0) {
    oss << comment_prefix << line << std::endl;
  }

  return oss.str();
}

std::string pyne::capitalize(std::string s) {
  unsigned int slen = s.length();
  if (slen == 0)
    return s;
  // uppercase the first character
  s[0] = toupper(s[0]);
  // change each subsequent element of the string to lower case
  for (unsigned int i = 1; i < slen; i++)
    s[i] = tolower(s[i]);
  return s;
}


std::string pyne::get_flag(char line[], int max_l) {
  char tempflag [10];
  for (int i = 0; i < max_l; i++) {
    if (line[i] == '\t' || line[i] == '\n' || line[i] == ' ' || line[i] == '\0') {
      tempflag[i] = '\0';
      break;
    } else
      tempflag[i] = line[i];
  }
  return std::string(tempflag);
}



std::string pyne::remove_substring(std::string s, std::string substr) {
  // Removes a substring from the string s
  int n_found = s.find(substr);
  while (0 <= n_found) {
    s.erase(n_found, substr.length());
    n_found = s.find(substr);
  }
  return s;
}


std::string pyne::remove_characters(std::string s, std::string chars) {
  // Removes all characters in the string chars from the string s
  for (int i = 0; i < chars.length(); i++) {
    s = remove_substring(s, chars.substr(i, 1));
  }
  return s;
}


std::string pyne::replace_all_substrings(std::string s, std::string substr, std::string repstr) {
  // Replaces all instance of substr in s with the string repstr
  int n_found = s.find(substr);
  while (0 <= n_found) {
    s.replace(n_found, substr.length(), repstr);
    n_found = s.find(substr);
  }
  return s;
}



std::string pyne::last_char(std::string s) {
  // Returns the last character in a string.
  return s.substr(s.length() - 1, 1);
}


std::string pyne::slice_from_end(std::string s, int n, int l) {
  // Returns the slice of a string using negative indices.
  return s.substr(s.length() + n, l);
}


bool pyne::ternary_ge(int a, int b, int c) {
  // Returns true id a <= b <= c and flase otherwise.
  return (a <= b && b <= c);
}


bool pyne::contains_substring(std::string s, std::string substr) {
  // Returns a boolean based on if the sub is in s.
  int n = s.find(substr);
  return (0 <= n && n < s.length());
}


std::string pyne::natural_naming(std::string name) {
  // Calculates a version on the string name that is a valid
  // variable name, ie it uses only word characters.
  std::string nat_name(name);

  // Replace Whitespace characters with underscores
  nat_name = pyne::replace_all_substrings(nat_name, " ",  "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\t", "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\n", "_");

  // Remove non-word characters
  int n = 0;
  while (n < nat_name.length()) {
    if (pyne::words.find(nat_name[n]) == std::string::npos)
      nat_name.erase(n, 1);
    else
      n++;
  }

  // Make sure that the name in non-empty before continuing
  if (nat_name.length() == 0)
    return nat_name;

  // Make sure that the name doesn't begin with a number.
  if (pyne::digits.find(nat_name[0]) != std::string::npos)
    nat_name.insert(0, "_");

  return nat_name;
}


std::vector<std::string> pyne::split_string(std::string particles_list, std::string delimiter) {
  std::vector<std::string> output_vector;
  size_t prev_pos = 0; //item start position
  size_t pos = 0; //item end position

  while ((pos = particles_list.find(delimiter, prev_pos)) != std::string::npos) {
    output_vector.push_back(particles_list.substr(prev_pos, pos));
    prev_pos = pos + delimiter.length();
  }
  // catch list with a single particle
  if (pos == std::string::npos && prev_pos == 0 && particles_list.length() > 0)
    output_vector.push_back(particles_list);

  return output_vector;
}


//
// Math Helpers
//

double pyne::slope(double x2, double y2, double x1, double y1) {
  // Finds the slope of a line.
  return (y2 - y1) / (x2 - x1);
}


double pyne::solve_line(double x, double x2, double y2, double x1, double y1) {
  return (slope(x2, y2, x1, y1) * (x - x2)) + y2;
}


double pyne::tanh(double x) {
  return std::tanh(x);
}

double pyne::coth(double x) {
  return 1.0 / std::tanh(x);
}



// File Helpers

bool pyne::file_exists(std::string strfilename) {
  // Thank you intarwebz for this function!
  // Sepcifically: http://www.techbytes.ca/techbyte103.html
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strfilename.c_str(), &stFileInfo);

  if (intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return (blnReturn);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(char* filename) {
  std::string file(filename);
  return pyne::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(std::string filename) {
  // remove all extra whitespace
  filename = pyne::remove_characters(" ", filename);
  // use stdlib call
  const char* full_filepath = realpath(filename.c_str(), NULL);
  return std::string(full_filepath);
}

// Message Helpers

bool pyne::USE_WARNINGS = true;

bool pyne::toggle_warnings() {
  USE_WARNINGS = !USE_WARNINGS;
  return USE_WARNINGS;
}

void pyne::warning(std::string s) {
  // Prints a warning message
  if (USE_WARNINGS) {
    std::cout << "\033[1;33m WARNING: \033[0m" << s << "\n";
  }
}




//
// end of src/utils.cpp
//


//
// start of src/state_map.cpp
//
//Mapping file for state ids to nuc ids
//This File was autogenerated!!
#ifndef PYNE_4HFU6PUEQJB3ZJ4UIFLVU4SPCM
#define PYNE_4HFU6PUEQJB3ZJ4UIFLVU4SPCM
namespace pyne {
namespace nucname {
#define TOTAL_STATE_MAPS 1016
std::map<int, int> state_id_map;
int map_nuc_ids [TOTAL_STATE_MAPS] = {110240001,
                                      130240001,
                                      90260001,
                                      130260001,
                                      130320002,
                                      170340001,
                                      170380001,
                                      190380001,
                                      190380015,
                                      210420002,
                                      210430001,
                                      210440004,
                                      230440001,
                                      210450001,
                                      210460002,
                                      230460001,
                                      210500001,
                                      250500001,
                                      250520001,
                                      260520042,
                                      260530022,
                                      270530003,
                                      270540001,
                                      210560001,
                                      210560004,
                                      250580001,
                                      270580001,
                                      270580002,
                                      230600000,
                                      230600001,
                                      250600001,
                                      270600001,
                                      260610004,
                                      250620001,
                                      270620001,
                                      230640001,
                                      250640002,
                                      260650003,
                                      260670002,
                                      290670023,
                                      290680003,
                                      280690001,
                                      280690008,
                                      300690001,
                                      340690004,
                                      270700001,
                                      290700001,
                                      290700003,
                                      350700006,
                                      280710002,
                                      300710001,
                                      320710002,
                                      310720002,
                                      350720001,
                                      300730001,
                                      300730002,
                                      320730002,
                                      340730001,
                                      360730004,
                                      310740002,
                                      350740002,
                                      290750001,
                                      290750002,
                                      300750001,
                                      320750002,
                                      330750004,
                                      280760004,
                                      290760001,
                                      350760002,
                                      300770002,
                                      320770001,
                                      330770004,
                                      340770001,
                                      350770001,
                                      300780004,
                                      310780004,
                                      350780004,
                                      370780003,
                                      390780001,
                                      300790002,
                                      320790001,
                                      330790007,
                                      340790001,
                                      350790001,
                                      360790001,
                                      310800001,
                                      350800002,
                                      390800001,
                                      390800003,
                                      320810001,
                                      340810001,
                                      360810002,
                                      370810001,
                                      330820001,
                                      340820015,
                                      350820001,
                                      370820001,
                                      410820003,
                                      340830001,
                                      360830002,
                                      370830002,
                                      380830002,
                                      390830001,
                                      310840001,
                                      350840001,
                                      360840019,
                                      360840061,
                                      370840002,
                                      390840002,
                                      410840007,
                                      360850001,
                                      370850003,
                                      380850002,
                                      390850001,
                                      400850002,
                                      410850003,
                                      410850005,
                                      370860002,
                                      380860014,
                                      390860002,
                                      410860001,
                                      410860002,
                                      380870001,
                                      390870001,
                                      400870002,
                                      410870001,
                                      430870005,
                                      350880003,
                                      410880001,
                                      430880000,
                                      430880001,
                                      390890001,
                                      400890001,
                                      410890001,
                                      420890002,
                                      430890001,
                                      370900001,
                                      390900002,
                                      400900003,
                                      410900002,
                                      410900007,
                                      430900001,
                                      430900006,
                                      390910001,
                                      400910040,
                                      410910001,
                                      420910001,
                                      430910001,
                                      440910001,
                                      450910001,
                                      410920001,
                                      450920001,
                                      390930002,
                                      410930001,
                                      420930016,
                                      430930001,
                                      440930001,
                                      410940001,
                                      430940001,
                                      470940001,
                                      470940002,
                                      410950001,
                                      430950001,
                                      450950001,
                                      460950005,
                                      470950002,
                                      390960005,
                                      430960001,
                                      450960001,
                                      370970002,
                                      390970001,
                                      390970029,
                                      410970001,
                                      430970001,
                                      450970001,
                                      370980001,
                                      390980005,
                                      410980001,
                                      450980001,
                                      410990001,
                                      430990002,
                                      450990001,
                                      470990002,
                                      371000001,
                                      391000004,
                                      411000001,
                                      411000009,
                                      411000012,
                                      431000002,
                                      431000004,
                                      451000004,
                                      471000001,
                                      451010001,
                                      471010002,
                                      411020001,
                                      431020001,
                                      451020005,
                                      471020001,
                                      441030005,
                                      451030001,
                                      471030002,
                                      491030001,
                                      411040004,
                                      451040003,
                                      471040001,
                                      491040003,
                                      451050001,
                                      471050001,
                                      491050001,
                                      411060003,
                                      451060001,
                                      471060001,
                                      491060001,
                                      431070000,
                                      461070002,
                                      471070001,
                                      491070001,
                                      401080009,
                                      411080003,
                                      451080004,
                                      471080002,
                                      491080001,
                                      461090002,
                                      471090001,
                                      491090001,
                                      491090026,
                                      451100001,
                                      471100002,
                                      491100001,
                                      421110001,
                                      461110002,
                                      471110001,
                                      481110003,
                                      491110001,
                                      491120001,
                                      491120004,
                                      491120009,
                                      461130002,
                                      471130001,
                                      481130001,
                                      491130001,
                                      501130001,
                                      451140005,
                                      491140001,
                                      491140005,
                                      531140005,
                                      461150001,
                                      471150001,
                                      481150001,
                                      491150001,
                                      521150001,
                                      451160001,
                                      471160001,
                                      471160004,
                                      491160001,
                                      491160004,
                                      511160003,
                                      551160001,
                                      441170003,
                                      471170001,
                                      481170002,
                                      491170001,
                                      501170002,
                                      521170003,
                                      471180004,
                                      491180001,
                                      491180003,
                                      511180007,
                                      531180002,
                                      551180001,
                                      441190002,
                                      471190000,
                                      471190001,
                                      481190002,
                                      491190001,
                                      501190002,
                                      511190072,
                                      521190002,
                                      551190001,
                                      451200002,
                                      471200002,
                                      491200001,
                                      491200002,
                                      511200001,
                                      531200013,
                                      551200001,
                                      571200000,
                                      461210001,
                                      481210002,
                                      491210001,
                                      501210001,
                                      521210002,
                                      551210001,
                                      451220002,
                                      471220001,
                                      471220002,
                                      491220001,
                                      491220005,
                                      511220005,
                                      511220006,
                                      551220007,
                                      551220008,
                                      481230003,
                                      491230001,
                                      501230001,
                                      521230002,
                                      551230005,
                                      461240004,
                                      471240001,
                                      471240002,
                                      471240003,
                                      491240002,
                                      501240016,
                                      511240001,
                                      511240002,
                                      551240025,
                                      471250007,
                                      471250010,
                                      481250001,
                                      491250001,
                                      501250001,
                                      521250002,
                                      541250002,
                                      571250005,
                                      461260003,
                                      461260004,
                                      461260005,
                                      471260001,
                                      471260004,
                                      491260001,
                                      511260001,
                                      511260002,
                                      481270006,
                                      491270001,
                                      491270009,
                                      501270001,
                                      521270002,
                                      541270002,
                                      561270002,
                                      571270001,
                                      581270001,
                                      461280004,
                                      511280001,
                                      571280001,
                                      471290001,
                                      481290001,
                                      481290004,
                                      491290001,
                                      491290010,
                                      491290012,
                                      491290013,
                                      501290001,
                                      501290017,
                                      501290018,
                                      501290025,
                                      511290011,
                                      511290012,
                                      511290023,
                                      521290001,
                                      541290002,
                                      551290010,
                                      561290001,
                                      571290002,
                                      601290001,
                                      601290003,
                                      491300001,
                                      491300002,
                                      491300003,
                                      501300002,
                                      511300001,
                                      531300001,
                                      551300004,
                                      561300030,
                                      591300002,
                                      491310001,
                                      491310004,
                                      501310001,
                                      521310001,
                                      521310033,
                                      541310002,
                                      561310002,
                                      571310006,
                                      581310001,
                                      591310002,
                                      501320006,
                                      511320001,
                                      521320006,
                                      521320022,
                                      531320003,
                                      541320030,
                                      571320004,
                                      581320030,
                                      491330001,
                                      521330002,
                                      531330016,
                                      531330059,
                                      531330065,
                                      541330001,
                                      561330002,
                                      581330001,
                                      591330003,
                                      601330001,
                                      611330005,
                                      621330000,
                                      511340002,
                                      521340003,
                                      531340005,
                                      541340007,
                                      551340003,
                                      601340017,
                                      611340000,
                                      611340001,
                                      521350010,
                                      541350002,
                                      551350010,
                                      561350002,
                                      581350004,
                                      591350004,
                                      601350001,
                                      611350000,
                                      611350003,
                                      501360003,
                                      531360006,
                                      551360001,
                                      561360005,
                                      571360051,
                                      611360000,
                                      611360001,
                                      631360001,
                                      561370002,
                                      581370002,
                                      601370004,
                                      501380003,
                                      531380002,
                                      551380003,
                                      581380005,
                                      591380005,
                                      611380001,
                                      571390021,
                                      581390002,
                                      601390002,
                                      601390035,
                                      611390001,
                                      621390004,
                                      631390003,
                                      641390001,
                                      511400003,
                                      591400003,
                                      591400015,
                                      601400009,
                                      611400008,
                                      631400004,
                                      601410002,
                                      621410002,
                                      631410001,
                                      641410004,
                                      651410001,
                                      591420001,
                                      591420024,
                                      601420004,
                                      611420012,
                                      631420031,
                                      641420019,
                                      641420020,
                                      651420003,
                                      621430002,
                                      621430043,
                                      641430002,
                                      651430001,
                                      661430003,
                                      551440004,
                                      591440001,
                                      651440004,
                                      651440006,
                                      651440007,
                                      671440003,
                                      641450002,
                                      651450004,
                                      661450002,
                                      681450002,
                                      571460001,
                                      651460022,
                                      651460026,
                                      661460011,
                                      691460001,
                                      651470001,
                                      661470002,
                                      681470002,
                                      691470001,
                                      591480000,
                                      591480001,
                                      611480003,
                                      651480001,
                                      671480001,
                                      671480012,
                                      681480008,
                                      651490001,
                                      661490027,
                                      671490001,
                                      681490002,
                                      631500001,
                                      651500002,
                                      671500001,
                                      691500005,
                                      581510001,
                                      621510012,
                                      631510002,
                                      651510003,
                                      671510001,
                                      681510021,
                                      691510001,
                                      691510012,
                                      701510001,
                                      701510005,
                                      701510010,
                                      611520004,
                                      611520014,
                                      631520001,
                                      631520016,
                                      651520006,
                                      671520001,
                                      691520006,
                                      691520018,
                                      691520019,
                                      701520006,
                                      621530006,
                                      641530003,
                                      641530008,
                                      651530003,
                                      671530001,
                                      691530001,
                                      601540003,
                                      611540001,
                                      631540013,
                                      651540001,
                                      651540002,
                                      711540015,
                                      721540006,
                                      641550006,
                                      661550009,
                                      671550002,
                                      691550001,
                                      711550001,
                                      711550004,
                                      611560002,
                                      651560002,
                                      651560004,
                                      671560001,
                                      671560012,
                                      711560001,
                                      721560004,
                                      731560001,
                                      641570012,
                                      661570005,
                                      711570001,
                                      731570001,
                                      731570004,
                                      601580004,
                                      611580002,
                                      651580003,
                                      651580023,
                                      671580001,
                                      671580009,
                                      691580002,
                                      691580003,
                                      711580000,
                                      731580001,
                                      731580015,
                                      621590006,
                                      641590002,
                                      661590009,
                                      671590003,
                                      731590001,
                                      601600003,
                                      671600001,
                                      671600006,
                                      691600002,
                                      711600001,
                                      611610005,
                                      671610002,
                                      681610014,
                                      691610001,
                                      711610004,
                                      721610002,
                                      671620003,
                                      691620020,
                                      711620008,
                                      711620009,
                                      751620001,
                                      671630003,
                                      751630001,
                                      621640005,
                                      671640003,
                                      691640001,
                                      771640001,
                                      661650002,
                                      751650001,
                                      771650001,
                                      641660009,
                                      671660001,
                                      691660006,
                                      711660001,
                                      711660002,
                                      681670003,
                                      711670001,
                                      751670001,
                                      671680001,
                                      711680013,
                                      771680001,
                                      701690001,
                                      711690001,
                                      751690001,
                                      771690001,
                                      671700001,
                                      711700008,
                                      771700001,
                                      711710001,
                                      721710001,
                                      771710001,
                                      781710002,
                                      691720011,
                                      711720001,
                                      711720005,
                                      751720001,
                                      771720001,
                                      791720001,
                                      771730000,
                                      771730003,
                                      791730001,
                                      711740003,
                                      771740001,
                                      701750007,
                                      711750053,
                                      791750001,
                                      701760005,
                                      711760001,
                                      731760012,
                                      731760090,
                                      791760001,
                                      791760002,
                                      691770000,
                                      701770006,
                                      711770029,
                                      711770203,
                                      721770048,
                                      721770107,
                                      791770002,
                                      711780003,
                                      721780005,
                                      721780109,
                                      731780000,
                                      731780059,
                                      731780094,
                                      731780139,
                                      711790006,
                                      721790005,
                                      721790046,
                                      731790117,
                                      741790002,
                                      751790137,
                                      791790007,
                                      811790001,
                                      721800006,
                                      731800002,
                                      721810025,
                                      721810078,
                                      761810001,
                                      811810002,
                                      721820009,
                                      721820026,
                                      731820001,
                                      731820029,
                                      741820062,
                                      751820001,
                                      761820029,
                                      721830007,
                                      731830032,
                                      741830007,
                                      751830004,
                                      751830005,
                                      751830070,
                                      751830074,
                                      761830002,
                                      781830001,
                                      811830002,
                                      811830005,
                                      821830001,
                                      721840005,
                                      751840005,
                                      771840007,
                                      781840034,
                                      791840003,
                                      741850006,
                                      781850002,
                                      791850001,
                                      801850004,
                                      811850003,
                                      751860004,
                                      771860001,
                                      811860000,
                                      811860005,
                                      831860001,
                                      791870002,
                                      801870001,
                                      811870002,
                                      821870001,
                                      831870002,
                                      751880007,
                                      811880001,
                                      731890001,
                                      751890033,
                                      751890034,
                                      761890001,
                                      771890006,
                                      771890085,
                                      781890004,
                                      781890005,
                                      791890003,
                                      791890006,
                                      791890200,
                                      801890002,
                                      811890001,
                                      821890001,
                                      821890014,
                                      831890002,
                                      831890003,
                                      731900002,
                                      741900006,
                                      751900003,
                                      761900032,
                                      771900002,
                                      771900037,
                                      791900014,
                                      811900000,
                                      811900001,
                                      811900006,
                                      831900000,
                                      831900001,
                                      761910001,
                                      771910003,
                                      771910071,
                                      791910004,
                                      801910035,
                                      811910002,
                                      821910002,
                                      831910002,
                                      831910005,
                                      831910028,
                                      751920002,
                                      751920003,
                                      761920047,
                                      761920112,
                                      771920003,
                                      771920015,
                                      791920004,
                                      791920015,
                                      811920002,
                                      811920008,
                                      821920011,
                                      821920014,
                                      821920017,
                                      821920020,
                                      821920021,
                                      831920001,
                                      841920006,
                                      851920000,
                                      851920001,
                                      771930002,
                                      781930005,
                                      791930004,
                                      801930003,
                                      811930002,
                                      821930001,
                                      831930001,
                                      841930001,
                                      851930001,
                                      851930002,
                                      751940001,
                                      751940002,
                                      751940003,
                                      771940007,
                                      771940012,
                                      791940003,
                                      791940008,
                                      811940001,
                                      831940001,
                                      831940002,
                                      851940000,
                                      851940001,
                                      761950002,
                                      761950004,
                                      771950002,
                                      781950007,
                                      791950004,
                                      791950055,
                                      801950003,
                                      811950002,
                                      821950002,
                                      831950001,
                                      841950002,
                                      851950001,
                                      861950001,
                                      751960001,
                                      771960004,
                                      791960003,
                                      791960054,
                                      811960006,
                                      831960002,
                                      831960003,
                                      841960015,
                                      761970001,
                                      771970002,
                                      781970009,
                                      791970004,
                                      801970004,
                                      811970002,
                                      821970002,
                                      831970001,
                                      841970002,
                                      851970001,
                                      861970001,
                                      761980006,
                                      761980010,
                                      771980001,
                                      791980051,
                                      811980007,
                                      811980012,
                                      831980001,
                                      831980003,
                                      851980001,
                                      871980001,
                                      781990008,
                                      791990006,
                                      801990007,
                                      811990003,
                                      821990003,
                                      831990001,
                                      841990002,
                                      861990001,
                                      792000011,
                                      812000010,
                                      832000001,
                                      832000003,
                                      852000001,
                                      852000003,
                                      802010013,
                                      812010003,
                                      822010004,
                                      832010001,
                                      842010003,
                                      862010001,
                                      872010001,
                                      882010000,
                                      782020003,
                                      822020014,
                                      852020001,
                                      852020002,
                                      872020001,
                                      822030006,
                                      822030053,
                                      832030006,
                                      842030005,
                                      862030001,
                                      882030001,
                                      812040029,
                                      822040021,
                                      832040008,
                                      832040038,
                                      852040001,
                                      872040001,
                                      872040002,
                                      802050008,
                                      822050009,
                                      842050010,
                                      842050017,
                                      882050001,
                                      812060045,
                                      832060016,
                                      872060001,
                                      872060002,
                                      892060001,
                                      812070002,
                                      822070003,
                                      832070036,
                                      842070014,
                                      862070007,
                                      882070001,
                                      802080004,
                                      832080018,
                                      802100002,
                                      802100005,
                                      832100002,
                                      822110014,
                                      832110021,
                                      842110015,
                                      852110076,
                                      872110013,
                                      872110019,
                                      832120005,
                                      832120012,
                                      842120030,
                                      852120004,
                                      882130005,
                                      822140004,
                                      852140006,
                                      862140004,
                                      862140005,
                                      872140001,
                                      902140004,
                                      832150009,
                                      862150013,
                                      902150003,
                                      822160004,
                                      872160001,
                                      922160001,
                                      832170005,
                                      892170010,
                                      902170001,
                                      912170001,
                                      872180002,
                                      922180001,
                                      892220001,
                                      912340002,
                                      922350001,
                                      922350174,
                                      942350010,
                                      932360001,
                                      952360001,
                                      942370003,
                                      922380119,
                                      932380128,
                                      942380047,
                                      942380052,
                                      952380001,
                                      942390106,
                                      942390111,
                                      952390011,
                                      932400001,
                                      942400102,
                                      952400057,
                                      962400002,
                                      962400003,
                                      942410106,
                                      942410107,
                                      952410075,
                                      962410007,
                                      932420007,
                                      942420044,
                                      942420045,
                                      952420002,
                                      952420141,
                                      962420004,
                                      962420005,
                                      972420002,
                                      972420003,
                                      942440032,
                                      952440001,
                                      952440112,
                                      952440113,
                                      962440009,
                                      962440013,
                                      962440014,
                                      972440004,
                                      982440002,
                                      942450024,
                                      952450021,
                                      962450061,
                                      972450003,
                                      1012450001,
                                      952460001,
                                      952460008,
                                      972460000,
                                      982460002,
                                      992460000,
                                      1012460000,
                                      1012460001,
                                      1002470001,
                                      1012470001,
                                      972480001,
                                      1002480006,
                                      992500001,
                                      1002500001,
                                      1002500002,
                                      1022500001,
                                      1022510002,
                                      1002530008,
                                      1022530003,
                                      1022530030,
                                      1022530031,
                                      1022530032,
                                      1032530000,
                                      1032530001,
                                      992540002,
                                      1012540000,
                                      1012540001,
                                      1022540011,
                                      1032550001,
                                      1032550027,
                                      992560001,
                                      1002560026,
                                      1042560007,
                                      1042560009,
                                      1042560012,
                                      1042570002,
                                      1052570002,
                                      1012580001,
                                      1052580001,
                                      1042610001,
                                      1072620001,
                                      1062630003,
                                      1062650001,
                                      1082650001,
                                      1082670002,
                                      1102700001,
                                      1102710001,
                                      1082770001,
                                     };
int map_metastable [TOTAL_STATE_MAPS] = {1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         4,
                                         1,
                                         2,
                                         3,
                                         4,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         4,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         4,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         4,
                                         5,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         4,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                        };
}
}
#endif
//
// end of src/state_map.cpp
//


//
// start of src/nucname.cpp
//
// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// name is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#ifndef PYNE_IS_AMALGAMATED
#include "nucname.h"
#include "state_map.cpp"
#endif


/*** Constructs the LL to zz Dictionary ***/
pyne::nucname::name_zz_t pyne::nucname::get_name_zz() {
  pyne::nucname::name_zz_t lzd;

  lzd["Be"] = 04;
  lzd["Ba"] = 56;
  lzd["Bh"] = 107;
  lzd["Bi"] = 83;
  lzd["Bk"] = 97;
  lzd["Br"] = 35;
  lzd["Ru"] = 44;
  lzd["Re"] = 75;
  lzd["Rf"] = 104;
  lzd["Rg"] = 111;
  lzd["Ra"] = 88;
  lzd["Rb"] = 37;
  lzd["Rn"] = 86;
  lzd["Rh"] = 45;
  lzd["Tm"] = 69;
  lzd["H"] = 01;
  lzd["P"] = 15;
  lzd["Ge"] = 32;
  lzd["Gd"] = 64;
  lzd["Ga"] = 31;
  lzd["Os"] = 76;
  lzd["Hs"] = 108;
  lzd["Zn"] = 30;
  lzd["Ho"] = 67;
  lzd["Hf"] = 72;
  lzd["Hg"] = 80;
  lzd["He"] = 02;
  lzd["Pr"] = 59;
  lzd["Pt"] = 78;
  lzd["Pu"] = 94;
  lzd["Pb"] = 82;
  lzd["Pa"] = 91;
  lzd["Pd"] = 46;
  lzd["Po"] = 84;
  lzd["Pm"] = 61;
  lzd["C"] = 6;
  lzd["K"] = 19;
  lzd["O"] = 8;
  lzd["S"] = 16;
  lzd["W"] = 74;
  lzd["Eu"] = 63;
  lzd["Es"] = 99;
  lzd["Er"] = 68;
  lzd["Md"] = 101;
  lzd["Mg"] = 12;
  lzd["Mo"] = 42;
  lzd["Mn"] = 25;
  lzd["Mt"] = 109;
  lzd["U"] = 92;
  lzd["Fr"] = 87;
  lzd["Fe"] = 26;
  lzd["Fm"] = 100;
  lzd["Ni"] = 28;
  lzd["No"] = 102;
  lzd["Na"] = 11;
  lzd["Nb"] = 41;
  lzd["Nd"] = 60;
  lzd["Ne"] = 10;
  lzd["Zr"] = 40;
  lzd["Np"] = 93;
  lzd["B"] = 05;
  lzd["Co"] = 27;
  lzd["Cm"] = 96;
  lzd["F"] = 9;
  lzd["Ca"] = 20;
  lzd["Cf"] = 98;
  lzd["Ce"] = 58;
  lzd["Cd"] = 48;
  lzd["V"] = 23;
  lzd["Cs"] = 55;
  lzd["Cr"] = 24;
  lzd["Cu"] = 29;
  lzd["Sr"] = 38;
  lzd["Kr"] = 36;
  lzd["Si"] = 14;
  lzd["Sn"] = 50;
  lzd["Sm"] = 62;
  lzd["Sc"] = 21;
  lzd["Sb"] = 51;
  lzd["Sg"] = 106;
  lzd["Se"] = 34;
  lzd["Yb"] = 70;
  lzd["Db"] = 105;
  lzd["Dy"] = 66;
  lzd["Ds"] = 110;
  lzd["La"] = 57;
  lzd["Cl"] = 17;
  lzd["Li"] = 03;
  lzd["Tl"] = 81;
  lzd["Lu"] = 71;
  lzd["Lr"] = 103;
  lzd["Th"] = 90;
  lzd["Ti"] = 22;
  lzd["Te"] = 52;
  lzd["Tb"] = 65;
  lzd["Tc"] = 43;
  lzd["Ta"] = 73;
  lzd["Ac"] = 89;
  lzd["Ag"] = 47;
  lzd["I"] = 53;
  lzd["Ir"] = 77;
  lzd["Am"] = 95;
  lzd["Al"] = 13;
  lzd["As"] = 33;
  lzd["Ar"] = 18;
  lzd["Au"] = 79;
  lzd["At"] = 85;
  lzd["In"] = 49;
  lzd["Y"] = 39;
  lzd["N"] = 07;
  lzd["Xe"] = 54;
  lzd["Cn"] = 112;
  lzd["Fl"] = 114;
  lzd["Lv"] = 116;
  lzd["Nh"] = 113;
  lzd["Mc"] = 115;
  lzd["Ts"] = 117;
  lzd["Og"] = 118;

  return lzd;
}
pyne::nucname::name_zz_t pyne::nucname::name_zz = pyne::nucname::get_name_zz();


/*** Constructs zz to LL dictionary **/
pyne::nucname::zzname_t pyne::nucname::get_zz_name() {
  zzname_t zld;
  for (name_zz_iter i = name_zz.begin(); i != name_zz.end(); i++) {
    zld[i->second] = i->first;
  }
  return zld;
}
pyne::nucname::zzname_t pyne::nucname::zz_name = pyne::nucname::get_zz_name();



/*** Constructs the fluka to zz Dictionary ***/
pyne::nucname::name_zz_t pyne::nucname::get_fluka_zz() {
  pyne::nucname::name_zz_t fzd;

  fzd["BERYLLIU"] = 40000000;
  fzd["BARIUM"]   = 560000000;
  fzd["BOHRIUM"]  = 1070000000;   // No fluka
  fzd["BISMUTH"]  = 830000000;
  fzd["BERKELIU"] = 970000000;    // No fluka
  fzd["BROMINE"]  = 350000000;
  fzd["RUTHENIU"] = 440000000;    // No fluka
  fzd["RHENIUM"]  = 750000000;
  fzd["RUTHERFO"] = 1040000000;
  fzd["ROENTGEN"] = 1110000000;
  fzd["RADIUM"]   = 880000000;    // No fluka
  fzd["RUBIDIUM"] = 370000000;    // No fluka
  fzd["RADON"]    = 860000000;    // no fluka
  fzd["RHODIUM"]  = 450000000;    // no fluka
  fzd["THULIUM"]  = 690000000;    // no fluka
  fzd["HYDROGEN"] = 10000000;
  fzd["PHOSPHO"]  = 150000000;
  fzd["GERMANIU"] = 320000000;
  fzd["GADOLINI"] = 640000000;
  fzd["GALLIUM"]  = 310000000;
  fzd["OSMIUM"]   = 760000000;    // no fluka
  fzd["HASSIUM"]  = 1080000000;
  fzd["ZINC"]     = 300000000;
  fzd["HOLMIUM"]  = 670000000;    // no fluka
  fzd["HAFNIUM"]  = 720000000;
  fzd["MERCURY"]  = 800000000;
  fzd["HELIUM"]   = 20000000;
  fzd["PRASEODY"] = 590000000;   // no fluka
  fzd["PLATINUM"] = 780000000;
  fzd["239-PU"]   = 940000000;   // "239-PU"
  fzd["LEAD"]     = 820000000;
  fzd["PROTACTI"] = 910000000;   // no fluka
  fzd["PALLADIU"] = 460000000;   // no fluka
  fzd["POLONIUM"] = 840000000;   // no fluka
  fzd["PROMETHI"] = 610000000;   // no fluka
  fzd["CARBON"]   = 60000000;
  fzd["POTASSIU"] = 190000000;
  fzd["OXYGEN"]   = 80000000;
  fzd["SULFUR"]   = 160000000;
  fzd["TUNGSTEN"] = 740000000;
  fzd["EUROPIUM"] = 630000000;
  fzd["EINSTEIN"] = 990000000;   // no fluka
  fzd["ERBIUM"]   = 680000000;   // no fluka
  fzd["MENDELEV"] = 1010000000;  // no fluka
  fzd["MAGNESIU"] = 120000000;
  fzd["MOLYBDEN"] = 420000000;
  fzd["MANGANES"] = 250000000;
  fzd["MEITNERI"] = 1090000000;  // no fluka
  fzd["URANIUM"]  = 920000000;
  fzd["FRANCIUM"] = 870000000;   // no fluka
  fzd["IRON"]     = 260000000;
  fzd["FERMIUM"]  = 1000000000;  // no fluka
  fzd["NICKEL"]   = 280000000;
  fzd["NITROGEN"] = 70000000;
  fzd["NOBELIUM"] = 1020000000;  // no fluka
  fzd["SODIUM"]   = 110000000;
  fzd["NIOBIUM"]  = 410000000;
  fzd["NEODYMIU"] = 600000000;
  fzd["NEON"]     = 100000000;
  fzd["ZIRCONIU"] = 400000000;
  fzd["NEPTUNIU"] = 930000000;   // no fluka
  fzd["BORON"]    = 50000000;
  fzd["COBALT"]   = 270000000;
  fzd["CURIUM"]   = 960000000;   // no fluka
  fzd["FLUORINE"] = 90000000;
  fzd["CALCIUM"]  = 200000000;
  fzd["CALIFORN"] = 980000000;   // no fluka
  fzd["CERIUM"]   = 580000000;
  fzd["CADMIUM"]  = 480000000;
  fzd["VANADIUM"] = 230000000;
  fzd["CESIUM"]   = 550000000;
  fzd["CHROMIUM"] = 240000000;
  fzd["COPPER"]   = 290000000;
  fzd["STRONTIU"] = 380000000;
  fzd["KRYPTON"]  = 360000000;
  fzd["SILICON"]  = 140000000;
  fzd["TIN"]      = 500000000;
  fzd["SAMARIUM"] = 620000000;
  fzd["SCANDIUM"] = 210000000;
  fzd["ANTIMONY"] = 510000000;
  fzd["SEABORGI"] = 1060000000;  // no fluka
  fzd["SELENIUM"] = 340000000;   // no fluka
  fzd["YTTERBIU"] = 700000000;   // no fluka
  fzd["DUBNIUM"]  = 1050000000;  // no fluka
  fzd["DYSPROSI"] = 660000000;   // no fluka
  fzd["DARMSTAD"] = 1100000000;  // no fluka
  fzd["LANTHANU"] = 570000000;
  fzd["CHLORINE"] = 170000000;
  fzd["LITHIUM"]  = 030000000;
  fzd["THALLIUM"] = 810000000;   // no fluka
  fzd["LUTETIUM"] = 710000000;   // no fluka
  fzd["LAWRENCI"] = 1030000000;  // no fluka
  fzd["THORIUM"]  = 900000000;   // no fluka
  fzd["TITANIUM"] = 220000000;
  fzd["TELLURIU"] = 520000000;   // no fluka
  fzd["TERBIUM"]  = 650000000;
  fzd["99-TC"]    = 430000000;   // "99-TC"
  fzd["TANTALUM"] = 730000000;
  fzd["ACTINIUM"] = 890000000;   // no fluka
  fzd["SILVER"]   = 470000000;
  fzd["IODINE"]   = 530000000;
  fzd["IRIDIUM"]  = 770000000;
  fzd["241-AM"]   = 950000000;   // "241-AM"
  fzd["ALUMINUM"] = 130000000;
  fzd["ARSENIC"]  = 330000000;
  fzd["ARGON"]    = 180000000;
  fzd["GOLD"]     = 790000000;
  fzd["ASTATINE"] = 850000000;   // no fluka
  fzd["INDIUM"]   = 490000000;
  fzd["YTTRIUM"]  = 390000000;
  fzd["XENON"]    = 540000000;
  fzd["COPERNIC"] = 1120000000;  // no fluka
  fzd["UNUNQUAD"] = 1140000000;  // no fluka:  UNUNQUADIUM,  "Flerovium"
  fzd["UNUNHEXI"] = 1160000000;  // no fluka:  UNUNHEXIUM , "Livermorium"
  fzd["HYDROG-1"] = 10010000;
  fzd["DEUTERIU"] = 10020000;
  fzd["TRITIUM"]  = 10040000;
  fzd["HELIUM-3"] = 20030000;
  fzd["HELIUM-4"] = 20040000;
  fzd["LITHIU-6"] = 30060000;
  fzd["LITHIU-7"] = 30070000;
  fzd["BORON-10"] = 50100000;
  fzd["BORON-11"] = 50110000;
  fzd["90-SR"]    = 380900000;   // fluka "90-SR"
  fzd["129-I"]    = 531290000;   // fluka "129-I"
  fzd["124-XE"]   = 541240000;   // fluka "124-XE"
  fzd["126-XE"]   = 541260000;   // fluka "126-XE"
  fzd["128-XE"]   = 541280000;   // fluka "128-XE"
  fzd["130-XE"]   = 541300000;   // fluka "130-XE"
  fzd["131-XE"]   = 541310000;   // fluka "131-XE"
  fzd["132-XE"]   = 541320000;   // fluka "132-XE"
  fzd["134-XE"]   = 541340000;   // fluka "134-XE"
  fzd["135-XE"]   = 541350000;   // fluka "135-XE"
  fzd["136-XE"]   = 541360000;   // fluka "136-XE"
  fzd["135-CS"]   = 551350000;   // fluka "135-CS"
  fzd["137-CS"]   = 551370000;   // fluka "137-CS"
  fzd["230-TH"]   = 902300000;   // fluka "230-TH"
  fzd["232-TH"]   = 902320000;   // fluka "232-TH"
  fzd["233-U"]    = 922330000;   // fluka "233-U"
  fzd["234-U"]    = 922340000;   // fluka "234-U"
  fzd["235-U"]    = 922350000;   // fluka "235-U"
  fzd["238-U"]    = 922380000;   // fluka "238-U"

  return fzd;
}
pyne::nucname::name_zz_t pyne::nucname::fluka_zz = pyne::nucname::get_fluka_zz();


/*** Constructs zz to fluka dictionary **/
pyne::nucname::zzname_t pyne::nucname::get_zz_fluka() {
  zzname_t zfd;
  for (name_zz_iter i = fluka_zz.begin(); i != fluka_zz.end(); i++) {
    zfd[i->second] = i->first;
  }
  return zfd;
}
pyne::nucname::zzname_t pyne::nucname::zz_fluka = pyne::nucname::get_zz_fluka();



/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

pyne::nucname::zz_group pyne::nucname::name_to_zz_group(pyne::nucname::name_group eg) {
  zz_group zg;
  for (name_group_iter i = eg.begin(); i != eg.end(); i++)
    zg.insert(name_zz[*i]);
  return zg;
}

// Lanthanides
pyne::nucname::name_t pyne::nucname::LAN_array[15] = {"La", "Ce", "Pr", "Nd",
                                                      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"
                                                     };
pyne::nucname::name_group pyne::nucname::LAN(pyne::nucname::LAN_array,
                                             pyne::nucname::LAN_array + 15);
pyne::nucname::zz_group pyne::nucname::lan = \
                                             pyne::nucname::name_to_zz_group(pyne::nucname::LAN);

// Actinides
pyne::nucname::name_t pyne::nucname::ACT_array[15] = {"Ac", "Th", "Pa", "U",
                                                      "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
                                                     };
pyne::nucname::name_group pyne::nucname::ACT(pyne::nucname::ACT_array, pyne::nucname::ACT_array + 15);
pyne::nucname::zz_group pyne::nucname::act = pyne::nucname::name_to_zz_group(pyne::nucname::ACT);

// Transuarnics
pyne::nucname::name_t pyne::nucname::TRU_array[22] = {"Np", "Pu", "Am", "Cm",
                                                      "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                                                      "Ds", "Rg", "Cn", "Fl", "Lv"
                                                     };
pyne::nucname::name_group pyne::nucname::TRU(pyne::nucname::TRU_array,
                                             pyne::nucname::TRU_array + 22);
pyne::nucname::zz_group pyne::nucname::tru = \
                                             pyne::nucname::name_to_zz_group(pyne::nucname::TRU);

//Minor Actinides
pyne::nucname::name_t pyne::nucname::MA_array[10] = {"Np", "Am", "Cm", "Bk",
                                                     "Cf", "Es", "Fm", "Md", "No", "Lr"
                                                    };
pyne::nucname::name_group pyne::nucname::MA(pyne::nucname::MA_array,
                                            pyne::nucname::MA_array + 10);
pyne::nucname::zz_group pyne::nucname::ma = \
                                            pyne::nucname::name_to_zz_group(pyne::nucname::MA);

//Fission Products
pyne::nucname::name_t pyne::nucname::FP_array[88] = {"Ag", "Al", "Ar", "As",
                                                     "At", "Au", "B",  "Ba", "Be", "Bi", "Br", "C",  "Ca", "Cd", "Ce", "Cl", "Co",
                                                     "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F",  "Fe", "Fr", "Ga", "Gd", "Ge", "H",
                                                     "He", "Hf", "Hg", "Ho", "I",  "In", "Ir", "K",  "Kr", "La", "Li", "Lu", "Mg",
                                                     "Mn", "Mo", "N",  "Na", "Nb", "Nd", "Ne", "Ni", "O",  "Os", "P",  "Pb", "Pd",
                                                     "Pm", "Po", "Pr", "Pt", "Ra", "Rb", "Re", "Rh", "Rn", "Ru", "S",  "Sb", "Sc",
                                                     "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Ti", "Tl", "Tm", "V",
                                                     "W",  "Xe", "Y",  "Yb", "Zn", "Zr"
                                                    };
pyne::nucname::name_group pyne::nucname::FP(pyne::nucname::FP_array,
                                            pyne::nucname::FP_array + 88);
pyne::nucname::zz_group pyne::nucname::fp = \
                                            pyne::nucname::name_to_zz_group(pyne::nucname::FP);


/***************************/
/*** isnuclide functions ***/
/***************************/

bool pyne::nucname::isnuclide(std::string nuc) {
  int n;
  try {
    n = id(nuc);
  } catch (NotANuclide) {
    return false;
  } catch (IndeterminateNuclideForm) {
    return false;
  }
  return isnuclide(n);
}

bool pyne::nucname::isnuclide(const char* nuc) {
  return isnuclide(std::string(nuc));
}

bool pyne::nucname::isnuclide(int nuc) {
  int n;
  try {
    n = id(nuc);
  } catch (NotANuclide) {
    return false;
  } catch (IndeterminateNuclideForm) {
    return false;
  }
  if (n <= 10000000)
    return false;
  int zzz = n / 10000000;
  int aaa = (n % 10000000) / 10000;
  if (aaa == 0)
    return false;  // is element
  else if (aaa < zzz)
    return false;
  return true;
}



/********************/
/*** id functions ***/
/********************/
int pyne::nucname::id(int nuc) {
  if (nuc < 0)
    throw NotANuclide(nuc, "");

  int newnuc;
  int zzz = nuc / 10000000;     // ZZZ ?
  int aaassss = nuc % 10000000; // AAA-SSSS ?
  int aaa = aaassss / 10000;    // AAA ?
  int ssss = aaassss % 10000;   // SSSS ?
  // Nuclide must already be in id form
  if (0 < zzz && zzz <= aaa && aaa <= zzz * 7) {
    // Normal nuclide
    if (5 < ssss) {
      // Unphysical metastable state warning
      warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return nuc;
  } else if (aaassss == 0 && 0 < zz_name.count(zzz)) {
    // Natural elemental nuclide:  ie for Uranium = 920000000
    return nuc;
  } else if (nuc < 1000 && 0 < zz_name.count(nuc))
    //  Gave Z-number
    return nuc * 10000000;

  // Not in id form, try  ZZZAAAM form.
  zzz = nuc / 10000;     // ZZZ ?
  aaassss = nuc % 10000; // AAA-SSSS ?
  aaa = aaassss / 10;    // AAA ?
  ssss = nuc % 10;       // SSSS ?
  if (zzz <= aaa && aaa <= zzz * 7) {
    // ZZZAAAM nuclide
    if (5 < ssss) {
      // Unphysical metastable state warning
      warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return (zzz * 10000000) + (aaa * 10000) + (nuc % 10);
  } else if (aaa <= zzz && zzz <= aaa * 7 && 0 < zz_name.count(aaa)) {
    // Cinder-form (aaazzzm), ie 2350920
    if (5 < ssss) {
      // Unphysical metastable state warning
      warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return (aaa * 10000000) + (zzz * 10000) + (nuc % 10);
  }
  //else if (aaassss == 0 && 0 == zz_name.count(nuc/1000) && 0 < zz_name.count(zzz))
  else if (aaassss == 0 && 0 < zz_name.count(zzz)) {
    // zzaaam form natural nuclide
    return zzz * 10000000;
  }

  if (nuc >= 1000000) {
    // From now we assume no metastable info has been given.
    throw IndeterminateNuclideForm(nuc, "");
  }

  // Nuclide is not in zzaaam form,
  // Try MCNP form, ie zzaaa
  // This is the same form as SZA for the 0th state.
  zzz = nuc / 1000;
  aaa = nuc % 1000;
  if (zzz <= aaa) {
    if (aaa - 400 < 0) {
      if (nuc == 95242)
        return nuc * 10000 + 1;  // special case MCNP Am-242m
      else
        return nuc * 10000;  // Nuclide in normal MCNP form
    } else {
      // Nuclide in MCNP metastable form
      if (nuc == 95642)
        return (95642 - 400) * 10000; // special case MCNP Am-242
      nuc = ((nuc - 400) * 10000) + 1;
      while (3.0 < (float ((nuc / 10000) % 1000) / float (nuc / 10000000)))
        nuc -= 999999;
      return nuc;
    }
  } else if (aaa == 0 && 0 < zz_name.count(zzz)) {
    // MCNP form natural nuclide
    return zzz * 10000000;
  }

  // Not a normal nuclide, might be a
  // Natural elemental nuclide.
  // ie 92 for Uranium = 920000
  if (0 < zz_name.count(nuc))
    return nuc * 10000000;
  throw IndeterminateNuclideForm(nuc, "");
}

int pyne::nucname::id(const char* nuc) {
  std::string newnuc(nuc);
  return id(newnuc);
}

int pyne::nucname::id(std::string nuc) {
  size_t npos = std::string::npos;
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int newnuc;
  std::string elem_name;
  int dash1 = nuc.find("-");
  int dash2;
  if (dash1 == npos)
    dash2 = npos;
  else
    dash2 = nuc.find("-", dash1 + 1);

  // nuc must be at least 4 characters or greater if it is in ZZLLAAAM form.
  if (nuc.length() >= 5 && dash1 != npos && dash2 != npos) {
    // Nuclide most likely in ZZLLAAAM Form, only form that contains two "-"'s.
    std::string zz = nuc.substr(0, dash1);
    std::string ll = nuc.substr(dash1 + 1, dash2);
    int zz_int = to_int(zz);
    // Verifying that the LL and ZZ point to the same element as secondary
    if (znum(ll) != zz_int)
      throw NotANuclide(nuc, "mismatched znum and chemical symbol");
    return zzllaaam_to_id(nuc);
  }

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  if (pyne::contains_substring(pyne::digits, nucstr.substr(0, 1))) {
    if (pyne::contains_substring(pyne::digits, nucstr.substr(nuclen - 1, nuclen))) {
      // Nuclide must actually be an integer that
      // just happens to be living in string form.
      newnuc = pyne::to_int(nucstr);
      newnuc = id(newnuc);
    } else {
      // probably in NIST-like form (242Am)
      // Here we know we have both digits and letters
      std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);
      newnuc = pyne::to_int(anum_str) * 10000;

      // Add the Z-number
      elem_name = pyne::remove_characters(nucstr, pyne::digits);
      elem_name = pyne::capitalize(elem_name);
      if (0 < name_zz.count(elem_name))
        newnuc = (10000000 * name_zz[elem_name]) + newnuc;
      else
        throw NotANuclide(nucstr, newnuc);
    }
  } else if (pyne::contains_substring(pyne::alphabet, nucstr.substr(0, 1))) {
    // Nuclide is probably in name form, or some variation therein
    std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

    // natural element form, a la 'U' -> 920000000
    if (anum_str.empty()) {
      elem_name = pyne::capitalize(nucstr);
      if (0 < name_zz.count(elem_name))
        return 10000000 * name_zz[elem_name];
    }

    int anum = pyne::to_int(anum_str);

    // bad form
    if (anum < 0)
      throw NotANuclide(nucstr, anum);

    // Figure out if we are meta-stable or not
    std::string end_char = pyne::last_char(nucstr);
    if (end_char == "M")
      newnuc = (10000 * anum) + 1;
    else if (pyne::contains_substring(pyne::digits, end_char))
      newnuc = (10000 * anum);
    else
      throw NotANuclide(nucstr, newnuc);

    // Add the Z-number
    elem_name = pyne::remove_characters(nucstr.substr(0, nuclen - 1), pyne::digits);
    elem_name = pyne::capitalize(elem_name);
    if (0 < name_zz.count(elem_name))
      newnuc = (10000000 * name_zz[elem_name]) + newnuc;
    else
      throw NotANuclide(nucstr, newnuc);
  } else {
    // Clearly not a nuclide
    throw NotANuclide(nuc, nucstr);
  }
  return newnuc;
}


/***************************/
/*** iselement functions ***/
/***************************/

bool pyne::nucname::iselement(std::string nuc) {
  int n;
  try {
    n = id(nuc);
  } catch (NotANuclide) {
    return false;
  }
  return iselement(n);
}

bool pyne::nucname::iselement(const char* nuc) {
  return iselement(std::string(nuc));
}

bool pyne::nucname::iselement(int nuc) {
  int n;
  try {
    n = id(nuc);
  } catch (NotANuclide) {
    return false;
  }

  if (n < 10000000)
    return false;
  int zzz = znum(n);
  int aaa = anum(n);
  if (zzz > 0 && aaa == 0)
    return true;  // is element
  return false;
}

/**********************/
/*** name functions ***/
/**********************/
std::string pyne::nucname::name(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  newnuc += zz_name[zzz];

  // Add A-number
  if (0 < aaa)
    newnuc += pyne::to_str(aaa);

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "M";

  return newnuc;
}



std::string pyne::nucname::name(const char* nuc) {
  std::string newnuc(nuc);
  return name(newnuc);
}


std::string pyne::nucname::name(std::string nuc) {
  return name(id(nuc));
}


/**********************/
/*** znum functions ***/
/**********************/
int pyne::nucname::znum(int nuc) {
  return id(nuc) / 10000000;
}

int pyne::nucname::znum(const char* nuc) {
  return id(nuc) / 10000000;
}

int pyne::nucname::znum(std::string nuc) {
  return id(nuc) / 10000000;
}

/**********************/
/*** anum functions ***/
/**********************/
int pyne::nucname::anum(int nuc) {
  return (id(nuc) / 10000) % 1000;
}

int pyne::nucname::anum(const char* nuc) {
  return (id(nuc) / 10000) % 1000;
}

int pyne::nucname::anum(std::string nuc) {
  return (id(nuc) / 10000) % 1000;
}

/**********************/
/*** snum functions ***/
/**********************/
int pyne::nucname::snum(int nuc) {
  return id(nuc) % 10000;
}

int pyne::nucname::snum(const char* nuc) {
  return id(nuc) % 10000;
}

int pyne::nucname::snum(std::string nuc) {
  return id(nuc) % 10000;
}

/************************/
/*** zzaaam functions ***/
/************************/
int pyne::nucname::zzaaam(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;
  int ssss = nucid % 10000;
  if (10 <= ssss)
    ssss = 9;
  return zzzaaa * 10 + ssss;
}


int pyne::nucname::zzaaam(const char* nuc) {
  std::string newnuc(nuc);
  return zzaaam(newnuc);
}


int pyne::nucname::zzaaam(std::string nuc) {
  return zzaaam(id(nuc));
}


int pyne::nucname::zzaaam_to_id(int nuc) {
  return (nuc / 10) * 10000 + (nuc % 10);
}


int pyne::nucname::zzaaam_to_id(const char* nuc) {
  return zzaaam_to_id(std::string(nuc));
}


int pyne::nucname::zzaaam_to_id(std::string nuc) {
  return zzaaam_to_id(pyne::to_int(nuc));
}

/************************/
/*** zzzaaa functions ***/
/************************/
int pyne::nucname::zzzaaa(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;

  return zzzaaa;
}


int pyne::nucname::zzzaaa(const char* nuc) {
  std::string newnuc(nuc);
  return zzzaaa(newnuc);
}


int pyne::nucname::zzzaaa(std::string nuc) {
  return zzzaaa(id(nuc));
}


int pyne::nucname::zzzaaa_to_id(int nuc) {
  return (nuc) * 10000;
}


int pyne::nucname::zzzaaa_to_id(const char* nuc) {
  return zzzaaa_to_id(std::string(nuc));
}


int pyne::nucname::zzzaaa_to_id(std::string nuc) {
  return zzzaaa_to_id(pyne::to_int(nuc));
}

/*************************/
/*** zzllaaam functions ***/
/*************************/
std::string pyne::nucname::zzllaaam(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int zzz = nucid / 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);
  //Adding ZZ
  newnuc += pyne::to_str(zzz);
  newnuc += "-";
  // Add LL
  newnuc += zz_name[zzz];
  // Add required dash
  newnuc += "-";
  // Add AAA
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);
  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "m";
  return newnuc;
}


std::string pyne::nucname::zzllaaam(const char* nuc) {
  std::string newnuc(nuc);
  return zzllaaam(newnuc);
}


std::string pyne::nucname::zzllaaam(std::string nuc) {
  return zzllaaam(id(nuc));
}


int pyne::nucname::zzllaaam_to_id(const char* nuc) {
  return zzllaaam_to_id(std::string(nuc));
}


int pyne::nucname::zzllaaam_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  std::string elem_name;

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  // Removing first two characters (redundant), for 1 digit nuclides, such
  // as 2-He-4, the first slash will be removed, and the second attempt to
  // remove the second slash will do nothing.
  nucstr.erase(0, 2);
  nucstr = pyne::remove_substring(nucstr, "-");
  // Does nothing if nuclide is short, otherwise removes the second "-" instance
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty() || pyne::contains_substring(nucstr, "NAT")) {
    elem_name = pyne::capitalize(pyne::remove_substring(nucstr, "NAT"));
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  int anum = pyne::to_int(anum_str);

  // Figure out if we are meta-stable or not
  std::string end_char = pyne::last_char(nucstr);
  if (end_char == "M")
    nucid = (10000 * anum) + 1;
  else if (pyne::contains_substring(pyne::digits, end_char))
    nucid = (10000 * anum);
  else
    throw NotANuclide(nucstr, nucid);

  // Add the Z-number
  elem_name = pyne::remove_characters(nucstr.substr(0, nuclen - 1), pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nucstr, nucid);
  return nucid;
}

/**********************/
/*** mcnp functions ***/
/**********************/
int pyne::nucname::mcnp(int nuc) {
  nuc = id(nuc);
  int ssss = nuc % 10000;
  int newnuc = nuc / 10000;

  // special case Am242(m)
  if (newnuc == 95242 && ssss < 2)
    ssss = (ssss + 1) % 2;

  // Handle the crazy MCNP meta-stable format
  if (0 != ssss && ssss < 10)
    newnuc += 300 + (ssss * 100);

  return newnuc;
}



int pyne::nucname::mcnp(const char* nuc) {
  std::string newnuc(nuc);
  return mcnp(newnuc);
}



int pyne::nucname::mcnp(std::string nuc) {
  return mcnp(id(nuc));
}

//
// MCNP -> id
//
int pyne::nucname::mcnp_to_id(int nuc) {
  int zzz = nuc / 1000;
  int aaa = nuc % 1000;
  if (zzz == 0)
    throw NotANuclide(nuc, "not in the MCNP format");
  else if (zzz <= aaa) {
    if (aaa - 400 < 0) {
      if (nuc == 95242)
        return nuc * 10000 + 1;  // special case MCNP Am-242m
      else
        return nuc * 10000;  // Nuclide in normal MCNP form
    } else {
      // Nuclide in MCNP metastable form
      if (nuc == 95642)
        return (95642 - 400) * 10000; // special case MCNP Am-242
      nuc = ((nuc - 400) * 10000) + 1;
      while (3.0 < (float ((nuc / 10000) % 1000) / float (nuc / 10000000)))
        nuc -= 999999;
      return nuc;
    }
  } else if (aaa == 0)
    // MCNP form natural nuclide
    return zzz * 10000000;
  throw IndeterminateNuclideForm(nuc, "");
}


int pyne::nucname::mcnp_to_id(const char* nuc) {
  return mcnp_to_id(std::string(nuc));
}


int pyne::nucname::mcnp_to_id(std::string nuc) {
  return mcnp_to_id(pyne::to_int(nuc));
}

/************************/
/*** openmc functions ***/
/************************/
std::string pyne::nucname::openmc(int nuc) {
  std::string nucname = name(nuc);

  // check aaa value
  if (iselement(nuc)) {
    nucname.append("0");
  }

  // format metadata
  if ('M' == nucname.back()) {
    nucname.back() = '_';
    nucname.append("m");
    int meta_id = snum(nuc);
    std::string meta_str = std::to_string(meta_id);
    nucname.append(meta_str);
  }
  return nucname;
}

std::string pyne::nucname::openmc(const char* nuc) {
  std::string newnuc(nuc);
  return openmc(newnuc);
}

std::string pyne::nucname::openmc(std::string nuc) {
  return openmc(id(nuc));
}

//
// OPENMC -> id
//
int pyne::nucname::openmc_to_id(const char* nuc) {
  return openmc_to_id(std::string(nuc));
}

int pyne::nucname::openmc_to_id(std::string nuc) {
  std::string nucname;
  name_zz_t zznames = get_name_zz();

  // first two characters
  std::string::iterator aaa_start;
  int zzz = 0;
  if (zznames.count(nuc.substr(0, 2)) == 1) {
    aaa_start = nuc.begin() + 2;
    zzz = zznames[nuc.substr(0, 2)];
  }
  // then try only the first
  else if (zznames.count(nuc.substr(0, 1)) == 1) {
    aaa_start = nuc.begin() + 1;
    zzz = zznames[nuc.substr(0, 1)];
  } else {
    throw NotANuclide(nuc, "Not in the OpenMC format");
  }

  // set aaa - stop on "-" if the character exists
  std::string::iterator aaa_end = std::find(nuc.begin(), nuc.end(), '_');
  int aaa = pyne::to_int(nuc.substr(aaa_start - nuc.begin(), aaa_end - aaa_start));

  // check for metastable state
  int m = 0;
  if (aaa_end != nuc.end()) {
    std::string::iterator m_start = aaa_end + 2; // move forward once to skip "_m" characters
    m = pyne::to_int(nuc.substr(m_start - nuc.begin(), nuc.end() - m_start));
  }

  // form integer id and return
  return (zzz * 10000000) + (aaa * 10000) + m;

}


/**********************/
/*** fluka functions ***/
/**********************/
std::string pyne::nucname::fluka(int nuc) {
  int x = id(nuc);
  if (zz_fluka.count(x) == 0) {
    throw NotANuclide(nuc, "fluka name could not be found");
  }
  return zz_fluka[x];
}


//
// FLUKA name -> id
//
int pyne::nucname::fluka_to_id(std::string name) {
  if (fluka_zz.count(name) == 0) {
    throw NotANuclide(-1, "No nuclide: fluka name could not be found");
  }
  return fluka_zz[name];
}

int pyne::nucname::fluka_to_id(char* name) {
  return fluka_to_id(std::string(name));
}


/*************************/
/*** serpent functions ***/
/*************************/
std::string pyne::nucname::serpent(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int zzz = nucid / 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  std::string llupper = pyne::to_upper(zz_name[zzz]);
  std::string lllower = pyne::to_lower(zz_name[zzz]);
  newnuc += llupper[0];
  for (int l = 1; l < lllower.size(); l++)
    newnuc += lllower[l];

  // Add required dash
  newnuc += "-";

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);
  else if (0 == aaassss)
    newnuc += "nat";

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "m";

  return newnuc;
}


std::string pyne::nucname::serpent(const char* nuc) {
  std::string newnuc(nuc);
  return serpent(newnuc);
}


std::string pyne::nucname::serpent(std::string nuc) {
  return serpent(id(nuc));
}

//
// Serpent -> id
//
//int pyne::nucname::serpent_to_id(int nuc)
//{
// Should be ZAID
//}


int pyne::nucname::serpent_to_id(const char* nuc) {
  return serpent_to_id(std::string(nuc));
}


int pyne::nucname::serpent_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  std::string elem_name;

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty() || pyne::contains_substring(nucstr, "NAT")) {
    elem_name = pyne::capitalize(pyne::remove_substring(nucstr, "NAT"));
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  int anum = pyne::to_int(anum_str);

  // Figure out if we are meta-stable or not
  std::string end_char = pyne::last_char(nucstr);
  if (end_char == "M")
    nucid = (10000 * anum) + 1;
  else if (pyne::contains_substring(pyne::digits, end_char))
    nucid = (10000 * anum);
  else
    throw NotANuclide(nucstr, nucid);

  // Add the Z-number
  elem_name = pyne::remove_characters(nucstr.substr(0, nuclen - 1), pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nucstr, nucid);
  return nucid;
}


/**********************/
/*** nist functions ***/
/**********************/
std::string pyne::nucname::nist(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);

  // Add name
  std::string name_upper = pyne::to_upper(zz_name[zzz]);
  std::string name_lower = pyne::to_lower(zz_name[zzz]);
  newnuc += name_upper[0];
  for (int l = 1; l < name_lower.size(); l++)
    newnuc += name_lower[l];

  // Add meta-stable flag
  // No metastable flag for NIST,
  // but could add star, by uncommenting below
  //if (0 < mod_10)
  //  newnuc += "*";

  return newnuc;
}


std::string pyne::nucname::nist(const char* nuc) {
  std::string newnuc(nuc);
  return nist(newnuc);
}


std::string pyne::nucname::nist(std::string nuc) {
  return nist(id(nuc));
}


//
// NIST -> id
//
//int pyne::nucname::nist_to_id(int nuc)
//{
// NON-EXISTANT
//};

int pyne::nucname::nist_to_id(const char* nuc) {
  return nist_to_id(std::string(nuc));
}

int pyne::nucname::nist_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  nuc = pyne::to_upper(nuc);
  std::string elem_name;
  int nuclen = nuc.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nuc, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty()) {
    elem_name = pyne::capitalize(nuc);
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  nucid = pyne::to_int(anum_str) * 10000;

  // Add the Z-number
  elem_name = pyne::remove_characters(nuc, pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nuc, nucid);
  return nucid;
}


/************************/
/*** cinder functions ***/
/************************/
int pyne::nucname::cinder(int nuc) {
  // cinder nuclides of form aaazzzm
  int nucid = id(nuc);
  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;
  if (10 <= ssss)
    ssss = 9;
  return (aaa * 10000) + (zzz * 10) + ssss;
}



int pyne::nucname::cinder(const char* nuc) {
  std::string newnuc(nuc);
  return cinder(newnuc);
}



int pyne::nucname::cinder(std::string nuc) {
  return cinder(id(nuc));
}

//
// Cinder -> Id
//
int pyne::nucname::cinder_to_id(int nuc) {
  int ssss = nuc % 10;
  int aaazzz = nuc / 10;
  int zzz = aaazzz % 1000;
  int aaa = aaazzz / 1000;
  return (zzz * 10000000) + (aaa * 10000) + ssss;
}


int pyne::nucname::cinder_to_id(const char* nuc) {
  return cinder_to_id(std::string(nuc));
}


int pyne::nucname::cinder_to_id(std::string nuc) {
  return cinder_to_id(pyne::to_int(nuc));
}




/**********************/
/*** ALARA functions ***/
/**********************/
std::string pyne::nucname::alara(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";
  std::string ll = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL, in lower case
  ll += zz_name[zzz];

  for (int i = 0; ll[i] != '\0'; i++)
    ll[i] = tolower(ll[i]);
  newnuc += ll;

  // Add A-number
  if (0 < aaassss) {
    newnuc += ":";
    newnuc += pyne::to_str(aaa);
  }

  // Note, ALARA input format does not use metastable flag
  return newnuc;
}


std::string pyne::nucname::alara(const char* nuc) {
  std::string newnuc(nuc);
  return alara(newnuc);
}


std::string pyne::nucname::alara(std::string nuc) {
  return alara(id(nuc));
}


//
// Cinder -> Id
//
//int pyne::nucname::alara_to_id(int nuc)
//{
// Not Possible
//}


int pyne::nucname::alara_to_id(const char* nuc) {
  return alara_to_id(std::string(nuc));
}


int pyne::nucname::alara_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  nuc = pyne::to_upper(pyne::remove_characters(nuc, ":"));
  std::string elem_name;
  int nuclen = nuc.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nuc, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty()) {
    elem_name = pyne::capitalize(nuc);
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  nucid = pyne::to_int(anum_str) * 10000;

  // Add the Z-number
  elem_name = pyne::remove_characters(nuc, pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nuc, nucid);
  return nucid;
}




/***********************/
/***  SZA functions  ***/
/***********************/
int pyne::nucname::sza(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;
  int sss = nucid % 10000;
  return sss * 1000000 + zzzaaa;
}


int pyne::nucname::sza(const char* nuc) {
  std::string newnuc(nuc);
  return sza(newnuc);
}


int pyne::nucname::sza(std::string nuc) {
  return sza(id(nuc));
}


int pyne::nucname::sza_to_id(int nuc) {
  int sss = nuc / 1000000;
  int zzzaaa = nuc % 1000000;
  if (5 < sss) {
    // Unphysical metastable state warning
    warning("You have indicated a metastable state of " + pyne::to_str(sss) + ". Metastable state above 5, possibly unphysical. ");
  }
  return zzzaaa * 10000 + sss;
}


int pyne::nucname::sza_to_id(const char* nuc) {
  std::string newnuc(nuc);
  return sza_to_id(newnuc);
}


int pyne::nucname::sza_to_id(std::string nuc) {
  return sza_to_id(pyne::to_int(nuc));
}


void pyne::nucname::_load_state_map() {
  for (int i = 0; i < TOTAL_STATE_MAPS; ++i) {
    state_id_map[map_nuc_ids[i]] = map_metastable[i];
  }
}

int pyne::nucname::state_id_to_id(int state) {
  int zzzaaa = (state / 10000) * 10000;
  int state_number = state % 10000;
  if (state_number == 0)
    return state;
  std::map<int, int>::iterator nuc_iter, nuc_end;

  nuc_iter = state_id_map.find(state);
  nuc_end = state_id_map.end();
  if (nuc_iter != nuc_end) {
    int m = (*nuc_iter).second;
    return zzzaaa + m;
  }

  if (state_id_map.empty())  {
    _load_state_map();
    return state_id_to_id(state);
  }
  return -1;
}


int pyne::nucname::id_to_state_id(int nuc_id) {
  int zzzaaa = (nuc_id / 10000) * 10000;
  int state = nuc_id % 10000;
  if (state == 0)
    return nuc_id;
  std::map<int, int>::iterator nuc_iter, nuc_end, it;

  nuc_iter = state_id_map.lower_bound(zzzaaa);
  nuc_end = state_id_map.upper_bound(zzzaaa + 9999);
  for (it = nuc_iter; it != nuc_end; ++it) {
    if (state == it->second) {
      return it->first;
    }
  }

  if (state_id_map.empty())  {
    _load_state_map();
    return id_to_state_id(nuc_id);
  }
  return -1;
}


/************************/
/*** ENSDF functions ***/
/************************/
//
// ENSDF  -> Id
//

int pyne::nucname::ensdf_to_id(const char* nuc) {
  return ensdf_to_id(std::string(nuc));
}

int pyne::nucname::ensdf_to_id(std::string nuc) {
  if (nuc.size() < 4) {
    return nucname::id(nuc);
  } else if (std::isdigit(nuc[3])) {
    int aaa = to_int(nuc.substr(0, 3));
    int zzz;
    std::string xx_str = nuc.substr(3, 2);
    zzz = to_int(xx_str) + 100;
    int nid = 10000 * aaa + 10000000 * zzz;
    return nid;
  } else {
    return nucname::id(nuc);
  }

}

//
// end of src/nucname.cpp
//


//
// start of src/rxname.cpp
//
#ifndef PYNE_IS_AMALGAMATED
#include "rxname.h"
#endif

std::string pyne::rxname::_names[NUM_RX_NAMES] = {
  "total",
  "scattering",
  "elastic",
  "nonelastic",
  "n",
  "misc",
  "continuum",
  "z_2nd",
  "z_2n",
  "z_2n_0",
  "z_2n_1",
  "z_2n_2",
  "z_3n",
  "z_3n_0",
  "z_3n_1",
  "z_3n_2",
  "fission",
  "fission_first",
  "fission_second",
  "fission_third",
  "na",
  "na_0",
  "na_1",
  "na_2",
  "n3a",
  "z_2na",
  "z_3na",
  "absorption",
  "np",
  "np_0",
  "np_1",
  "np_2",
  "npd",
  "n2a",
  "z_2n2a",
  "nd",
  "nd_0",
  "nd_1",
  "nd_2",
  "nt",
  "nt_0",
  "nt_1",
  "nt_2",
  "nHe3",
  "nHe3_0",
  "nHe3_1",
  "nHe3_2",
  "nd3a",
  "nt2a",
  "z_4n",
  "z_4n_0",
  "z_4n_1",
  "fission_fourth",
  "z_2np",
  "z_3np",
  "n2p",
  "npa",
  "n_0",
  "n_1",
  "n_2",
  "n_3",
  "n_4",
  "n_5",
  "n_6",
  "n_7",
  "n_8",
  "n_9",
  "n_10",
  "n_11",
  "n_12",
  "n_13",
  "n_14",
  "n_15",
  "n_16",
  "n_17",
  "n_18",
  "n_19",
  "n_20",
  "n_21",
  "n_22",
  "n_23",
  "n_24",
  "n_25",
  "n_26",
  "n_27",
  "n_28",
  "n_29",
  "n_30",
  "n_31",
  "n_32",
  "n_33",
  "n_34",
  "n_35",
  "n_36",
  "n_37",
  "n_38",
  "n_39",
  "n_40",
  "n_continuum",
  "disappearance",
  "gamma",
  "gamma_0",
  "gamma_1",
  "gamma_2",
  "p",
  "d",
  "t",
  "He3",
  "a",
  "z_2a",
  "z_3a",
  "z_2p",
  "z_2p_0",
  "z_2p_1",
  "z_2p_2",
  "pa",
  "t2a",
  "d2a",
  "pd",
  "pt",
  "da",
  "resonance_parameters",
  "n_total",
  "gamma_total",
  "p_total",
  "d_total",
  "t_total",
  "He3_total",
  "a_total",
  "pionp",
  "pion0",
  "pionm",
  "muonp",
  "muonm",
  "kaonp",
  "kaon0_long",
  "kaon0_short",
  "kaonm",
  "antip",
  "antin",
  "mubar",
  "epsilon",
  "y",
  "erel_total",
  "erel_elastic",
  "erel_nonelastic",
  "erel_n",
  "erel_misc",
  "erel_continuum",
  "erel_2nd",
  "erel_2n",
  "erel_3n",
  "erel_fission",
  "erel_fission_first",
  "erel_fission_second",
  "erel_fission_third",
  "erel_na",
  "erel_n3a",
  "erel_2na",
  "erel_3na",
  "erel_absorption",
  "erel_np",
  "erel_n2a",
  "erel_2n2a",
  "erel_nd",
  "erel_nt",
  "erel_nHe3",
  "erel_nd3a",
  "erel_nt2a",
  "erel_4n",
  "erel_fission_fourth",
  "erel_2np",
  "erel_3np",
  "erel_n2p",
  "erel_npa",
  "erel_n_0",
  "erel_n_1",
  "erel_n_2",
  "erel_n_3",
  "erel_n_4",
  "erel_n_5",
  "erel_n_6",
  "erel_n_7",
  "erel_n_8",
  "erel_n_9",
  "erel_n_10",
  "erel_n_11",
  "erel_n_12",
  "erel_n_13",
  "erel_n_14",
  "erel_n_15",
  "erel_n_16",
  "erel_n_17",
  "erel_n_18",
  "erel_n_19",
  "erel_n_20",
  "erel_n_21",
  "erel_n_22",
  "erel_n_23",
  "erel_n_24",
  "erel_n_25",
  "erel_n_26",
  "erel_n_27",
  "erel_n_28",
  "erel_n_29",
  "erel_n_30",
  "erel_n_31",
  "erel_n_32",
  "erel_n_33",
  "erel_n_34",
  "erel_n_35",
  "erel_n_36",
  "erel_n_37",
  "erel_n_38",
  "erel_n_39",
  "erel_n_40",
  "erel_n_continuum",
  "erel_disappearance",
  "erel_gamma",
  "erel_p",
  "erel_d",
  "erel_t",
  "erel_He3",
  "erel_a",
  "erel_2a",
  "erel_3a",
  "erel_2p",
  "erel_pa",
  "erel_t2a",
  "erel_d2a",
  "erel_pd",
  "erel_pt",
  "erel_da",
  "damage",
  "heading",
  "nubar",
  "fission_product_yield_independent",
  "nubar_delayed",
  "nubar_prompt",
  "decay",
  "energy_per_fission",
  "fission_product_yield_cumulative",
  "gamma_delayed",
  "stopping_power",
  "photon_total",
  "photon_coherent",
  "photon_incoherent",
  "scattering_factor_imag",
  "scattering_factor_real",
  "pair_prod_elec",
  "pair_prod",
  "pair_prod_nuc",
  "absorption_photoelectric",
  "photoexcitation",
  "scattering_electroatomic",
  "bremsstrahlung",
  "excitation_electroatomic",
  "atomic_relaxation",
  "k_photoelectric",
  "l1_photoelectric",
  "l2_photoelectric",
  "l3_photoelectric",
  "m1_photoelectric",
  "m2_photoelectric",
  "m3_photoelectric",
  "m4_photoelectric",
  "m5_photoelectric",
  "n1_photoelectric",
  "n2_photoelectric",
  "n3_photoelectric",
  "n4_photoelectric",
  "n5_photoelectric",
  "n6_photoelectric",
  "n7_photoelectric",
  "o1_photoelectric",
  "o2_photoelectric",
  "o3_photoelectric",
  "o4_photoelectric",
  "o5_photoelectric",
  "o6_photoelectric",
  "o7_photoelectric",
  "o8_photoelectric",
  "o9_photoelectric",
  "p1_photoelectric",
  "p2_photoelectric",
  "p3_photoelectric",
  "p4_photoelectric",
  "p5_photoelectric",
  "p6_photoelectric",
  "p7_photoelectric",
  "p8_photoelectric",
  "p9_photoelectric",
  "p10_photoelectric",
  "p11_photoelectric",
  "q1_photoelectric",
  "q2_photoelectric",
  "q3_photoelectric",
  "p_0",
  "p_1",
  "p_2",
  "p_3",
  "p_4",
  "p_5",
  "p_6",
  "p_7",
  "p_8",
  "p_9",
  "p_10",
  "p_11",
  "p_12",
  "p_13",
  "p_14",
  "p_15",
  "p_16",
  "p_17",
  "p_18",
  "p_19",
  "p_20",
  "p_21",
  "p_22",
  "p_23",
  "p_24",
  "p_25",
  "p_26",
  "p_27",
  "p_28",
  "p_29",
  "p_30",
  "p_31",
  "p_32",
  "p_33",
  "p_34",
  "p_35",
  "p_36",
  "p_37",
  "p_38",
  "p_39",
  "p_40",
  "p_41",
  "p_42",
  "p_43",
  "p_44",
  "p_45",
  "p_46",
  "p_47",
  "p_48",
  "p_continuum",
  "d_0",
  "d_1",
  "d_2",
  "d_3",
  "d_4",
  "d_5",
  "d_6",
  "d_7",
  "d_8",
  "d_9",
  "d_10",
  "d_11",
  "d_12",
  "d_13",
  "d_14",
  "d_15",
  "d_16",
  "d_17",
  "d_18",
  "d_19",
  "d_20",
  "d_21",
  "d_22",
  "d_23",
  "d_24",
  "d_25",
  "d_26",
  "d_27",
  "d_28",
  "d_29",
  "d_30",
  "d_31",
  "d_32",
  "d_33",
  "d_34",
  "d_35",
  "d_36",
  "d_37",
  "d_38",
  "d_39",
  "d_40",
  "d_41",
  "d_42",
  "d_43",
  "d_44",
  "d_45",
  "d_46",
  "d_47",
  "d_48",
  "d_continuum",
  "t_0",
  "t_1",
  "t_2",
  "t_3",
  "t_4",
  "t_5",
  "t_6",
  "t_7",
  "t_8",
  "t_9",
  "t_10",
  "t_11",
  "t_12",
  "t_13",
  "t_14",
  "t_15",
  "t_16",
  "t_17",
  "t_18",
  "t_19",
  "t_20",
  "t_21",
  "t_22",
  "t_23",
  "t_24",
  "t_25",
  "t_26",
  "t_27",
  "t_28",
  "t_29",
  "t_30",
  "t_31",
  "t_32",
  "t_33",
  "t_34",
  "t_35",
  "t_36",
  "t_37",
  "t_38",
  "t_39",
  "t_40",
  "t_41",
  "t_42",
  "t_43",
  "t_44",
  "t_45",
  "t_46",
  "t_47",
  "t_48",
  "t_continuum",
  "He3_0",
  "He3_1",
  "He3_2",
  "He3_3",
  "He3_4",
  "He3_5",
  "He3_6",
  "He3_7",
  "He3_8",
  "He3_9",
  "He3_10",
  "He3_11",
  "He3_12",
  "He3_13",
  "He3_14",
  "He3_15",
  "He3_16",
  "He3_17",
  "He3_18",
  "He3_19",
  "He3_20",
  "He3_21",
  "He3_22",
  "He3_23",
  "He3_24",
  "He3_25",
  "He3_26",
  "He3_27",
  "He3_28",
  "He3_29",
  "He3_30",
  "He3_31",
  "He3_32",
  "He3_33",
  "He3_34",
  "He3_35",
  "He3_36",
  "He3_37",
  "He3_38",
  "He3_39",
  "He3_40",
  "He3_41",
  "He3_42",
  "He3_43",
  "He3_44",
  "He3_45",
  "He3_46",
  "He3_47",
  "He3_48",
  "He3_continuum",
  "a_0",
  "a_1",
  "a_2",
  "a_3",
  "a_4",
  "a_5",
  "a_6",
  "a_7",
  "a_8",
  "a_9",
  "a_10",
  "a_11",
  "a_12",
  "a_13",
  "a_14",
  "a_15",
  "a_16",
  "a_17",
  "a_18",
  "a_19",
  "a_20",
  "a_21",
  "a_22",
  "a_23",
  "a_24",
  "a_25",
  "a_26",
  "a_27",
  "a_28",
  "a_29",
  "a_30",
  "a_31",
  "a_32",
  "a_33",
  "a_34",
  "a_35",
  "a_36",
  "a_37",
  "a_38",
  "a_39",
  "a_40",
  "a_41",
  "a_42",
  "a_43",
  "a_44",
  "a_45",
  "a_46",
  "a_47",
  "a_48",
  "a_continuum",
  "lumped_covar",
  "excited",
  "bminus",
  "bplus",
  "ec",
  "bminus_n",
  "bminus_a",
  "it",
  "bplus_a",
  "ec_bplus",
  "bplus_p",
  "bminus_2n",
  "bminus_3n",
  "bminus_4n",
  "ecp",
  "eca",
  "bplus_2p",
  "ec_2p",
  "decay_2bminus",
  "bminus_p",
  "decay_14c",
  "bplus_3p",
  "sf",
  "decay_2bplus",
  "decay_2ec",
  "ec_3p",
  "bminus_sf"
};
std::set<std::string> pyne::rxname::names(pyne::rxname::_names,
                                          pyne::rxname::_names + NUM_RX_NAMES);


std::map<std::string, unsigned int> pyne::rxname::altnames;
std::map<unsigned int, std::string> pyne::rxname::id_name;
std::map<std::string, unsigned int> pyne::rxname::name_id;
std::map<unsigned int, unsigned int> pyne::rxname::id_mt;
std::map<unsigned int, unsigned int> pyne::rxname::mt_id;
std::map<unsigned int, std::string> pyne::rxname::labels;
std::map<unsigned int, std::string> pyne::rxname::docs;
std::map<std::pair<std::string, int>, unsigned int> pyne::rxname::offset_id;
std::map<std::pair<std::string, unsigned int>, int> pyne::rxname::id_offset;

void* pyne::rxname::_fill_maps() {
  using std::make_pair;
  std::string rx;
  unsigned int rxid;
  unsigned int _mts [NUM_RX_NAMES] = {
    1,
    0,
    2,
    3,
    4,
    5,
    10,
    11,
    16,
    0,
    0,
    0,
    17,
    0,
    0,
    0,
    18,
    19,
    20,
    21,
    22,
    0,
    0,
    0,
    23,
    24,
    25,
    27,
    28,
    0,
    0,
    0,
    0,
    29,
    30,
    32,
    0,
    0,
    0,
    33,
    0,
    0,
    0,
    34,
    0,
    0,
    0,
    35,
    36,
    37,
    0,
    0,
    38,
    41,
    42,
    44,
    45,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    91,
    101,
    102,
    0,
    0,
    0,
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    111,
    0,
    0,
    0,
    112,
    113,
    114,
    115,
    116,
    117,
    151,
    201,
    202,
    203,
    204,
    205,
    206,
    207,
    208,
    209,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    218,
    251,
    252,
    253,
    301,
    302,
    303,
    304,
    305,
    310,
    311,
    316,
    317,
    318,
    319,
    320,
    321,
    322,
    323,
    324,
    325,
    327,
    328,
    329,
    330,
    332,
    333,
    334,
    335,
    336,
    337,
    338,
    341,
    342,
    344,
    345,
    350,
    351,
    352,
    353,
    354,
    355,
    356,
    357,
    358,
    359,
    360,
    361,
    362,
    363,
    364,
    365,
    366,
    367,
    368,
    369,
    370,
    371,
    372,
    373,
    374,
    375,
    376,
    377,
    378,
    379,
    380,
    381,
    382,
    383,
    384,
    385,
    386,
    387,
    388,
    389,
    390,
    391,
    401,
    402,
    403,
    404,
    405,
    406,
    407,
    408,
    409,
    411,
    412,
    413,
    414,
    415,
    416,
    417,
    444,
    451,
    452,
    454,
    455,
    456,
    457,
    458,
    459,
    460,
    500,
    501,
    502,
    504,
    505,
    506,
    515,
    516,
    517,
    522,
    523,
    526,
    527,
    528,
    533,
    534,
    535,
    536,
    537,
    538,
    539,
    540,
    541,
    542,
    543,
    544,
    545,
    546,
    547,
    548,
    549,
    550,
    551,
    552,
    553,
    554,
    555,
    556,
    557,
    558,
    559,
    560,
    561,
    562,
    563,
    564,
    565,
    566,
    567,
    568,
    569,
    570,
    571,
    572,
    600,
    601,
    602,
    603,
    604,
    605,
    606,
    607,
    608,
    609,
    610,
    611,
    612,
    613,
    614,
    615,
    616,
    617,
    618,
    619,
    620,
    621,
    622,
    623,
    624,
    625,
    626,
    627,
    628,
    629,
    630,
    631,
    632,
    633,
    634,
    635,
    636,
    637,
    638,
    639,
    640,
    641,
    642,
    643,
    644,
    645,
    646,
    647,
    648,
    649,
    650,
    651,
    652,
    653,
    654,
    655,
    656,
    657,
    658,
    659,
    660,
    661,
    662,
    663,
    664,
    665,
    666,
    667,
    668,
    669,
    670,
    671,
    672,
    673,
    674,
    675,
    676,
    677,
    678,
    679,
    680,
    681,
    682,
    683,
    684,
    685,
    686,
    687,
    688,
    689,
    690,
    691,
    692,
    693,
    694,
    695,
    696,
    697,
    698,
    699,
    700,
    701,
    702,
    703,
    704,
    705,
    706,
    707,
    708,
    709,
    710,
    711,
    712,
    713,
    714,
    715,
    716,
    717,
    718,
    719,
    720,
    721,
    722,
    723,
    724,
    725,
    726,
    727,
    728,
    729,
    730,
    731,
    732,
    733,
    734,
    735,
    736,
    737,
    738,
    739,
    740,
    741,
    742,
    743,
    744,
    745,
    746,
    747,
    748,
    749,
    750,
    751,
    752,
    753,
    754,
    755,
    756,
    757,
    758,
    759,
    760,
    761,
    762,
    763,
    764,
    765,
    766,
    767,
    768,
    769,
    770,
    771,
    772,
    773,
    774,
    775,
    776,
    777,
    778,
    779,
    780,
    781,
    782,
    783,
    784,
    785,
    786,
    787,
    788,
    789,
    790,
    791,
    792,
    793,
    794,
    795,
    796,
    797,
    798,
    799,
    800,
    801,
    802,
    803,
    804,
    805,
    806,
    807,
    808,
    809,
    810,
    811,
    812,
    813,
    814,
    815,
    816,
    817,
    818,
    819,
    820,
    821,
    822,
    823,
    824,
    825,
    826,
    827,
    828,
    829,
    830,
    831,
    832,
    833,
    834,
    835,
    836,
    837,
    838,
    839,
    840,
    841,
    842,
    843,
    844,
    845,
    846,
    847,
    848,
    849,
    851,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
  std::string _labels[NUM_RX_NAMES] = {
    "(z,total)",
    "(z,scattering)",
    "(z,elastic)",
    "(z,nonelastic)",
    "(z,n)",
    "(misc)",
    "(z,continuum)",
    "(z,2nd)",
    "(z,2n)",
    "(z,2n0)",
    "(z,2n1)",
    "(z,2n2)",
    "(z,3n)",
    "(z,3n0)",
    "(z,3n1)",
    "(z,3n2)",
    "(z,fission)",
    "(z,f)",
    "(z,nf)",
    "(z,2nf)",
    "(z,n+a)",
    "(z,n+a0)",
    "(z,n+a1)",
    "(z,n+a2)",
    "(z,n+3a)",
    "(z,2n+a)",
    "(z,3n+a)",
    "(z,abs) Absorption",
    "(z,n+p)",
    "(z,n+p0)",
    "(z,n+p1)",
    "(z,n+p2)",
    "(z,n+p+d)",
    "(z,n+2a)",
    "(z,2n+2a)",
    "(z,nd)",
    "(z,nd0)",
    "(z,nd1)",
    "(z,nd2)",
    "(z,nt)",
    "(z,nt0)",
    "(z,nt1)",
    "(z,nt2)",
    "(z,n+He3)",
    "(z,n+He3-0)",
    "(z,n+He3-1)",
    "(z,n+He3-2)",
    "(z,n+d+3a)",
    "(z,n+t+2a)",
    "(z,4n)",
    "(z,4n0)",
    "(z,4n1)",
    "(z,3nf)",
    "(z,2n+p)",
    "(z,3n+p)",
    "(z,n+2p)",
    "(z,npa)",
    "(z,n0)",
    "(z,n1)",
    "(z,n2)",
    "(z,n3)",
    "(z,n4)",
    "(z,n5)",
    "(z,n6)",
    "(z,n7)",
    "(z,n8)",
    "(z,n9)",
    "(z,n10)",
    "(z,n11)",
    "(z,n12)",
    "(z,n13)",
    "(z,n14)",
    "(z,n15)",
    "(z,n16)",
    "(z,n17)",
    "(z,n18)",
    "(z,n19)",
    "(z,n20)",
    "(z,n21)",
    "(z,n22)",
    "(z,n23)",
    "(z,n24)",
    "(z,n25)",
    "(z,n26)",
    "(z,n27)",
    "(z,n28)",
    "(z,n29)",
    "(z,n30)",
    "(z,n31)",
    "(z,n32)",
    "(z,n33)",
    "(z,n34)",
    "(z,n35)",
    "(z,n36)",
    "(z,n37)",
    "(z,n38)",
    "(z,n39)",
    "(z,n40)",
    "(z,nc)",
    "(z,disap) Neutron disappearance",
    "(z,gamma)",
    "(z,gamma0)",
    "(z,gamma1)",
    "(z,gamma2)",
    "(z,p)",
    "(z,d)",
    "(z,t)",
    "(z,3He)",
    "(z,a)",
    "(z,2a)",
    "(z,3a)",
    "(z,2p)",
    "(z,2p0)",
    "(z,2p1)",
    "(z,2p2)",
    "(z,pa)",
    "(z,t2a)",
    "(z,d2a)",
    "(z,pd)",
    "(z,pt)",
    "(z,da)",
    "Resonance Parameters",
    "(z,Xn)",
    "(z,Xgamma)",
    "(z,Xp)",
    "(z,Xd)",
    "(z,Xt)",
    "(z,X3He)",
    "(z,Xa)",
    "(z,Xpi+) Total pi+ meson production",
    "(z,Xpi0) Total pi0 meson production",
    "(z,Xpi-) Total pi- meson production",
    "(z,Xmu+) Total anti-muon production",
    "(z,Xmu-) Total muon production",
    "(z,Xk+) Total positive kaon production",
    "(z,Xk0long) Total long-lived neutral kaon production",
    "(z,Xk0short) Total short-lived neutral kaon production",
    "(z,Xk-) Total negative kaon production",
    "(z,Xp-) Total anti-proton production",
    "(z,Xn-) Total anti-neutron production",
    "Average cosine of scattering angle",
    "Average logarithmic energy decrement",
    "Average xi^2/(2*xi)",
    "Energy Release from (z,total)",
    "Energy Release from (z,elastic)",
    "Energy Release from (z,nonelastic)",
    "Energy Release from (z,inelastic)",
    "Energy Release from (misc)",
    "Energy Release from (z,continuum)",
    "Energy Release from (z,2nd)",
    "Energy Release from (z,2n)",
    "Energy Release from (z,3n)",
    "Energy Release from (z,fission)",
    "Energy Release from (z,f)",
    "Energy Release from (z,nf)",
    "Energy Release from (z,2nf)",
    "Energy Release from (z,n+a)",
    "Energy Release from (z,n+3a)",
    "Energy Release from (z,2n+a)",
    "Energy Release from (z,3n+a)",
    "Energy Release from (z,abs) Absorption",
    "Energy Release from (z,n+p)",
    "Energy Release from (z,n+2a)",
    "Energy Release from (z,2n+2a)",
    "Energy Release from (z,nd)",
    "Energy Release from (z,nt)",
    "Energy Release from (z,n+He3)",
    "Energy Release from (z,n+d+3a)",
    "Energy Release from (z,n+t+2a)",
    "Energy Release from (z,4n)",
    "Energy Release from (z,3nf)",
    "Energy Release from (z,2n+p)",
    "Energy Release from (z,3n+p)",
    "Energy Release from (z,n+2p)",
    "Energy Release from (z,npa)",
    "Energy Release from (z,n0)",
    "Energy Release from (z,n1)",
    "Energy Release from (z,n2)",
    "Energy Release from (z,n3)",
    "Energy Release from (z,n4)",
    "Energy Release from (z,n5)",
    "Energy Release from (z,n6)",
    "Energy Release from (z,n7)",
    "Energy Release from (z,n8)",
    "Energy Release from (z,n9)",
    "Energy Release from (z,n10)",
    "Energy Release from (z,n11)",
    "Energy Release from (z,n12)",
    "Energy Release from (z,n13)",
    "Energy Release from (z,n14)",
    "Energy Release from (z,n15)",
    "Energy Release from (z,n16)",
    "Energy Release from (z,n17)",
    "Energy Release from (z,n18)",
    "Energy Release from (z,n19)",
    "Energy Release from (z,n20)",
    "Energy Release from (z,n21)",
    "Energy Release from (z,n22)",
    "Energy Release from (z,n23)",
    "Energy Release from (z,n24)",
    "Energy Release from (z,n25)",
    "Energy Release from (z,n26)",
    "Energy Release from (z,n27)",
    "Energy Release from (z,n28)",
    "Energy Release from (z,n29)",
    "Energy Release from (z,n30)",
    "Energy Release from (z,n31)",
    "Energy Release from (z,n32)",
    "Energy Release from (z,n33)",
    "Energy Release from (z,n34)",
    "Energy Release from (z,n35)",
    "Energy Release from (z,n36)",
    "Energy Release from (z,n37)",
    "Energy Release from (z,n38)",
    "Energy Release from (z,n39)",
    "Energy Release from (z,n40)",
    "Energy Release from (z,nc)",
    "Energy Release from (z,disap) Neutron disappearance",
    "Energy Release from (z,gamma)",
    "Energy Release from (z,p)",
    "Energy Release from (z,d)",
    "Energy Release from (z,t)",
    "Energy Release from (z,3He)",
    "Energy Release from (z,a)",
    "Energy Release from (z,2a)",
    "Energy Release from (z,3a)",
    "Energy Release from (z,2p)",
    "Energy Release from (z,pa)",
    "Energy Release from (z,t2a)",
    "Energy Release from (z,d2a)",
    "Energy Release from (z,pd)",
    "Energy Release from (z,pt)",
    "Energy Release from (z,da)",
    "(damage)",
    "Descriptive Data",
    "Total Neutrons per Fission",
    "Independent fission product yield",
    "Delayed Neutron Data",
    "Prompt Neutrons per Fission",
    "Radioactive Decay Data",
    "Energy Release Due to Fission",
    "Cumulative Fission Product Yield",
    "Delayed Photon Data",
    "Total charged-particle stopping power",
    "Total photon interaction",
    "Photon coherent scattering",
    "Photon incoherent scattering",
    "Imaginary scattering factor",
    "Real scattering factor",
    "Pair production, electron field",
    "Total pair production",
    "Pair production, nuclear field",
    "Photoelectric absorption",
    "Photo-excitation cross section",
    "Electro-atomic scattering",
    "Electro-atomic bremsstrahlung",
    "Electro-atomic excitation cross section",
    "Atomic relaxation data",
    "K (1s1/2) subshell",
    "L1 (2s1/2) subshell",
    "L2 (2p1/2) subshell",
    "L3 (2p3/2) subshell",
    "M1 (3s1/2) subshell",
    "M2 (3p1/2) subshell",
    "M3 (3p3/2) subshell",
    "M4 (3d1/2) subshell",
    "M5 (3d1/2) subshell",
    "N1 (4s1/2) subshell",
    "N2 (4p1/2) subshell",
    "N3 (4p3/2) subshell",
    "N4 (4d3/2) subshell",
    "N5 (4d5/2) subshell",
    "N6 (4f5/2) subshell",
    "N7 (4f7/2) subshell",
    "O1 (5s1/2) subshell",
    "O2 (5p1/2) subshell",
    "O3 (5p3/2) subshell",
    "O4 (5d3/2) subshell",
    "O5 (5d5/2) subshell",
    "O6 (5f5/2) subshell",
    "O7 (5f7/2) subshell",
    "O8 (5g7/2) subshell",
    "O9 (5g9/2) subshell",
    "P1 (6s1/2) subshell",
    "P2 (6p1/2) subshell",
    "P3 (6p3/2) subshell",
    "P4 (6d3/2) subshell",
    "P5 (6d5/2) subshell",
    "P6 (6f5/2) subshell",
    "P7 (6f7/2) subshell",
    "P8 (6g7/2) subshell",
    "P9 (6g9/2) subshell",
    "P10 (6h9/2) subshell",
    "P11 (6h11/2) subshell",
    "Q1 (7s1/2) subshell",
    "Q2 (7p1/2) subshell",
    "Q3 (7p3/2) subshell",
    "(z,p0)",
    "(z,p1)",
    "(z,p2)",
    "(z,p3)",
    "(z,p4)",
    "(z,p5)",
    "(z,p6)",
    "(z,p7)",
    "(z,p8)",
    "(z,p9)",
    "(z,p10)",
    "(z,p11)",
    "(z,p12)",
    "(z,p13)",
    "(z,p14)",
    "(z,p15)",
    "(z,p16)",
    "(z,p17)",
    "(z,p18)",
    "(z,p19)",
    "(z,p20)",
    "(z,p21)",
    "(z,p22)",
    "(z,p23)",
    "(z,p24)",
    "(z,p25)",
    "(z,p26)",
    "(z,p27)",
    "(z,p28)",
    "(z,p29)",
    "(z,p30)",
    "(z,p31)",
    "(z,p32)",
    "(z,p33)",
    "(z,p34)",
    "(z,p35)",
    "(z,p36)",
    "(z,p37)",
    "(z,p38)",
    "(z,p39)",
    "(z,p40)",
    "(z,p41)",
    "(z,p42)",
    "(z,p43)",
    "(z,p44)",
    "(z,p45)",
    "(z,p46)",
    "(z,p47)",
    "(z,p48)",
    "(z,pc)",
    "(z,d0)",
    "(z,d1)",
    "(z,d2)",
    "(z,d3)",
    "(z,d4)",
    "(z,d5)",
    "(z,d6)",
    "(z,d7)",
    "(z,d8)",
    "(z,d9)",
    "(z,d10)",
    "(z,d11)",
    "(z,d12)",
    "(z,d13)",
    "(z,d14)",
    "(z,d15)",
    "(z,d16)",
    "(z,d17)",
    "(z,d18)",
    "(z,d19)",
    "(z,d20)",
    "(z,d21)",
    "(z,d22)",
    "(z,d23)",
    "(z,d24)",
    "(z,d25)",
    "(z,d26)",
    "(z,d27)",
    "(z,d28)",
    "(z,d29)",
    "(z,d30)",
    "(z,d31)",
    "(z,d32)",
    "(z,d33)",
    "(z,d34)",
    "(z,d35)",
    "(z,d36)",
    "(z,d37)",
    "(z,d38)",
    "(z,d39)",
    "(z,d40)",
    "(z,d41)",
    "(z,d42)",
    "(z,d43)",
    "(z,d44)",
    "(z,d45)",
    "(z,d46)",
    "(z,d47)",
    "(z,d48)",
    "(z,dc)",
    "(z,t0)",
    "(z,t1)",
    "(z,t2)",
    "(z,t3)",
    "(z,t4)",
    "(z,t5)",
    "(z,t6)",
    "(z,t7)",
    "(z,t8)",
    "(z,t9)",
    "(z,t10)",
    "(z,t11)",
    "(z,t12)",
    "(z,t13)",
    "(z,t14)",
    "(z,t15)",
    "(z,t16)",
    "(z,t17)",
    "(z,t18)",
    "(z,t19)",
    "(z,t20)",
    "(z,t21)",
    "(z,t22)",
    "(z,t23)",
    "(z,t24)",
    "(z,t25)",
    "(z,t26)",
    "(z,t27)",
    "(z,t28)",
    "(z,t29)",
    "(z,t30)",
    "(z,t31)",
    "(z,t32)",
    "(z,t33)",
    "(z,t34)",
    "(z,t35)",
    "(z,t36)",
    "(z,t37)",
    "(z,t38)",
    "(z,t39)",
    "(z,t40)",
    "(z,t41)",
    "(z,t42)",
    "(z,t43)",
    "(z,t44)",
    "(z,t45)",
    "(z,t46)",
    "(z,t47)",
    "(z,t48)",
    "(z,tc)",
    "(z,3He0)",
    "(z,3He1)",
    "(z,3He2)",
    "(z,3He3)",
    "(z,3He4)",
    "(z,3He5)",
    "(z,3He6)",
    "(z,3He7)",
    "(z,3He8)",
    "(z,3He9)",
    "(z,3He10)",
    "(z,3He11)",
    "(z,3He12)",
    "(z,3He13)",
    "(z,3He14)",
    "(z,3He15)",
    "(z,3He16)",
    "(z,3He17)",
    "(z,3He18)",
    "(z,3He19)",
    "(z,3He20)",
    "(z,3He21)",
    "(z,3He22)",
    "(z,3He23)",
    "(z,3He24)",
    "(z,3He25)",
    "(z,3He26)",
    "(z,3He27)",
    "(z,3He28)",
    "(z,3He29)",
    "(z,3He30)",
    "(z,3He31)",
    "(z,3He32)",
    "(z,3He33)",
    "(z,3He34)",
    "(z,3He35)",
    "(z,3He36)",
    "(z,3He37)",
    "(z,3He38)",
    "(z,3He39)",
    "(z,3He40)",
    "(z,3He41)",
    "(z,3He42)",
    "(z,3He43)",
    "(z,3He44)",
    "(z,3He45)",
    "(z,3He46)",
    "(z,3He47)",
    "(z,3He48)",
    "(z,3Hec)",
    "(z,a0)",
    "(z,a1)",
    "(z,a2)",
    "(z,a3)",
    "(z,a4)",
    "(z,a5)",
    "(z,a6)",
    "(z,a7)",
    "(z,a8)",
    "(z,a9)",
    "(z,a10)",
    "(z,a11)",
    "(z,a12)",
    "(z,a13)",
    "(z,a14)",
    "(z,a15)",
    "(z,a16)",
    "(z,a17)",
    "(z,a18)",
    "(z,a19)",
    "(z,a20)",
    "(z,a21)",
    "(z,a22)",
    "(z,a23)",
    "(z,a24)",
    "(z,a25)",
    "(z,a26)",
    "(z,a27)",
    "(z,a28)",
    "(z,a29)",
    "(z,a30)",
    "(z,a31)",
    "(z,a32)",
    "(z,a33)",
    "(z,a34)",
    "(z,a35)",
    "(z,a36)",
    "(z,a37)",
    "(z,a38)",
    "(z,a39)",
    "(z,a40)",
    "(z,a41)",
    "(z,a42)",
    "(z,a43)",
    "(z,a44)",
    "(z,a45)",
    "(z,a46)",
    "(z,a47)",
    "(z,a48)",
    "(z,ac)",
    "Lumped Covariances",
    "Any Excited State",
    "(z,b-)",
    "(z,b+)",
    "(z,ec)",
    "(z,b-n)",
    "(z,b-a)",
    "(z,it)",
    "(z,b+a)",
    "(z,ec+b+)",
    "(z,b+p)",
    "(z,b-2n)",
    "(z,b-3n)",
    "(z,b-4n)",
    "(z,ecp)",
    "(z,eca)",
    "(z,b+2p)",
    "(z,ec2p)",
    "(z,2b-)",
    "(z,b-p)",
    "(z,14c)",
    "(z,b+3p)",
    "(z,sf)",
    "(z,2b+)",
    "(z,2ec)",
    "(z,ec3p)",
    "(z,b-sf)"
  };
  std::string _docs[NUM_RX_NAMES] = {
    "(n,total) Neutron total",
    "Total scattering",
    "(z,z0) Elastic scattering",
    "(z,nonelas) Nonelastic neutron",
    "(z,n) One neutron in exit channel",
    "(z,anything) Miscellaneous",
    "(z,contin) Total continuum reaction",
    "(z,2nd) Production of 2n and d",
    "(z,2n) Production of 2n",
    "(z,2n0) Production of 2n, ground state",
    "(z,2n1) Production of 2n, 1st excited state",
    "(z,2n2) Production of 2n, 2nd excited state",
    "(z,3n) Production of 3n",
    "(z,3n0) Production of 3n, ground state",
    "(z,3n1) Production of 3n, 1st excited state",
    "(z,3n2) Production of 3n, 2nd excited state",
    "(z,fiss) Particle-induced fission",
    "(z,f) First-chance fission",
    "(z,nf) Second chance fission",
    "(z,2nf) Third-chance fission",
    "(z,na) Production of n and alpha",
    "(z,na0) Production of n and alpha, ground state",
    "(z,na1) Production of n and alpha, 1st excited state",
    "(z,na2) Production of n and alpha, 2nd excited state",
    "(z,n3a) Production of n and 3 alphas",
    "(z,2na) Production of 2n and alpha",
    "(z,3na) Production of 3n and alpha",
    "(n,abs) Absorption",
    "(z,np) Production of n and p",
    "(z,np0) Production of n and p, ground state",
    "(z,np1) Production of n and p, 1st excited state",
    "(z,np2) Production of n and p, 2nd excited state",
    "(z,npd) Production of n, p, and d",
    "(z,n2a) Production of n and 2 alphas",
    "(z,2n2a) Production of 2n and 2 alphas",
    "(z,nd) Production of n and d",
    "(z,nd0) Production of n and d, ground state",
    "(z,nd1) Production of n and d, 1st excited state",
    "(z,nd2) Production of n and d, 2nd excited state",
    "(z,nt) Production of n and t",
    "(z,nt0) Production of n and t, ground state",
    "(z,nt1) Production of n and t, 1st excited state",
    "(z,nt2) Production of n and t, 2nd excited state",
    "(z,n3He) Production of n and He-3",
    "(z,n3He-0) Production of n and He-3, ground state",
    "(z,n3He-1) Production of n and He-3, 1st excited state",
    "(z,n3He-2) Production of n and He-3, 2nd excited state",
    "(z,nd2a) Production of n, d, and alpha",
    "(z,nt2a) Production of n, t, and 2 alphas",
    "(z,4n) Production of 4n",
    "(z,4n0) Production of 4n, ground state",
    "(z,4n1) Production of 4n, 1st excited state",
    "(z,3nf) Fourth-chance fission",
    "(z,2np) Production of 2n and p",
    "(z,3np) Production of 3n and p",
    "(z,n2p) Production of n and 2p",
    "(z,npa) Production of n, p, and alpha",
    "(z,n0) Production of n, ground state",
    "(z,n1) Production of n, 1st excited state",
    "(z,n2) Production of n, 2nd excited state",
    "(z,n3) Production of n, 3rd excited state",
    "(z,n4) Production of n, 4th excited state",
    "(z,n5) Production of n, 5th excited state",
    "(z,n6) Production of n, 6th excited state",
    "(z,n7) Production of n, 7th excited state",
    "(z,n8) Production of n, 8th excited state",
    "(z,n9) Production of n, 9th excited state",
    "(z,n10) Production of n, 10th excited state",
    "(z,n11) Production of n, 11th excited state",
    "(z,n12) Production of n, 12th excited state",
    "(z,n13) Production of n, 13th excited state",
    "(z,n14) Production of n, 14th excited state",
    "(z,n15) Production of n, 15th excited state",
    "(z,n16) Production of n, 16th excited state",
    "(z,n17) Production of n, 17th excited state",
    "(z,n18) Production of n, 18th excited state",
    "(z,n19) Production of n, 19th excited state",
    "(z,n20) Production of n, 20th excited state",
    "(z,n21) Production of n, 21st excited state",
    "(z,n22) Production of n, 22nd excited state",
    "(z,n23) Production of n, 23rd excited state",
    "(z,n24) Production of n, 24th excited state",
    "(z,n25) Production of n, 25th excited state",
    "(z,n26) Production of n, 26th excited state",
    "(z,n27) Production of n, 27th excited state",
    "(z,n28) Production of n, 28th excited state",
    "(z,n29) Production of n, 29th excited state",
    "(z,n30) Production of n, 30th excited state",
    "(z,n31) Production of n, 31st excited state",
    "(z,n32) Production of n, 32nd excited state",
    "(z,n33) Production of n, 33rd excited state",
    "(z,n34) Production of n, 34th excited state",
    "(z,n35) Production of n, 35th excited state",
    "(z,n36) Production of n, 36th excited state",
    "(z,n37) Production of n, 37th excited state",
    "(z,n38) Production of n, 38th excited state",
    "(z,n39) Production of n, 39th excited state",
    "(z,n40) Production of n, 40th excited state",
    "(z,nc) Production of n in continuum",
    "(n,disap) Neutron disappearance",
    "(z,gamma) Radiative capture",
    "(z,gamma0) Radiative capture, ground state",
    "(z,gamma1) Radiative capture, 1st excited state",
    "(z,gamma2) Radiative capture, 2st excited state",
    "(z,p) Production of p",
    "(z,d) Production of d",
    "(z,t) Production of t",
    "(z,3He) Production of He-3",
    "(z,a) Production of alpha",
    "(z,2a) Production of 2 alphas",
    "(z,3a) Production of 3 alphas",
    "(z,2p) Production of 2p",
    "(z,2p0) Production of 2p, ground state",
    "(z,2p1) Production of 2p, 1st excited state",
    "(z,2p2) Production of 2p, 2nd excited state",
    "(z,pa) Production of p and alpha",
    "(z,t2a) Production of t and 2 alphas",
    "(z,d2a) Production of d and 2 alphas",
    "(z,pd) Production of p and d",
    "(z,pt) Production of p and t",
    "(z,da) Production of d and a",
    "Resonance Parameters",
    "(z,Xn) Total neutron production",
    "(z,Xgamma) Total gamma production",
    "(z,Xp) Total proton production",
    "(z,Xd) Total deuteron production",
    "(z,Xt) Total triton production",
    "(z,X3He) Total He-3 production",
    "(z,Xa) Total alpha production",
    "(z,Xpi+) Total pi+ meson production",
    "(z,Xpi0) Total pi0 meson production",
    "(z,Xpi-) Total pi- meson production",
    "(z,Xmu+) Total anti-muon production",
    "(z,Xmu-) Total muon production",
    "(z,Xk+) Total positive kaon production",
    "(z,Xk0long) Total long-lived neutral kaon production",
    "(z,Xk0short) Total short-lived neutral kaon production",
    "(z,Xk-) Total negative kaon production",
    "(z,Xp-) Total anti-proton production",
    "(z,Xn-) Total anti-neutron production",
    "Average cosine of scattering angle",
    "Average logarithmic energy decrement",
    "Average xi^2/(2*xi)",
    "Energy Release from (n,total) Neutron total",
    "Energy Release from (z,z0) Elastic scattering",
    "Energy Release from (z,nonelas) Nonelastic neutron",
    "Energy Release from (z,n) One neutron in exit channel",
    "Energy Release from (z,anything) Miscellaneous",
    "Energy Release from (z,contin) Total continuum reaction",
    "Energy Release from (z,2nd) Production of 2n and d",
    "Energy Release from (z,2n) Production of 2n",
    "Energy Release from (z,3n) Production of 3n",
    "Energy Release from (z,fiss) Particle-induced fission",
    "Energy Release from (z,f) First-chance fission",
    "Energy Release from (z,nf) Second chance fission",
    "Energy Release from (z,2nf) Third-chance fission",
    "Energy Release from (z,na) Production of n and alpha",
    "Energy Release from (z,n3a) Production of n and 3 alphas",
    "Energy Release from (z,2na) Production of 2n and alpha",
    "Energy Release from (z,3na) Production of 3n and alpha",
    "Energy Release from (n,abs) Absorption",
    "Energy Release from (z,np) Production of n and p",
    "Energy Release from (z,n2a) Production of n and 2 alphas",
    "Energy Release from (z,2n2a) Production of 2n and 2 alphas",
    "Energy Release from (z,nd) Production of n and d",
    "Energy Release from (z,nt) Production of n and t",
    "Energy Release from (z,n3He) Production of n and He-3",
    "Energy Release from (z,nd2a) Production of n, d, and alpha",
    "Energy Release from (z,nt2a) Production of n, t, and 2 alphas",
    "Energy Release from (z,4n) Production of 4n",
    "Energy Release from (z,3nf) Fourth-chance fission",
    "Energy Release from (z,2np) Production of 2n and p",
    "Energy Release from (z,3np) Production of 3n and p",
    "Energy Release from (z,n2p) Production of n and 2p",
    "Energy Release from (z,npa) Production of n, p, and alpha",
    "Energy Release from (z,n0) Production of n, ground state",
    "Energy Release from (z,n1) Production of n, 1st excited state",
    "Energy Release from (z,n2) Production of n, 2nd excited state",
    "Energy Release from (z,n3) Production of n, 3rd excited state",
    "Energy Release from (z,n4) Production of n, 4th excited state",
    "Energy Release from (z,n5) Production of n, 5th excited state",
    "Energy Release from (z,n6) Production of n, 6th excited state",
    "Energy Release from (z,n7) Production of n, 7th excited state",
    "Energy Release from (z,n8) Production of n, 8th excited state",
    "Energy Release from (z,n9) Production of n, 9th excited state",
    "Energy Release from (z,n10) Production of n, 10th excited state",
    "Energy Release from (z,n11) Production of n, 11th excited state",
    "Energy Release from (z,n12) Production of n, 12th excited state",
    "Energy Release from (z,n13) Production of n, 13th excited state",
    "Energy Release from (z,n14) Production of n, 14th excited state",
    "Energy Release from (z,n15) Production of n, 15th excited state",
    "Energy Release from (z,n16) Production of n, 16th excited state",
    "Energy Release from (z,n17) Production of n, 17th excited state",
    "Energy Release from (z,n18) Production of n, 18th excited state",
    "Energy Release from (z,n19) Production of n, 19th excited state",
    "Energy Release from (z,n20) Production of n, 20th excited state",
    "Energy Release from (z,n21) Production of n, 21st excited state",
    "Energy Release from (z,n22) Production of n, 22nd excited state",
    "Energy Release from (z,n23) Production of n, 23rd excited state",
    "Energy Release from (z,n24) Production of n, 24th excited state",
    "Energy Release from (z,n25) Production of n, 25th excited state",
    "Energy Release from (z,n26) Production of n, 26th excited state",
    "Energy Release from (z,n27) Production of n, 27th excited state",
    "Energy Release from (z,n28) Production of n, 28th excited state",
    "Energy Release from (z,n29) Production of n, 29th excited state",
    "Energy Release from (z,n30) Production of n, 30th excited state",
    "Energy Release from (z,n31) Production of n, 31st excited state",
    "Energy Release from (z,n32) Production of n, 32nd excited state",
    "Energy Release from (z,n33) Production of n, 33rd excited state",
    "Energy Release from (z,n34) Production of n, 34th excited state",
    "Energy Release from (z,n35) Production of n, 35th excited state",
    "Energy Release from (z,n36) Production of n, 36th excited state",
    "Energy Release from (z,n37) Production of n, 37th excited state",
    "Energy Release from (z,n38) Production of n, 38th excited state",
    "Energy Release from (z,n39) Production of n, 39th excited state",
    "Energy Release from (z,n40) Production of n, 40th excited state",
    "Energy Release from (z,nc) Production of n in continuum",
    "Energy Release from (n,disap) Neutron disappearance",
    "Energy Release from (z,gamma) Radiative capture",
    "Energy Release from (z,p) Production of p",
    "Energy Release from (z,d) Production of d",
    "Energy Release from (z,t) Production of t",
    "Energy Release from (z,3He) Production of He-3",
    "Energy Release from (z,a) Production of alpha",
    "Energy Release from (z,2a) Production of 2 alphas",
    "Energy Release from (z,3a) Production of 3 alphas",
    "Energy Release from (z,2p) Production of 2p",
    "Energy Release from (z,pa) Production of p and alpha",
    "Energy Release from (z,t2a) Production of t and 2 alphas",
    "Energy Release from (z,d2a) Production of d and 2 alphas",
    "Energy Release from (z,pd) Production of p and d",
    "Energy Release from (z,pt) Production of p and t",
    "Energy Release from (z,da) Production of d and a",
    "(damage)",
    "Descriptive Data",
    "Total Neutrons per Fission",
    "Independent fission product yield",
    "Delayed Neutron Data",
    "Prompt Neutrons per Fission",
    "Radioactive Decay Data",
    "Energy Release Due to Fission",
    "Cumulative Fission Product Yield",
    "Delayed Photon Data",
    "Total charged-particle stopping power",
    "Total photon interaction",
    "Photon coherent scattering",
    "Photon incoherent scattering",
    "Imaginary scattering factor",
    "Real scattering factor",
    "Pair production, electron field",
    "Total pair production",
    "Pair production, nuclear field",
    "Photoelectric absorption",
    "Photo-excitation cross section",
    "Electro-atomic scattering",
    "Electro-atomic bremsstrahlung",
    "Electro-atomic excitation cross section",
    "Atomic relaxation data",
    "K (1s1/2) subshell",
    "L1 (2s1/2) subshell",
    "L2 (2p1/2) subshell",
    "L3 (2p3/2) subshell",
    "M1 (3s1/2) subshell",
    "M2 (3p1/2) subshell",
    "M3 (3p3/2) subshell",
    "M4 (3d1/2) subshell",
    "M5 (3d1/2) subshell",
    "N1 (4s1/2) subshell",
    "N2 (4p1/2) subshell",
    "N3 (4p3/2) subshell",
    "N4 (4d3/2) subshell",
    "N5 (4d5/2) subshell",
    "N6 (4f5/2) subshell",
    "N7 (4f7/2) subshell",
    "O1 (5s1/2) subshell",
    "O2 (5p1/2) subshell",
    "O3 (5p3/2) subshell",
    "O4 (5d3/2) subshell",
    "O5 (5d5/2) subshell",
    "O6 (5f5/2) subshell",
    "O7 (5f7/2) subshell",
    "O8 (5g7/2) subshell",
    "O9 (5g9/2) subshell",
    "P1 (6s1/2) subshell",
    "P2 (6p1/2) subshell",
    "P3 (6p3/2) subshell",
    "P4 (6d3/2) subshell",
    "P5 (6d5/2) subshell",
    "P6 (6f5/2) subshell",
    "P7 (6f7/2) subshell",
    "P8 (6g7/2) subshell",
    "P9 (6g9/2) subshell",
    "P10 (6h9/2) subshell",
    "P11 (6h11/2) subshell",
    "Q1 (7s1/2) subshell",
    "Q2 (7p1/2) subshell",
    "Q3 (7p3/2) subshell",
    "(n,p0)",
    "(n,p1)",
    "(n,p2)",
    "(n,p3)",
    "(n,p4)",
    "(n,p5)",
    "(n,p6)",
    "(n,p7)",
    "(n,p8)",
    "(n,p9)",
    "(n,p10)",
    "(n,p11)",
    "(n,p12)",
    "(n,p13)",
    "(n,p14)",
    "(n,p15)",
    "(n,p16)",
    "(n,p17)",
    "(n,p18)",
    "(n,p19)",
    "(n,p20)",
    "(n,p21)",
    "(n,p22)",
    "(n,p23)",
    "(n,p24)",
    "(n,p25)",
    "(n,p26)",
    "(n,p27)",
    "(n,p28)",
    "(n,p29)",
    "(n,p30)",
    "(n,p31)",
    "(n,p32)",
    "(n,p33)",
    "(n,p34)",
    "(n,p35)",
    "(n,p36)",
    "(n,p37)",
    "(n,p38)",
    "(n,p39)",
    "(n,p40)",
    "(n,p41)",
    "(n,p42)",
    "(n,p43)",
    "(n,p44)",
    "(n,p45)",
    "(n,p46)",
    "(n,p47)",
    "(n,p48)",
    "(n,pc)",
    "(n,d0)",
    "(n,d1)",
    "(n,d2)",
    "(n,d3)",
    "(n,d4)",
    "(n,d5)",
    "(n,d6)",
    "(n,d7)",
    "(n,d8)",
    "(n,d9)",
    "(n,d10)",
    "(n,d11)",
    "(n,d12)",
    "(n,d13)",
    "(n,d14)",
    "(n,d15)",
    "(n,d16)",
    "(n,d17)",
    "(n,d18)",
    "(n,d19)",
    "(n,d20)",
    "(n,d21)",
    "(n,d22)",
    "(n,d23)",
    "(n,d24)",
    "(n,d25)",
    "(n,d26)",
    "(n,d27)",
    "(n,d28)",
    "(n,d29)",
    "(n,d30)",
    "(n,d31)",
    "(n,d32)",
    "(n,d33)",
    "(n,d34)",
    "(n,d35)",
    "(n,d36)",
    "(n,d37)",
    "(n,d38)",
    "(n,d39)",
    "(n,d40)",
    "(n,d41)",
    "(n,d42)",
    "(n,d43)",
    "(n,d44)",
    "(n,d45)",
    "(n,d46)",
    "(n,d47)",
    "(n,d48)",
    "(n,dc)",
    "(z,t0)",
    "(z,t1)",
    "(z,t2)",
    "(z,t3)",
    "(z,t4)",
    "(z,t5)",
    "(z,t6)",
    "(z,t7)",
    "(z,t8)",
    "(z,t9)",
    "(z,t10)",
    "(z,t11)",
    "(z,t12)",
    "(z,t13)",
    "(z,t14)",
    "(z,t15)",
    "(z,t16)",
    "(z,t17)",
    "(z,t18)",
    "(z,t19)",
    "(z,t20)",
    "(z,t21)",
    "(z,t22)",
    "(z,t23)",
    "(z,t24)",
    "(z,t25)",
    "(z,t26)",
    "(z,t27)",
    "(z,t28)",
    "(z,t29)",
    "(z,t30)",
    "(z,t31)",
    "(z,t32)",
    "(z,t33)",
    "(z,t34)",
    "(z,t35)",
    "(z,t36)",
    "(z,t37)",
    "(z,t38)",
    "(z,t39)",
    "(z,t40)",
    "(z,t41)",
    "(z,t42)",
    "(z,t43)",
    "(z,t44)",
    "(z,t45)",
    "(z,t46)",
    "(z,t47)",
    "(z,t48)",
    "(n,tc)",
    "(n,3He0)",
    "(n,3He1)",
    "(n,3He2)",
    "(n,3He3)",
    "(n,3He4)",
    "(n,3He5)",
    "(n,3He6)",
    "(n,3He7)",
    "(n,3He8)",
    "(n,3He9)",
    "(n,3He10)",
    "(n,3He11)",
    "(n,3He12)",
    "(n,3He13)",
    "(n,3He14)",
    "(n,3He15)",
    "(n,3He16)",
    "(n,3He17)",
    "(n,3He18)",
    "(n,3He19)",
    "(n,3He20)",
    "(n,3He21)",
    "(n,3He22)",
    "(n,3He23)",
    "(n,3He24)",
    "(n,3He25)",
    "(n,3He26)",
    "(n,3He27)",
    "(n,3He28)",
    "(n,3He29)",
    "(n,3He30)",
    "(n,3He31)",
    "(n,3He32)",
    "(n,3He33)",
    "(n,3He34)",
    "(n,3He35)",
    "(n,3He36)",
    "(n,3He37)",
    "(n,3He38)",
    "(n,3He39)",
    "(n,3He40)",
    "(n,3He41)",
    "(n,3He42)",
    "(n,3He43)",
    "(n,3He44)",
    "(n,3He45)",
    "(n,3He46)",
    "(n,3He47)",
    "(n,3He48)",
    "(n,3Hec)",
    "(z,a0)",
    "(z,a1)",
    "(z,a2)",
    "(z,a3)",
    "(z,a4)",
    "(z,a5)",
    "(z,a6)",
    "(z,a7)",
    "(z,a8)",
    "(z,a9)",
    "(z,a10)",
    "(z,a11)",
    "(z,a12)",
    "(z,a13)",
    "(z,a14)",
    "(z,a15)",
    "(z,a16)",
    "(z,a17)",
    "(z,a18)",
    "(z,a19)",
    "(z,a20)",
    "(z,a21)",
    "(z,a22)",
    "(z,a23)",
    "(z,a24)",
    "(z,a25)",
    "(z,a26)",
    "(z,a27)",
    "(z,a28)",
    "(z,a29)",
    "(z,a30)",
    "(z,a31)",
    "(z,a32)",
    "(z,a33)",
    "(z,a34)",
    "(z,a35)",
    "(z,a36)",
    "(z,a37)",
    "(z,a38)",
    "(z,a39)",
    "(z,a40)",
    "(z,a41)",
    "(z,a42)",
    "(z,a43)",
    "(z,a44)",
    "(z,a45)",
    "(z,a46)",
    "(z,a47)",
    "(z,a48)",
    "(n,ac)",
    "Lumped-Reaction Covariances",
    "production of any excited state nucleus",
    "(z,b-)",
    "(z,b+)",
    "(z,ec)",
    "(z,b-n)",
    "(z,b-a)",
    "(z,it)",
    "(z,b+a)",
    "(z,ec+b+)",
    "(z,b+p)",
    "(z,b-2n)",
    "(z,b-3n)",
    "(z,b-4n)",
    "(z,ecp)",
    "(z,eca)",
    "(z,b+2p)",
    "(z,ec2p)",
    "(z,2b-)",
    "(z,b-p)",
    "(z,14c)",
    "(z,b+3p)",
    "(z,sf)",
    "(z,2b+)",
    "(z,2ec)",
    "(z,ec3p)",
    "(z,b-sf)"
  };

  // fill the maps
  for (int i = 0; i < NUM_RX_NAMES; i++) {
    rx = _names[i];
    rxid = pyne::rxname::hash(rx);
    id_name[rxid] = rx;
    name_id[rx] = rxid;
    if (0 < _mts[i]) {
      id_mt[rxid] = _mts[i];
      mt_id[_mts[i]] = rxid;
    }
    labels[rxid] = _labels[i];
    docs[rxid] = _docs[i];
  }

  // set alternative names
  altnames["tot"] = name_id["total"];
  altnames["s"] = name_id["scattering"];
  altnames["scat"] = name_id["scattering"];
  altnames["e"] = name_id["elastic"];
  altnames["elas"] = name_id["elastic"];
  altnames["i"] = name_id["n"];
  altnames["inel"] = name_id["n"];
  altnames["inelastic"] = name_id["n"];
  altnames["abs"] = name_id["absorption"];
  altnames["fis"] = name_id["fission"];
  altnames["fiss"] = name_id["fission"];
  altnames["alpha"] = name_id["a"];
  altnames["deut"] = name_id["d"];
  altnames["deuteron"] = name_id["d"];
  altnames["deuterium"] = name_id["d"];
  altnames["trit"] = name_id["t"];
  altnames["triton"] = name_id["t"];
  altnames["tritium"] = name_id["t"];
  altnames["proton"] = name_id["p"];
  altnames["h"] = name_id["He3"];  // 'h' stands for helion
  altnames["he3"] = name_id["He3"];
  altnames["HE3"] = name_id["He3"];
  altnames["3HE"] = name_id["He3"];
  altnames["3He"] = name_id["He3"];
  altnames["3he"] = name_id["He3"];
  altnames["he-3"] = name_id["He3"];
  altnames["HE-3"] = name_id["He3"];
  altnames["*"] = name_id["excited"];
  altnames["2n"] = name_id["z_2n"];
  altnames["2p"] = name_id["z_2p"];
  altnames["3h"] = name_id["t"];
  altnames["g"] = name_id["it"];
  altnames["b-"] = name_id["bminus"];
  altnames["b+"] = name_id["bplus"];
  altnames["b-n"] = name_id["bminus_n"];
  altnames["b-a"] = name_id["bminus_a"];
  altnames["b+a"] = name_id["bplus_a"];
  altnames["ec+b+"] = name_id["ec_bplus"];
  altnames["b+p"] = name_id["bplus_p"];
  altnames["b-2n"] = name_id["bminus_2n"];
  altnames["b-3n"] = name_id["bminus_3n"];
  altnames["b-4n"] = name_id["bminus_4n"];
  altnames["b+2p"] = name_id["bplus_2p"];
  altnames["ec2p"] = name_id["ec_2p"];
  altnames["ec3p"] = name_id["ec_3p"];
  altnames["2b-"] = name_id["decay_2bminus"];
  altnames["b-p"] = name_id["bminus_p"];
  altnames["14c"] = name_id["decay_14c"];
  altnames["b+3p"] = name_id["bplus_3p"];
  altnames["2b+"] = name_id["decay_2bplus"];
  altnames["2ec"] = name_id["decay_2ec"];
  altnames["b-sf"] = name_id["bminus_sf"];


  // set the nuclide difference mappings, offset_id
  // offset_id[incident particle type "n", "p", ...][delta Z num][delta A num][rxid]
  // offset_id mapping may be ambiquious so they must come before the id_offsets!
  // the following should be sorted by (dz, da, ds)
  // neutrons:
  offset_id[make_pair("n", offset(-4, -8))] = name_id["n2a"];
  offset_id[make_pair("n", offset(-4, -7))] = name_id["z_2a"];
  offset_id[make_pair("n", offset(-2, -5))] = name_id["z_2na"];
  offset_id[make_pair("n", offset(-2, -4))] = name_id["na"];
  offset_id[make_pair("n", offset(-2, -4, 1))] = name_id["na_1"];
  offset_id[make_pair("n", offset(-2, -4, 2))] = name_id["na_2"];
  offset_id[make_pair("n", offset(-2, -3))] = name_id["a"];
  offset_id[make_pair("n", offset(-2, -3, 1))] = name_id["a_1"];
  offset_id[make_pair("n", offset(-2, -3, 2))] = name_id["a_2"];
  offset_id[make_pair("n", offset(-2, -2))] = name_id["He3"];
  offset_id[make_pair("n", offset(-2, -2, 1))] = name_id["He3_1"];
  offset_id[make_pair("n", offset(-2, -2, 2))] = name_id["He3_2"];
  offset_id[make_pair("n", offset(-2, -1))] = name_id["z_2p"];
  offset_id[make_pair("n", offset(-2, -1, 1))] = name_id["z_2p_1"];
  offset_id[make_pair("n", offset(-2, -1, 2))] = name_id["z_2p_2"];
  offset_id[make_pair("n", offset(-1, -3))] = name_id["nt"];
  offset_id[make_pair("n", offset(-1, -3, 1))] = name_id["nt_1"];
  offset_id[make_pair("n", offset(-1, -3, 2))] = name_id["nt_2"];
  offset_id[make_pair("n", offset(-1, -2))] = name_id["t"];
  offset_id[make_pair("n", offset(-1, -2, 1))] = name_id["t_1"];
  offset_id[make_pair("n", offset(-1, -2, 2))] = name_id["t_2"];
  offset_id[make_pair("n", offset(-1, -1))] = name_id["d"];
  offset_id[make_pair("n", offset(-1, -1, 1))] = name_id["d_1"];
  offset_id[make_pair("n", offset(-1, -1, 2))] = name_id["d_2"];
  offset_id[make_pair("n", offset(-1, 0))] = name_id["p"];
  offset_id[make_pair("n", offset(-1, 0, 1))] = name_id["p_1"];
  offset_id[make_pair("n", offset(-1, 0, 2))] = name_id["p_2"];
  offset_id[make_pair("n", offset(0, -3))] = name_id["z_4n"];
  offset_id[make_pair("n", offset(0, -3, 1))] = name_id["z_4n_1"];
  offset_id[make_pair("n", offset(0, -2))] = name_id["z_3n"];
  offset_id[make_pair("n", offset(0, -2, 1))] = name_id["z_3n_1"];
  offset_id[make_pair("n", offset(0, -2, 2))] = name_id["z_3n_2"];
  offset_id[make_pair("n", offset(0, -1))] = name_id["z_2n"];
  offset_id[make_pair("n", offset(0, -1, 1))] = name_id["z_2n_1"];
  offset_id[make_pair("n", offset(0, -1, 2))] = name_id["z_2n_2"];
  offset_id[make_pair("n", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("n", offset(0, 0, 1))] = name_id["n_1"];
  offset_id[make_pair("n", offset(0, 0, 2))] = name_id["n_2"];
  offset_id[make_pair("n", offset(0, 1))] = name_id["absorption"];
  offset_id[make_pair("n", offset(0, 1, 1))] = name_id["gamma_1"];
  offset_id[make_pair("n", offset(0, 1, 2))] = name_id["gamma_2"];
  // proton:
  offset_id[make_pair("p", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("p", offset(1, 1))] = name_id["absorption"];
  offset_id[make_pair("p", offset(1, 0))] = name_id["n"];
  offset_id[make_pair("p", offset(1, -1))] = name_id["z_2n"];
  offset_id[make_pair("p", offset(1, -2))] = name_id["z_3n"];
  offset_id[make_pair("p", offset(1, -3))] = name_id["z_4n"];
  offset_id[make_pair("p", offset(-1, -1))] = name_id["z_2p"];
  offset_id[make_pair("p", offset(0, -1))] = name_id["d"];
  offset_id[make_pair("p", offset(0, -2))] = name_id["t"];
  offset_id[make_pair("p", offset(-1, -2))] = name_id["He3"];
  offset_id[make_pair("p", offset(-1, -3))] = name_id["a"];
  // deuterium:
  offset_id[make_pair("d", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("d", offset(1, 2))] = name_id["absorption"];
  offset_id[make_pair("d", offset(1, 1))] = name_id["n"];
  offset_id[make_pair("d", offset(1, 0))] = name_id["z_2n"];
  offset_id[make_pair("d", offset(1, -1))] = name_id["z_3n"];
  offset_id[make_pair("d", offset(1, -2))] = name_id["z_4n"];
  offset_id[make_pair("d", offset(0, 1))] = name_id["p"];
  offset_id[make_pair("d", offset(-1, 0))] = name_id["z_2p"];
  offset_id[make_pair("d", offset(0, -1))] = name_id["t"];
  offset_id[make_pair("d", offset(-1, -1))] = name_id["He3"];
  offset_id[make_pair("d", offset(-1, -2))] = name_id["a"];
  // tritium:
  offset_id[make_pair("t", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("t", offset(1, 3))] = name_id["absorption"];
  offset_id[make_pair("t", offset(1, 2))] = name_id["n"];
  offset_id[make_pair("t", offset(1, 1))] = name_id["z_2n"];
  offset_id[make_pair("t", offset(1, 0))] = name_id["z_3n"];
  offset_id[make_pair("t", offset(1, -1))] = name_id["z_4n"];
  offset_id[make_pair("t", offset(0, 2))] = name_id["p"];
  offset_id[make_pair("t", offset(-1, 1))] = name_id["z_2p"];
  offset_id[make_pair("t", offset(0, 1))] = name_id["d"];
  offset_id[make_pair("t", offset(-1, 0))] = name_id["He3"];
  offset_id[make_pair("t", offset(-1, -1))] = name_id["a"];
  // He3:
  offset_id[make_pair("He3", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("He3", offset(2, 3))] = name_id["absorption"];
  offset_id[make_pair("He3", offset(2, 2))] = name_id["n"];
  offset_id[make_pair("He3", offset(2, 1))] = name_id["z_2n"];
  offset_id[make_pair("He3", offset(2, 0))] = name_id["z_3n"];
  offset_id[make_pair("He3", offset(2, -1))] = name_id["z_4n"];
  offset_id[make_pair("He3", offset(1, 2))] = name_id["p"];
  offset_id[make_pair("He3", offset(0, 1))] = name_id["z_2p"];
  offset_id[make_pair("He3", offset(1, 1))] = name_id["d"];
  offset_id[make_pair("He3", offset(1, 0))] = name_id["t"];
  offset_id[make_pair("He3", offset(0, -1))] = name_id["a"];
  // alpha:
  offset_id[make_pair("a", offset(0, 0))] = name_id["scattering"];
  offset_id[make_pair("a", offset(2, 4))] = name_id["absorption"];
  offset_id[make_pair("a", offset(2, 3))] = name_id["n"];
  offset_id[make_pair("a", offset(2, 2))] = name_id["z_2n"];
  offset_id[make_pair("a", offset(2, 1))] = name_id["z_3n"];
  offset_id[make_pair("a", offset(2, 0))] = name_id["z_4n"];
  offset_id[make_pair("a", offset(1, 3))] = name_id["p"];
  offset_id[make_pair("a", offset(0, 2))] = name_id["z_2p"];
  offset_id[make_pair("a", offset(1, 2))] = name_id["d"];
  offset_id[make_pair("a", offset(1, 1))] = name_id["t"];
  offset_id[make_pair("a", offset(0, 1))] = name_id["He3"];
  // gamma:
  offset_id[make_pair("gamma", offset(0, -1))] = name_id["n"];
  offset_id[make_pair("gamma", offset(0, -2))] = name_id["z_2n"];
  offset_id[make_pair("gamma", offset(0, -3))] = name_id["z_3n"];
  offset_id[make_pair("gamma", offset(0, -4))] = name_id["z_4n"];
  offset_id[make_pair("gamma", offset(-1, -1))] = name_id["p"];
  offset_id[make_pair("gamma", offset(-2, -2))] = name_id["z_2p"];
  offset_id[make_pair("gamma", offset(-1, -2))] = name_id["d"];
  offset_id[make_pair("gamma", offset(-1, -3))] = name_id["t"];
  offset_id[make_pair("gamma", offset(-2, -3))] = name_id["He3"];
  offset_id[make_pair("gamma", offset(-2, -4))] = name_id["a"];
  // decay:
  offset_id[make_pair("decay", offset(0, -1))] = name_id["n"];
  offset_id[make_pair("decay", offset(0, -2))] = name_id["z_2n"];
  offset_id[make_pair("decay", offset(0, -3))] = name_id["z_3n"];
  offset_id[make_pair("decay", offset(0, -4))] = name_id["z_4n"];
  offset_id[make_pair("decay", offset(-1, -1))] = name_id["p"];
  offset_id[make_pair("decay", offset(-2, -2))] = name_id["z_2p"];
  offset_id[make_pair("decay", offset(-1, -2))] = name_id["d"];
  offset_id[make_pair("decay", offset(-1, -3))] = name_id["t"];
  offset_id[make_pair("decay", offset(-2, -3))] = name_id["He3"];
  offset_id[make_pair("decay", offset(-2, -4))] = name_id["a"];
  offset_id[make_pair("decay", offset(1, 0))] = name_id["bminus"];
  offset_id[make_pair("decay", offset(-1, 0))] = name_id["bplus"];
  offset_id[make_pair("decay", offset(1, -1))] = name_id["bminus_n"];
  offset_id[make_pair("decay", offset(-1, -4))] = name_id["bminus_a"];
  offset_id[make_pair("decay", offset(0, 0))] = name_id["it"];
  offset_id[make_pair("decay", offset(-3, -4))] = name_id["bplus_a"];
  offset_id[make_pair("decay", offset(-2, -1))] = name_id["bplus_p"];
  offset_id[make_pair("decay", offset(1, -2))] = name_id["bminus_2n"];
  offset_id[make_pair("decay", offset(1, -3))] = name_id["bminus_3n"];
  offset_id[make_pair("decay", offset(1, -4))] = name_id["bminus_4n"];
  offset_id[make_pair("decay", offset(-3, -2))] = name_id["bplus_2p"];
  offset_id[make_pair("decay", offset(-4, -3))] = name_id["bplus_3p"];
  offset_id[make_pair("decay", offset(2, 0))] = name_id["decay_2bminus"];
  offset_id[make_pair("decay", offset(-2, 0))] = name_id["decay_2bplus"];
  offset_id[make_pair("decay", offset(-6, -14))] = name_id["decay_14c"];

  // pre-loaded child offsets
  std::map<std::pair<std::string, int>, unsigned int>::iterator ioffid;
  for (ioffid = offset_id.begin(); ioffid != offset_id.end(); ioffid++) {
    id_offset[make_pair(ioffid->first.first, ioffid->second)] = ioffid->first.second;
  }
  // neutrons:
  id_offset[make_pair("n", name_id["nHe3"])] = offset(-2, -3);
  id_offset[make_pair("n", name_id["nHe3_1"])] = offset(-2, -3, 2);
  id_offset[make_pair("n", name_id["nHe3_2"])] = offset(-2, -3, 2);
  id_offset[make_pair("n", name_id["z_3np"])] = offset(-1, -3);
  id_offset[make_pair("n", name_id["nd"])] = offset(-1, -2);
  id_offset[make_pair("n", name_id["nd_1"])] = offset(-1, -2, 1);
  id_offset[make_pair("n", name_id["nd_2"])] = offset(-1, -2, 2);
  id_offset[make_pair("n", name_id["np"])] = offset(-1, -1);
  id_offset[make_pair("n", name_id["np_1"])] = offset(-1, -1, 1);
  id_offset[make_pair("n", name_id["np_2"])] = offset(-1, -1, 2);
  id_offset[make_pair("n", name_id["n"])] = offset(0, 0);
  id_offset[make_pair("n", name_id["gamma"])] = offset(0, 1);
  // decay:
  id_offset[make_pair("decay", name_id["bminus_p"])] = offset(0, -1);
  id_offset[make_pair("decay", name_id["ec_2p"])] = offset(-3, -2);
  id_offset[make_pair("decay", name_id["ec_3p"])] = offset(-4, -3);
  id_offset[make_pair("decay", name_id["ec"])] = offset(-1, 0);
  id_offset[make_pair("decay", name_id["ec_bplus"])] = offset(-1, 0);
  id_offset[make_pair("decay", name_id["ecp"])] = offset(-2, -1);
  id_offset[make_pair("decay", name_id["eca"])] = offset(-3, -4);
  id_offset[make_pair("decay", name_id["decay_2ec"])] = offset(-2, 0);
  return NULL;
}
void* pyne::rxname::_ = pyne::rxname::_fill_maps();


unsigned int pyne::rxname::hash(std::string s) {
  return pyne::rxname::hash(s.c_str());
}

unsigned int pyne::rxname::hash(const char* s) {
  // Modified from http://cboard.cprogramming.com/tech-board/114650-string-hashing-algorithm.html#post853145
  // starting from h = 32*2^5 > 1000, rather than 0, to reserve space for MT numbers
  int c;
  unsigned int h = 32;
  while ((c = *s++)) {
    h = ((h << 5) + h) ^ c;
  }
  return h;
}


// ************************
// *** name functions *****
// ************************

std::string pyne::rxname::name(char* s) {
  return pyne::rxname::name(std::string(s));
}

std::string pyne::rxname::name(std::string s) {
  if (0 < names.count(s))
    return s;
  if (0 < altnames.count(s))
    return id_name[altnames[s]];
  // see if id in string form
  int i = 0;
  int I = s.length();
  int found = 0;
  while (0 <= found && i < I) {
    found = pyne::digits.find(s[i]);
    i++;
  }
  if (0 <= found)
    return pyne::rxname::name(atoi(s.c_str()));
  // dead...
  throw NotAReaction(s, "???");
}


std::string pyne::rxname::name(int n) {
  return pyne::rxname::name((unsigned int) n);
}

std::string pyne::rxname::name(unsigned int n) {
  if (0 < id_name.count(n))
    return id_name[n];
  if (0 < mt_id.count(n))
    return id_name[mt_id[n]];
  throw NotAReaction(n, "???");
}


std::string pyne::rxname::name(int from_nuc, int to_nuc, std::string z) {
  // This assumes nuclides are in id form
  std::pair<std::string, int> key = std::make_pair(z, to_nuc - from_nuc);
  if (0 == offset_id.count(key))
    throw IndeterminateReactionForm("z=" + z + ", " + pyne::to_str(from_nuc) + \
                                    ", " + pyne::to_str(to_nuc), "???");
  return id_name[offset_id[key]];
}

std::string pyne::rxname::name(std::string from_nuc, int to_nuc, std::string z) {
  return pyne::rxname::name(pyne::nucname::id(from_nuc),
                            pyne::nucname::id(to_nuc), z);
}

std::string pyne::rxname::name(int from_nuc, std::string to_nuc, std::string z) {
  return pyne::rxname::name(pyne::nucname::id(from_nuc),
                            pyne::nucname::id(to_nuc), z);
}

std::string pyne::rxname::name(std::string from_nuc, std::string to_nuc, std::string z) {
  return pyne::rxname::name(pyne::nucname::id(from_nuc),
                            pyne::nucname::id(to_nuc), z);
}



// **********************
// *** id functions *****
// **********************
unsigned int pyne::rxname::id(int x) {
  return name_id[pyne::rxname::name(x)];
}

unsigned int pyne::rxname::id(unsigned int x) {
  if (0 < id_name.count(x))
    return x;
  if (0 < mt_id.count(x))
    return mt_id[x];
  return name_id[pyne::rxname::name(x)];
}

unsigned int pyne::rxname::id(const char* x) {
  return name_id[pyne::rxname::name(x)];
}

unsigned int pyne::rxname::id(std::string x) {
  if (0 < names.count(x))
    return name_id[x];
  if (0 < altnames.count(x))
    return altnames[x];
  return name_id[pyne::rxname::name(x)];
}

unsigned int pyne::rxname::id(int from_nuc, int to_nuc, std::string z) {
  // This assumes nuclides are in id form
  std::pair<std::string, int> key = std::make_pair(z, to_nuc - from_nuc);
  if (0 == offset_id.count(key))
    throw IndeterminateReactionForm("z=" + z + ", " + pyne::to_str(from_nuc) + \
                                    ", " + pyne::to_str(to_nuc), "???");
  return offset_id[key];
}

unsigned int pyne::rxname::id(int from_nuc, std::string to_nuc, std::string z) {
  return pyne::rxname::id(pyne::nucname::id(from_nuc),
                          pyne::nucname::id(to_nuc), z);
}

unsigned int pyne::rxname::id(std::string from_nuc, int to_nuc, std::string z) {
  return pyne::rxname::id(pyne::nucname::id(from_nuc),
                          pyne::nucname::id(to_nuc), z);
}

unsigned int pyne::rxname::id(std::string from_nuc, std::string to_nuc, std::string z) {
  return pyne::rxname::id(pyne::nucname::id(from_nuc),
                          pyne::nucname::id(to_nuc), z);
}


// **********************
// *** MT functions *****
// **********************
unsigned int pyne::rxname::mt(int x) {
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(unsigned int x) {
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(char* x) {
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(std::string x) {
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(int from_nuc, int to_nuc, std::string z) {
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(int from_nuc, std::string to_nuc, std::string z) {
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(std::string from_nuc, int to_nuc, std::string z) {
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}

unsigned int pyne::rxname::mt(std::string from_nuc, std::string to_nuc, std::string z) {
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
}


// ***********************
// *** label functions ***
// ***********************
std::string pyne::rxname::label(int x) {
  return labels[pyne::rxname::id(x)];
}

std::string pyne::rxname::label(unsigned int x) {
  return labels[pyne::rxname::id(x)];
}

std::string pyne::rxname::label(char* x) {
  return labels[pyne::rxname::id(x)];
}

std::string pyne::rxname::label(std::string x) {
  return labels[pyne::rxname::id(x)];
}

std::string pyne::rxname::label(int from_nuc, int to_nuc, std::string z) {
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::label(int from_nuc, std::string to_nuc, std::string z) {
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::label(std::string from_nuc, int to_nuc, std::string z) {
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::label(std::string from_nuc, std::string to_nuc, std::string z) {
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
}


// *********************
// *** doc functions ***
// *********************
std::string pyne::rxname::doc(int x) {
  return docs[pyne::rxname::id(x)];
}

std::string pyne::rxname::doc(unsigned int x) {
  return docs[pyne::rxname::id(x)];
}

std::string pyne::rxname::doc(char* x) {
  return docs[pyne::rxname::id(x)];
}

std::string pyne::rxname::doc(std::string x) {
  return docs[pyne::rxname::id(x)];
}

std::string pyne::rxname::doc(int from_nuc, int to_nuc, std::string z) {
  return docs[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::doc(int from_nuc, std::string to_nuc, std::string z) {
  return docs[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::doc(std::string from_nuc, int to_nuc, std::string z) {
  return docs[pyne::rxname::id(from_nuc, to_nuc, z)];
}

std::string pyne::rxname::doc(std::string from_nuc, std::string to_nuc, std::string z) {
  return docs[pyne::rxname::id(from_nuc, to_nuc, z)];
}


// ***********************
// *** child functions ***
// ***********************

int pyne::rxname::child(int nuc, unsigned int rx, std::string z) {
  // This assumes nuclides are in id form
  std::pair<std::string, unsigned int> key = std::make_pair(z, rx);
  if (0 == id_offset.count(key))
    throw IndeterminateReactionForm("z=" + z + ", rx=" + pyne::to_str(rx), "???");
  int to_nuc = pyne::nucname::groundstate(nuc) + id_offset[key];
  if (!pyne::nucname::isnuclide(to_nuc))
    throw pyne::nucname::NotANuclide(nuc, to_nuc);
  return to_nuc;
}

int pyne::rxname::child(int nuc, std::string rx, std::string z) {
  return child(nuc, id(rx), z);
}

int pyne::rxname::child(std::string nuc, unsigned int rx, std::string z) {
  return child(pyne::nucname::id(nuc), rx, z);
}

int pyne::rxname::child(std::string nuc, std::string rx, std::string z) {
  return child(pyne::nucname::id(nuc), id(rx), z);
}

// ************************
// *** parent functions ***
// ************************

int pyne::rxname::parent(int nuc, unsigned int rx, std::string z) {
  // This assumes nuclides are in id form
  std::pair<std::string, unsigned int> key = std::make_pair(z, rx);
  if (0 == id_offset.count(key))
    throw IndeterminateReactionForm("z=" + z + ", rx=" + pyne::to_str(rx), "???");
  int from_nuc = nuc - id_offset[key];
  if (!pyne::nucname::isnuclide(from_nuc))
    throw pyne::nucname::NotANuclide(from_nuc, nuc);
  return from_nuc;
}

int pyne::rxname::parent(int nuc, std::string rx, std::string z) {
  return parent(nuc, id(rx), z);
}

int pyne::rxname::parent(std::string nuc, unsigned int rx, std::string z) {
  return parent(pyne::nucname::id(nuc), rx, z);
}

int pyne::rxname::parent(std::string nuc, std::string rx, std::string z) {
  return parent(pyne::nucname::id(nuc), id(rx), z);
}

//
// end of src/rxname.cpp
//


//
// start of src/particle.cpp
//
#ifndef PYNE_IS_AMALGAMATED
#include "particle.h"
#endif

std::string pyne::particle::_names[NUM_PARTICLES] = {
  // leptons
  "Electron",
  "Positron",
  "ElectronNeutrino",
  "ElectronAntiNeutrino",
  "Muon",
  "AntiMuon",
  "MuonNeutrino",
  "MuonAntiNeutrino",
  "Tauon",
  "AntiTauon",
  "TauNeutrino",
  "TauAntiNeutrino",
  // gauge bosons
  "Photon",
  // light mesons
  "Pion",
  "AntiPion",
  // strange mesons
  "Kaon",
  "AntiKaon",
  "KaonZeroShort",
  "KaonZero",
  "AntiKaonZero",
  // light baryons
  "Neutron",
  "AntiNeutron",
  "Proton",
  "AntiProton",
  // strange baryons
  "Lambda",
  "AntiLambda",
  "Sigma-",
  "AntiSigma-",
  "Sigma+",
  "AntiSigma+",
  "Sigma",
  "AntiSigmaZero"
  // Charmed baryons
};

int pyne::particle::_pdcids[NUM_PARTICLES] = {
  11,
  -11,
  12,
  -12,
  13,
  -13,
  14,
  -14,
  15,
  -15,
  16,
  -16,
  // gauge bosons
  22,
  // light mesons
  211,
  -211,
  // strange mesons
  321,
  -321,
  310,
  311,
  -311,
  // light baryons
  2112,
  -2112,
  2212,
  -2212,
  // strange Baryons
  3122,
  -3122,
  3112,
  3112,
  3222,
  -3222,
  3212,
  -3212
  // charmed baryons
};

std::set<std::string> pyne::particle::names(pyne::particle::_names,
                                            pyne::particle::_names + NUM_PARTICLES);

std::set<int> pyne::particle::pdc_nums(pyne::particle::_pdcids,
                                       pyne::particle::_pdcids + NUM_PARTICLES);

std::map<std::string, int> pyne::particle::altnames;
std::map<int, std::string> pyne::particle::id_name;
std::map<std::string, int> pyne::particle::name_id;
std::map<std::string, std::string> pyne::particle::docs;

std::map<std::string, std::string> pyne::particle::part_to_fluka;
std::map<std::string, std::string> pyne::particle::part_to_mcnp;
std::map<std::string, std::string> pyne::particle::part_to_mcnp6;
std::map<std::string, std::string> pyne::particle::part_to_geant4;


void* pyne::particle::_fill_maps() {
  using std::make_pair;

  std::string _docs[NUM_PARTICLES] = {
    // leptons
    "Electron",
    "Positron",
    "Electron Neutrino",
    "Electron Anti Neutrino",
    "Muon Neutrino",
    "Anti Muon",
    "Muon Neutrino",
    "Muon Anti Neutrino",
    "Tauon",
    "Anti Tauon",
    "Tau Neutrino",
    "Tau Anti Neutrino",
    // gauge bosons
    "Photon",
    // light mesons
    "Pion",
    "Anti Pion",
    // strange mesons
    "Kaon",
    "Anti Kaon",
    "Kaon Zero Short",
    "Kaon Zero",
    "Anti Kaon Zero",
    // light baryons
    "Neutron",
    "Anti Neutron",
    "Proton",
    "Anti Proton",
    // strange baryons
    "Lambda",
    "Anti Lambda",
    "Sigma-",
    "Anti Sigma-",
    "Sigma+",
    "Anti Sigma+",
    "Sigma",
    "Anti Sigma Zero"
    // Charmed baryons
  };

  int pid;  // particle id
  for (int i = 0 ; i < NUM_PARTICLES ; i++) {
    pid = _pdcids[i];
    // make id to name map
    id_name[pid] = _names[i];
    // make name to id map
    name_id[_names[i]] = pid;
    // make doc correspondence
    docs[_names[i]] = _docs[i];
  }

  // make the alternates
  altnames["Hydrogen"] = name_id["Proton"];
  altnames["Protium"] = name_id["Proton"];
  altnames["Beta"] = name_id["Electron"];
  altnames["Beta-"] = name_id["Electron"];
  altnames["Beta+"] = name_id["Positron"];
  altnames["Gamma"] = name_id["Photon"];
  altnames["X-Ray"] = name_id["Photon"];

  part_to_mcnp["Neutron"] = "n";
  part_to_mcnp["Photon"] = "p";
  part_to_mcnp["Electron"] = "e";

  part_to_mcnp6["Neutron"] = "n";
  part_to_mcnp6["Photon"] = "p";
  part_to_mcnp6["Electron"] = "e";
  part_to_mcnp6["Proton"] = "h";

  part_to_fluka["Electron"] = "ELECTRON";
  part_to_fluka["Positron"] = "POSITRON";
  part_to_fluka["ElectronNeutrino"] = "NEUTRIE";
  part_to_fluka["ElectronAntiNeutrino"] = "ANEUTRIE";
  part_to_fluka["Muon"] = "MUON+";
  part_to_fluka["AntiMuon"] = "MUON-";
  part_to_fluka["MuonNeutrino"] = "NEUTRIM";
  part_to_fluka["MuonAntiNeutrino"] = "ANEUTRIM";
  part_to_fluka["Tauon"] = "TAU+";
  part_to_fluka["Anti Tauon"] = "TAU-";
  part_to_fluka["TauNeutrino"] = "NEUTRIT";
  part_to_fluka["TauAntiNeutrino"] = "ANEUTRIT";
  // gauge bosons
  part_to_fluka["Photon"] = "PHOTON";
  // light mesons
  part_to_fluka["Pion"] = "PION-";
  part_to_fluka["Anti Pion"] = "PION+";
  // strange mesons
  part_to_fluka["Kaon"] = "KAON+";
  part_to_fluka["AntiKaon"] = "KAON-";
  part_to_fluka["KaonZero Short"] = "KAONSHRT";
  part_to_fluka["KaonZero"] = "KAONZERO";
  part_to_fluka["AntiKaonZero"] = "AKAONZER";
  // light baryons
  part_to_fluka["Neutron"] = "NEUTRON";
  part_to_fluka["AntiNeutron"] = "ANEUTRON";
  part_to_fluka["Proton"] = "PROTON";
  part_to_fluka["AntiProton"] = "APROTON";
  // strange baryons
  part_to_fluka["Lambda"] = "LAMBDA";
  part_to_fluka["AntiLambda"] = "ALAMBDA";
  part_to_fluka["Sigma-"] = "SIGMA-";
  part_to_fluka["Anti Sigma-"] = "ASIGMA-";
  part_to_fluka["Sigma+"] = "SIGMA+";
  part_to_fluka["Anti Sigma+"] = "ASIGMA+";
  part_to_fluka["Sigma"] = "SIGMAZER";
  part_to_fluka["AntiSigmaZero"] = "ASIGMAZE";

  part_to_geant4["Electron"] = "e-";
  part_to_geant4["Positron"] = "e+";
  part_to_geant4["ElectronNeutrino"] = "nu_e";
  part_to_geant4["ElectronAntiNeutrino"] = "anti_nu_e";
  part_to_geant4["Muon"] = "mu+";
  part_to_geant4["AntiMuon"] = "mu-";
  part_to_geant4["MuonNeutrino"] = "nu_mu";
  part_to_geant4["MuonAntiNeutrino"] = "anti_nu_mu";
  part_to_geant4["Tauon"] = "tau+";
  part_to_geant4["Anti Tauon"] = "tau-";
  part_to_geant4["TauNeutrino"] = "nu_tau";
  part_to_geant4["TauAntiNeutrino"] = "anti_nu_tau";
  // gauge bosons
  part_to_geant4["Photon"] = "gamma";
  // light mesons
  part_to_geant4["Pion"] = "pi-";
  part_to_geant4["Anti Pion"] = "pi+";
  // strange mesons
  part_to_geant4["Kaon"] = "kaon+";
  part_to_geant4["AntiKaon"] = "kaon-";
  part_to_geant4["KaonZero Short"] = "kaon0S";
  part_to_geant4["KaonZero"] = "kaon0";
  // light baryons
  part_to_geant4["Neutron"] = "neutron";
  part_to_geant4["AntiNeutron"] = "anti_neutron";
  part_to_geant4["Proton"] = "proton";
  part_to_geant4["AntiProton"] = "anti_proton";
  // strange baryons
  part_to_geant4["Lambda"] = "lambda";
  part_to_geant4["AntiLambda"] = "anti_lambda";
  part_to_geant4["Sigma-"] = "sigma-";
  part_to_geant4["Anti Sigma-"] = "anti_sigma-";
  part_to_geant4["Sigma+"] = "sigma+";
  part_to_geant4["Anti Sigma+"] = "anti_sigma+";
  part_to_geant4["Sigma"] = "sigma0";
  part_to_geant4["AntiSigmaZero"] = "anti_sigma0";
  return NULL;
}

void* pyne::particle::filler = pyne::particle::_fill_maps();

// is hydrogen
bool pyne::particle::is_hydrogen(int s) {
  if (s == name_id["Proton"])
    return true;
  if (pyne::particle::is_hydrogen(pyne::nucname::name(s)))
    return true;
  return false;
}

bool pyne::particle::is_hydrogen(char* s) {
  return pyne::particle::is_hydrogen(std::string(s));
}

bool pyne::particle::is_hydrogen(std::string s) {
  // check std name
  if (name_id[s] == name_id["Proton"])
    return true;
  if (altnames[s] == name_id["Proton"])
    return true;
  if (pyne::nucname::name(s).find("H1") != std::string::npos)
    return true;
  return false;
}
// heavy ion
bool pyne::particle::is_heavy_ion(int s) {
  return pyne::particle::is_heavy_ion(std::string(id_name[s]));
}

bool pyne::particle::is_heavy_ion(char* s) {
  return pyne::particle::is_heavy_ion(std::string(s));
}

bool pyne::particle::is_heavy_ion(std::string s) {
  if (pyne::nucname::isnuclide(s)) {
    if (pyne::particle::is_hydrogen(s))
      return false;
    else
      return true;
  }
  return false;
}

// is valid functions
bool pyne::particle::is_valid(int s) {
  if (pyne::nucname::isnuclide(s))
    return true;
  else
    return pyne::particle::is_valid(std::string(id_name[s]));
}

bool pyne::particle::is_valid(char* s) {
  return pyne::particle::is_valid(std::string(s));
}

bool pyne::particle::is_valid(std::string s) {
  // check std name
  if (0 < names.count(s))
    return true;
  // check alternative name
  if (0 < altnames.count(s))
    return true;
  // check if is a heavy ion
  if (pyne::nucname::isnuclide(s))
    return true;
  else
    return false;
}

// pdc functions
int pyne::particle::id(int s) {
  if (0 < pdc_nums.count(s))
    return s;
  else
    return 0;
}

int pyne::particle::id(char* s) {
  return pyne::particle::id(std::string(s));
}

int pyne::particle::id(std::string s) {
  if (pyne::nucname::isnuclide(s)) {
    if (pyne::particle::is_hydrogen(s))
      return name_id["Proton"];
    if (pyne::particle::is_heavy_ion(s))
      return 0;
  }

  if (0 < pdc_nums.count(name_id[s]))
    return name_id[s];
  if (0 < pdc_nums.count(altnames[s]))
    return altnames[s];
  return 0;
}

// name functions
std::string pyne::particle::name(int s) {
  if (s < 9999999)
    return pyne::particle::name(id_name[s]);
  if (pyne::nucname::isnuclide(s))
    return pyne::particle::name(pyne::nucname::name(s));
  return pyne::particle::name(id_name[s]);
}

std::string pyne::particle::name(char* s) {
  return pyne::particle::name(std::string(s));
}

std::string pyne::particle::name(std::string s) {
  // check if is a hydrogen
  if (pyne::nucname::isnuclide(s)) {
    if (pyne::particle::is_hydrogen(s))
      return "Proton";
    if (pyne::particle::is_heavy_ion(s))
      return s;
  }
  // check std name
  if (0 < names.count(s))
    return s;
  // check alternative name
  if (0 < altnames.count(s))
    return id_name[altnames[s]];
  // check for heavy ion
  else
    throw NotAParticle(s);
}

// convert name to mcnp id
std::string pyne::particle::mcnp(int s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(char* s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(std::string s) {
  if (0 < part_to_mcnp.count(pyne::particle::name(s)))
    return part_to_mcnp[pyne::particle::name(s)];
  else {
    std::cout << "Not a valid MCNP5 particle" << std::endl;
    return "?";
  }
}

// convert name to mcnp6 id
std::string pyne::particle::mcnp6(int s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(char* s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(std::string s) {
  if (0 < part_to_mcnp6.count(pyne::particle::name(s)))
    return part_to_mcnp6[pyne::particle::name(s)];
  else
    return "?";
}

// convert name to fluka id
std::string pyne::particle::fluka(int s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(char* s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "HEAVYION";
  else if (0 < part_to_fluka.count(pyne::particle::name(s)))
    return part_to_fluka[pyne::particle::name(s)];
  else {
    std::cout << "Not a valid Fluka particle" << std::endl;
    return "???????";
  }
}

// convert name to geant4 id
std::string pyne::particle::geant4(int s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(char* s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "GenericIon";
  else if (0 < part_to_geant4.count(pyne::particle::name(s)))
    return part_to_geant4[pyne::particle::name(s)];
  else {
    std::cout << "Not a valid Geant4 particle" << std::endl;
    return "???????";
  }
}


// describe functions
std::string pyne::particle::describe(int s) {
  if (pyne::nucname::isnuclide(s))
    return pyne::particle::describe(pyne::nucname::name(s));
  return pyne::particle::describe(id_name[s]);
}

std::string pyne::particle::describe(char* s) {
  return pyne::particle::describe(std::string(s));
}

std::string pyne::particle::describe(std::string s) {
  // check if is a hydrogen
  if (pyne::nucname::isnuclide(s)) {
    if (pyne::particle::is_hydrogen(s))
      return docs[pyne::particle::name(s)];
    if (pyne::particle::is_heavy_ion(s))
      return "Is a heavy ion";
  }
  // check std name
  if (0 < names.count(s))
    return docs[s];
  // check alternative name
  if (0 < altnames.count(s))
    return docs[id_name[altnames[s]]];
  // check if is a heavy ion
  else
    throw NotAParticle(s);
}
//
// end of src/particle.cpp
//


//
// start of src/data.cpp
//
// Implements basic nuclear data functions.
#ifndef PYNE_IS_AMALGAMATED
#include "data.h"
#include "atomic_data.h"
#endif

//
// Math Helpers
//

const double pyne::pi = 3.14159265359;
const double pyne::N_A = 6.0221415e+23;
const double pyne::barns_per_cm2 = 1e24;
const double pyne::cm2_per_barn = 1e-24;
const double pyne::sec_per_day = 24.0 * 3600.0;
const double pyne::MeV_per_K = 8.617343e-11;
const double pyne::MeV_per_MJ = 6.2415096471204E+18;
const double pyne::Bq_per_Ci = 3.7e10;
const double pyne::Ci_per_Bq = 2.7027027e-11;


/********************************/
/*** data_checksums Functions ***/
/********************************/

std::map<std::string, std::string> pyne::get_data_checksums() {
  std::map<std::string, std::string> temp_map;
  // Initialization of dataset hashes
  temp_map["/atomic_mass"] = "10edfdc662e35bdfab91beb89285efff";
  temp_map["/material_library"] = "8b10864378fbd88538434679acf908cc";
  temp_map["/neutron/eaf_xs"] = "29622c636c4a3a46802207b934f9516c";
  temp_map["/neutron/scattering_lengths"] = "a24d391cc9dc0fc146392740bb97ead4";
  temp_map["/neutron/simple_xs"] = "3d6e086977783dcdf07e5c6b0c2416be";
  temp_map["/decay"] = "4f41f3e46f4306cc44449f08a20922e0";
  temp_map["/dose_factors"] = "dafa32c24b2303850a0bebdf3e6b122e";
  return temp_map;
}

std::map<std::string, std::string> pyne::data_checksums =
    pyne::get_data_checksums();

/*****************************/
/*** atomic_mass Functions ***/
/*****************************/
std::map<int, double> pyne::atomic_mass_map = std::map<int, double>();

void pyne::_load_atomic_mass_map() {
  // Loads the important parts of atomic_wight table into atomic_mass_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH)) {
    pyne::_load_atomic_mass_map_memory();
    return;
  }

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic_mass_data));
  H5Tinsert(desc, "nuc",   HOFFSET(atomic_mass_data, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "mass",  HOFFSET(atomic_mass_data, mass),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "error", HOFFSET(atomic_mass_data, error), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "abund", HOFFSET(atomic_mass_data, abund), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atomic_mass_set = H5Dopen2(nuc_data_h5, "/atomic_mass", H5P_DEFAULT);
  hid_t atomic_mass_space = H5Dget_space(atomic_mass_set);
  int atomic_mass_length = H5Sget_simple_extent_npoints(atomic_mass_space);

  // Read in the data
  atomic_mass_data* atomic_mass_array = new atomic_mass_data[atomic_mass_length];
  H5Dread(atomic_mass_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_mass_array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(atomic_mass_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of structs, put it in the map
  for (int n = 0; n < atomic_mass_length; n++) {
    atomic_mass_map.insert(std::pair<int, double>(atomic_mass_array[n].nuc, \
                                                  atomic_mass_array[n].mass));
    natural_abund_map.insert(std::pair<int, double>(atomic_mass_array[n].nuc, \
                                                    atomic_mass_array[n].abund));
  }

  delete[] atomic_mass_array;
}


double pyne::atomic_mass(int nuc) {
  // Find the nuclide's mass in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = atomic_mass_map.find(nuc);
  nuc_end = atomic_mass_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end) {
    return (*nuc_iter).second;
  }

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (atomic_mass_map.empty()) {
    // Don't fail if we can't load the library
    _load_atomic_mass_map();
    return atomic_mass(nuc);
  }

  double aw;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state mass...not strictly true, but good guess.
  if (0 < nucid % 10000) {
    aw = atomic_mass((nucid / 10000) * 10000);
    if (atomic_mass_map.count(nuc) != 1) {
      atomic_mass_map.insert(std::pair<int, double>(nuc, aw));
    }
    return aw;
  };

  // Finally, if none of these work,
  // take a best guess based on the
  // aaa number.
  aw = (double)((nucid / 10000) % 1000);
  if (atomic_mass_map.count(nuc) != 1) {
    atomic_mass_map.insert(std::pair<int, double>(nuc, aw));
  }
  return aw;
}


double pyne::atomic_mass(char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
}


double pyne::atomic_mass(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
}


/*******************************/
/*** natural_abund functions ***/
/*******************************/

std::map<int, double> pyne::natural_abund_map = std::map<int, double>();

double pyne::natural_abund(int nuc) {
  // Find the nuclide's natural abundance
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = natural_abund_map.find(nuc);
  nuc_end = natural_abund_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (natural_abund_map.empty()) {
    // Don't fail if we can't load the library
    _load_atomic_mass_map();
    return natural_abund(nuc);
  }

  double na;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state abundance...not strictly true, but good guess.
  if (0 < nucid % 10000) {
    na = natural_abund((nucid / 10000) * 10000);
    natural_abund_map[nuc] = na;
    return na;
  }

  // Finally, if none of these work,
  // take a best guess based on the
  // aaa number.
  na = 0.0;
  natural_abund_map[nuc] = na;
  return na;
}


double pyne::natural_abund(char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
}


double pyne::natural_abund(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
}



/*************************/
/*** Q_value Functions ***/
/*************************/

void pyne::_load_q_val_map() {
  // Loads the important parts of q_value table into q_value_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(q_val_data));
  H5Tinsert(desc, "nuc", HOFFSET(q_val_data, nuc),  H5T_NATIVE_INT);
  H5Tinsert(desc, "q_val", HOFFSET(q_val_data, q_val), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "gamma_frac", HOFFSET(q_val_data, gamma_frac), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t q_val_set = H5Dopen2(nuc_data_h5, "/decay/q_values", H5P_DEFAULT);
  hid_t q_val_space = H5Dget_space(q_val_set);
  int q_val_length = H5Sget_simple_extent_npoints(q_val_space);

  // Read in the data
  q_val_data* q_val_array = new q_val_data[q_val_length];
  H5Dread(q_val_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, q_val_array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(q_val_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of structs, put it in the map
  for (int n = 0; n < q_val_length; n++) {
    q_val_map[q_val_array[n].nuc] = q_val_array[n].q_val;
    gamma_frac_map[q_val_array[n].nuc] = q_val_array[n].gamma_frac;
  }

  delete[] q_val_array;
}

std::map<int, double> pyne::q_val_map = std::map<int, double>();

double pyne::q_val(int nuc) {
  // Find the nuclide's q_val in MeV/fission
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = q_val_map.find(nuc);
  nuc_end = q_val_map.end();

  // First check if we already have the nuc q_val in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the nuc_data.h5 if the map is empty.
  if (q_val_map.empty()) {
    _load_q_val_map();
    return q_val(nuc);
  };

  double qv;
  int nucid = nucname::id(nuc);
  if (nucid != nuc)
    return q_val(nucid);

  // If nuclide is not found, return 0
  qv = 0.0;
  q_val_map[nuc] = qv;
  return qv;
}


double pyne::q_val(const char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return q_val(nuc_zz);
}


double pyne::q_val(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return q_val(nuc_zz);
}


/****************************/
/*** gamma_frac functions ***/
/****************************/

std::map<int, double> pyne::gamma_frac_map = std::map<int, double>();

double pyne::gamma_frac(int nuc) {
  // Find the nuclide's fraction of Q that comes from gammas
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = gamma_frac_map.find(nuc);
  nuc_end = gamma_frac_map.end();

  // First check if we already have the gamma_frac in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from nuc_data.h5 if the map is empty.
  if (gamma_frac_map.empty()) {
    _load_q_val_map();
    return gamma_frac(nuc);
  }

  double gf;
  int nucid = nucname::id(nuc);
  if (nucid != nuc)
    return gamma_frac(nucid);

  // If nuclide is not found, return 0
  gf = 0.0;
  gamma_frac_map[nucid] = gf;
  return gf;
}


double pyne::gamma_frac(const char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return gamma_frac(nuc_zz);
}


double pyne::gamma_frac(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return gamma_frac(nuc_zz);
}


/*****************************/
/*** Dose Factor Functions ***/
/*****************************/
/***************************************************************************
This data is from: [Exposure Scenarios and Unit Dose Factors for the Hanford
Immobilized Low-Activity Tank Waste Performance Assessment, ref.
HNF-SD-WM-TI-707 Rev. 1 December 1999] Appendix O of HNF-5636 [DATA PACKAGES
FOR THE HANFORD IMMOBILIZED LOW-ACTIVITY TANK WASTE PERFORMANCE ASSESSMENT:
2001 VERSION]

Liability Disclaimer: The PyNE Development Team shall not be liable for any
loss or injury resulting from decisions made with this data.
**************************************************************************/

void pyne::_load_dose_map(std::map<int, dose>& dm, std::string source_path) {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Defining string type for lung model data
  hid_t string_type_;
  string_type_ = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_, 1);
  H5Tset_strpad(string_type_, H5T_STR_NULLPAD);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(dose));
  status = H5Tinsert(desc, "nuc", HOFFSET(dose, nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "ext_air_dose", HOFFSET(dose, ext_air_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ratio", HOFFSET(dose, ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ext_soil_dose", HOFFSET(dose, ext_soil_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ingest_dose", HOFFSET(dose, ingest_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "fluid_frac", HOFFSET(dose, fluid_frac), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "inhale_dose", HOFFSET(dose, inhale_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "lung_mod", HOFFSET(dose, lung_mod), string_type_);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Convert source_path to proper format for HD5open
  const char* c = source_path.c_str();

  // Open the data set
  hid_t dose_set = H5Dopen2(nuc_data_h5, c, H5P_DEFAULT);
  hid_t dose_space = H5Dget_space(dose_set);
  int dose_length = H5Sget_simple_extent_npoints(dose_space);

  // Read in the data
  dose* dose_array = new dose[dose_length];
  H5Dread(dose_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, dose_array);

  // Put array of structs in the map
  for (int n = 0; n < dose_length; n++) {
    dm[dose_array[n].nuc] = dose_array[n];
  }

  // Close the nuc_data library
  H5Dclose(dose_set);
  H5Tclose(string_type_);
  H5Fclose(nuc_data_h5);

  delete[] dose_array;
}

///
/// Functions for Source Location in nuc_data.h5
/// and related Map Pointers
///

std::string source_string(int source) {
  std::string source_location;
  if (source == 1) {
    source_location = "/dose_factors/DOE";
  } else if (source == 2) {
    source_location = "/dose_factors/GENII";
  } else {
    source_location = "/dose_factors/EPA";
  }
  return source_location;
}

std::map<int, pyne::dose>& dose_source_map(int source) {
  std::map<int, pyne::dose>* dm;
  if (source == 1) {
    dm = &pyne::doe_dose_map;
  } else if (source == 2) {
    dm = &pyne::genii_dose_map;
  } else {
    dm = &pyne::epa_dose_map;
  }
  if (dm->empty()) {
    std::string source_path = source_string(source);
    _load_dose_map(*dm, source_path);
  }
  return *dm;
}

std::map<int, pyne::dose> pyne::epa_dose_map;
std::map<int, pyne::dose> pyne::doe_dose_map;
std::map<int, pyne::dose> pyne::genii_dose_map;

///
/// Functions for External Air and
/// Ratio of External Air to Inhalation Dose Factors
///

/// External Air
double pyne::ext_air_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].ext_air_dose;
  } else {
    return -1;
  }
}

double pyne::ext_air_dose(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_air_dose(nuc_zz, source);
}


double pyne::ext_air_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_air_dose(nuc_zz, source);
}

/// Dose Ratio
double pyne::dose_ratio(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].ratio;
  } else {
    return -1;
  }
}

double pyne::dose_ratio(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_ratio(nuc_zz, source);
}


double pyne::dose_ratio(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_ratio(nuc_zz, source);
}

///
/// Function for External Soil Dose Factors
///

double pyne::ext_soil_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].ext_soil_dose;
  } else {
    return -1;
  }
}

double pyne::ext_soil_dose(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_soil_dose(nuc_zz, source);
}


double pyne::ext_soil_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_soil_dose(nuc_zz, source);
}

///
/// Functions for Ingestion Dose Factors and
/// Fraction of activity that is absorbed by body fluids
///

/// Ingestion
double pyne::ingest_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].ingest_dose;
  } else {
    return -1;
  }
}

double pyne::ingest_dose(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ingest_dose(nuc_zz, source);
}

double pyne::ingest_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ingest_dose(nuc_zz, source);
}

/// Fluid Fraction
double pyne::dose_fluid_frac(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].fluid_frac;
  } else {
    return -1;
  }
}

double pyne::dose_fluid_frac(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_fluid_frac(nuc_zz, source);
}

double pyne::dose_fluid_frac(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_fluid_frac(nuc_zz, source);
}

///
/// Functions for Inhalation Dose Factors and
/// Lung Model used to obtain dose factors
///

/// Inhalation
double pyne::inhale_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return dm[nucid].inhale_dose;
  } else {
    return -1;
  }
}

double pyne::inhale_dose(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return inhale_dose(nuc_zz, source);
}

double pyne::inhale_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return inhale_dose(nuc_zz, source);
}

/// Lung Model
std::string pyne::dose_lung_model(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid) == 1) {
    return std::string(1, dm[nucid].lung_mod);
  } else {
    return "Nada";
  }
}

std::string pyne::dose_lung_model(const char* nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_lung_model(nuc_zz, source);
}


std::string pyne::dose_lung_model(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_lung_model(nuc_zz, source);
}


/***********************************/
/*** scattering length functions ***/
/***********************************/
std::map<int, xd_complex_t> pyne::b_coherent_map = std::map<int, xd_complex_t>();
std::map<int, xd_complex_t> pyne::b_incoherent_map = std::map<int, xd_complex_t>();
std::map<int, double> pyne::b_map = std::map<int, double>();


void pyne::_load_scattering_lengths() {
  // Loads the important parts of atomic_wight table into atomic_mass_map
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(scattering_lengths));
  status = H5Tinsert(desc, "nuc", HOFFSET(scattering_lengths, nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "b_coherent", HOFFSET(scattering_lengths, b_coherent),
                     h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "b_incoherent", HOFFSET(scattering_lengths, b_incoherent),
                     h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "xs_coherent", HOFFSET(scattering_lengths, xs_coherent),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs_incoherent", HOFFSET(scattering_lengths, xs_incoherent),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs", HOFFSET(scattering_lengths, xs), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t scat_len_set = H5Dopen2(nuc_data_h5, "/neutron/scattering_lengths", H5P_DEFAULT);
  hid_t scat_len_space = H5Dget_space(scat_len_set);
  int scat_len_length = H5Sget_simple_extent_npoints(scat_len_space);

  // Read in the data
  scattering_lengths* scat_len_array = new scattering_lengths[scat_len_length];
  status = H5Dread(scat_len_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, scat_len_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(scat_len_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for (int n = 0; n < scat_len_length; n++) {
    b_coherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_coherent;
    b_incoherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_incoherent;
  }

  delete[] scat_len_array;
}



//
// Coherent functions
//


xd_complex_t pyne::b_coherent(int nuc) {
  // Find the nuclide's bound scattering length in cm
  std::map<int, xd_complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_coherent_map.find(nuc);
  nuc_end = b_coherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (b_coherent_map.empty()) {
    _load_scattering_lengths();
    return b_coherent(nuc);
  }

  xd_complex_t bc;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (anum == nucname::anum((*nuc_iter).first)) {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    }
    nuc_iter++;
  }

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (znum == nucname::znum((*nuc_iter).first)) {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    }
    nuc_iter++;
  }

  // Finally, if none of these work,
  // just return zero...
  bc.re = 0.0;
  bc.im = 0.0;
  b_coherent_map[nuc] = bc;
  return bc;
}


xd_complex_t pyne::b_coherent(char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
}


xd_complex_t pyne::b_coherent(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
}



//
// Incoherent functions
//


xd_complex_t pyne::b_incoherent(int nuc) {
  // Find the nuclide's bound inchoherent scattering length in cm
  std::map<int, xd_complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_incoherent_map.find(nuc);
  nuc_end = b_incoherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (b_incoherent_map.empty()) {
    _load_scattering_lengths();
    return b_incoherent(nuc);
  }

  xd_complex_t bi;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (anum == nucname::anum((*nuc_iter).first)) {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    }
    nuc_iter++;
  }

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (znum == nucname::znum((*nuc_iter).first)) {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    }
    nuc_iter++;
  }

  // Finally, if none of these work,
  // just return zero...
  bi.re = 0.0;
  bi.im = 0.0;
  b_incoherent_map[nuc] = bi;
  return bi;
}


xd_complex_t pyne::b_incoherent(char* nuc) {
  return b_incoherent(nucname::id(nuc));
}


xd_complex_t pyne::b_incoherent(std::string nuc) {
  return b_incoherent(nucname::id(nuc));
}



//
// b functions
//

double pyne::b(int nuc) {
  // Find the nuclide's bound scattering length in cm
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = b_map.find(nuc);
  nuc_end = b_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, calculate the value from coherent and incoherent lengths
  xd_complex_t bc = b_coherent(nuc);
  xd_complex_t bi = b_incoherent(nuc);

  double b_val = sqrt(bc.re * bc.re + bc.im * bc.im + bi.re * bi.re + bi.im * bi.im);

  return b_val;
}


double pyne::b(char* nuc) {
  int nucid = nucname::id(nuc);
  return b(nucid);
}


double pyne::b(std::string nuc) {
  int nucid = nucname::id(nuc);
  return b(nucid);
}



//
// Fission Product Yield Data
//
std::map<std::pair<int, int>, double> pyne::wimsdfpy_data = \
                                                            std::map<std::pair<int, int>, double>();

void pyne::_load_wimsdfpy() {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(wimsdfpy));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(wimsdfpy, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(wimsdfpy, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yields", HOFFSET(wimsdfpy, yields),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t wimsdfpy_set = H5Dopen2(nuc_data_h5, "/neutron/wimsd_fission_products",
                                H5P_DEFAULT);
  hid_t wimsdfpy_space = H5Dget_space(wimsdfpy_set);
  int wimsdfpy_length = H5Sget_simple_extent_npoints(wimsdfpy_space);

  // Read in the data
  wimsdfpy* wimsdfpy_array = new wimsdfpy[wimsdfpy_length];
  status = H5Dread(wimsdfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, wimsdfpy_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(wimsdfpy_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for (int n = 0; n < wimsdfpy_length; n++) {
    wimsdfpy_data[std::make_pair(wimsdfpy_array[n].from_nuc,
                                 wimsdfpy_array[n].to_nuc)] = wimsdfpy_array[n].yields;
  }

  delete[] wimsdfpy_array;
}


std::map<std::pair<int, int>, pyne::ndsfpysub> pyne::ndsfpy_data = \
                                                                   std::map<std::pair<int, int>, pyne::ndsfpysub>();

void pyne::_load_ndsfpy() {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ndsfpy));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ndsfpy, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ndsfpy, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yield_thermal", HOFFSET(ndsfpy, yield_thermal),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_thermal_err", HOFFSET(ndsfpy, yield_thermal_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast", HOFFSET(ndsfpy, yield_fast),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast_err", HOFFSET(ndsfpy, yield_fast_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV", HOFFSET(ndsfpy, yield_14MeV),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV_err", HOFFSET(ndsfpy, yield_14MeV_err),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t ndsfpy_set = H5Dopen2(nuc_data_h5, "/neutron/nds_fission_products",
                              H5P_DEFAULT);
  hid_t ndsfpy_space = H5Dget_space(ndsfpy_set);
  int ndsfpy_length = H5Sget_simple_extent_npoints(ndsfpy_space);

  // Read in the data
  ndsfpy* ndsfpy_array = new ndsfpy[ndsfpy_length];
  status = H5Dread(ndsfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ndsfpy_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(ndsfpy_set);
  status = H5Fclose(nuc_data_h5);

  ndsfpysub ndsfpysub_temp;

  // Ok now that we have the array of structs, put it in the maps
  for (int n = 0; n < ndsfpy_length; n++) {
    ndsfpysub_temp.yield_thermal = ndsfpy_array[n].yield_thermal;
    ndsfpysub_temp.yield_thermal_err = ndsfpy_array[n].yield_thermal_err;
    ndsfpysub_temp.yield_fast = ndsfpy_array[n].yield_fast;
    ndsfpysub_temp.yield_fast_err = ndsfpy_array[n].yield_fast_err;
    ndsfpysub_temp.yield_14MeV = ndsfpy_array[n].yield_14MeV;
    ndsfpysub_temp.yield_14MeV_err = ndsfpy_array[n].yield_14MeV_err;
    ndsfpy_data[std::make_pair(ndsfpy_array[n].from_nuc,
                               ndsfpy_array[n].to_nuc)] = ndsfpysub_temp;
  }



  delete[] ndsfpy_array;
}

double pyne::fpyield(std::pair<int, int> from_to, int source, bool get_error) {
  // Note that this may be expanded eventually to include other
  // sources of fission product data.

  // Find the parent/child pair branch ratio as a fraction
  if (source == 0) {
    std::map<std::pair<int, int>, double>::iterator fpy_iter, fpy_end;
    fpy_iter = wimsdfpy_data.find(from_to);
    fpy_end = wimsdfpy_data.end();
    if (fpy_iter != fpy_end)
      //if (get_error == true) return 0;
      return (*fpy_iter).second;
  } else {
    std::map<std::pair<int, int>, ndsfpysub>::iterator fpy_iter, fpy_end;
    fpy_iter = ndsfpy_data.find(from_to);
    fpy_end = ndsfpy_data.end();
    if (fpy_iter != fpy_end) {
      switch (source) {
        case 1:
          if (get_error)
            return (*fpy_iter).second.yield_thermal_err;
          return (*fpy_iter).second.yield_thermal;
          break;
        case 2:
          if (get_error)
            return (*fpy_iter).second.yield_fast_err;
          return (*fpy_iter).second.yield_fast;
          break;
        case 3:
          if (get_error)
            return (*fpy_iter).second.yield_14MeV_err;
          return (*fpy_iter).second.yield_14MeV;
          break;
      }
    }
  }


  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if ((source == 0) && (wimsdfpy_data.empty())) {
    _load_wimsdfpy();
    return fpyield(from_to, 0, get_error);
  } else if (ndsfpy_data.empty()) {
    _load_ndsfpy();
    return fpyield(from_to, source, get_error);
  }

  // Finally, if none of these work,
  // assume the value is stable
  double fpy = 0.0;
  wimsdfpy_data[from_to] = fpy;
  return fpy;
}

double pyne::fpyield(int from_nuc, int to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}

double pyne::fpyield(char* from_nuc, char* to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}

double pyne::fpyield(std::string from_nuc, std::string to_nuc, int source,
                     bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}


/***********************/
/*** decay functions ***/
/***********************/

//
// Data access tools
//

bool pyne::swapmapcompare::operator()(const std::pair<int, double>& lhs,
                                      const std::pair<int, double>& rhs) const {
  return lhs.second < rhs.second || (!(rhs.second < lhs.second) &&
                                     lhs.first < rhs.first);
}

template<typename T, typename U> std::vector<T> pyne::data_access(
    double energy_min, double energy_max, size_t valoffset, std::map<std::pair<int,
    double>, U>&  data) {
  typename std::map<std::pair<int, double>, U, swapmapcompare>::iterator
  nuc_iter, nuc_end, it;
  std::map<std::pair<int, double>, U, swapmapcompare> dc(data.begin(),
                                                         data.end());
  std::vector<T> result;
  if (energy_max < energy_min) {
    double temp = energy_max;
    energy_max = energy_min;
    energy_min = temp;
  }
  nuc_iter = dc.lower_bound(std::make_pair(0, energy_min));
  nuc_end = dc.upper_bound(std::make_pair(9999999999, energy_max));
  T* ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it != nuc_end; ++it) {
    ret = (T*)((char*) & (it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<T, U>(energy_min, energy_max, valoffset, data);
  }
  return result;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
                                                                  double min, double max, size_t valoffset,
                                                                  std::map<std::pair<int, double>, U>&  data) {
  typename std::map<std::pair<int, double>, U>::iterator nuc_iter, nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent, min));
  nuc_end = data.upper_bound(std::make_pair(parent, max));
  T* ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it != nuc_end; ++it) {
    ret = (T*)((char*) & (it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<T, U>(parent, min, max, valoffset, data);
  }
  return result;
}

template<typename T, typename U> T pyne::data_access(std::pair<int, int>
                                                     from_to, size_t valoffset, std::map<std::pair<int, int>, U>& data) {
  typename std::map<std::pair<int, int>, U>::iterator nuc_iter, nuc_end;

  nuc_iter = data.find(from_to);
  nuc_end = data.end();
  T* ret;
  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end) {
    ret = (T*)((char*) & (nuc_iter->second) + valoffset);
    return *ret;
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<T, U>(from_to, valoffset, data);
  }
  // This is okay for now because we only return ints and doubles
  return 0;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
                                                                  size_t valoffset, std::map<std::pair<int, int>, U>& data) {
  typename std::map<std::pair<int, int>, U>::iterator nuc_iter, nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent, 0));
  nuc_end = data.upper_bound(std::make_pair(parent, 9999999999));
  T* ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it != nuc_end; ++it) {
    ret = (T*)((char*) & (it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<T, U>(parent, valoffset, data);
  }
  return result;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
                                                                  size_t valoffset, std::map<std::pair<int, unsigned int>, U>& data) {
  typename std::map<std::pair<int, unsigned int>, U>::iterator nuc_iter,
           nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent, 0));
  nuc_end = data.upper_bound(std::make_pair(parent, UINT_MAX));
  T* ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it != nuc_end; ++it) {
    ret = (T*)((char*) & (it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<T, U>(parent, valoffset, data);
  }
  return result;
}

template<typename U> double pyne::data_access(int nuc,
                                              size_t valoffset, std::map<int, U>& data) {
  typename std::map<int, U>::iterator nuc_iter,
           nuc_end;
  nuc_iter = data.find(nuc);
  nuc_end = data.end();
  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end) {
    return *(double*)((char*) & (nuc_iter->second) + valoffset);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty()) {
    _load_data<U>();
    return data_access<U>(nuc, valoffset, data);
  }
  throw pyne::nucname::NotANuclide(nuc, "");
}


//
// Load atomic data
//

std::map<int, pyne::atomic> pyne::atomic_data_map;

template<> void pyne::_load_data<pyne::atomic>() {
  // Loads the atomic table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic));
  status = H5Tinsert(desc, "z", HOFFSET(atomic, z),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "k_shell_fluor", HOFFSET(atomic, k_shell_fluor),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_fluor_error", HOFFSET(atomic, k_shell_fluor_error),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_shell_fluor", HOFFSET(atomic, l_shell_fluor),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_shell_fluor_error", HOFFSET(atomic, l_shell_fluor_error),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "prob", HOFFSET(atomic, prob),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_be", HOFFSET(atomic, k_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_be_err", HOFFSET(atomic, k_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "li_shell_be", HOFFSET(atomic, li_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "li_shell_be_err", HOFFSET(atomic, li_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "mi_shell_be", HOFFSET(atomic, mi_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "mi_shell_be_err", HOFFSET(atomic, mi_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ni_shell_be", HOFFSET(atomic, ni_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ni_shell_be_err", HOFFSET(atomic, ni_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_to_ka", HOFFSET(atomic, kb_to_ka),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_to_ka_err", HOFFSET(atomic, kb_to_ka_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_to_ka1", HOFFSET(atomic, ka2_to_ka1),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_to_ka1_err", HOFFSET(atomic, ka2_to_ka1_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_auger", HOFFSET(atomic, l_auger),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_auger", HOFFSET(atomic, k_auger),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka1_x_ray_en", HOFFSET(atomic, ka1_x_ray_en),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka1_x_ray_en_err", HOFFSET(atomic, ka1_x_ray_en_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_x_ray_en", HOFFSET(atomic, ka2_x_ray_en),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_x_ray_en_err", HOFFSET(atomic, ka2_x_ray_en_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_x_ray_en", HOFFSET(atomic, kb_x_ray_en),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_x_ray_en", HOFFSET(atomic, l_x_ray_en),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);
  // Open the data set
  hid_t atomic_set = H5Dopen2(nuc_data_h5, "/decay/atomic", H5P_DEFAULT);
  hid_t atomic_space = H5Dget_space(atomic_set);
  int atomic_length = H5Sget_simple_extent_npoints(atomic_space);

  // Read in the data
  atomic* atomic_array = new atomic[atomic_length];
  status = H5Dread(atomic_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   atomic_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(atomic_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < atomic_length; ++i) {
    atomic_data_map[atomic_array[i].z] = atomic_array[i];
  }

  delete[] atomic_array;

}

std::vector<std::pair<double, double> >
pyne::calculate_xray_data(int z, double k_conv, double l_conv) {
  double xk = 0;
  double xka = 0;
  double xka1 = 0;
  double xka2 = 0;
  double xkb = 0;
  double xl = 0;
  if (!isnan(k_conv)) {
    xk = data_access<atomic> (z, offsetof(atomic, k_shell_fluor),
                              atomic_data_map) * k_conv;
    xka = xk / (1.0 + data_access<atomic> (z, offsetof(atomic,
                                                       kb_to_ka), atomic_data_map));
    xka1 = xka / (1.0 + data_access<atomic> (z, offsetof(atomic,
                                                         ka2_to_ka1), atomic_data_map));
    xka2 = xka - xka1;
    xkb = xk - xka;
    if (!isnan(l_conv)) {
      xl = (l_conv + k_conv * data_access<atomic> (z, offsetof(atomic,
                                                               prob), atomic_data_map)) * data_access<atomic> (z, offsetof(atomic,
                                                                   l_shell_fluor), atomic_data_map);
    }
  } else if (!isnan(l_conv)) {
    xl = l_conv * data_access<atomic> (z, offsetof(atomic,
                                                   l_shell_fluor), atomic_data_map);
  }
  std::vector<std::pair<double, double> > result;
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
                                                                   ka1_x_ray_en), atomic_data_map), xka1));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
                                                                   ka2_x_ray_en), atomic_data_map), xka2));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
                                                                   kb_x_ray_en), atomic_data_map), xkb));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
                                                                   l_x_ray_en), atomic_data_map), xl));

  return result;
}


//
// Load level data
//

std::map<std::pair<int, double>, pyne::level_data> pyne::level_data_lvl_map;
std::map<std::pair<int, unsigned int>,
    pyne::level_data> pyne::level_data_rx_map;


template<> void pyne::_load_data<pyne::level_data>() {

  // Loads the level table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(level_data));
  status = H5Tinsert(desc, "nuc_id", HOFFSET(level_data, nuc_id),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "rx_id", HOFFSET(level_data, rx_id),
                     H5T_NATIVE_UINT);
  status = H5Tinsert(desc, "half_life", HOFFSET(level_data, half_life),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "level", HOFFSET(level_data, level),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(level_data, branch_ratio),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "metastable", HOFFSET(level_data, metastable),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "special", HOFFSET(level_data, special),
                     H5T_C_S1);
  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);
  // Open the data set
  hid_t level_set = H5Dopen2(nuc_data_h5, "/decay/level_list", H5P_DEFAULT);
  hid_t level_space = H5Dget_space(level_set);
  int level_length = H5Sget_simple_extent_npoints(level_space);

  // Read in the data
  level_data* level_array = new level_data[level_length];
  status = H5Dread(level_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   level_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(level_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < level_length; ++i) {
    if (level_array[i].rx_id == 0)
      level_data_lvl_map[std::make_pair(level_array[i].nuc_id,
                                        level_array[i].level)] = level_array[i];
    else
      level_data_rx_map[std::make_pair(level_array[i].nuc_id,
                                       level_array[i].rx_id)] = level_array[i];
  }

  delete[] level_array;
}

//
// level id
//
int pyne::id_from_level(int nuc, double level, std::string special) {
  int nostate = (nuc / 10000) * 10000;
  if (level_data_lvl_map.empty()) {
    _load_data<level_data>();
  }

  std::map<std::pair<int, double>, level_data>::iterator nuc_lower, nuc_upper;

  nuc_lower = level_data_lvl_map.lower_bound(std::make_pair(nostate, 0.0));
  nuc_upper = level_data_lvl_map.upper_bound(std::make_pair(nostate + 9999,
                                                            DBL_MAX));
  double minv = DBL_MAX;
  //by default return input nuc_id with level stripped
  int ret_id = nuc;
  for (std::map<std::pair<int, double>, level_data>::iterator it = nuc_lower;
       it != nuc_upper; ++it) {
    if ((std::abs(level - it->second.level) < minv) &&
        ((char)it->second.special == special.c_str()[0]) &&
        !isnan(it->second.level)) {
      minv = std::abs(level - it->second.level);
      ret_id = it->second.nuc_id;
    }
  }
  // This value was chosen so important transitions in U-235 are not missed
  if (minv > 3.0)
    return -nuc;
  else
    return ret_id;
}

int pyne::id_from_level(int nuc, double level) {
  return id_from_level(nuc, level, " ");
}
//
// Metastable id data
//

int pyne::metastable_id(int nuc, int m) {
  int nostate = (nuc / 10000) * 10000;
  if (m == 0)
    return nostate;
  if (level_data_lvl_map.empty()) {
    _load_data<level_data>();
  }

  std::map<std::pair<int, double>, level_data>::iterator nuc_lower, nuc_upper;

  nuc_lower = level_data_lvl_map.lower_bound(std::make_pair(nostate, 0.0));
  nuc_upper = level_data_lvl_map.upper_bound(std::make_pair(nostate + 9999,
                                                            DBL_MAX));
  for (std::map<std::pair<int, double>, level_data>::iterator it = nuc_lower;
       it != nuc_upper; ++it) {
    if (it->second.metastable == m)
      return it->second.nuc_id;
  }

  return -1;
}

int pyne::metastable_id(int nuc) {
  return metastable_id(nuc, 1);
}

//
// Decay children data
//


std::set<int> pyne::decay_children(int nuc) {
  // make sure spontaneous fission data is loaded
  if (wimsdfpy_data.empty())
    _load_wimsdfpy();

  std::vector<unsigned int> part = data_access<unsigned int, level_data>(nuc,
                                                                         offsetof(level_data, rx_id), level_data_rx_map);
  std::vector<unsigned int>::iterator it = part.begin();
  std::set<int> result;
  for (; it != part.end(); ++it) {
    switch (*it) {
      case 36125: {
        // internal transition, rx == 'it'
        result.insert((nuc / 10000) * 10000);
        break;
      }
      case 36565:
      case 1794828612: {
        // spontaneous fission, rx == 'sf', 36565
        // beta- & spontaneous fission, rx == 'b-sf', 1794828612
        std::map<std::pair<int, int>, double>::iterator sf = wimsdfpy_data.begin();
        for (; sf != wimsdfpy_data.end(); ++sf)
          if (sf->first.first == nuc)
            result.insert(sf->first.second);
        break;
      }
      default: {
        result.insert((rxname::child(nuc, *it, "decay") / 10000) * 10000);
      }
    }
  }
  return result;
}

std::set<int> pyne::decay_children(char* nuc) {
  return decay_children(nucname::id(nuc));
}

std::set<int> pyne::decay_children(std::string nuc) {
  return decay_children(nucname::id(nuc));
}

//
// Excitation state energy data
//

double pyne::state_energy(int nuc) {
  std::vector<double> result = data_access<double, level_data>(nuc, 0.0,
                                                               DBL_MAX, offsetof(level_data, level), level_data_lvl_map);
  if (result.size() == 1)
    return result[0] / 1000.0;
  return 0.0;
}

double pyne::state_energy(char* nuc) {
  return state_energy(nucname::id(nuc));
}


double pyne::state_energy(std::string nuc) {
  return state_energy(nucname::id(nuc));
}


//
// Decay constant data
//

double pyne::decay_const(int nuc) {
  std::vector<double> result = data_access<double, level_data>(nuc, 0.0,
                                                               DBL_MAX, offsetof(level_data, half_life), level_data_lvl_map);
  if (result.size() == 1) {
    return log(2.0) / result[0];
  }
  return 0.0;
}


double pyne::decay_const(char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
}


double pyne::decay_const(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
}


//
// Half-life data
//
double pyne::half_life(int nuc) {
  std::vector<double> result = data_access<double, level_data>(nuc, 0.0,
                                                               DBL_MAX, offsetof(level_data, half_life), level_data_lvl_map);
  if (result.size() == 1) {
    return result[0];
  }
  return 1.0 / 0.0;
}


double pyne::half_life(char* nuc) {
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
}

double pyne::half_life(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
}


//
// Branch ratio data
//
double pyne::branch_ratio(std::pair<int, int> from_to) {
  using std::vector;
  using pyne::nucname::groundstate;
  // make sure spontaneous fission data is loaded
  if (wimsdfpy_data.empty())
    _load_wimsdfpy();

  vector<unsigned int> part1 = \
                               data_access<unsigned int, level_data>(from_to.first, offsetof(level_data, rx_id),
                                                                     level_data_rx_map);
  vector<double> part2 = \
                         data_access<double, level_data>(from_to.first, offsetof(level_data, branch_ratio),
                                                         level_data_rx_map);
  double result = 0.0;
  if ((from_to.first == from_to.second) && (decay_const(from_to.first) == 0.0))
    return 1.0;
  for (vector<unsigned int>::size_type i = 0; i < part1.size(); ++i) {
    if ((part1[i] == 36125) &&
        (groundstate(from_to.first) == groundstate(from_to.second)) &&
        (from_to.second % 10000 == 0)) {
      // internal transition, rx == 'it', 36125
      result = part2[i] * 0.01;
      break;
    } else if (part1[i] == 36565 || part1[i] == 1794828612) {
      // spontaneous fission, rx == 'sf', 36565
      // beta- & spontaneous fission, rx == 'b-sf', 1794828612
      result += part2[i] * 0.01 * wimsdfpy_data[from_to];
    } else if ((part1[i] != 0) && (groundstate(rxname::child(from_to.first,
                                                             part1[i], "decay")) == from_to.second)) {
      result += part2[i] * 0.01;
    }
  }
  return result;
}

double pyne::branch_ratio(int from_nuc, int to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

double pyne::branch_ratio(char* from_nuc, char* to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

double pyne::branch_ratio(std::string from_nuc, std::string to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

std::map<std::pair<int, int>, pyne::decay> pyne::decay_data = \
                                                              std::map<std::pair<int, int>, pyne::decay>();

template<> void pyne::_load_data<pyne::decay>() {

  // Loads the decay table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(decay));
  status = H5Tinsert(desc, "parent", HOFFSET(decay, parent),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "child", HOFFSET(decay, child),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "decay", HOFFSET(decay, decay),
                     H5T_NATIVE_UINT);
  status = H5Tinsert(desc, "half_life", HOFFSET(decay, half_life),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "half_life_error", HOFFSET(decay,
                                                      half_life_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(decay, branch_ratio),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio_error", HOFFSET(decay, branch_ratio_error),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio", HOFFSET(decay,
                                                          photon_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio_err", HOFFSET(decay,
                                                              photon_branch_ratio_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio", HOFFSET(decay,
                                                        beta_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio_err", HOFFSET(decay,
                                                            beta_branch_ratio_error), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t decay_set = H5Dopen2(nuc_data_h5, "/decay/decays", H5P_DEFAULT);
  hid_t decay_space = H5Dget_space(decay_set);
  int decay_length = H5Sget_simple_extent_npoints(decay_space);

  // Read in the data
  decay* decay_array = new decay[decay_length];
  status = H5Dread(decay_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   decay_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(decay_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < decay_length; ++i) {
    decay_data[std::make_pair(decay_array[i].parent, decay_array[i].child)] = \
                                                                              decay_array[i];
  }
  delete[] decay_array;
}


std::vector<int> pyne::decay_data_children(int parent) {
  std::vector<int> result = data_access<int, decay>(parent,
                                                    offsetof(decay, child), decay_data);
  return result;
}

std::pair<double, double> pyne::decay_half_life(std::pair<int, int> from_to) {
  return std::make_pair(data_access<double, decay>(from_to, offsetof(
                                                       decay, half_life), decay_data), data_access<double, decay>(
                            from_to, offsetof(decay, half_life_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_half_lifes(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
                                                         offsetof(decay, half_life), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
                                                         offsetof(decay, half_life_error), decay_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::pair<double, double> pyne::decay_branch_ratio(std::pair<int, int> from_to) {
  return std::make_pair(data_access<double, decay>(from_to, offsetof(decay,
                                                                     branch_ratio), decay_data), data_access<double, decay>(from_to, offsetof(decay,
                                                                         branch_ratio_error), decay_data));
}

std::vector<double> pyne::decay_branch_ratios(int parent) {
  return data_access<double, decay>(parent, offsetof(decay,
                                                     branch_ratio), decay_data);
}

std::pair<double, double> pyne::decay_photon_branch_ratio(std::pair<int, int>
                                                          from_to) {
  return std::make_pair(data_access<double, decay>(from_to,
                                                   offsetof(decay, photon_branch_ratio), decay_data),
                        data_access<double, decay>(from_to, offsetof(decay,
                                                                     photon_branch_ratio_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_photon_branch_ratios(
    int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
                                                         offsetof(decay, photon_branch_ratio), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
                                                         offsetof(decay, photon_branch_ratio_error), decay_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::pair<double, double> pyne::decay_beta_branch_ratio(std::pair<int, int>
                                                        from_to) {
  return std::make_pair(data_access<double, decay>(from_to,
                                                   offsetof(decay, beta_branch_ratio), decay_data),
                        data_access<double, decay>(from_to, offsetof(decay,
                                                                     beta_branch_ratio_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_beta_branch_ratios(
    int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
                                                         offsetof(decay, beta_branch_ratio), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
                                                         offsetof(decay, beta_branch_ratio_error), decay_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::map<std::pair<int, double>, pyne::gamma> pyne::gamma_data;

template<> void pyne::_load_data<pyne::gamma>() {

  // Loads the gamma table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(gamma));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(gamma, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(gamma, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "parent_nuc", HOFFSET(gamma, parent_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "child_nuc", HOFFSET(gamma, child_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "energy", HOFFSET(gamma, energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "energy_err", HOFFSET(gamma, energy_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity", HOFFSET(gamma,
                                                       photon_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity_err", HOFFSET(gamma,
                                                           photon_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity", HOFFSET(gamma,
                                                     conv_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity_err", HOFFSET(gamma,
                                                         conv_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity", HOFFSET(gamma,
                                                      total_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity_err", HOFFSET(gamma,
                                                          total_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(gamma, k_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(gamma, l_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(gamma, m_conv_e),
                     H5T_NATIVE_DOUBLE);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t gamma_set = H5Dopen2(nuc_data_h5, "/decay/gammas", H5P_DEFAULT);
  hid_t gamma_space = H5Dget_space(gamma_set);
  int gamma_length = H5Sget_simple_extent_npoints(gamma_space);

  // Read in the data
  gamma* gamma_array = new gamma[gamma_length];
  status = H5Dread(gamma_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   gamma_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(gamma_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < gamma_length; ++i) {
    if ((gamma_array[i].parent_nuc != 0) && !isnan(gamma_array[i].energy))
      gamma_data[std::make_pair(gamma_array[i].parent_nuc,
                                gamma_array[i].energy)] = gamma_array[i];
  }
  delete[] gamma_array;
}

std::vector<std::pair<double, double> > pyne::gamma_energy(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, energy), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, energy_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_energy(double energy,
                                                           double error) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(energy + error,
                                                         energy - error, offsetof(gamma, energy), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(energy + error,
                                                         energy - error, offsetof(gamma, energy_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_photon_intensity(
    int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, photon_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, photon_intensity_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_photon_intensity(
    double energy, double error) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(energy + error,
                                                         energy - error, offsetof(gamma, photon_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(energy + error,
                                                         energy - error, offsetof(gamma, photon_intensity_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_conversion_intensity(
    int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, conv_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, conv_intensity_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_total_intensity(
    int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, total_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
                                                         DBL_MAX, offsetof(gamma, total_intensity_err), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<int, int> > pyne::gamma_from_to(int parent) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(parent, 0.0, DBL_MAX,
                                                   offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(parent, 0.0, DBL_MAX,
                                                   offsetof(gamma, to_nuc), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<std::pair<int, int> > pyne::gamma_from_to(double energy,
                                                      double error) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(energy + error,
                                                   energy - error, offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(energy + error,
                                                   energy - error, offsetof(gamma, to_nuc), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}


std::vector<std::pair<int, int> > pyne::gamma_parent_child(double energy,
                                                           double error) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(energy + error,
                                                   energy - error, offsetof(gamma, parent_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(energy + error,
                                                   energy - error, offsetof(gamma, child_nuc), gamma_data);
  for (int i = 0; i < part1.size(); ++i) {
    result.push_back(std::make_pair(part1[i], part2[i]));
  }
  return result;
}

std::vector<int> pyne::gamma_parent(double energy, double error) {
  return data_access<int, gamma>(energy + error, energy - error,
                                 offsetof(gamma, parent_nuc), gamma_data);
}

std::vector<int> pyne::gamma_child(double energy, double error) {
  return data_access<int, gamma>(energy + error, energy - error,
                                 offsetof(gamma, child_nuc), gamma_data);
}

std::vector<int> pyne::gamma_child(int parent) {
  return data_access<int, gamma>(parent, 0.0, DBL_MAX,
                                 offsetof(gamma, child_nuc), gamma_data);
}

std::vector<std::pair<double, double> > pyne::gamma_xrays(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
                                                          offsetof(gamma, k_conv_e), gamma_data);
  std::vector<double> l_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
                                                          offsetof(gamma, l_conv_e), gamma_data);
  std::vector<int> children = data_access<int, gamma>(parent, 0.0,
                                                      DBL_MAX, offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
      decay_photon_branch_ratios(parent);
  for (int i = 0; i < k_list.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
                                   k_list[i] * decay_br[j].first, l_list[i] * decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  return result;
}

std::map<std::pair<int, double>, pyne::alpha> pyne::alpha_data;

template<> void pyne::_load_data<pyne::alpha>() {

  // Loads the alpha table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(alpha));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(alpha, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(alpha, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "energy", HOFFSET(alpha, energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(alpha, intensity),
                     H5T_NATIVE_DOUBLE);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t alpha_set = H5Dopen2(nuc_data_h5, "/decay/alphas", H5P_DEFAULT);
  hid_t alpha_space = H5Dget_space(alpha_set);
  int alpha_length = H5Sget_simple_extent_npoints(alpha_space);

  // Read in the data
  alpha* alpha_array = new alpha[alpha_length];
  status = H5Dread(alpha_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   alpha_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(alpha_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < alpha_length; ++i) {
    if ((alpha_array[i].from_nuc != 0) && !isnan(alpha_array[i].energy))
      alpha_data[std::make_pair(alpha_array[i].from_nuc, alpha_array[i].energy)]
        = alpha_array[i];
  }
  delete[] alpha_array;
}

std::vector<double > pyne::alpha_energy(int parent) {
  return data_access<double, alpha>(parent, 0.0, DBL_MAX,
                                    offsetof(alpha, energy), alpha_data);
}
std::vector<double> pyne::alpha_intensity(int parent) {
  return data_access<double, alpha>(parent, 0.0, DBL_MAX,
                                    offsetof(alpha, intensity), alpha_data);
}

std::vector<int> pyne::alpha_parent(double energy, double error) {
  return data_access<int, alpha>(energy + error, energy - error,
                                 offsetof(alpha, from_nuc), alpha_data);
}

std::vector<int> pyne::alpha_child(double energy, double error) {
  return data_access<int, alpha>(energy + error, energy - error,
                                 offsetof(alpha, to_nuc), alpha_data);
}

std::vector<int> pyne::alpha_child(int parent) {
  return data_access<int, alpha>(parent, 0.0, DBL_MAX,
                                 offsetof(alpha, to_nuc), alpha_data);
}

std::map<std::pair<int, double>, pyne::beta> pyne::beta_data;

template<> void pyne::_load_data<pyne::beta>() {

  // Loads the beta table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(beta));
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(beta,
                                                      endpoint_energy), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(beta, avg_energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(beta, intensity),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(beta, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(beta, to_nuc),
                     H5T_NATIVE_INT);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t beta_set = H5Dopen2(nuc_data_h5, "/decay/betas", H5P_DEFAULT);
  hid_t beta_space = H5Dget_space(beta_set);
  int beta_length = H5Sget_simple_extent_npoints(beta_space);

  // Read in the data
  beta* beta_array = new beta[beta_length];
  status = H5Dread(beta_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, beta_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(beta_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < beta_length; ++i) {
    if ((beta_array[i].from_nuc != 0) && !isnan(beta_array[i].avg_energy))
      beta_data[std::make_pair(beta_array[i].from_nuc, beta_array[i].avg_energy)]
        = beta_array[i];
  }
  delete[] beta_array;
}

std::vector<double > pyne::beta_endpoint_energy(int parent) {
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                                   offsetof(beta, endpoint_energy), beta_data);
}

std::vector<double > pyne::beta_average_energy(int parent) {
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                                   offsetof(beta, avg_energy), beta_data);
}

std::vector<double> pyne::beta_intensity(int parent) {
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                                   offsetof(beta, intensity), beta_data);
}

std::vector<int> pyne::beta_parent(double energy, double error) {
  return data_access<int, beta>(energy + error, energy - error,
                                offsetof(beta, from_nuc), beta_data);
}

std::vector<int> pyne::beta_child(double energy, double error) {
  return data_access<int, beta>(energy + error, energy - error,
                                offsetof(beta, to_nuc), beta_data);
}

std::vector<int> pyne::beta_child(int parent) {
  return data_access<int, beta>(parent, 0.0, DBL_MAX,
                                offsetof(beta, to_nuc), beta_data);
}


std::map<std::pair<int, double>, pyne::ecbp> pyne::ecbp_data;

template<> void pyne::_load_data<pyne::ecbp>() {

  // Loads the ecbp table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ecbp));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ecbp, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ecbp, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(ecbp,
                                                      endpoint_energy), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(ecbp, avg_energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_plus_intensity", HOFFSET(ecbp,
                                                          beta_plus_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ec_intensity", HOFFSET(ecbp, ec_intensity),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(ecbp, k_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(ecbp, l_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(ecbp, m_conv_e),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t ecbp_set = H5Dopen2(nuc_data_h5, "/decay/ecbp", H5P_DEFAULT);
  hid_t ecbp_space = H5Dget_space(ecbp_set);
  int ecbp_length = H5Sget_simple_extent_npoints(ecbp_space);

  // Read in the data
  ecbp* ecbp_array = new ecbp[ecbp_length];
  status = H5Dread(ecbp_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecbp_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(ecbp_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < ecbp_length; ++i) {
    if ((ecbp_array[i].from_nuc != 0) && !isnan(ecbp_array[i].avg_energy))
      ecbp_data[std::make_pair(ecbp_array[i].from_nuc, ecbp_array[i].avg_energy)]
        = ecbp_array[i];
  }
  delete[] ecbp_array;
}

std::vector<double > pyne::ecbp_endpoint_energy(int parent) {
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                   offsetof(ecbp, endpoint_energy), ecbp_data);
}

std::vector<double > pyne::ecbp_average_energy(int parent) {
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                   offsetof(ecbp, avg_energy), ecbp_data);
}

std::vector<double> pyne::ec_intensity(int parent) {
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                   offsetof(ecbp, ec_intensity), ecbp_data);
}

std::vector<double> pyne::bp_intensity(int parent) {
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                   offsetof(ecbp, beta_plus_intensity), ecbp_data);
}

std::vector<int> pyne::ecbp_parent(double energy, double error) {
  return data_access<int, ecbp>(energy + error, energy - error,
                                offsetof(ecbp, from_nuc), ecbp_data);
}

std::vector<int> pyne::ecbp_child(double energy, double error) {
  return data_access<int, ecbp>(energy + error, energy - error,
                                offsetof(ecbp, to_nuc), ecbp_data);
}

std::vector<int> pyne::ecbp_child(int parent) {
  return data_access<int, ecbp>(parent, 0.0, DBL_MAX,
                                offsetof(ecbp, to_nuc), ecbp_data);
}

std::vector<std::pair<double, double> > pyne::ecbp_xrays(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                                         offsetof(ecbp, k_conv_e), ecbp_data);
  std::vector<double> l_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                                         offsetof(ecbp, l_conv_e), ecbp_data);
  std::vector<int> children = data_access<int, ecbp>(parent, 0.0, DBL_MAX,
                                                     offsetof(ecbp, to_nuc), ecbp_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
      decay_beta_branch_ratios(parent);
  for (int i = 0; i < k_list.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
                                   k_list[i] * decay_br[j].first, l_list[i] * decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  return result;
}

///////////////////////////
/// Combined decay data ///
///////////////////////////

std::vector<std::pair<double, double> > pyne::gammas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<std::pair<double, double> > energies = gamma_energy(parent_state_id);
  std::vector<std::pair<double, double> > intensities =
      gamma_photon_intensity(parent_state_id);
  std::vector<std::pair<int, int> > children = gamma_from_to(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<std::pair<double, double> > decay_br =
      decay_photon_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i].first) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i].first,
                                        decay_c * intensities[i].first * decay_br[j].first));
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::alphas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<double> energies = alpha_energy(parent_state_id);
  std::vector<double> intensities = alpha_intensity(parent_state_id);
  std::vector<int> children = alpha_child(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<double> decay_br = decay_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i],
                                        decay_c * decay_br[j]*intensities[i]));
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::betas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<double> energies = beta_average_energy(parent_state_id);
  std::vector<double> intensities = beta_intensity(parent_state_id);
  std::vector<int> children = beta_child(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<std::pair<double, double> > decay_br =
      decay_beta_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i],
                                        decay_c * decay_br[j].first * intensities[i]));
        break;
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::xrays(int parent) {
  double decay_c = decay_const(parent);
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                                         offsetof(ecbp, k_conv_e), ecbp_data);
  std::vector<double> l_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                                                         offsetof(ecbp, l_conv_e), ecbp_data);
  std::vector<int> children = data_access<int, ecbp>(parent, 0.0, DBL_MAX,
                                                     offsetof(ecbp, to_nuc), ecbp_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
      decay_beta_branch_ratios(parent);
  for (int i = 0; i < k_list.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
                                   k_list[i] * decay_br[j].first, l_list[i] * decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  std::vector<double> gk_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
                                                           offsetof(gamma, k_conv_e), gamma_data);
  std::vector<double> gl_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
                                                           offsetof(gamma, l_conv_e), gamma_data);
  std::vector<int> gchildren = data_access<int, gamma>(parent, 0.0,
                                                       DBL_MAX, offsetof(gamma, from_nuc), gamma_data);
  std::vector<std::pair<double, double> > decay_nrbr =
      decay_photon_branch_ratios(parent);
  for (int i = 0; i < gk_list.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(gchildren[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(gchildren[i]),
                                   gk_list[i] * decay_nrbr[j].first, gl_list[i] * decay_nrbr[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }

  for (int i = 0; i < result.size(); ++i)
    result[i].second = result[i].second * decay_c;
  return result;
}

//////////////////////////////////////////
//////////// simple xs data //////////////
//////////////////////////////////////////

typedef struct simple_xs {
  int nuc;
  double sigma_t;
  double sigma_s;
  double sigma_e;
  double sigma_i;
  double sigma_a;
  double sigma_gamma;
  double sigma_f;
  double sigma_alpha;
  double sigma_proton;
  double sigma_deut;
  double sigma_trit;
  double sigma_2n;
  double sigma_3n;
  double sigma_4n;
} simple_xs;

std::map<std::string, std::map<int, std::map<int, double> > > pyne::simple_xs_map;

// loads the simple cross section data for the specified energy band from
// the nuc_data.h5 file into memory.
static void _load_simple_xs_map(std::string energy) {
  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  using pyne::rxname::id;
  std::map<unsigned int, size_t> rxns;
  rxns[id("tot")] = offsetof(simple_xs, sigma_t);
  rxns[id("scat")] = offsetof(simple_xs, sigma_s);
  rxns[id("elas")] = offsetof(simple_xs, sigma_e);
  rxns[id("inel")] = offsetof(simple_xs, sigma_i);
  rxns[id("abs")] = offsetof(simple_xs, sigma_a);
  rxns[id("gamma")] = offsetof(simple_xs, sigma_gamma);
  rxns[id("fiss")] = offsetof(simple_xs, sigma_f);
  rxns[id("alpha")] = offsetof(simple_xs, sigma_alpha);
  rxns[id("proton")] = offsetof(simple_xs, sigma_proton);
  rxns[id("deut")] = offsetof(simple_xs, sigma_deut);
  rxns[id("trit")] = offsetof(simple_xs, sigma_trit);
  rxns[id("z_2n")] = offsetof(simple_xs, sigma_2n);
  rxns[id("z_3n")] = offsetof(simple_xs, sigma_3n);
  rxns[id("z_4n")] = offsetof(simple_xs, sigma_4n);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(simple_xs));
  H5Tinsert(desc, "nuc",   HOFFSET(simple_xs, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "sigma_t",  HOFFSET(simple_xs, sigma_t),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_s", HOFFSET(simple_xs, sigma_s), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_e", HOFFSET(simple_xs, sigma_e), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_i", HOFFSET(simple_xs, sigma_i), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_a", HOFFSET(simple_xs, sigma_a), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_gamma", HOFFSET(simple_xs, sigma_gamma), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_f", HOFFSET(simple_xs, sigma_f), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_alpha", HOFFSET(simple_xs, sigma_alpha), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_proton", HOFFSET(simple_xs, sigma_proton), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_deut", HOFFSET(simple_xs, sigma_deut), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_trit", HOFFSET(simple_xs, sigma_trit), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_2n", HOFFSET(simple_xs, sigma_2n), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_3n", HOFFSET(simple_xs, sigma_3n), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_4n", HOFFSET(simple_xs, sigma_4n), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // build path to prober simple xs table
  std::string path = "/neutron/simple_xs/" + energy;

  // Open the data set
  hid_t simple_xs_set = H5Dopen2(nuc_data_h5, path.c_str(), H5P_DEFAULT);
  hid_t simple_xs_space = H5Dget_space(simple_xs_set);
  int n = H5Sget_simple_extent_npoints(simple_xs_space);

  // Read in the data
  simple_xs* array = new simple_xs[n];
  H5Dread(simple_xs_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(simple_xs_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the map
  for (int i = 0; i < n; i++) {
    std::map<unsigned int, size_t>::iterator it;
    for (it = rxns.begin(); it != rxns.end(); ++it) {
      double xs = *(double*)((char*)&array[i] + it->second);
      pyne::simple_xs_map[energy][array[i].nuc][it->first] = xs;
    }
  }
  delete[] array;
}

double pyne::simple_xs(int nuc, int rx_id, std::string energy) {
  std::set<std::string> energies;
  energies.insert("thermal");
  energies.insert("thermal_maxwell_ave");
  energies.insert("resonance_integral");
  energies.insert("fourteen_MeV");
  energies.insert("fission_spectrum_ave");

  if (energies.count(energy) == 0) {
    throw InvalidSimpleXS("Energy '" + energy +
                          "' is not a valid simple_xs group");
  } else if (simple_xs_map.count(energy) == 0) {
    _load_simple_xs_map(energy);
  }

  if (simple_xs_map[energy].count(nuc) == 0) {
    throw InvalidSimpleXS(rxname::name(rx_id) +
                          " is not a valid simple_xs nuclide");
  } else if (simple_xs_map[energy][nuc].count(rx_id) == 0) {
    throw InvalidSimpleXS(rxname::name(rx_id) +
                          " is not a valid simple_xs reaction");
  }

  return simple_xs_map[energy][nuc][rx_id];
}

double pyne::simple_xs(int nuc, std::string rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
double pyne::simple_xs(std::string nuc, int rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
double pyne::simple_xs(std::string nuc, std::string rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
//
// end of src/data.cpp
//


//
// start of src/jsoncpp.cpp
//
/// Json-cpp amalgated source (https://github.com/open-source-parsers/jsoncpp/).
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
#else
#include "json.h"
#endif


// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_tool.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

#ifndef LIB_JSONCPP_JSON_TOOL_H_INCLUDED
# define LIB_JSONCPP_JSON_TOOL_H_INCLUDED

/* This header provides common string manipulation support, such as UTF-8,
 * portable conversion from/to string...
 *
 * It is an internal header that must not be exposed.
 */

namespace Json {

/// Converts a unicode code-point to UTF-8.
static inline std::string
codePointToUTF8(unsigned int cp) {
  std::string result;

  // based on description from http://en.wikipedia.org/wiki/UTF-8

  if (cp <= 0x7f) {
    result.resize(1);
    result[0] = static_cast<char>(cp);
  } else if (cp <= 0x7FF) {
    result.resize(2);
    result[1] = static_cast<char>(0x80 | (0x3f & cp));
    result[0] = static_cast<char>(0xC0 | (0x1f & (cp >> 6)));
  } else if (cp <= 0xFFFF) {
    result.resize(3);
    result[2] = static_cast<char>(0x80 | (0x3f & cp));
    result[1] = 0x80 | static_cast<char>((0x3f & (cp >> 6)));
    result[0] = 0xE0 | static_cast<char>((0xf & (cp >> 12)));
  } else if (cp <= 0x10FFFF) {
    result.resize(4);
    result[3] = static_cast<char>(0x80 | (0x3f & cp));
    result[2] = static_cast<char>(0x80 | (0x3f & (cp >> 6)));
    result[1] = static_cast<char>(0x80 | (0x3f & (cp >> 12)));
    result[0] = static_cast<char>(0xF0 | (0x7 & (cp >> 18)));
  }

  return result;
}


/// Returns true if ch is a control character (in range [0,32[).
static inline bool
isControlCharacter(char ch) {
  return ch > 0 && ch <= 0x1F;
}


enum {
  /// Constant that specify the size of the buffer that must be passed to uintToString.
  uintToStringBufferSize = 3 * sizeof(LargestUInt) + 1
};

// Defines a char buffer for use with uintToString().
typedef char UIntToStringBuffer[uintToStringBufferSize];


/** Converts an unsigned integer to string.
 * @param value Unsigned interger to convert to string
 * @param current Input/Output string buffer.
 *        Must have at least uintToStringBufferSize chars free.
 */
static inline void
uintToString(LargestUInt value,
             char*& current) {
  *--current = 0;
  do {
    *--current = char(value % 10) + '0';
    value /= 10;
  } while (value != 0);
}

} // namespace Json {

#endif // LIB_JSONCPP_JSON_TOOL_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_tool.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_reader.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/reader.h>
# include <json/value.h>
# include "json_tool.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <utility>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>

#if _MSC_VER >= 1400 // VC++ 8.0
#pragma warning( disable : 4996 )   // disable warning about strdup being deprecated.
#endif

namespace Json {

// Implementation of class Features
// ////////////////////////////////

Features::Features()
  : allowComments_(true)
  , strictRoot_(false) {
}


Features
Features::all() {
  return Features();
}


Features
Features::strictMode() {
  Features features;
  features.allowComments_ = false;
  features.strictRoot_ = true;
  return features;
}

// Implementation of class Reader
// ////////////////////////////////


static inline bool
in(Reader::Char c, Reader::Char c1, Reader::Char c2, Reader::Char c3, Reader::Char c4) {
  return c == c1  ||  c == c2  ||  c == c3  ||  c == c4;
}

static inline bool
in(Reader::Char c, Reader::Char c1, Reader::Char c2, Reader::Char c3, Reader::Char c4, Reader::Char c5) {
  return c == c1  ||  c == c2  ||  c == c3  ||  c == c4  ||  c == c5;
}


static bool
containsNewLine(Reader::Location begin,
                Reader::Location end) {
  for (; begin < end; ++begin)
    if (*begin == '\n'  ||  *begin == '\r')
      return true;
  return false;
}


// Class Reader
// //////////////////////////////////////////////////////////////////

Reader::Reader()
  : features_(Features::all()) {
}


Reader::Reader(const Features& features)
  : features_(features) {
}


bool
Reader::parse(const std::string& document,
              Value& root,
              bool collectComments) {
  document_ = document;
  const char* begin = document_.c_str();
  const char* end = begin + document_.length();
  return parse(begin, end, root, collectComments);
}


bool
Reader::parse(std::istream& sin,
              Value& root,
              bool collectComments) {
  //std::istream_iterator<char> begin(sin);
  //std::istream_iterator<char> end;
  // Those would allow streamed input from a file, if parse() were a
  // template function.

  // Since std::string is reference-counted, this at least does not
  // create an extra copy.
  std::string doc;
  std::getline(sin, doc, (char)EOF);
  return parse(doc, root, collectComments);
}

bool
Reader::parse(const char* beginDoc, const char* endDoc,
              Value& root,
              bool collectComments) {
  if (!features_.allowComments_) {
    collectComments = false;
  }

  begin_ = beginDoc;
  end_ = endDoc;
  collectComments_ = collectComments;
  current_ = begin_;
  lastValueEnd_ = 0;
  lastValue_ = 0;
  commentsBefore_ = "";
  errors_.clear();
  while (!nodes_.empty())
    nodes_.pop();
  nodes_.push(&root);

  bool successful = readValue();
  Token token;
  skipCommentTokens(token);
  if (collectComments_  &&  !commentsBefore_.empty())
    root.setComment(commentsBefore_, commentAfter);
  if (features_.strictRoot_) {
    if (!root.isArray()  &&  !root.isObject()) {
      // Set error location to start of doc, ideally should be first token found in doc
      token.type_ = tokenError;
      token.start_ = beginDoc;
      token.end_ = endDoc;
      addError("A valid JSON document must be either an array or an object value.",
               token);
      return false;
    }
  }
  return successful;
}


bool
Reader::readValue() {
  Token token;
  skipCommentTokens(token);
  bool successful = true;

  if (collectComments_  &&  !commentsBefore_.empty()) {
    currentValue().setComment(commentsBefore_, commentBefore);
    commentsBefore_ = "";
  }


  switch (token.type_) {
    case tokenObjectBegin:
      successful = readObject(token);
      break;
    case tokenArrayBegin:
      successful = readArray(token);
      break;
    case tokenNumber:
      successful = decodeNumber(token);
      break;
    case tokenString:
      successful = decodeString(token);
      break;
    case tokenTrue:
      currentValue() = true;
      break;
    case tokenFalse:
      currentValue() = false;
      break;
    case tokenNull:
      currentValue() = Value();
      break;
    default:
      return addError("Syntax error: value, object or array expected.", token);
  }

  if (collectComments_) {
    lastValueEnd_ = current_;
    lastValue_ = &currentValue();
  }

  return successful;
}


void
Reader::skipCommentTokens(Token& token) {
  if (features_.allowComments_) {
    do {
      readToken(token);
    } while (token.type_ == tokenComment);
  } else {
    readToken(token);
  }
}


bool
Reader::expectToken(TokenType type, Token& token, const char* message) {
  readToken(token);
  if (token.type_ != type)
    return addError(message, token);
  return true;
}


bool
Reader::readToken(Token& token) {
  skipSpaces();
  token.start_ = current_;
  Char c = getNextChar();
  bool ok = true;
  switch (c) {
    case '{':
      token.type_ = tokenObjectBegin;
      break;
    case '}':
      token.type_ = tokenObjectEnd;
      break;
    case '[':
      token.type_ = tokenArrayBegin;
      break;
    case ']':
      token.type_ = tokenArrayEnd;
      break;
    case '"':
      token.type_ = tokenString;
      ok = readString();
      break;
    case '/':
      token.type_ = tokenComment;
      ok = readComment();
      break;
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case '-':
      token.type_ = tokenNumber;
      readNumber();
      break;
    case 't':
      token.type_ = tokenTrue;
      ok = match("rue", 3);
      break;
    case 'f':
      token.type_ = tokenFalse;
      ok = match("alse", 4);
      break;
    case 'n':
      token.type_ = tokenNull;
      ok = match("ull", 3);
      break;
    case ',':
      token.type_ = tokenArraySeparator;
      break;
    case ':':
      token.type_ = tokenMemberSeparator;
      break;
    case 0:
      token.type_ = tokenEndOfStream;
      break;
    default:
      ok = false;
      break;
  }
  if (!ok)
    token.type_ = tokenError;
  token.end_ = current_;
  return true;
}


void
Reader::skipSpaces() {
  while (current_ != end_) {
    Char c = *current_;
    if (c == ' '  ||  c == '\t'  ||  c == '\r'  ||  c == '\n')
      ++current_;
    else
      break;
  }
}


bool
Reader::match(Location pattern,
              int patternLength) {
  if (end_ - current_ < patternLength)
    return false;
  int index = patternLength;
  while (index--)
    if (current_[index] != pattern[index])
      return false;
  current_ += patternLength;
  return true;
}


bool
Reader::readComment() {
  Location commentBegin = current_ - 1;
  Char c = getNextChar();
  bool successful = false;
  if (c == '*')
    successful = readCStyleComment();
  else if (c == '/')
    successful = readCppStyleComment();
  if (!successful)
    return false;

  if (collectComments_) {
    CommentPlacement placement = commentBefore;
    if (lastValueEnd_  &&  !containsNewLine(lastValueEnd_, commentBegin)) {
      if (c != '*'  ||  !containsNewLine(commentBegin, current_))
        placement = commentAfterOnSameLine;
    }

    addComment(commentBegin, current_, placement);
  }
  return true;
}


void
Reader::addComment(Location begin,
                   Location end,
                   CommentPlacement placement) {
  assert(collectComments_);
  if (placement == commentAfterOnSameLine) {
    assert(lastValue_ != 0);
    lastValue_->setComment(std::string(begin, end), placement);
  } else {
    if (!commentsBefore_.empty())
      commentsBefore_ += "\n";
    commentsBefore_ += std::string(begin, end);
  }
}


bool
Reader::readCStyleComment() {
  while (current_ != end_) {
    Char c = getNextChar();
    if (c == '*'  &&  *current_ == '/')
      break;
  }
  return getNextChar() == '/';
}


bool
Reader::readCppStyleComment() {
  while (current_ != end_) {
    Char c = getNextChar();
    if (c == '\r'  ||  c == '\n')
      break;
  }
  return true;
}


void
Reader::readNumber() {
  while (current_ != end_) {
    if (!(*current_ >= '0'  &&  *current_ <= '9')  &&
        !in(*current_, '.', 'e', 'E', '+', '-'))
      break;
    ++current_;
  }
}

bool
Reader::readString() {
  Char c = 0;
  while (current_ != end_) {
    c = getNextChar();
    if (c == '\\')
      getNextChar();
    else if (c == '"')
      break;
  }
  return c == '"';
}


bool
Reader::readObject(Token& /*tokenStart*/) {
  Token tokenName;
  std::string name;
  currentValue() = Value(objectValue);
  while (readToken(tokenName)) {
    bool initialTokenOk = true;
    while (tokenName.type_ == tokenComment  &&  initialTokenOk)
      initialTokenOk = readToken(tokenName);
    if (!initialTokenOk)
      break;
    if (tokenName.type_ == tokenObjectEnd  &&  name.empty())    // empty object
      return true;
    if (tokenName.type_ != tokenString)
      break;

    name = "";
    if (!decodeString(tokenName, name))
      return recoverFromError(tokenObjectEnd);

    Token colon;
    if (!readToken(colon) ||  colon.type_ != tokenMemberSeparator) {
      return addErrorAndRecover("Missing ':' after object member name",
                                colon,
                                tokenObjectEnd);
    }
    Value& value = currentValue()[ name ];
    nodes_.push(&value);
    bool ok = readValue();
    nodes_.pop();
    if (!ok)   // error already set
      return recoverFromError(tokenObjectEnd);

    Token comma;
    if (!readToken(comma)
        || (comma.type_ != tokenObjectEnd  &&
            comma.type_ != tokenArraySeparator &&
            comma.type_ != tokenComment)) {
      return addErrorAndRecover("Missing ',' or '}' in object declaration",
                                comma,
                                tokenObjectEnd);
    }
    bool finalizeTokenOk = true;
    while (comma.type_ == tokenComment &&
           finalizeTokenOk)
      finalizeTokenOk = readToken(comma);
    if (comma.type_ == tokenObjectEnd)
      return true;
  }
  return addErrorAndRecover("Missing '}' or object member name",
                            tokenName,
                            tokenObjectEnd);
}


bool
Reader::readArray(Token& /*tokenStart*/) {
  currentValue() = Value(arrayValue);
  skipSpaces();
  if (*current_ == ']') {   // empty array
    Token endArray;
    readToken(endArray);
    return true;
  }
  int index = 0;
  for (;;) {
    Value& value = currentValue()[ index++ ];
    nodes_.push(&value);
    bool ok = readValue();
    nodes_.pop();
    if (!ok)   // error already set
      return recoverFromError(tokenArrayEnd);

    Token token;
    // Accept Comment after last item in the array.
    ok = readToken(token);
    while (token.type_ == tokenComment  &&  ok) {
      ok = readToken(token);
    }
    bool badTokenType = (token.type_ != tokenArraySeparator  &&
                         token.type_ != tokenArrayEnd);
    if (!ok  ||  badTokenType) {
      return addErrorAndRecover("Missing ',' or ']' in array declaration",
                                token,
                                tokenArrayEnd);
    }
    if (token.type_ == tokenArrayEnd)
      break;
  }
  return true;
}


bool
Reader::decodeNumber(Token& token) {
  bool isDouble = false;
  for (Location inspect = token.start_; inspect != token.end_; ++inspect) {
    isDouble = isDouble
               ||  in(*inspect, '.', 'e', 'E', '+')
               || (*inspect == '-'  &&  inspect != token.start_);
  }
  if (isDouble)
    return decodeDouble(token);
  // Attempts to parse the number as an integer. If the number is
  // larger than the maximum supported value of an integer then
  // we decode the number as a double.
  Location current = token.start_;
  bool isNegative = *current == '-';
  if (isNegative)
    ++current;
  Value::LargestUInt maxIntegerValue = isNegative ? Value::LargestUInt(-Value::minLargestInt)
                                       : Value::maxLargestUInt;
  Value::LargestUInt threshold = maxIntegerValue / 10;
  Value::UInt lastDigitThreshold = Value::UInt(maxIntegerValue % 10);
  assert(lastDigitThreshold >= 0  &&  lastDigitThreshold <= 9);
  Value::LargestUInt value = 0;
  while (current < token.end_) {
    Char c = *current++;
    if (c < '0'  ||  c > '9')
      return addError("'" + std::string(token.start_, token.end_) + "' is not a number.", token);
    Value::UInt digit(c - '0');
    if (value >= threshold) {
      // If the current digit is not the last one, or if it is
      // greater than the last digit of the maximum integer value,
      // the parse the number as a double.
      if (current != token.end_  ||  digit > lastDigitThreshold) {
        return decodeDouble(token);
      }
    }
    value = value * 10 + digit;
  }
  if (isNegative)
    currentValue() = -Value::LargestInt(value);
  else if (value <= Value::LargestUInt(Value::maxInt))
    currentValue() = Value::LargestInt(value);
  else
    currentValue() = value;
  return true;
}


bool
Reader::decodeDouble(Token& token) {
  double value = 0;
  const int bufferSize = 32;
  int count;
  int length = int(token.end_ - token.start_);
  if (length <= bufferSize) {
    Char buffer[bufferSize + 1];
    memcpy(buffer, token.start_, length);
    buffer[length] = 0;
    count = sscanf(buffer, "%lf", &value);
  } else {
    std::string buffer(token.start_, token.end_);
    count = sscanf(buffer.c_str(), "%lf", &value);
  }

  if (count != 1)
    return addError("'" + std::string(token.start_, token.end_) + "' is not a number.", token);
  currentValue() = value;
  return true;
}


bool
Reader::decodeString(Token& token) {
  std::string decoded;
  if (!decodeString(token, decoded))
    return false;
  currentValue() = decoded;
  return true;
}


bool
Reader::decodeString(Token& token, std::string& decoded) {
  decoded.reserve(token.end_ - token.start_ - 2);
  Location current = token.start_ + 1; // skip '"'
  Location end = token.end_ - 1;      // do not include '"'
  while (current != end) {
    Char c = *current++;
    if (c == '"')
      break;
    else if (c == '\\') {
      if (current == end)
        return addError("Empty escape sequence in string", token, current);
      Char escape = *current++;
      switch (escape) {
        case '"': decoded += '"'; break;
        case '/': decoded += '/'; break;
        case '\\': decoded += '\\'; break;
        case 'b': decoded += '\b'; break;
        case 'f': decoded += '\f'; break;
        case 'n': decoded += '\n'; break;
        case 'r': decoded += '\r'; break;
        case 't': decoded += '\t'; break;
        case 'u': {
          unsigned int unicode;
          if (!decodeUnicodeCodePoint(token, current, end, unicode))
            return false;
          decoded += codePointToUTF8(unicode);
        }
        break;
        default:
          return addError("Bad escape sequence in string", token, current);
      }
    } else {
      decoded += c;
    }
  }
  return true;
}

bool
Reader::decodeUnicodeCodePoint(Token& token,
                               Location& current,
                               Location end,
                               unsigned int& unicode) {

  if (!decodeUnicodeEscapeSequence(token, current, end, unicode))
    return false;
  if (unicode >= 0xD800 && unicode <= 0xDBFF) {
    // surrogate pairs
    if (end - current < 6)
      return addError("additional six characters expected to parse unicode surrogate pair.", token, current);
    unsigned int surrogatePair;
    if (*(current++) == '\\' && *(current++) == 'u') {
      if (decodeUnicodeEscapeSequence(token, current, end, surrogatePair)) {
        unicode = 0x10000 + ((unicode & 0x3FF) << 10) + (surrogatePair & 0x3FF);
      } else
        return false;
    } else
      return addError("expecting another \\u token to begin the second half of a unicode surrogate pair", token, current);
  }
  return true;
}

bool
Reader::decodeUnicodeEscapeSequence(Token& token,
                                    Location& current,
                                    Location end,
                                    unsigned int& unicode) {
  if (end - current < 4)
    return addError("Bad unicode escape sequence in string: four digits expected.", token, current);
  unicode = 0;
  for (int index = 0; index < 4; ++index) {
    Char c = *current++;
    unicode *= 16;
    if (c >= '0'  &&  c <= '9')
      unicode += c - '0';
    else if (c >= 'a'  &&  c <= 'f')
      unicode += c - 'a' + 10;
    else if (c >= 'A'  &&  c <= 'F')
      unicode += c - 'A' + 10;
    else
      return addError("Bad unicode escape sequence in string: hexadecimal digit expected.", token, current);
  }
  return true;
}


bool
Reader::addError(const std::string& message,
                 Token& token,
                 Location extra) {
  ErrorInfo info;
  info.token_ = token;
  info.message_ = message;
  info.extra_ = extra;
  errors_.push_back(info);
  return false;
}


bool
Reader::recoverFromError(TokenType skipUntilToken) {
  int errorCount = int(errors_.size());
  Token skip;
  for (;;) {
    if (!readToken(skip))
      errors_.resize(errorCount);   // discard errors caused by recovery
    if (skip.type_ == skipUntilToken  ||  skip.type_ == tokenEndOfStream)
      break;
  }
  errors_.resize(errorCount);
  return false;
}


bool
Reader::addErrorAndRecover(const std::string& message,
                           Token& token,
                           TokenType skipUntilToken) {
  addError(message, token);
  return recoverFromError(skipUntilToken);
}


Value&
Reader::currentValue() {
  return *(nodes_.top());
}


Reader::Char
Reader::getNextChar() {
  if (current_ == end_)
    return 0;
  return *current_++;
}


void
Reader::getLocationLineAndColumn(Location location,
                                 int& line,
                                 int& column) const {
  Location current = begin_;
  Location lastLineStart = current;
  line = 0;
  while (current < location  &&  current != end_) {
    Char c = *current++;
    if (c == '\r') {
      if (*current == '\n')
        ++current;
      lastLineStart = current;
      ++line;
    } else if (c == '\n') {
      lastLineStart = current;
      ++line;
    }
  }
  // column & line start at 1
  column = int(location - lastLineStart) + 1;
  ++line;
}


std::string
Reader::getLocationLineAndColumn(Location location) const {
  int line, column;
  getLocationLineAndColumn(location, line, column);
  char buffer[18 + 16 + 16 + 1];
  sprintf(buffer, "Line %d, Column %d", line, column);
  return buffer;
}


// Deprecated. Preserved for backward compatibility
std::string
Reader::getFormatedErrorMessages() const {
  return getFormattedErrorMessages();
}


std::string
Reader::getFormattedErrorMessages() const {
  std::string formattedMessage;
  for (Errors::const_iterator itError = errors_.begin();
       itError != errors_.end();
       ++itError) {
    const ErrorInfo& error = *itError;
    formattedMessage += "* " + getLocationLineAndColumn(error.token_.start_) + "\n";
    formattedMessage += "  " + error.message_ + "\n";
    if (error.extra_)
      formattedMessage += "See " + getLocationLineAndColumn(error.extra_) + " for detail.\n";
  }
  return formattedMessage;
}


std::istream& operator>>(std::istream& sin, Value& root) {
  Json::Reader reader;
  bool ok = reader.parse(sin, root, true);
  //JSON_ASSERT( ok );
  if (!ok)
    throw std::runtime_error(reader.getFormattedErrorMessages());
  return sin;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_reader.cpp
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_batchallocator.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

#ifndef JSONCPP_BATCHALLOCATOR_H_INCLUDED
# define JSONCPP_BATCHALLOCATOR_H_INCLUDED

# include <stdlib.h>
# include <assert.h>

# ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION

namespace Json {

/* Fast memory allocator.
 *
 * This memory allocator allocates memory for a batch of object (specified by
 * the page size, the number of object in each page).
 *
 * It does not allow the destruction of a single object. All the allocated objects
 * can be destroyed at once. The memory can be either released or reused for future
 * allocation.
 *
 * The in-place new operator must be used to construct the object using the pointer
 * returned by allocate.
 */
template<typename AllocatedType
         , const unsigned int objectPerAllocation>
class BatchAllocator {
 public:
  typedef AllocatedType Type;

  BatchAllocator(unsigned int objectsPerPage = 255)
    : freeHead_(0)
    , objectsPerPage_(objectsPerPage) {
//      printf( "Size: %d => %s\n", sizeof(AllocatedType), typeid(AllocatedType).name() );
    assert(sizeof(AllocatedType) * objectPerAllocation >= sizeof(AllocatedType*));    // We must be able to store a slist in the object free space.
    assert(objectsPerPage >= 16);
    batches_ = allocateBatch(0);     // allocated a dummy page
    currentBatch_ = batches_;
  }

  ~BatchAllocator() {
    for (BatchInfo* batch = batches_; batch;) {
      BatchInfo* nextBatch = batch->next_;
      free(batch);
      batch = nextBatch;
    }
  }

  /// allocate space for an array of objectPerAllocation object.
  /// @warning it is the responsability of the caller to call objects constructors.
  AllocatedType* allocate() {
    if (freeHead_) {   // returns node from free list.
      AllocatedType* object = freeHead_;
      freeHead_ = *(AllocatedType**)object;
      return object;
    }
    if (currentBatch_->used_ == currentBatch_->end_) {
      currentBatch_ = currentBatch_->next_;
      while (currentBatch_  &&  currentBatch_->used_ == currentBatch_->end_)
        currentBatch_ = currentBatch_->next_;

      if (!currentBatch_) {    // no free batch found, allocate a new one
        currentBatch_ = allocateBatch(objectsPerPage_);
        currentBatch_->next_ = batches_; // insert at the head of the list
        batches_ = currentBatch_;
      }
    }
    AllocatedType* allocated = currentBatch_->used_;
    currentBatch_->used_ += objectPerAllocation;
    return allocated;
  }

  /// Release the object.
  /// @warning it is the responsability of the caller to actually destruct the object.
  void release(AllocatedType* object) {
    assert(object != 0);
    *(AllocatedType**)object = freeHead_;
    freeHead_ = object;
  }

 private:
  struct BatchInfo {
    BatchInfo* next_;
    AllocatedType* used_;
    AllocatedType* end_;
    AllocatedType buffer_[objectPerAllocation];
  };

  // disabled copy constructor and assignement operator.
  BatchAllocator(const BatchAllocator&);
  void operator =(const BatchAllocator&);

  static BatchInfo* allocateBatch(unsigned int objectsPerPage) {
    const unsigned int mallocSize = sizeof(BatchInfo) - sizeof(AllocatedType) * objectPerAllocation
                                    + sizeof(AllocatedType) * objectPerAllocation * objectsPerPage;
    BatchInfo* batch = static_cast<BatchInfo*>(malloc(mallocSize));
    batch->next_ = 0;
    batch->used_ = batch->buffer_;
    batch->end_ = batch->buffer_ + objectsPerPage;
    return batch;
  }

  BatchInfo* batches_;
  BatchInfo* currentBatch_;
  /// Head of a single linked list within the allocated space of freeed object
  AllocatedType* freeHead_;
  unsigned int objectsPerPage_;
};


} // namespace Json

# endif // ifndef JSONCPP_DOC_INCLUDE_IMPLEMENTATION

#endif // JSONCPP_BATCHALLOCATOR_H_INCLUDED


// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_batchallocator.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_valueiterator.inl
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

// included by json_value.cpp

namespace Json {

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueIteratorBase
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueIteratorBase::ValueIteratorBase()
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  : current_()
  , isNull_(true) {
}
#else
  : isArray_(true)
  , isNull_(true) {
  iterator_.array_ = ValueInternalArray::IteratorState();
}
#endif


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueIteratorBase::ValueIteratorBase(const Value::ObjectValues::iterator& current)
  : current_(current)
  , isNull_(false) {
}
#else
ValueIteratorBase::ValueIteratorBase(const ValueInternalArray::IteratorState& state)
  : isArray_(true) {
  iterator_.array_ = state;
}


ValueIteratorBase::ValueIteratorBase(const ValueInternalMap::IteratorState& state)
  : isArray_(false) {
  iterator_.map_ = state;
}
#endif

Value&
ValueIteratorBase::deref() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  return current_->second;
#else
  if (isArray_)
    return ValueInternalArray::dereference(iterator_.array_);
  return ValueInternalMap::value(iterator_.map_);
#endif
}


void
ValueIteratorBase::increment() {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  ++current_;
#else
  if (isArray_)
    ValueInternalArray::increment(iterator_.array_);
  ValueInternalMap::increment(iterator_.map_);
#endif
}


void
ValueIteratorBase::decrement() {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  --current_;
#else
  if (isArray_)
    ValueInternalArray::decrement(iterator_.array_);
  ValueInternalMap::decrement(iterator_.map_);
#endif
}


ValueIteratorBase::difference_type
ValueIteratorBase::computeDistance(const SelfType& other) const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
# ifdef JSON_USE_CPPTL_SMALLMAP
  return current_ - other.current_;
# else
  // Iterator for null value are initialized using the default
  // constructor, which initialize current_ to the default
  // std::map::iterator. As begin() and end() are two instance
  // of the default std::map::iterator, they can not be compared.
  // To allow this, we handle this comparison specifically.
  if (isNull_  &&  other.isNull_) {
    return 0;
  }


  // Usage of std::distance is not portable (does not compile with Sun Studio 12 RogueWave STL,
  // which is the one used by default).
  // Using a portable hand-made version for non random iterator instead:
  //   return difference_type( std::distance( current_, other.current_ ) );
  difference_type myDistance = 0;
  for (Value::ObjectValues::iterator it = current_; it != other.current_; ++it) {
    ++myDistance;
  }
  return myDistance;
# endif
#else
  if (isArray_)
    return ValueInternalArray::distance(iterator_.array_, other.iterator_.array_);
  return ValueInternalMap::distance(iterator_.map_, other.iterator_.map_);
#endif
}


bool
ValueIteratorBase::isEqual(const SelfType& other) const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  if (isNull_) {
    return other.isNull_;
  }
  return current_ == other.current_;
#else
  if (isArray_)
    return ValueInternalArray::equals(iterator_.array_, other.iterator_.array_);
  return ValueInternalMap::equals(iterator_.map_, other.iterator_.map_);
#endif
}


void
ValueIteratorBase::copy(const SelfType& other) {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  current_ = other.current_;
#else
  if (isArray_)
    iterator_.array_ = other.iterator_.array_;
  iterator_.map_ = other.iterator_.map_;
#endif
}


Value
ValueIteratorBase::key() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  const Value::CZString czstring = (*current_).first;
  if (czstring.c_str()) {
    if (czstring.isStaticString())
      return Value(StaticString(czstring.c_str()));
    return Value(czstring.c_str());
  }
  return Value(czstring.index());
#else
  if (isArray_)
    return Value(ValueInternalArray::indexOf(iterator_.array_));
  bool isStatic;
  const char* memberName = ValueInternalMap::key(iterator_.map_, isStatic);
  if (isStatic)
    return Value(StaticString(memberName));
  return Value(memberName);
#endif
}


UInt
ValueIteratorBase::index() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  const Value::CZString czstring = (*current_).first;
  if (!czstring.c_str())
    return czstring.index();
  return Value::UInt(-1);
#else
  if (isArray_)
    return Value::UInt(ValueInternalArray::indexOf(iterator_.array_));
  return Value::UInt(-1);
#endif
}


const char*
ValueIteratorBase::memberName() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  const char* name = (*current_).first.c_str();
  return name ? name : "";
#else
  if (!isArray_)
    return ValueInternalMap::key(iterator_.map_);
  return "";
#endif
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueConstIterator
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueConstIterator::ValueConstIterator() {
}


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueConstIterator::ValueConstIterator(const Value::ObjectValues::iterator& current)
  : ValueIteratorBase(current) {
}
#else
ValueConstIterator::ValueConstIterator(const ValueInternalArray::IteratorState& state)
  : ValueIteratorBase(state) {
}

ValueConstIterator::ValueConstIterator(const ValueInternalMap::IteratorState& state)
  : ValueIteratorBase(state) {
}
#endif

ValueConstIterator&
ValueConstIterator::operator =(const ValueIteratorBase& other) {
  copy(other);
  return *this;
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueIterator
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueIterator::ValueIterator() {
}


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueIterator::ValueIterator(const Value::ObjectValues::iterator& current)
  : ValueIteratorBase(current) {
}
#else
ValueIterator::ValueIterator(const ValueInternalArray::IteratorState& state)
  : ValueIteratorBase(state) {
}

ValueIterator::ValueIterator(const ValueInternalMap::IteratorState& state)
  : ValueIteratorBase(state) {
}
#endif

ValueIterator::ValueIterator(const ValueConstIterator& other)
  : ValueIteratorBase(other) {
}

ValueIterator::ValueIterator(const ValueIterator& other)
  : ValueIteratorBase(other) {
}

ValueIterator&
ValueIterator::operator =(const SelfType& other) {
  copy(other);
  return *this;
}

} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_valueiterator.inl
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_value.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/value.h>
# include <json/writer.h>
# ifndef JSON_USE_SIMPLE_INTERNAL_ALLOCATOR
#  include "json_batchallocator.h"
# endif // #ifndef JSON_USE_SIMPLE_INTERNAL_ALLOCATOR
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <iostream>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <cassert>
#ifdef JSON_USE_CPPTL
# include <cpptl/conststring.h>
#endif
#include <cstddef>    // size_t

#define JSON_ASSERT_UNREACHABLE assert( false )
#define JSON_ASSERT( condition ) assert( condition );  // @todo <= change this into an exception throw
#define JSON_FAIL_MESSAGE( message ) throw std::runtime_error( message );
#define JSON_ASSERT_MESSAGE( condition, message ) if (!( condition )) JSON_FAIL_MESSAGE( message )

namespace Json {

const Value Value::null;
const Int Value::minInt = Int(~(UInt(-1) / 2));
const Int Value::maxInt = Int(UInt(-1) / 2);
const UInt Value::maxUInt = UInt(-1);
const Int64 Value::minInt64 = Int64(~(UInt64(-1) / 2));
const Int64 Value::maxInt64 = Int64(UInt64(-1) / 2);
const UInt64 Value::maxUInt64 = UInt64(-1);
const LargestInt Value::minLargestInt = LargestInt(~(LargestUInt(-1) / 2));
const LargestInt Value::maxLargestInt = LargestInt(LargestUInt(-1) / 2);
const LargestUInt Value::maxLargestUInt = LargestUInt(-1);


/// Unknown size marker
static const unsigned int unknown = (unsigned) - 1;


/** Duplicates the specified string value.
 * @param value Pointer to the string to duplicate. Must be zero-terminated if
 *              length is "unknown".
 * @param length Length of the value. if equals to unknown, then it will be
 *               computed using strlen(value).
 * @return Pointer on the duplicate instance of string.
 */
static inline char*
duplicateStringValue(const char* value,
                     unsigned int length = unknown) {
  if (length == unknown)
    length = (unsigned int)strlen(value);
  char* newString = static_cast<char*>(malloc(length + 1));
  JSON_ASSERT_MESSAGE(newString != 0, "Failed to allocate string value buffer");
  memcpy(newString, value, length);
  newString[length] = 0;
  return newString;
}


/** Free the string duplicated by duplicateStringValue().
 */
static inline void
releaseStringValue(char* value) {
  if (value)
    free(value);
}

} // namespace Json


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// ValueInternals...
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
#if !defined(JSON_IS_AMALGAMATION)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
#  include "json_internalarray.inl"
#  include "json_internalmap.inl"
# endif // JSON_VALUE_USE_INTERNAL_MAP

# include "json_valueiterator.inl"
#endif // if !defined(JSON_IS_AMALGAMATION)

namespace Json {

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::CommentInfo
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////


Value::CommentInfo::CommentInfo()
  : comment_(0) {
}

Value::CommentInfo::~CommentInfo() {
  if (comment_)
    releaseStringValue(comment_);
}


void
Value::CommentInfo::setComment(const char* text) {
  if (comment_)
    releaseStringValue(comment_);
  JSON_ASSERT(text != 0);
  JSON_ASSERT_MESSAGE(text[0] == '\0' || text[0] == '/', "Comments must start with /");
  // It seems that /**/ style comments are acceptable as well.
  comment_ = duplicateStringValue(text);
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::CZString
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
# ifndef JSON_VALUE_USE_INTERNAL_MAP

// Notes: index_ indicates if the string was allocated when
// a string is stored.

Value::CZString::CZString(ArrayIndex index)
  : cstr_(0)
  , index_(index) {
}

Value::CZString::CZString(const char* cstr, DuplicationPolicy allocate)
  : cstr_(allocate == duplicate ? duplicateStringValue(cstr)
          : cstr)
  , index_(allocate) {
}

Value::CZString::CZString(const CZString& other)
  : cstr_(other.index_ != noDuplication &&  other.cstr_ != 0
          ?  duplicateStringValue(other.cstr_)
          : other.cstr_)
  , index_(other.cstr_ ? (other.index_ == noDuplication ? noDuplication : duplicate)
           : other.index_) {
}

Value::CZString::~CZString() {
  if (cstr_  &&  index_ == duplicate)
    releaseStringValue(const_cast<char*>(cstr_));
}

void
Value::CZString::swap(CZString& other) {
  std::swap(cstr_, other.cstr_);
  std::swap(index_, other.index_);
}

Value::CZString&
Value::CZString::operator =(const CZString& other) {
  CZString temp(other);
  swap(temp);
  return *this;
}

bool
Value::CZString::operator<(const CZString& other) const  {
  if (cstr_)
    return strcmp(cstr_, other.cstr_) < 0;
  return index_ < other.index_;
}

bool
Value::CZString::operator==(const CZString& other) const  {
  if (cstr_)
    return strcmp(cstr_, other.cstr_) == 0;
  return index_ == other.index_;
}


ArrayIndex
Value::CZString::index() const {
  return index_;
}


const char*
Value::CZString::c_str() const {
  return cstr_;
}

bool
Value::CZString::isStaticString() const {
  return index_ == noDuplication;
}

#endif // ifndef JSON_VALUE_USE_INTERNAL_MAP


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::Value
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

/*! \internal Default constructor initialization must be equivalent to:
 * memset( this, 0, sizeof(Value) )
 * This optimization is used in ValueInternalMap fast allocator.
 */
Value::Value(ValueType type)
  : type_(type)
  , allocated_(0)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  switch (type) {
    case nullValue:
      break;
    case intValue:
    case uintValue:
      value_.int_ = 0;
      break;
    case realValue:
      value_.real_ = 0.0;
      break;
    case stringValue:
      value_.string_ = 0;
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue:
      value_.map_ = new ObjectValues();
      break;
#else
    case arrayValue:
      value_.array_ = arrayAllocator()->newArray();
      break;
    case objectValue:
      value_.map_ = mapAllocator()->newMap();
      break;
#endif
    case booleanValue:
      value_.bool_ = false;
      break;
    default:
      JSON_ASSERT_UNREACHABLE;
  }
}


#if defined(JSON_HAS_INT64)
Value::Value(UInt value)
  : type_(uintValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.uint_ = value;
}

Value::Value(Int value)
  : type_(intValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.int_ = value;
}

#endif // if defined(JSON_HAS_INT64)


Value::Value(Int64 value)
  : type_(intValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.int_ = value;
}


Value::Value(UInt64 value)
  : type_(uintValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.uint_ = value;
}

Value::Value(double value)
  : type_(realValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.real_ = value;
}

Value::Value(const char* value)
  : type_(stringValue)
  , allocated_(true)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.string_ = duplicateStringValue(value);
}


Value::Value(const char* beginValue,
             const char* endValue)
  : type_(stringValue)
  , allocated_(true)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.string_ = duplicateStringValue(beginValue,
                                        (unsigned int)(endValue - beginValue));
}


Value::Value(const std::string& value)
  : type_(stringValue)
  , allocated_(true)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.string_ = duplicateStringValue(value.c_str(),
                                        (unsigned int)value.length());

}

Value::Value(const StaticString& value)
  : type_(stringValue)
  , allocated_(false)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.string_ = const_cast<char*>(value.c_str());
}


# ifdef JSON_USE_CPPTL
Value::Value(const CppTL::ConstString& value)
  : type_(stringValue)
  , allocated_(true)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.string_ = duplicateStringValue(value, value.length());
}
# endif

Value::Value(bool value)
  : type_(booleanValue)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  value_.bool_ = value;
}


Value::Value(const Value& other)
  : type_(other.type_)
  , comments_(0)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
  , itemIsUsed_(0)
#endif
{
  switch (type_) {
    case nullValue:
    case intValue:
    case uintValue:
    case realValue:
    case booleanValue:
      value_ = other.value_;
      break;
    case stringValue:
      if (other.value_.string_) {
        value_.string_ = duplicateStringValue(other.value_.string_);
        allocated_ = true;
      } else
        value_.string_ = 0;
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue:
      value_.map_ = new ObjectValues(*other.value_.map_);
      break;
#else
    case arrayValue:
      value_.array_ = arrayAllocator()->newArrayCopy(*other.value_.array_);
      break;
    case objectValue:
      value_.map_ = mapAllocator()->newMapCopy(*other.value_.map_);
      break;
#endif
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  if (other.comments_) {
    comments_ = new CommentInfo[numberOfCommentPlacement];
    for (int comment = 0; comment < numberOfCommentPlacement; ++comment) {
      const CommentInfo& otherComment = other.comments_[comment];
      if (otherComment.comment_)
        comments_[comment].setComment(otherComment.comment_);
    }
  }
}


Value::~Value() {
  switch (type_) {
    case nullValue:
    case intValue:
    case uintValue:
    case realValue:
    case booleanValue:
      break;
    case stringValue:
      if (allocated_)
        releaseStringValue(value_.string_);
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue:
      delete value_.map_;
      break;
#else
    case arrayValue:
      arrayAllocator()->destructArray(value_.array_);
      break;
    case objectValue:
      mapAllocator()->destructMap(value_.map_);
      break;
#endif
    default:
      JSON_ASSERT_UNREACHABLE;
  }

  if (comments_)
    delete[] comments_;
}

Value&
Value::operator=(const Value& other) {
  Value temp(other);
  swap(temp);
  return *this;
}

void
Value::swap(Value& other) {
  ValueType temp = type_;
  type_ = other.type_;
  other.type_ = temp;
  std::swap(value_, other.value_);
  int temp2 = allocated_;
  allocated_ = other.allocated_;
  other.allocated_ = temp2;
}

ValueType
Value::type() const {
  return type_;
}


int
Value::compare(const Value& other) const {
  if (*this < other)
    return -1;
  if (*this > other)
    return 1;
  return 0;
}


bool
Value::operator <(const Value& other) const {
  int typeDelta = type_ - other.type_;
  if (typeDelta)
    return typeDelta < 0 ? true : false;
  switch (type_) {
    case nullValue:
      return false;
    case intValue:
      return value_.int_ < other.value_.int_;
    case uintValue:
      return value_.uint_ < other.value_.uint_;
    case realValue:
      return value_.real_ < other.value_.real_;
    case booleanValue:
      return value_.bool_ < other.value_.bool_;
    case stringValue:
      return (value_.string_ == 0  &&  other.value_.string_)
             || (other.value_.string_
                 &&  value_.string_
                 && strcmp(value_.string_, other.value_.string_) < 0);
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue: {
      int delta = int(value_.map_->size() - other.value_.map_->size());
      if (delta)
        return delta < 0;
      return (*value_.map_) < (*other.value_.map_);
    }
#else
    case arrayValue:
      return value_.array_->compare(*(other.value_.array_)) < 0;
    case objectValue:
      return value_.map_->compare(*(other.value_.map_)) < 0;
#endif
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return false;  // unreachable
}

bool
Value::operator <=(const Value& other) const {
  return !(other < *this);
}

bool
Value::operator >=(const Value& other) const {
  return !(*this < other);
}

bool
Value::operator >(const Value& other) const {
  return other < *this;
}

bool
Value::operator ==(const Value& other) const {
  //if ( type_ != other.type_ )
  // GCC 2.95.3 says:
  // attempt to take address of bit-field structure member `Json::Value::type_'
  // Beats me, but a temp solves the problem.
  int temp = other.type_;
  if (type_ != temp)
    return false;
  switch (type_) {
    case nullValue:
      return true;
    case intValue:
      return value_.int_ == other.value_.int_;
    case uintValue:
      return value_.uint_ == other.value_.uint_;
    case realValue:
      return value_.real_ == other.value_.real_;
    case booleanValue:
      return value_.bool_ == other.value_.bool_;
    case stringValue:
      return (value_.string_ == other.value_.string_)
             || (other.value_.string_
                 &&  value_.string_
                 && strcmp(value_.string_, other.value_.string_) == 0);
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue:
      return value_.map_->size() == other.value_.map_->size()
             && (*value_.map_) == (*other.value_.map_);
#else
    case arrayValue:
      return value_.array_->compare(*(other.value_.array_)) == 0;
    case objectValue:
      return value_.map_->compare(*(other.value_.map_)) == 0;
#endif
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return false;  // unreachable
}

bool
Value::operator !=(const Value& other) const {
  return !(*this == other);
}

const char*
Value::asCString() const {
  JSON_ASSERT(type_ == stringValue);
  return value_.string_;
}


std::string
Value::asString() const {
  switch (type_) {
    case nullValue:
      return "";
    case stringValue:
      return value_.string_ ? value_.string_ : "";
    case booleanValue:
      return value_.bool_ ? "true" : "false";
    case intValue:
    case uintValue:
    case realValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to string");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return ""; // unreachable
}

# ifdef JSON_USE_CPPTL
CppTL::ConstString
Value::asConstString() const {
  return CppTL::ConstString(asString().c_str());
}
# endif


Value::Int
Value::asInt() const {
  switch (type_) {
    case nullValue:
      return 0;
    case intValue:
      JSON_ASSERT_MESSAGE(value_.int_ >= minInt  &&  value_.int_ <= maxInt, "unsigned integer out of signed int range");
      return Int(value_.int_);
    case uintValue:
      JSON_ASSERT_MESSAGE(value_.uint_ <= UInt(maxInt), "unsigned integer out of signed int range");
      return Int(value_.uint_);
    case realValue:
      JSON_ASSERT_MESSAGE(value_.real_ >= minInt  &&  value_.real_ <= maxInt, "Real out of signed integer range");
      return Int(value_.real_);
    case booleanValue:
      return value_.bool_ ? 1 : 0;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to int");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}


Value::UInt
Value::asUInt() const {
  switch (type_) {
    case nullValue:
      return 0;
    case intValue:
      JSON_ASSERT_MESSAGE(value_.int_ >= 0, "Negative integer can not be converted to unsigned integer");
      JSON_ASSERT_MESSAGE(value_.int_ <= maxUInt, "signed integer out of UInt range");
      return UInt(value_.int_);
    case uintValue:
      JSON_ASSERT_MESSAGE(value_.uint_ <= maxUInt, "unsigned integer out of UInt range");
      return UInt(value_.uint_);
    case realValue:
      JSON_ASSERT_MESSAGE(value_.real_ >= 0  &&  value_.real_ <= maxUInt,  "Real out of unsigned integer range");
      return UInt(value_.real_);
    case booleanValue:
      return value_.bool_ ? 1 : 0;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to uint");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}


# if defined(JSON_HAS_INT64)

Value::Int64
Value::asInt64() const {
  switch (type_) {
    case nullValue:
      return 0;
    case intValue:
      return value_.int_;
    case uintValue:
      JSON_ASSERT_MESSAGE(value_.uint_ <= UInt64(maxInt64), "unsigned integer out of Int64 range");
      return value_.uint_;
    case realValue:
      JSON_ASSERT_MESSAGE(value_.real_ >= minInt64  &&  value_.real_ <= maxInt64, "Real out of Int64 range");
      return Int(value_.real_);
    case booleanValue:
      return value_.bool_ ? 1 : 0;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to Int64");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}


Value::UInt64
Value::asUInt64() const {
  switch (type_) {
    case nullValue:
      return 0;
    case intValue:
      JSON_ASSERT_MESSAGE(value_.int_ >= 0, "Negative integer can not be converted to UInt64");
      return value_.int_;
    case uintValue:
      return value_.uint_;
    case realValue:
      JSON_ASSERT_MESSAGE(value_.real_ >= 0  &&  value_.real_ <= maxUInt64,  "Real out of UInt64 range");
      return UInt(value_.real_);
    case booleanValue:
      return value_.bool_ ? 1 : 0;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to UInt64");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}
# endif // if defined(JSON_HAS_INT64)


LargestInt
Value::asLargestInt() const {
#if defined(JSON_NO_INT64)
  return asInt();
#else
  return asInt64();
#endif
}


LargestUInt
Value::asLargestUInt() const {
#if defined(JSON_NO_INT64)
  return asUInt();
#else
  return asUInt64();
#endif
}


double
Value::asDouble() const {
  switch (type_) {
    case nullValue:
      return 0.0;
    case intValue:
      return static_cast<double>(value_.int_);
    case uintValue:
#if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<double>(value_.uint_);
#else // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<double>(Int(value_.uint_ / 2)) * 2 + Int(value_.uint_ & 1);
#endif // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
    case realValue:
      return value_.real_;
    case booleanValue:
      return value_.bool_ ? 1.0 : 0.0;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to double");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}

float
Value::asFloat() const {
  switch (type_) {
    case nullValue:
      return 0.0f;
    case intValue:
      return static_cast<float>(value_.int_);
    case uintValue:
#if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<float>(value_.uint_);
#else // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<float>(Int(value_.uint_ / 2)) * 2 + Int(value_.uint_ & 1);
#endif // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
    case realValue:
      return static_cast<float>(value_.real_);
    case booleanValue:
      return value_.bool_ ? 1.0f : 0.0f;
    case stringValue:
    case arrayValue:
    case objectValue:
      JSON_FAIL_MESSAGE("Type is not convertible to float");
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0.0f; // unreachable;
}

bool
Value::asBool() const {
  switch (type_) {
    case nullValue:
      return false;
    case intValue:
    case uintValue:
      return value_.int_ != 0;
    case realValue:
      return value_.real_ != 0.0;
    case booleanValue:
      return value_.bool_;
    case stringValue:
      return value_.string_  &&  value_.string_[0] != 0;
    case arrayValue:
    case objectValue:
      return value_.map_->size() != 0;
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return false; // unreachable;
}


bool
Value::isConvertibleTo(ValueType other) const {
  switch (type_) {
    case nullValue:
      return true;
    case intValue:
      return (other == nullValue  &&  value_.int_ == 0)
             || other == intValue
             || (other == uintValue  && value_.int_ >= 0)
             || other == realValue
             || other == stringValue
             || other == booleanValue;
    case uintValue:
      return (other == nullValue  &&  value_.uint_ == 0)
             || (other == intValue  && value_.uint_ <= (unsigned)maxInt)
             || other == uintValue
             || other == realValue
             || other == stringValue
             || other == booleanValue;
    case realValue:
      return (other == nullValue  &&  value_.real_ == 0.0)
             || (other == intValue  &&  value_.real_ >= minInt  &&  value_.real_ <= maxInt)
             || (other == uintValue  &&  value_.real_ >= 0  &&  value_.real_ <= maxUInt)
             || other == realValue
             || other == stringValue
             || other == booleanValue;
    case booleanValue:
      return (other == nullValue  &&  value_.bool_ == false)
             || other == intValue
             || other == uintValue
             || other == realValue
             || other == stringValue
             || other == booleanValue;
    case stringValue:
      return other == stringValue
             || (other == nullValue  && (!value_.string_  ||  value_.string_[0] == 0));
    case arrayValue:
      return other == arrayValue
             || (other == nullValue  &&  value_.map_->size() == 0);
    case objectValue:
      return other == objectValue
             || (other == nullValue  &&  value_.map_->size() == 0);
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return false; // unreachable;
}


/// Number of values in array or object
ArrayIndex
Value::size() const {
  switch (type_) {
    case nullValue:
    case intValue:
    case uintValue:
    case realValue:
    case booleanValue:
    case stringValue:
      return 0;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:  // size of the array is highest index + 1
      if (!value_.map_->empty()) {
        ObjectValues::const_iterator itLast = value_.map_->end();
        --itLast;
        return (*itLast).first.index() + 1;
      }
      return 0;
    case objectValue:
      return ArrayIndex(value_.map_->size());
#else
    case arrayValue:
      return Int(value_.array_->size());
    case objectValue:
      return Int(value_.map_->size());
#endif
    default:
      JSON_ASSERT_UNREACHABLE;
  }
  return 0; // unreachable;
}


bool
Value::empty() const {
  if (isNull() || isArray() || isObject())
    return size() == 0u;
  else
    return false;
}


bool
Value::operator!() const {
  return isNull();
}


void
Value::clear() {
  JSON_ASSERT(type_ == nullValue  ||  type_ == arrayValue  || type_ == objectValue);

  switch (type_) {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
    case objectValue:
      value_.map_->clear();
      break;
#else
    case arrayValue:
      value_.array_->clear();
      break;
    case objectValue:
      value_.map_->clear();
      break;
#endif
    default:
      break;
  }
}

void
Value::resize(ArrayIndex newSize) {
  JSON_ASSERT(type_ == nullValue  ||  type_ == arrayValue);
  if (type_ == nullValue)
    *this = Value(arrayValue);
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  ArrayIndex oldSize = size();
  if (newSize == 0)
    clear();
  else if (newSize > oldSize)
    (*this)[ newSize - 1 ];
  else {
    for (ArrayIndex index = newSize; index < oldSize; ++index) {
      value_.map_->erase(index);
    }
    assert(size() == newSize);
  }
#else
  value_.array_->resize(newSize);
#endif
}


Value&
Value::operator[](ArrayIndex index) {
  JSON_ASSERT(type_ == nullValue  ||  type_ == arrayValue);
  if (type_ == nullValue)
    *this = Value(arrayValue);
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  CZString key(index);
  ObjectValues::iterator it = value_.map_->lower_bound(key);
  if (it != value_.map_->end()  && (*it).first == key)
    return (*it).second;

  ObjectValues::value_type defaultValue(key, null);
  it = value_.map_->insert(it, defaultValue);
  return (*it).second;
#else
  return value_.array_->resolveReference(index);
#endif
}


Value&
Value::operator[](int index) {
  JSON_ASSERT(index >= 0);
  return (*this)[ ArrayIndex(index) ];
}


const Value&
Value::operator[](ArrayIndex index) const {
  JSON_ASSERT(type_ == nullValue  ||  type_ == arrayValue);
  if (type_ == nullValue)
    return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  CZString key(index);
  ObjectValues::const_iterator it = value_.map_->find(key);
  if (it == value_.map_->end())
    return null;
  return (*it).second;
#else
  Value* value = value_.array_->find(index);
  return value ? *value : null;
#endif
}


const Value&
Value::operator[](int index) const {
  JSON_ASSERT(index >= 0);
  return (*this)[ ArrayIndex(index) ];
}


Value&
Value::operator[](const char* key) {
  return resolveReference(key, false);
}


Value&
Value::resolveReference(const char* key,
                        bool isStatic) {
  JSON_ASSERT(type_ == nullValue  ||  type_ == objectValue);
  if (type_ == nullValue)
    *this = Value(objectValue);
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  CZString actualKey(key, isStatic ? CZString::noDuplication
                     : CZString::duplicateOnCopy);
  ObjectValues::iterator it = value_.map_->lower_bound(actualKey);
  if (it != value_.map_->end()  && (*it).first == actualKey)
    return (*it).second;

  ObjectValues::value_type defaultValue(actualKey, null);
  it = value_.map_->insert(it, defaultValue);
  Value& value = (*it).second;
  return value;
#else
  return value_.map_->resolveReference(key, isStatic);
#endif
}


Value
Value::get(ArrayIndex index,
           const Value& defaultValue) const {
  const Value* value = &((*this)[index]);
  return value == &null ? defaultValue : *value;
}


bool
Value::isValidIndex(ArrayIndex index) const {
  return index < size();
}



const Value&
Value::operator[](const char* key) const {
  JSON_ASSERT(type_ == nullValue  ||  type_ == objectValue);
  if (type_ == nullValue)
    return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  CZString actualKey(key, CZString::noDuplication);
  ObjectValues::const_iterator it = value_.map_->find(actualKey);
  if (it == value_.map_->end())
    return null;
  return (*it).second;
#else
  const Value* value = value_.map_->find(key);
  return value ? *value : null;
#endif
}


Value&
Value::operator[](const std::string& key) {
  return (*this)[ key.c_str() ];
}


const Value&
Value::operator[](const std::string& key) const {
  return (*this)[ key.c_str() ];
}

Value&
Value::operator[](const StaticString& key) {
  return resolveReference(key, true);
}


# ifdef JSON_USE_CPPTL
Value&
Value::operator[](const CppTL::ConstString& key) {
  return (*this)[ key.c_str() ];
}


const Value&
Value::operator[](const CppTL::ConstString& key) const {
  return (*this)[ key.c_str() ];
}
# endif


Value&
Value::append(const Value& value) {
  return (*this)[size()] = value;
}


Value
Value::get(const char* key,
           const Value& defaultValue) const {
  const Value* value = &((*this)[key]);
  return value == &null ? defaultValue : *value;
}


Value
Value::get(const std::string& key,
           const Value& defaultValue) const {
  return get(key.c_str(), defaultValue);
}

Value
Value::removeMember(const char* key) {
  JSON_ASSERT(type_ == nullValue  ||  type_ == objectValue);
  if (type_ == nullValue)
    return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  CZString actualKey(key, CZString::noDuplication);
  ObjectValues::iterator it = value_.map_->find(actualKey);
  if (it == value_.map_->end())
    return null;
  Value old(it->second);
  value_.map_->erase(it);
  return old;
#else
  Value* value = value_.map_->find(key);
  if (value) {
    Value old(*value);
    value_.map_.remove(key);
    return old;
  } else {
    return null;
  }
#endif
}

Value
Value::removeMember(const std::string& key) {
  return removeMember(key.c_str());
}

# ifdef JSON_USE_CPPTL
Value
Value::get(const CppTL::ConstString& key,
           const Value& defaultValue) const {
  return get(key.c_str(), defaultValue);
}
# endif

bool
Value::isMember(const char* key) const {
  const Value* value = &((*this)[key]);
  return value != &null;
}


bool
Value::isMember(const std::string& key) const {
  return isMember(key.c_str());
}


# ifdef JSON_USE_CPPTL
bool
Value::isMember(const CppTL::ConstString& key) const {
  return isMember(key.c_str());
}
#endif

Value::Members
Value::getMemberNames() const {
  JSON_ASSERT(type_ == nullValue  ||  type_ == objectValue);
  if (type_ == nullValue)
    return Value::Members();
  Members members;
  members.reserve(value_.map_->size());
#ifndef JSON_VALUE_USE_INTERNAL_MAP
  ObjectValues::const_iterator it = value_.map_->begin();
  ObjectValues::const_iterator itEnd = value_.map_->end();
  for (; it != itEnd; ++it)
    members.push_back(std::string((*it).first.c_str()));
#else
  ValueInternalMap::IteratorState it;
  ValueInternalMap::IteratorState itEnd;
  value_.map_->makeBeginIterator(it);
  value_.map_->makeEndIterator(itEnd);
  for (; !ValueInternalMap::equals(it, itEnd); ValueInternalMap::increment(it))
    members.push_back(std::string(ValueInternalMap::key(it)));
#endif
  return members;
}
//
//# ifdef JSON_USE_CPPTL
//EnumMemberNames
//Value::enumMemberNames() const
//{
//   if ( type_ == objectValue )
//   {
//      return CppTL::Enum::any(  CppTL::Enum::transform(
//         CppTL::Enum::keys( *(value_.map_), CppTL::Type<const CZString &>() ),
//         MemberNamesTransform() ) );
//   }
//   return EnumMemberNames();
//}
//
//
//EnumValues
//Value::enumValues() const
//{
//   if ( type_ == objectValue  ||  type_ == arrayValue )
//      return CppTL::Enum::anyValues( *(value_.map_),
//                                     CppTL::Type<const Value &>() );
//   return EnumValues();
//}
//
//# endif


bool
Value::isNull() const {
  return type_ == nullValue;
}


bool
Value::isBool() const {
  return type_ == booleanValue;
}


bool
Value::isInt() const {
  return type_ == intValue;
}


bool
Value::isUInt() const {
  return type_ == uintValue;
}


bool
Value::isIntegral() const {
  return type_ == intValue
         ||  type_ == uintValue
         ||  type_ == booleanValue;
}


bool
Value::isDouble() const {
  return type_ == realValue;
}


bool
Value::isNumeric() const {
  return isIntegral() || isDouble();
}


bool
Value::isString() const {
  return type_ == stringValue;
}


bool
Value::isArray() const {
  return type_ == nullValue  ||  type_ == arrayValue;
}


bool
Value::isObject() const {
  return type_ == nullValue  ||  type_ == objectValue;
}


void
Value::setComment(const char* comment,
                  CommentPlacement placement) {
  if (!comments_)
    comments_ = new CommentInfo[numberOfCommentPlacement];
  comments_[placement].setComment(comment);
}


void
Value::setComment(const std::string& comment,
                  CommentPlacement placement) {
  setComment(comment.c_str(), placement);
}


bool
Value::hasComment(CommentPlacement placement) const {
  return comments_ != 0  &&  comments_[placement].comment_ != 0;
}

std::string
Value::getComment(CommentPlacement placement) const {
  if (hasComment(placement))
    return comments_[placement].comment_;
  return "";
}


std::string
Value::toStyledString() const {
  StyledWriter writer;
  return writer.write(*this);
}


Value::const_iterator
Value::begin() const {
  switch (type_) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
      if (value_.array_) {
        ValueInternalArray::IteratorState it;
        value_.array_->makeBeginIterator(it);
        return const_iterator(it);
      }
      break;
    case objectValue:
      if (value_.map_) {
        ValueInternalMap::IteratorState it;
        value_.map_->makeBeginIterator(it);
        return const_iterator(it);
      }
      break;
#else
    case arrayValue:
    case objectValue:
      if (value_.map_)
        return const_iterator(value_.map_->begin());
      break;
#endif
    default:
      break;
  }
  return const_iterator();
}

Value::const_iterator
Value::end() const {
  switch (type_) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
      if (value_.array_) {
        ValueInternalArray::IteratorState it;
        value_.array_->makeEndIterator(it);
        return const_iterator(it);
      }
      break;
    case objectValue:
      if (value_.map_) {
        ValueInternalMap::IteratorState it;
        value_.map_->makeEndIterator(it);
        return const_iterator(it);
      }
      break;
#else
    case arrayValue:
    case objectValue:
      if (value_.map_)
        return const_iterator(value_.map_->end());
      break;
#endif
    default:
      break;
  }
  return const_iterator();
}


Value::iterator
Value::begin() {
  switch (type_) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
      if (value_.array_) {
        ValueInternalArray::IteratorState it;
        value_.array_->makeBeginIterator(it);
        return iterator(it);
      }
      break;
    case objectValue:
      if (value_.map_) {
        ValueInternalMap::IteratorState it;
        value_.map_->makeBeginIterator(it);
        return iterator(it);
      }
      break;
#else
    case arrayValue:
    case objectValue:
      if (value_.map_)
        return iterator(value_.map_->begin());
      break;
#endif
    default:
      break;
  }
  return iterator();
}

Value::iterator
Value::end() {
  switch (type_) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
    case arrayValue:
      if (value_.array_) {
        ValueInternalArray::IteratorState it;
        value_.array_->makeEndIterator(it);
        return iterator(it);
      }
      break;
    case objectValue:
      if (value_.map_) {
        ValueInternalMap::IteratorState it;
        value_.map_->makeEndIterator(it);
        return iterator(it);
      }
      break;
#else
    case arrayValue:
    case objectValue:
      if (value_.map_)
        return iterator(value_.map_->end());
      break;
#endif
    default:
      break;
  }
  return iterator();
}


// class PathArgument
// //////////////////////////////////////////////////////////////////

PathArgument::PathArgument()
  : kind_(kindNone) {
}


PathArgument::PathArgument(ArrayIndex index)
  : index_(index)
  , kind_(kindIndex) {
}


PathArgument::PathArgument(const char* key)
  : key_(key)
  , kind_(kindKey) {
}


PathArgument::PathArgument(const std::string& key)
  : key_(key.c_str())
  , kind_(kindKey) {
}

// class Path
// //////////////////////////////////////////////////////////////////

Path::Path(const std::string& path,
           const PathArgument& a1,
           const PathArgument& a2,
           const PathArgument& a3,
           const PathArgument& a4,
           const PathArgument& a5) {
  InArgs in;
  in.push_back(&a1);
  in.push_back(&a2);
  in.push_back(&a3);
  in.push_back(&a4);
  in.push_back(&a5);
  makePath(path, in);
}


void
Path::makePath(const std::string& path,
               const InArgs& in) {
  const char* current = path.c_str();
  const char* end = current + path.length();
  InArgs::const_iterator itInArg = in.begin();
  while (current != end) {
    if (*current == '[') {
      ++current;
      if (*current == '%')
        addPathInArg(path, in, itInArg, PathArgument::kindIndex);
      else {
        ArrayIndex index = 0;
        for (; current != end && *current >= '0'  &&  *current <= '9'; ++current)
          index = index * 10 + ArrayIndex(*current - '0');
        args_.push_back(index);
      }
      if (current == end  ||  *current++ != ']')
        invalidPath(path, int(current - path.c_str()));
    } else if (*current == '%') {
      addPathInArg(path, in, itInArg, PathArgument::kindKey);
      ++current;
    } else if (*current == '.') {
      ++current;
    } else {
      const char* beginName = current;
      while (current != end  &&  !strchr("[.", *current))
        ++current;
      args_.push_back(std::string(beginName, current));
    }
  }
}


void
Path::addPathInArg(const std::string& path,
                   const InArgs& in,
                   InArgs::const_iterator& itInArg,
                   PathArgument::Kind kind) {
  if (itInArg == in.end()) {
    // Error: missing argument %d
  } else if ((*itInArg)->kind_ != kind) {
    // Error: bad argument type
  } else {
    args_.push_back(**itInArg);
  }
}


void
Path::invalidPath(const std::string& path,
                  int location) {
  // Error: invalid path.
}


const Value&
Path::resolve(const Value& root) const {
  const Value* node = &root;
  for (Args::const_iterator it = args_.begin(); it != args_.end(); ++it) {
    const PathArgument& arg = *it;
    if (arg.kind_ == PathArgument::kindIndex) {
      if (!node->isArray()  ||  node->isValidIndex(arg.index_)) {
        // Error: unable to resolve path (array value expected at position...
      }
      node = &((*node)[arg.index_]);
    } else if (arg.kind_ == PathArgument::kindKey) {
      if (!node->isObject()) {
        // Error: unable to resolve path (object value expected at position...)
      }
      node = &((*node)[arg.key_]);
      if (node == &Value::null) {
        // Error: unable to resolve path (object has no member named '' at position...)
      }
    }
  }
  return *node;
}


Value
Path::resolve(const Value& root,
              const Value& defaultValue) const {
  const Value* node = &root;
  for (Args::const_iterator it = args_.begin(); it != args_.end(); ++it) {
    const PathArgument& arg = *it;
    if (arg.kind_ == PathArgument::kindIndex) {
      if (!node->isArray()  ||  node->isValidIndex(arg.index_))
        return defaultValue;
      node = &((*node)[arg.index_]);
    } else if (arg.kind_ == PathArgument::kindKey) {
      if (!node->isObject())
        return defaultValue;
      node = &((*node)[arg.key_]);
      if (node == &Value::null)
        return defaultValue;
    }
  }
  return *node;
}


Value&
Path::make(Value& root) const {
  Value* node = &root;
  for (Args::const_iterator it = args_.begin(); it != args_.end(); ++it) {
    const PathArgument& arg = *it;
    if (arg.kind_ == PathArgument::kindIndex) {
      if (!node->isArray()) {
        // Error: node is not an array at position ...
      }
      node = &((*node)[arg.index_]);
    } else if (arg.kind_ == PathArgument::kindKey) {
      if (!node->isObject()) {
        // Error: node is not an object at position...
      }
      node = &((*node)[arg.key_]);
    }
  }
  return *node;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_value.cpp
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_writer.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at https://github.com/open-source-parsers/jsoncpp/blob/master/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/writer.h>
# include "json_tool.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <utility>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#if _MSC_VER >= 1400 // VC++ 8.0
#pragma warning( disable : 4996 )   // disable warning about strdup being deprecated.
#endif

namespace Json {

static bool containsControlCharacter(const char* str) {
  while (*str) {
    if (isControlCharacter(*(str++)))
      return true;
  }
  return false;
}


std::string valueToString(LargestInt value) {
  UIntToStringBuffer buffer;
  char* current = buffer + sizeof(buffer);
  bool isNegative = value < 0;
  if (isNegative)
    value = -value;
  uintToString(LargestUInt(value), current);
  if (isNegative)
    *--current = '-';
  assert(current >= buffer);
  return current;
}


std::string valueToString(LargestUInt value) {
  UIntToStringBuffer buffer;
  char* current = buffer + sizeof(buffer);
  uintToString(value, current);
  assert(current >= buffer);
  return current;
}

#if defined(JSON_HAS_INT64)

std::string valueToString(Int value) {
  return valueToString(LargestInt(value));
}


std::string valueToString(UInt value) {
  return valueToString(LargestUInt(value));
}

#endif // # if defined(JSON_HAS_INT64)


std::string valueToString(double value) {
  char buffer[32];
#if defined(_MSC_VER) && defined(__STDC_SECURE_LIB__) // Use secure version with visual studio 2005 to avoid warning.
  sprintf_s(buffer, sizeof(buffer), "%#.16g", value);
#else
  sprintf(buffer, "%#.16g", value);
#endif
  char* ch = buffer + strlen(buffer) - 1;
  if (*ch != '0')
    return buffer; // nothing to truncate, so save time
  while (ch > buffer && *ch == '0') {
    --ch;
  }
  char* last_nonzero = ch;
  while (ch >= buffer) {
    switch (*ch) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        --ch;
        continue;
      case '.':
        // Truncate zeroes to save bytes in output, but keep one.
        *(last_nonzero + 2) = '\0';
        return buffer;
      default:
        return buffer;
    }
  }
  return buffer;
}


std::string valueToString(bool value) {
  return value ? "true" : "false";
}

std::string valueToQuotedString(const char* value) {
  // Not sure how to handle unicode...
  if (strpbrk(value, "\"\\\b\f\n\r\t") == NULL && !containsControlCharacter(value))
    return std::string("\"") + value + "\"";
  // We have to walk value and escape any special characters.
  // Appending to std::string is not efficient, but this should be rare.
  // (Note: forward slashes are *not* rare, but I am not escaping them.)
  std::string::size_type maxsize = strlen(value) * 2 + 3; // allescaped+quotes+NULL
  std::string result;
  result.reserve(maxsize); // to avoid lots of mallocs
  result += "\"";
  for (const char* c = value; *c != 0; ++c) {
    switch (*c) {
      case '\"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
        break;
      case '\b':
        result += "\\b";
        break;
      case '\f':
        result += "\\f";
        break;
      case '\n':
        result += "\\n";
        break;
      case '\r':
        result += "\\r";
        break;
      case '\t':
        result += "\\t";
        break;
      //case '/':
      // Even though \/ is considered a legal escape in JSON, a bare
      // slash is also legal, so I see no reason to escape it.
      // (I hope I am not misunderstanding something.
      // blep notes: actually escaping \/ may be useful in javascript to avoid </
      // sequence.
      // Should add a flag to allow this compatibility mode and prevent this
      // sequence from occurring.
      default:
        if (isControlCharacter(*c)) {
          std::ostringstream oss;
          oss << "\\u" << std::hex << std::uppercase << std::setfill('0') << std::setw(4) << static_cast<int>(*c);
          result += oss.str();
        } else {
          result += *c;
        }
        break;
    }
  }
  result += "\"";
  return result;
}

// Class Writer
// //////////////////////////////////////////////////////////////////
Writer::~Writer() {
}


// Class FastWriter
// //////////////////////////////////////////////////////////////////

FastWriter::FastWriter()
  : yamlCompatiblityEnabled_(false) {
}


void
FastWriter::enableYAMLCompatibility() {
  yamlCompatiblityEnabled_ = true;
}


std::string
FastWriter::write(const Value& root) {
  document_ = "";
  writeValue(root);
  document_ += "\n";
  return document_;
}


void
FastWriter::writeValue(const Value& value) {
  switch (value.type()) {
    case nullValue:
      document_ += "null";
      break;
    case intValue:
      document_ += valueToString(value.asLargestInt());
      break;
    case uintValue:
      document_ += valueToString(value.asLargestUInt());
      break;
    case realValue:
      document_ += valueToString(value.asDouble());
      break;
    case stringValue:
      document_ += valueToQuotedString(value.asCString());
      break;
    case booleanValue:
      document_ += valueToString(value.asBool());
      break;
    case arrayValue: {
      document_ += "[";
      int size = value.size();
      for (int index = 0; index < size; ++index) {
        if (index > 0)
          document_ += ",";
        writeValue(value[index]);
      }
      document_ += "]";
    }
    break;
    case objectValue: {
      Value::Members members(value.getMemberNames());
      document_ += "{";
      for (Value::Members::iterator it = members.begin();
           it != members.end();
           ++it) {
        const std::string& name = *it;
        if (it != members.begin())
          document_ += ",";
        document_ += valueToQuotedString(name.c_str());
        document_ += yamlCompatiblityEnabled_ ? ": "
                     : ":";
        writeValue(value[name]);
      }
      document_ += "}";
    }
    break;
  }
}


// Class StyledWriter
// //////////////////////////////////////////////////////////////////

StyledWriter::StyledWriter()
  : rightMargin_(74)
  , indentSize_(3) {
}


std::string
StyledWriter::write(const Value& root) {
  document_ = "";
  addChildValues_ = false;
  indentString_ = "";
  writeCommentBeforeValue(root);
  writeValue(root);
  writeCommentAfterValueOnSameLine(root);
  document_ += "\n";
  return document_;
}


void
StyledWriter::writeValue(const Value& value) {
  switch (value.type()) {
    case nullValue:
      pushValue("null");
      break;
    case intValue:
      pushValue(valueToString(value.asLargestInt()));
      break;
    case uintValue:
      pushValue(valueToString(value.asLargestUInt()));
      break;
    case realValue:
      pushValue(valueToString(value.asDouble()));
      break;
    case stringValue:
      pushValue(valueToQuotedString(value.asCString()));
      break;
    case booleanValue:
      pushValue(valueToString(value.asBool()));
      break;
    case arrayValue:
      writeArrayValue(value);
      break;
    case objectValue: {
      Value::Members members(value.getMemberNames());
      if (members.empty())
        pushValue("{}");
      else {
        writeWithIndent("{");
        indent();
        Value::Members::iterator it = members.begin();
        for (;;) {
          const std::string& name = *it;
          const Value& childValue = value[name];
          writeCommentBeforeValue(childValue);
          writeWithIndent(valueToQuotedString(name.c_str()));
          document_ += " : ";
          writeValue(childValue);
          if (++it == members.end()) {
            writeCommentAfterValueOnSameLine(childValue);
            break;
          }
          document_ += ",";
          writeCommentAfterValueOnSameLine(childValue);
        }
        unindent();
        writeWithIndent("}");
      }
    }
    break;
  }
}


void
StyledWriter::writeArrayValue(const Value& value) {
  unsigned size = value.size();
  if (size == 0)
    pushValue("[]");
  else {
    bool isArrayMultiLine = isMultineArray(value);
    if (isArrayMultiLine) {
      writeWithIndent("[");
      indent();
      bool hasChildValue = !childValues_.empty();
      unsigned index = 0;
      for (;;) {
        const Value& childValue = value[index];
        writeCommentBeforeValue(childValue);
        if (hasChildValue)
          writeWithIndent(childValues_[index]);
        else {
          writeIndent();
          writeValue(childValue);
        }
        if (++index == size) {
          writeCommentAfterValueOnSameLine(childValue);
          break;
        }
        document_ += ",";
        writeCommentAfterValueOnSameLine(childValue);
      }
      unindent();
      writeWithIndent("]");
    } else { // output on a single line
      assert(childValues_.size() == size);
      document_ += "[ ";
      for (unsigned index = 0; index < size; ++index) {
        if (index > 0)
          document_ += ", ";
        document_ += childValues_[index];
      }
      document_ += " ]";
    }
  }
}


bool
StyledWriter::isMultineArray(const Value& value) {
  int size = value.size();
  bool isMultiLine = size * 3 >= rightMargin_ ;
  childValues_.clear();
  for (int index = 0; index < size  &&  !isMultiLine; ++index) {
    const Value& childValue = value[index];
    isMultiLine = isMultiLine  ||
                  ((childValue.isArray()  ||  childValue.isObject())  &&
                   childValue.size() > 0);
  }
  if (!isMultiLine) {   // check if line length > max line length
    childValues_.reserve(size);
    addChildValues_ = true;
    int lineLength = 4 + (size - 1) * 2; // '[ ' + ', '*n + ' ]'
    for (int index = 0; index < size  &&  !isMultiLine; ++index) {
      writeValue(value[index]);
      lineLength += int(childValues_[index].length());
      isMultiLine = isMultiLine  &&  hasCommentForValue(value[index]);
    }
    addChildValues_ = false;
    isMultiLine = isMultiLine  ||  lineLength >= rightMargin_;
  }
  return isMultiLine;
}


void
StyledWriter::pushValue(const std::string& value) {
  if (addChildValues_)
    childValues_.push_back(value);
  else
    document_ += value;
}


void
StyledWriter::writeIndent() {
  if (!document_.empty()) {
    char last = document_[document_.length() - 1];
    if (last == ' ')       // already indented
      return;
    if (last != '\n')      // Comments may add new-line
      document_ += '\n';
  }
  document_ += indentString_;
}


void
StyledWriter::writeWithIndent(const std::string& value) {
  writeIndent();
  document_ += value;
}


void
StyledWriter::indent() {
  indentString_ += std::string(indentSize_, ' ');
}


void
StyledWriter::unindent() {
  assert(int(indentString_.size()) >= indentSize_);
  indentString_.resize(indentString_.size() - indentSize_);
}


void
StyledWriter::writeCommentBeforeValue(const Value& root) {
  if (!root.hasComment(commentBefore))
    return;
  document_ += normalizeEOL(root.getComment(commentBefore));
  document_ += "\n";
}


void
StyledWriter::writeCommentAfterValueOnSameLine(const Value& root) {
  if (root.hasComment(commentAfterOnSameLine))
    document_ += " " + normalizeEOL(root.getComment(commentAfterOnSameLine));

  if (root.hasComment(commentAfter)) {
    document_ += "\n";
    document_ += normalizeEOL(root.getComment(commentAfter));
    document_ += "\n";
  }
}


bool
StyledWriter::hasCommentForValue(const Value& value) {
  return value.hasComment(commentBefore)
         ||  value.hasComment(commentAfterOnSameLine)
         ||  value.hasComment(commentAfter);
}


std::string
StyledWriter::normalizeEOL(const std::string& text) {
  std::string normalized;
  normalized.reserve(text.length());
  const char* begin = text.c_str();
  const char* end = begin + text.length();
  const char* current = begin;
  while (current != end) {
    char c = *current++;
    if (c == '\r') {   // mac or dos EOL
      if (*current == '\n')   // convert dos EOL
        ++current;
      normalized += '\n';
    } else // handle unix EOL & other char
      normalized += c;
  }
  return normalized;
}


// Class StyledStreamWriter
// //////////////////////////////////////////////////////////////////

StyledStreamWriter::StyledStreamWriter(std::string indentation)
  : document_(NULL)
  , rightMargin_(74)
  , indentation_(indentation) {
}


void
StyledStreamWriter::write(std::ostream& out, const Value& root) {
  document_ = &out;
  addChildValues_ = false;
  indentString_ = "";
  writeCommentBeforeValue(root);
  writeValue(root);
  writeCommentAfterValueOnSameLine(root);
  *document_ << "\n";
  document_ = NULL; // Forget the stream, for safety.
}


void
StyledStreamWriter::writeValue(const Value& value) {
  switch (value.type()) {
    case nullValue:
      pushValue("null");
      break;
    case intValue:
      pushValue(valueToString(value.asLargestInt()));
      break;
    case uintValue:
      pushValue(valueToString(value.asLargestUInt()));
      break;
    case realValue:
      pushValue(valueToString(value.asDouble()));
      break;
    case stringValue:
      pushValue(valueToQuotedString(value.asCString()));
      break;
    case booleanValue:
      pushValue(valueToString(value.asBool()));
      break;
    case arrayValue:
      writeArrayValue(value);
      break;
    case objectValue: {
      Value::Members members(value.getMemberNames());
      if (members.empty())
        pushValue("{}");
      else {
        writeWithIndent("{");
        indent();
        Value::Members::iterator it = members.begin();
        for (;;) {
          const std::string& name = *it;
          const Value& childValue = value[name];
          writeCommentBeforeValue(childValue);
          writeWithIndent(valueToQuotedString(name.c_str()));
          *document_ << " : ";
          writeValue(childValue);
          if (++it == members.end()) {
            writeCommentAfterValueOnSameLine(childValue);
            break;
          }
          *document_ << ",";
          writeCommentAfterValueOnSameLine(childValue);
        }
        unindent();
        writeWithIndent("}");
      }
    }
    break;
  }
}


void
StyledStreamWriter::writeArrayValue(const Value& value) {
  unsigned size = value.size();
  if (size == 0)
    pushValue("[]");
  else {
    bool isArrayMultiLine = isMultineArray(value);
    if (isArrayMultiLine) {
      writeWithIndent("[");
      indent();
      bool hasChildValue = !childValues_.empty();
      unsigned index = 0;
      for (;;) {
        const Value& childValue = value[index];
        writeCommentBeforeValue(childValue);
        if (hasChildValue)
          writeWithIndent(childValues_[index]);
        else {
          writeIndent();
          writeValue(childValue);
        }
        if (++index == size) {
          writeCommentAfterValueOnSameLine(childValue);
          break;
        }
        *document_ << ",";
        writeCommentAfterValueOnSameLine(childValue);
      }
      unindent();
      writeWithIndent("]");
    } else { // output on a single line
      assert(childValues_.size() == size);
      *document_ << "[ ";
      for (unsigned index = 0; index < size; ++index) {
        if (index > 0)
          *document_ << ", ";
        *document_ << childValues_[index];
      }
      *document_ << " ]";
    }
  }
}


bool
StyledStreamWriter::isMultineArray(const Value& value) {
  int size = value.size();
  bool isMultiLine = size * 3 >= rightMargin_ ;
  childValues_.clear();
  for (int index = 0; index < size  &&  !isMultiLine; ++index) {
    const Value& childValue = value[index];
    isMultiLine = isMultiLine  ||
                  ((childValue.isArray()  ||  childValue.isObject())  &&
                   childValue.size() > 0);
  }
  if (!isMultiLine) {   // check if line length > max line length
    childValues_.reserve(size);
    addChildValues_ = true;
    int lineLength = 4 + (size - 1) * 2; // '[ ' + ', '*n + ' ]'
    for (int index = 0; index < size  &&  !isMultiLine; ++index) {
      writeValue(value[index]);
      lineLength += int(childValues_[index].length());
      isMultiLine = isMultiLine  &&  hasCommentForValue(value[index]);
    }
    addChildValues_ = false;
    isMultiLine = isMultiLine  ||  lineLength >= rightMargin_;
  }
  return isMultiLine;
}


void
StyledStreamWriter::pushValue(const std::string& value) {
  if (addChildValues_)
    childValues_.push_back(value);
  else
    *document_ << value;
}


void
StyledStreamWriter::writeIndent() {
  /*
    Some comments in this method would have been nice. ;-)

   if ( !document_.empty() )
   {
      char last = document_[document_.length()-1];
      if ( last == ' ' )     // already indented
         return;
      if ( last != '\n' )    // Comments may add new-line
         *document_ << '\n';
   }
  */
  *document_ << '\n' << indentString_;
}


void
StyledStreamWriter::writeWithIndent(const std::string& value) {
  writeIndent();
  *document_ << value;
}


void
StyledStreamWriter::indent() {
  indentString_ += indentation_;
}


void
StyledStreamWriter::unindent() {
  assert(indentString_.size() >= indentation_.size());
  indentString_.resize(indentString_.size() - indentation_.size());
}


void
StyledStreamWriter::writeCommentBeforeValue(const Value& root) {
  if (!root.hasComment(commentBefore))
    return;
  *document_ << normalizeEOL(root.getComment(commentBefore));
  *document_ << "\n";
}


void
StyledStreamWriter::writeCommentAfterValueOnSameLine(const Value& root) {
  if (root.hasComment(commentAfterOnSameLine))
    *document_ << " " + normalizeEOL(root.getComment(commentAfterOnSameLine));

  if (root.hasComment(commentAfter)) {
    *document_ << "\n";
    *document_ << normalizeEOL(root.getComment(commentAfter));
    *document_ << "\n";
  }
}


bool
StyledStreamWriter::hasCommentForValue(const Value& value) {
  return value.hasComment(commentBefore)
         ||  value.hasComment(commentAfterOnSameLine)
         ||  value.hasComment(commentAfter);
}


std::string
StyledStreamWriter::normalizeEOL(const std::string& text) {
  std::string normalized;
  normalized.reserve(text.length());
  const char* begin = text.c_str();
  const char* end = begin + text.length();
  const char* current = begin;
  while (current != end) {
    char c = *current++;
    if (c == '\r') {   // mac or dos EOL
      if (*current == '\n')   // convert dos EOL
        ++current;
      normalized += '\n';
    } else // handle unix EOL & other char
      normalized += c;
  }
  return normalized;
}


std::ostream& operator<<(std::ostream& sout, const Value& root) {
  Json::StyledStreamWriter writer;
  writer.write(sout, root);
  return sout;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_writer.cpp
// //////////////////////////////////////////////////////////////////////





//
// end of src/jsoncpp.cpp
//


//
// start of src/jsoncustomwriter.cpp
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

#ifndef PYNE_IS_AMALGAMATED
#include "json.h"
#include "jsoncustomwriter.h"
#endif

namespace Json {

CustomWriter::CustomWriter(std::string opencurly,
                           std::string closecurly,
                           std::string opensquare,
                           std::string closesquare,
                           std::string colon,
                           std::string comma,
                           std::string indent,
                           int maxWidth)
  : opencurly_(opencurly)
  , closecurly_(closecurly)
  , opensquare_(opensquare)
  , closesquare_(closesquare)
  , colon_(colon)
  , comma_(comma)
  , indent_(indent)
  , maxWidth_(maxWidth) {
}


std::string
CustomWriter::write(const Value& root) {
  document_ = "";
  indentString_ = "";
  writeValue(root, document_, false);
  document_ += "\n";
  return document_;
}


void
CustomWriter::writeValue(const Value& value, std::string& doc, bool forceSingleLine) {
  switch (value.type()) {
    case nullValue:
      doc += "null";
      break;
    case intValue:
      doc += valueToString(value.asLargestInt());
      break;
    case uintValue:
      doc += valueToString(value.asLargestUInt());
      break;
    case realValue:
      doc += valueToString(value.asDouble());
      break;
    case stringValue:
      doc += valueToQuotedString(value.asCString());
      break;
    case booleanValue:
      doc += valueToString(value.asBool());
      break;
    case arrayValue: {
      bool isMulti = false;
      if (!forceSingleLine) {
        std::string valLine = "";
        writeValue(value, valLine, true);
        if (valLine.length() > maxWidth_) {
          isMulti = true;
        } else {
          doc += valLine;
          break;
        }
      }
      doc += opensquare_;
      if (isMulti)
        indent();
      for (int index = 0; index < value.size(); ++index) {
        if (isMulti) {
          doc += "\n";
          doc += indentString_;
        }
        writeValue(value[index], doc, false);
        if (index < value.size() - 1)
          doc += comma_;
      }
      if (isMulti) {
        unindent();
        doc += "\n";
        doc += indentString_;
      }
      doc += closesquare_;
    }
    break;
    case objectValue: {
      bool isMulti = false;
      if (!forceSingleLine) {
        std::string valLine = "";
        writeValue(value, valLine, true);
        if (valLine.length() > maxWidth_) {
          isMulti = true;
        } else {
          doc += valLine;
          break;
        }
      }
      Value::Members members(value.getMemberNames());
      doc += opencurly_;
      if (isMulti)
        indent();
      for (Value::Members::iterator it = members.begin();
           it != members.end();
           ++it) {
        if (isMulti) {
          doc += "\n";
          doc += indentString_;

        }
        const std::string& name = *it;
        doc += valueToQuotedString(name.c_str());
        doc += colon_;
        writeValue(value[name], doc, forceSingleLine);
        if (!(it + 1 == members.end()))
          doc += comma_;
      }
      if (isMulti) {
        unindent();
        doc += "\n";
        doc += indentString_;
      }
      doc += closecurly_;
    }
    break;
  }
}


void
CustomWriter::indent() {
  indentString_ += indent_;
}


void
CustomWriter::unindent() {
  int idSize = int(indent_.size());
  int idsSize = int(indentString_.size());
  if (idsSize >= idSize)
    indentString_.resize(idsSize - idSize);
}

}
//
// end of src/jsoncustomwriter.cpp
//


//
// start of src/material.cpp
//
// Material.cpp
// The very central Material class
// -- Anthony Scopatz

#include <string>
#include <vector>
#include <iomanip>  // std::setprecision
#include <math.h>   // modf
#include <stdexcept>

#ifndef PYNE_IS_AMALGAMATED
#include "transmuters.h"
#include "material.h"
#endif

// h5wrap template
template double h5wrap::get_array_index(hid_t, int, hid_t);
const int mcnp_line_length = 79;


/***************************/
/*** Protected Functions ***/
/***************************/

double pyne::Material::get_comp_sum() {
  // Sums the weights in the composition dictionary
  double sum = 0.0;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    sum = sum + i->second;
  }
  return sum;
}


void pyne::Material::norm_comp() {
  double sum = get_comp_sum();
  if (sum != 1.0 && sum != 0.0) {
    for (comp_iter i = comp.begin(); i != comp.end(); i++)
      i->second = i->second / sum;
  }

  if (mass < 0.0)
    mass = sum;
}


void pyne::Material::_load_comp_protocol0(hid_t db, std::string datapath, int row) {
  // Clear current content
  comp.clear();

  hid_t matgroup = H5Gopen2(db, datapath.c_str(), H5P_DEFAULT);
  hid_t nucset;
  double nucvalue;
  ssize_t nuckeylen;
  std::string nuckey;

  // get the number of members in the material group
  H5G_info_t group_info;
  H5Gget_info(matgroup, &group_info);
  hsize_t matG = group_info.nlinks;

  // Iterate over datasets in the group.
  for (int matg = 0; matg < matG; matg++) {
    nuckeylen = 1 + H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                       NULL, 0, H5P_DEFAULT);
    char* nkey = new char[nuckeylen];
    nuckeylen = H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                   nkey, nuckeylen, H5P_DEFAULT);
    nuckey = nkey;
    nucset = H5Dopen2(matgroup, nkey, H5P_DEFAULT);
    nucvalue = h5wrap::get_array_index<double>(nucset, row);

    if (nuckey == "Mass" || nuckey == "MASS" || nuckey == "mass")
      mass = nucvalue;
    else
      comp[pyne::nucname::id(nuckey)] = nucvalue;

    H5Dclose(nucset);
    delete[] nkey;
  }

  // Set meta data
  atoms_per_molecule = -1.0;
}


void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath,
                                          int row) {
  std::string nucpath;
  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Grab the nucpath
  if (detect_nuclidelist(data_set, nucpath)) {
    H5Dclose(data_set);
    _load_comp_protocol1(db, datapath, nucpath, row);
  } else {
    H5Dclose(data_set);
    throw std::runtime_error("Can't find location of the nuclide list: nucpath attribute not found!");
  }
}


void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath,
                                          std::string nucpath, int row) {
  // Clear current content
  comp.clear();

  if (!h5wrap::path_exists(db, nucpath))
    throw std::runtime_error("No path found at the location: " + nucpath);

  if (!h5wrap::path_exists(db, datapath))
    throw std::runtime_error("No path found at the location: " + datapath);

  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  hsize_t data_offset[1] = {static_cast<hsize_t>(row)};
  if (row < 0) {
    // Handle negative row indices
    hid_t data_space = H5Dget_space(data_set);
    hsize_t data_dims[1];
    H5Sget_simple_extent_dims(data_space, data_dims, NULL);
    data_offset[0] += data_dims[0];
  }

  // Grab the nuclides
  std::vector<int> nuclides =
      h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
  int nuc_size = nuclides.size();
  hsize_t nuc_dims[1] = {static_cast<hsize_t>(nuc_size)};

  // Get the data hyperslab
  hid_t data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL,
                      data_count, NULL);

  // Get memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, NULL);

  // Get material type
  size_t material_data_size =
      sizeof(pyne::material_data) + sizeof(double) * (nuc_size - 1);
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type =
      H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule",
            HOFFSET(pyne::material_data, atoms_per_mol), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp),
            comp_values_array_type);

  // make the data array, have to over-allocate
  material_data* mat_data = new material_data[material_data_size];

  // Finally, get data and put in on this instance
  H5Dread(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  mass = (*mat_data).mass;
  density = (*mat_data).density;
  atoms_per_molecule = (*mat_data).atoms_per_mol;
  for (int i = 0; i < nuc_size; i++) {
    if ((double)(*mat_data).comp[i] != 0) {
      comp[nuclides[i]] = (double)(*mat_data).comp[i];
    }
  }
  delete[] mat_data;

  //
  // Get metadata from associated dataset, if available
  //
  std::string attrpath = datapath + "_metadata";
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (!attrpath_exists)
    return;

  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;
  hvl_t attrdata[1];

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // Get the metadata from the file
  metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL,
                      data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, NULL);
  H5Dread(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT,
          attrdata);

  // convert to in-memory JSON
  Json::Reader reader;
  reader.parse((char*)attrdata[0].p, (char*)attrdata[0].p + attrdata[0].len,
               metadata, false);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);
}


void pyne::Material::from_hdf5(char* filename, char* datapath, int row, int protocol) {
  std::string fname(filename);
  std::string dpath(datapath);
  from_hdf5(fname, dpath, row, protocol);
}


int pyne::Material::detect_hdf5_layout(hid_t db, std::string path) {
  // Check hdf5 material layout:
  // return options are:
  //     -"-1": path and "/material" do not exist
  //     - "0": path and/or "/material" exist but either as a group or a dataset
  //     - "1": path exists as a dataset -> old layout
  //     - "2": "/material" exists as a group-> new layout

  // Initialize test variables
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Test if datapath exist as a non-dataset
  H5O_info_t path_info;
  status = H5Oget_info_by_name(db, path.c_str(), &path_info, H5P_DEFAULT);
  bool path_exists = (status == 0);

  // Reset status and test "/material" path exist as a non-dataset
  H5O_info_t matpath_info;
  status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  status = H5Oget_info_by_name(db, "/material", &matpath_info, H5P_DEFAULT);
  bool matpath_exists = (status == 0);
  if (!matpath_exists && !path_exists) {
    return prot1_layout::path_donotexists;
  } else if (path_info.type == H5O_TYPE_DATASET || matpath_info.type == H5O_TYPE_DATASET) {
    return prot1_layout::old_layout;
  } else if (matpath_info.type == H5O_TYPE_GROUP) {
    return prot1_layout::new_layout;
  } else {
    return prot1_layout::unknown;
  }
}


void pyne::Material::from_hdf5(std::string filename, std::string datapath, int row, int protocol) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);
  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(filename);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  // Clear current content
  comp.clear();
  // Load via various protocols
  if (protocol == 0) {
    bool datapath_exists = h5wrap::path_exists(db, datapath);
    if (!datapath_exists)
      throw h5wrap::PathNotFound(filename, datapath);
    _load_comp_protocol0(db, datapath, row);
  } else if (protocol == 1) {

    int prot1_hdf5_layout = detect_hdf5_layout(db, datapath);
    switch (prot1_hdf5_layout) {
      case prot1_layout::path_donotexists:
        throw std::runtime_error("/material and " + datapath + " paths do not exist.");
        break;

      case prot1_layout::unknown:
        throw std::runtime_error(datapath + " is not a dataset and /material entity is not a group nor a dataset.");
        break;

      case old_layout:
        _load_comp_protocol1(db, datapath, row);
        break;

      case prot1_layout::new_layout:
        std::string full_datapath = "/material" + datapath + "/composition";
        std::string nucpath = "/material" + datapath + "/nuclidelist";
        _load_comp_protocol1(db, full_datapath, nucpath, row);
        break;
    }
  } else
    throw pyne::MaterialProtocolError();

  // Close the database
  status = H5Fclose(db);

  // Renormalize the composition, just to be safe.
  norm_comp();
}


void pyne::Material::deprecated_write_hdf5(char* filename, char* datapath, char* nucpath, float row, int chunksize) {
  std::string fname(filename);
  std::string groupname(datapath);
  std::string nuclist(nucpath);
  deprecated_write_hdf5(fname, groupname, nuclist, row, chunksize);
}


std::vector<int> pyne::Material::write_hdf5_nucpath(hid_t db, std::string nucpath) {
  //
  // Read in nuclist if available, write it out if not
  //
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  std::vector<int> nuclides;
  int nuc_size;
  hsize_t nuc_dims[1];

  // if nucpath exist: get it and check against the one we have!
  if (nucpath_exists) {
    nuclides =
        h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
    nuc_size = nuclides.size();
    nuc_dims[0] = nuc_size;
    bool missing_nucs = false;
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
      missing_nucs |= !std::binary_search(nuclides.begin(), nuclides.end(), i->first);
      if (missing_nucs) {
        break;
      }
    }
    if (missing_nucs)
      std::cout
          << "One or more nuclides are missing from the existing nuclides "
          "list, material will likely not be written correctly."
          << std::endl;

  } else {
    nuclides = std::vector<int>();
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
      nuclides.push_back(i->first);
    nuc_size = nuclides.size();

    // Create the data if it doesn't exist
    int nuc_data[nuc_size];
    for (int n = 0; n != nuc_size; n++)
      nuc_data[n] = nuclides[n];
    nuc_dims[0] = nuc_size;
    hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
    hid_t nuc_set = H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuc_space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nuc_data);
    H5Fflush(db, H5F_SCOPE_GLOBAL);

    H5Dclose(nuc_set);
  }
  return nuclides;
}


void pyne::Material::write_hdf5_datapath(hid_t db, std::string datapath,
                                         float row, int chunksize,
                                         std::vector<int> nuclides) {
  int row_num = (int)row;

  hsize_t nuc_dims[1];
  int nuc_size = nuc_dims[0] = nuclides.size();
  //
  // Write out the data itself to the file
  //
  hid_t data_set, data_space, data_hyperslab;
  int data_rank = 1;
  hsize_t data_dims[1] = {1};
  hsize_t data_max_dims[1] = {H5S_UNLIMITED};
  hsize_t data_offset[1] = {0};

  size_t material_data_size = sizeof(pyne::material_data) + sizeof(double) * (nuc_size - 1);
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule", HOFFSET(pyne::material_data, atoms_per_mol),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp),
            comp_values_array_type);

  material_data* mat_data  = new material_data[material_data_size];
  (*mat_data).mass = mass;
  (*mat_data).density = density;
  (*mat_data).atoms_per_mol = atoms_per_molecule;
  for (int n = 0; n != nuc_size; n++) {
    if (0 < comp.count(nuclides[n]))
      (*mat_data).comp[n] = comp[nuclides[n]];
    else
      (*mat_data).comp[n] = 0.0;
  }

  // get / make the data set
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (datapath_exists) {
    data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    data_space = H5Dget_space(data_set);
    data_rank = H5Sget_simple_extent_dims(data_space, data_dims, data_max_dims);

    // Determine the row size.
    if (std::signbit(row))
      row_num = data_dims[0] + row;  // careful, row is negative

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(data_set, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    // Get full space
    data_space = H5Screate_simple(1, data_dims, data_max_dims);

    // Make data set properties to enable chunking
    hid_t data_set_params = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] = {static_cast<hsize_t>(chunksize)};
    H5Pset_chunk(data_set_params, 1, chunk_dims);
    H5Pset_deflate(data_set_params, 1);

    // Create the data set
    data_set = H5Dcreate2(db, datapath.c_str(), desc, data_space, H5P_DEFAULT,
                          data_set_params, H5P_DEFAULT);
    H5Dset_extent(data_set, data_dims);

    H5Fflush(db, H5F_SCOPE_GLOBAL);
  }

  // Get the data hyperslab
  data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

  // Get a memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, data_max_dims);

  // Write the row...
  H5Dwrite(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  // Close out the Dataset
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(data_set);
  H5Sclose(data_space);
  H5Tclose(desc);

  //
  // Write out the metadata to the file
  //
  std::string attrpath = datapath + "_metadata";
  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // get / make the data set
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (attrpath_exists) {
    metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
    metadatapace = H5Dget_space(metadataet);
    attrrank = H5Sget_simple_extent_dims(metadatapace, data_dims, data_max_dims);

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(metadataet, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    hid_t metadataetparams;
    hsize_t attrchunkdims [1];

    // Make data set properties to enable chunking
    metadataetparams = H5Pcreate(H5P_DATASET_CREATE);
    attrchunkdims[0] = chunksize;
    H5Pset_chunk(metadataetparams, 1, attrchunkdims);
    H5Pset_deflate(metadataetparams, 1);

    hvl_t attrfillvalue [1];
    attrfillvalue[0].len = 3;
    attrfillvalue[0].p = (char*) "{}\n";
    H5Pset_fill_value(metadataetparams, attrtype, &attrfillvalue);

    // make dataset
    metadatapace = H5Screate_simple(1, data_dims, data_max_dims);
    metadataet = H5Dcreate2(db, attrpath.c_str(), attrtype, metadatapace,
                            H5P_DEFAULT, metadataetparams, H5P_DEFAULT);
    H5Dset_extent(metadataet, data_dims);
    H5Pclose(metadataetparams);
  }

  // set the attr string
  hvl_t attrdata [1];
  Json::FastWriter writer;
  std::string metadatatr = writer.write(metadata);
  attrdata[0].p = (char*) metadatatr.c_str();
  attrdata[0].len = metadatatr.length();

  // write the attr
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, data_max_dims);
  H5Dwrite(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT, attrdata);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);

  delete[] mat_data;
}


void pyne::Material::write_hdf5(std::string filename, std::string datapath,
                                float row, int chunksize) {
  if (datapath.front() != '/')
    datapath = '/' + datapath;

  hid_t material_grp_id;  // Holder of HDF5 Id of the "/material" group
  hid_t data_id;  // Holder of HDF5 Id of the data group to write the data
  // (located in "/material/datapath")

  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  // Create new/open datafile.

  // This complicated algorithm is required to allow backward compatibility with
  // previous version of write_hdf5 (where data were written in a hdf5 DATASET)
  // FILE EXIST ?
  //    NO -> create file with a "/material" group, write the data in it
  //    YES:
  //      - DATAPTH exist:
  //        - YES && != /material -> Detect NUCPATH + old write_hdf5
  //        - NO -> CONTINUE
  //          - "/material" EXIST:
  //            - NO -> create "/material" group and write everything in it
  //            - YES:
  //              - "/material" is a DATASET -> detect NUCPATH + old write_hdf5
  //              - "/material" is a GROUP -> add data in it
  //              - "/material" isn't a GROUP nor a DATASET -> fail and complain

  hid_t db;
  if (!pyne::file_exists(filename)) {
    // Create the file
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  } else {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  }
  int prot1_hdf5_layout = detect_hdf5_layout(db, datapath);
  switch (prot1_hdf5_layout) {
    case prot1_layout::unknown: {
      throw std::runtime_error(
          datapath +
          " is not a dataset and /material entity is neither a group nor a dataset.");
      break;
    }

    // old layout
    case prot1_layout::old_layout: {
      std::string nucpath = "/nucid";
      hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

      if (h5wrap::path_exists(db, datapath)) {
        // if path exists grab the nucpath
        bool nucpath_detetcted = detect_nuclidelist(data_set, nucpath);
        if (!nucpath_detetcted) {  // can't find a valid nuclide list path
          // from datapath... fail
          throw std::runtime_error(
              "Can't find the nuclide list path in the existing datapath. "
              "Can't add your material to the datapath.");
        }
      }
      deprecated_write_hdf5(db, datapath, nucpath, row, chunksize);
      H5Fclose(db);
      return;
      break;
    }

    // no "/material or matname in the hdf5 file -> build the new layout
    case prot1_layout::path_donotexists: {
      material_grp_id =
          H5Gcreate2(db, "/material", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      break;
    }

    // /material as a group -> new layout: open the "/material" group
    case prot1_layout::new_layout: {
      material_grp_id = H5Gopen2(db, "/material", H5P_DEFAULT);
      break;
    }
  }

  //datapath is provided with a "/" so need to open with fullpath
  // Group "/material/datapath" does not exist create it
  if (!h5wrap::path_exists(db, "/material" + datapath)) {
    data_id = H5Gcreate2(db, ("/material" + datapath).c_str(), H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
  } else {
    data_id = H5Gopen2(db, ("/material" + datapath).c_str(), H5P_DEFAULT);
  }
  // write nuclide list
  std::string nucpath = "nuclidelist";
  std::vector<int> nuclides = write_hdf5_nucpath(data_id, nucpath);
  // write data
  std::string full_datapath = "composition";
  write_hdf5_datapath(data_id, full_datapath, row, chunksize, nuclides);

  // close all groups and files
  H5Gclose(data_id);
  H5Gclose(material_grp_id);
  H5Fclose(db);
}


void pyne::Material::deprecated_write_hdf5(std::string filename, std::string datapath,
                                           std::string nucpath, float row, int chunksize) {
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  // Create new/open datafile.
  hid_t db;
  if (pyne::file_exists(filename)) {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  } else
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  deprecated_write_hdf5(db, datapath, nucpath, row, chunksize);

  H5Fclose(db);
}


void pyne::Material::deprecated_write_hdf5(hid_t db, std::string datapath,
                                           std::string nucpath, float row, int chunksize) {
  int row_num = (int)row;

  //
  // Read in nuclist if available, write it out if not
  //
  std::vector<int> nuclides = write_hdf5_nucpath(db, nucpath);

  // Check if datapath already exist
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  // write mat composition in datapath
  write_hdf5_datapath(db, datapath, row, chunksize, nuclides);

  // if datapath has just been create register location of nucpath
  if (!datapath_exists) {
    hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    // Add attribute pointing to nuc path
    hid_t nuc_attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(nuc_attr_type, nucpath.length());
    hid_t nuc_attr_space = H5Screate(H5S_SCALAR);
    hid_t nuc_attr = H5Acreate2(data_set, "nucpath", nuc_attr_type,
                                nuc_attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(nuc_attr, nuc_attr_type, nucpath.c_str());
    H5Fflush(db, H5F_SCOPE_GLOBAL);
    H5Dclose(data_set);
  }
}


std::string pyne::Material::openmc(std::string frac_type) {
  std::ostringstream oss;

  std::set<int> carbon_set; carbon_set.insert(nucname::id("C"));
  pyne::Material temp_mat = this->expand_elements(carbon_set);

  // vars for consistency
  std::string new_quote = "\"";
  std::string end_quote = "\" ";
  std::string indent = "  ";

  // open the material element
  oss << "<material id=" ;

  // add the mat number
  if (temp_mat.metadata.isMember("mat_number")) {
    int mat_num = temp_mat.metadata["mat_number"].asInt();
    oss << new_quote << mat_num << end_quote;
  }
  // mat numbers are required for openmc
  else {
    throw pyne::ValueError("No material number found in metadata. This is not valid for use in OpenMC.");
    oss << new_quote << "?" << end_quote;
  }

  // add name if specified
  if (temp_mat.metadata.isMember("name")) {
    oss << "name=" << new_quote << temp_mat.metadata["name"].asString() << end_quote;
  }
  // close the material tag
  oss << ">";
  // new line
  oss << std::endl;

  //indent
  oss << indent;

  // specify density
  oss << "<density ";
  // if density is negtaive, report to user
  if (temp_mat.density < 0.0) {
    throw pyne::ValueError("A density < 0.0 was found. This is not valid for use in OpenMC.");
  }
  std::string density_str = std::to_string(temp_mat.density);
  // remove trailing zeros
  density_str.erase(density_str.find_last_not_of('0') + 1, std::string::npos);
  oss << "value=" <<  std::fixed << new_quote << density_str << end_quote;
  oss << "units=" << new_quote << "g/cc" << end_quote << "/>";
  // new line
  oss << std::endl;

  std::map<int, double> fracs;
  std::string frac_attrib;
  if (frac_type == "atom") {
    fracs = temp_mat.to_atom_frac();
    frac_attrib = "ao=";
  } else {
    fracs = temp_mat.comp;
    frac_attrib = "wo=";
  }

  // add nuclides
  for (comp_map::iterator f = fracs.begin(); f != fracs.end(); f++) {
    if (f->second == 0.0) {
      continue;
    }
    //indent
    oss << "  ";
    // start a new nuclide element
    oss << "<nuclide name=" << new_quote;
    oss << pyne::nucname::openmc(f->first);
    oss << end_quote;
    oss << frac_attrib;
    oss << std::setprecision(4) << std::scientific << new_quote << f->second << end_quote;
    oss << "/>";
    // new line
    oss << std::endl;
  }

  // other OpenMC material properties
  if (temp_mat.metadata.isMember("sab")) {
    oss << indent;
    oss << "<sab name=";
    oss << new_quote << temp_mat.metadata["sab"].asString() << end_quote;
    oss << "/>";
    oss << std::endl;
  }

  if (temp_mat.metadata.isMember("temperature")) {
    oss << indent;
    oss << "<temperature>";
    oss << new_quote << temp_mat.metadata["temperature"].asString() << end_quote;
    oss << "</temperature>";
    oss << std::endl;
  }

  if (temp_mat.metadata.isMember("macroscopic")) {
    oss << indent;
    oss << "<macroscopic name=";
    oss << new_quote << temp_mat.metadata["macroscropic"].asString() << end_quote;
    oss << "/>";
    oss << std::endl;
  }

  if (temp_mat.metadata.isMember("isotropic")) {
    oss << indent;
    oss << "<isotropic>";
    oss << new_quote << temp_mat.metadata["isotropic"].asString() << end_quote;
    oss << "</isotropic>";
    oss << std::endl;
  }

  // close the material node
  oss << "</material>" << std::endl;

  return oss.str();
}


///---------------------------------------------------------------------------//
std::string pyne::Material::get_uwuw_name() {
  // standard uwuw material name is : "mat:<Name of Material>/rho:<density>"
  if (! metadata.isMember("name")) {
    pyne::warning("The material has no name");
    return "";
  }
  std::ostringstream uwuw_name;
  uwuw_name << "mat:";
  uwuw_name << metadata["name"].asString();
  if (density > 0) {
    uwuw_name << "/rho:" << std::setprecision(5) << density;
  } else {
    pyne::warning("No Density defined for this Material");
  }

  return uwuw_name.str();
}


///---------------------------------------------------------------------------//
std::string pyne::Material::mcnp(std::string frac_type, bool mult_den) {
  //////////////////// Begin card creation ///////////////////////
  std::ostringstream oss;

  std::string comment_prefix = "C ";

  // 'name'
  if (metadata.isMember("name")) {
    oss << "C name: " << metadata["name"].asString() << std::endl;
  }
  // 'density'
  if (density != -1.0) {
    std::stringstream ds;
    ds << std::setprecision(5) << std::fixed << "C density = " << density << std::endl;
    oss << ds.str();
  }
  // 'source'
  if (metadata.isMember("source")) {
    oss << "C source: " << metadata["source"].asString() << std::endl;
  }
  // Metadata comments
  if (metadata.isMember("comments")) {
    std::string comment_string = "comments: " + metadata["comments"].asString();
    oss << pyne::comment_line_wrapping(comment_string, comment_prefix, mcnp_line_length);
  }

  // Metadata mat_num
  oss << "m";
  if (metadata.isMember("mat_number")) {
    int mat_num = metadata["mat_number"].asInt();
    oss << mat_num << std::endl;
  } else {
    oss << "?" << std::endl;
  }

  // Set up atom or mass frac map
  std::map<int, double> fracs = get_density_frac(frac_type, mult_den);
  std::string frac_sign = "";

  // write the frac map
  oss << mcnp_frac(fracs, frac_type);

  return oss.str();
}


///---------------------------------------------------------------------------//
std::string pyne::Material::phits(std::string frac_type, bool mult_den) {
  //////////////////// Begin card creation ///////////////////////
  std::ostringstream oss;

  std::string comment_prefix = "C ";

  // 'name'
  if (metadata.isMember("name")) {
    oss << "C name: " << metadata["name"].asString() << std::endl;
  }
  // Metadata comments
  if (metadata.isMember("comments")) {
    std::string comment_string = "comments: " + metadata["comments"].asString();
    oss << pyne::comment_line_wrapping(comment_string, comment_prefix, mcnp_line_length);
  }

  // Metadata mat_num
  oss << "M[ ";
  if (metadata.isMember("mat_number")) {
    int mat_num = metadata["mat_number"].asInt();
    oss << mat_num;
  } else {
    oss << "?";
  }
  oss << " ]" << std::endl;

  // check for metadata
  std::string keyworkds[7] = {"GAS", "ESTEP", "NLIB", "PLIB", "PNLIB", "ELIB", "HLIB"};
  for (auto keyword : keyworkds) {
    if (metadata.isMember(keyword)) {
      oss << "     " << keyword << "=" << metadata[keyword].asInt() << std::endl;
    }
  }
  // COND should be "<" or "=" or ">" if present
  if (metadata.isMember("COND")) {
    oss << "     COND" << metadata["COND"].asString() << "0" << std::endl;
  }

  // Set up atom or mass frac map
  std::map<int, double> fracs = get_density_frac(frac_type, mult_den);
  std::string frac_sign = "";

  // write the frac map
  oss << mcnp_frac(fracs, frac_type);

  return oss.str();
}


std::string pyne::Material::mcnp_frac(std::map<int, double> fracs, std::string frac_type) {
  std::string frac_sign = "";
  if ("atom" != frac_type) {
    frac_sign = "-";
  }

  // iterate through frac map
  // This is an awkward pre-C++11 way to put an int to a string
  std::ostringstream oss;
  for (pyne::comp_iter i = fracs.begin(); i != fracs.end(); ++i) {
    if (i->second > 0.0) {
      // Clear first
      std::stringstream ss;
      std::string nucmcnp;
      std::string table_item;
      ss << pyne::nucname::mcnp(i->first);
      nucmcnp = ss.str();

      int mcnp_id;
      mcnp_id = pyne::nucname::mcnp(i->first);
      // Spaces are important for tests
      table_item = metadata["table_ids"][nucmcnp].asString();
      if (!table_item.empty()) {
        oss << "     " << mcnp_id << "." << table_item << " ";
      } else {
        oss << "     " << mcnp_id << " ";
      }
      // The int needs a little formatting
      std::stringstream fs;
      fs << std::setprecision(4) << std::scientific << frac_sign << i->second << std::endl;
      oss << fs.str();
    }
  }
  return oss.str();
}


///---------------------------------------------------------------------------//
/// Create a set out of the static string array.
std::set<std::string> fluka_builtin(pyne::fluka_mat_strings,
                                    pyne::fluka_mat_strings + pyne::FLUKA_MAT_NUM);

///---------------------------------------------------------------------------//
/// not_fluka_builtin
///---------------------------------------------------------------------------//
/// Convenience function
/// This is written as a negative because that is what we care about
bool pyne::Material::not_fluka_builtin(std::string fluka_name) {
  return (fluka_builtin.find(fluka_name) == fluka_builtin.end());
}

///---------------------------------------------------------------------------//
/// fluka
///---------------------------------------------------------------------------//
/// Main external call
std::string pyne::Material::fluka(int id, std::string frac_type) {
  std::stringstream rs;

  // Element, one nucid
  if (comp.size() == 1) {
    rs << fluka_material_str(id);
  } else if (comp.size() > 1) {
    // Compound
    rs << fluka_compound_str(id, frac_type);
  } else {
    rs << "There is no nuclide information in the Material Object" << std::endl;
  }
  return rs.str();
}

///---------------------------------------------------------------------------//
/// fluka_material_str
///---------------------------------------------------------------------------//
///
/// Requirement:  the material upon which this function is called has
///               exactly one nucid component, i.e. it is elemental
/// Do not assume fluka_name is defined in the metadata.  This function
/// may be called from a user-defined material, i.e. on that is not
/// read out of a UW^2-tagged geometry file, and thus does not have
/// certain metadata.
std::string pyne::Material::fluka_material_str(int id) {
  std::stringstream ms;
  std::string fluka_name; // needed to determine if built-in

  int nucid = comp.begin()->first;

  // NOTE:  first part of 'if' may never be called
  if (metadata.isMember("fluka_name")) {
    fluka_name = metadata["fluka_name"].asString();
  } else {  // Should be elemental
    if (comp.size() > 1) {
      std::cerr << "Error: this mix is a compound, there should be a fluka_name defined."
                << std::endl;
      return ms.str();
    }
    fluka_name = nucname::fluka(nucid);
  }

  if (not_fluka_builtin(fluka_name)) {
    ms << fluka_material_component(id, nucid, fluka_name);
  }

  // could be empty
  return ms.str();
}

///---------------------------------------------------------------------------//
/// fluka_material_component
///---------------------------------------------------------------------------//
/// Material has only one component,
/// Density is either object density or it is ignored ==> use object density
/// This function is not called for a compound, but it is called on the
/// material-ized components of compounds
std::string pyne::Material::fluka_material_component(int fid, int nucid,
                                                     std::string fluka_name) {
  int znum = pyne::nucname::znum(nucid);

  double atomic_mass;
  if (0 != pyne::NUC_DATA_PATH.length()) {
    // for compounds (i.e., unrecognized nucids), this will be 0
    atomic_mass = pyne::atomic_mass(nucid);
  } else {
    atomic_mass = 1.0;
  }

  return fluka_material_line(znum, atomic_mass, fid, fluka_name);
}

///---------------------------------------------------------------------------//
/// fluka_material_line
///---------------------------------------------------------------------------//
/// Given all the info, return the Material string
std::string pyne::Material::fluka_material_line(int znum, double atomic_mass,
                                                int fid, std::string fluka_name) {
  std::stringstream ls;

  if (metadata.isMember("comments")) {
    std::string comment = metadata["comments"].asString();
    ls << "* " << comment;
    ls << std::endl;
  }
  ls << std::setw(10) << std::left << "MATERIAL";
  ls << std::setprecision(0) << std::fixed << std::showpoint <<
     std::setw(10) << std::right << (float)znum;

  ls << fluka_format_field(atomic_mass);
  // Note this is the current object density, and may or may not be meaningful
  ls << fluka_format_field(std::sqrt(density * density));

  ls << std::setprecision(0) << std::fixed << std::showpoint <<
     std::setw(10) << std::right << (float)fid;
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::left << fluka_name << std::endl;

  return ls.str();
}

///---------------------------------------------------------------------------//
/// fluka_format_field
///---------------------------------------------------------------------------//
/// Convenience function that returns a 10-character formatted string
/// 999 -> 999.
/// 999.12 -> 999.12
/// 999.123 -> 999.123
/// 999.1234 -> 999.123
std::string pyne::Material::fluka_format_field(float field) {
  std::stringstream ls;
  double intpart;
  modf(field, &intpart);
  if (field == intpart) {
    ls << std::setprecision(0) << std::fixed << std::showpoint
       << std::setw(10) << std::right << field;
  } else {
    // This will print however many digits after the decimal, up to a max of six
    ls.unsetf(std::ios::showpoint);
    ls.unsetf(std::ios::floatfield);
    ls.precision(6);
    ls << std::setw(10) << std::right << field;
  }

  return ls.str();
}

///---------------------------------------------------------------------------//
/// fluka_compound_str
///---------------------------------------------------------------------------//
/// Returns
/// -- MATERIAL line for compound
/// -- COMPOUND lines
std::string pyne::Material::fluka_compound_str(int id, std::string frac_type) {
  std::stringstream ss;
  std::map<double, std::string> frac_name_map;
  std::string compound_string = "";
  std::vector<std::string> material_names;

  // The nucid doesn't make sense for a compound
  int znum = 1;
  double atomic_mass = 1.;
  // This better be true
  std::string compound_name;
  if (metadata.isMember("fluka_name")) {
    compound_name = metadata["fluka_name"].asString();
  } else {
    std::cerr << "Error:  metadata \"fluka_name\" expected." << std::endl;
    compound_name = "NotFound";
  }
  ss << fluka_material_line(znum, atomic_mass, id, compound_name);

  std::string frac_sign;
  if ("atom" == frac_type) {
    frac_sign = "";
  } else {
    frac_sign = "-";
  }

  std::stringstream temp_s;
  temp_s << std::scientific;
  temp_s << std::setprecision(3);

  int counter = comp.size();
  pyne::comp_iter nuc = comp.begin();
  // This will pick up multiples of 3 components
  while (counter >= 3) {
    ss << std::setw(10) << std::left  << "COMPOUND";

    temp_s << frac_sign << nuc->second;

    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");  // reset the stringstream for reuse

    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    ss << std::setw(10) << std::left << compound_name;
    ss << std::endl;

    counter -= 3;
  }

  // Get the last (or only, as the case may be) one or two fractions
  if (nuc != comp.end()) {
    ss << std::setw(10) << std::left  << "COMPOUND";
    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    if (nuc != comp.end()) {
      temp_s << frac_sign << nuc->second;
      ss << std::setw(10) << std::right << temp_s.str();
      ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
      nuc++;
      temp_s.str("");
    } else {
      ss << std::setw(10) << std::right << "";
      ss << std::setw(10) << std::right << "";
    }

    ss << std::setw(10) << std::right << "";
    ss << std::setw(10) << std::right << "";
    ss << std::setw(10) << std::left << compound_name;
    ss << std::endl;
  }

  return ss.str();
}


void pyne::Material::from_text(char* filename) {
  std::string fname(filename);
  from_text(fname);
}


void pyne::Material::from_text(std::string filename) {
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // New filestream
  std::ifstream f;
  f.open(filename.c_str());

  // Read in
  comp.clear();
  std::string keystr, valstr;

  f >> keystr;
  while (!f.eof()) {

    if (0 == keystr.length())
      continue;

    if (keystr == "Mass") {
      f >> valstr;
      mass = pyne::to_dbl(valstr);
    } else if (keystr == "Density") {
      f >> valstr;
      density = pyne::to_dbl(valstr);
    } else if (keystr == "APerM") {
      f >> valstr;
      atoms_per_molecule = pyne::to_dbl(valstr);
    } else if (pyne::nucname::isnuclide(keystr) ||
               pyne::nucname::iselement(keystr)) {
      f >> valstr;
      if (comp.count(pyne::nucname::id(keystr)) > 0) {
        comp[pyne::nucname::id(keystr)] += pyne::to_dbl(valstr);
      } else {
        comp[pyne::nucname::id(keystr)] = pyne::to_dbl(valstr);
      }
    } else {
      getline(f, valstr);
      valstr = valstr.substr(0, valstr.length() - 1);
      metadata[keystr] = valstr;

    }
    f >> keystr;
  }

  f.close();
  norm_comp();
}


void pyne::Material::write_text(char* filename) {
  std::string fname(filename);
  write_text(fname);
}


void pyne::Material::write_text(std::string filename) {
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);

  Json::Reader reader;
  std::vector<std::string> obj = metadata.getMemberNames();

  if (0 <= mass)
    f << "Mass    " << mass << "\n";

  if (0 <= density)
    f << "Density "  << density << "\n";

  if (0 <= atoms_per_molecule)
    f << "APerM   " << atoms_per_molecule << "\n";

  for (int i = 0; i < metadata.size(); i = i + 2) {
    f << metadata.get(obj.at(i), "") << metadata.get(obj.at(i + 1), "");
  }

  std::string nuc_name;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    nuc_name = pyne::nucname::name(i->first) + "  ";
    while (nuc_name.length() < 8)
      nuc_name += " ";
    f << nuc_name << i->second << "\n";
  }

  f.close();
}


void pyne::Material::load_json(Json::Value json) {
  Json::Value::Members keys = json["comp"].getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  comp.clear();
  for (; ikey != ikey_end; ++ikey)
    comp[nucname::id(*ikey)] = json["comp"][*ikey].asDouble();
  norm_comp();
  mass = json["mass"].asDouble();
  density = json["density"].asDouble();
  atoms_per_molecule = json["atoms_per_molecule"].asDouble();
  metadata = json["metadata"];
}


Json::Value pyne::Material::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);
  Json::Value jcomp = Json::Value(Json::objectValue);
  json["mass"] = mass;
  json["density"] = density;
  json["atoms_per_molecule"] = atoms_per_molecule;
  json["metadata"] = metadata;
  for (comp_iter i = comp.begin(); i != comp.end(); i++)
    jcomp[nucname::name(i->first)] = (i->second);
  json["comp"] = jcomp;
  return json;
}


void pyne::Material::from_json(char* filename) {
  std::string fname(filename);
  from_json(fname);
}


void pyne::Material::from_json(std::string filename) {
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
  f.seekg(0, std::ios::end);
  s.resize(f.tellg());
  f.seekg(0, std::ios::beg);
  f.read(&s[0], s.size());
  f.close();
  Json::Reader reader;
  Json::Value json;
  reader.parse(s, json);
  load_json(json);
}


void pyne::Material::write_json(char* filename) {
  std::string fname(filename);
  write_json(fname);
}


void pyne::Material::write_json(std::string filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}


/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Material::Material() {
  // Empty Material constructor
  mass = -1.0;
  density = -1.0;
  atoms_per_molecule = -1.0;
  metadata = Json::Value(Json::objectValue);
}


pyne::Material::Material(pyne::comp_map cm, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  density = d;
  atoms_per_molecule = apm;
  metadata = attributes;
  if (!comp.empty())
    norm_comp();
}


pyne::Material::Material(char* filename, double m, double d, double apm,
                         Json::Value attributes) {
  mass = m;
  density = d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  std::string fname(filename);
  if (!pyne::file_exists(fname))
    throw pyne::FileNotFound(fname);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(fname.c_str());
  if (ish5)
    from_hdf5(fname);
  else
    from_text(fname);
}


pyne::Material::Material(std::string filename, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  density = d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename);
  else
    from_text(filename);
}


pyne::Material::~Material() {
}


/*--- Method definitions ---*/
std::ostream& operator<<(std::ostream& os, pyne::Material mat) {
  //print the Mass Stream to stdout
  os << "\tMass: " << mat.mass << "\n";
  os << "\t---------\n";
  for (pyne::comp_iter i = mat.comp.begin(); i != mat.comp.end(); i++) {
    os << "\t" << pyne::nucname::name(i->first) << "\t" << i->second << "\n";
  }
  return os;
}


// Note this refines << for an inheritor of std::ostream.
std::ostringstream& operator<<(std::ostringstream& os, pyne::Material mat) {
  return os;
}


void pyne::Material::normalize() {
  // normalizes the mass
  mass = 1.0;
}


pyne::comp_map pyne::Material::mult_by_mass() {
  // bypass calculation if already normalized.
  if (mass == 1.0)
    return comp;

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    cm[i->first] = (i->second) * mass;
  }
  return cm;
}


pyne::comp_map pyne::Material::activity() {
  pyne::comp_map act;
  double masspermole = mass * pyne::N_A;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
    act[i->first] = masspermole * (i->second) * decay_const(i->first) / \
                    atomic_mass(i->first);
  }
  return act;
}


pyne::comp_map pyne::Material::decay_heat() {
  pyne::comp_map dh;
  double masspermole = mass * pyne::N_A;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
    dh[i->first] = masspermole * (i->second) * \
                   decay_const(metastable_id(i->first, nucname::snum(i->first))) * \
                   q_val(i->first) / atomic_mass(i->first) / pyne::MeV_per_MJ;
  }
  return dh;
}


pyne::comp_map pyne::Material::dose_per_g(std::string dose_type, int source) {
  pyne::comp_map dose;
  const double pCi_per_Bq = 27.027027;
  if (dose_type == "ext_air") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = Ci_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ext_air_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "ext_soil") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = Ci_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ext_soil_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "ingest") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = pCi_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ingest_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "inhale") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = pCi_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * inhale_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else {
    throw std::invalid_argument("Dose type must be one of: ext_air, ext_soil, ingest, inhale.");
  }
  return dose;
}


double pyne::Material::molecular_mass(double apm) {
  // Calculate the atomic weight of the Material
  double inverseA = 0.0;

  for (pyne::comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++)
    inverseA += (nuc->second) / pyne::atomic_mass(nuc->first);

  if (inverseA == 0.0)
    return inverseA;

  // select the atoms per mol
  double atsperm = 1.0; // default to 1.0
  if (0.0 <= apm) {
    atsperm = apm;            // take the function argument, if valid
    if (atoms_per_molecule < 0.0)
      atoms_per_molecule = apm;     // Store the function argument on class, if class has no value
  } else if (0.0 <= atoms_per_molecule)
    atsperm = atoms_per_molecule;  // select the class's value

  return atsperm / inverseA;
}


pyne::Material pyne::Material::expand_elements(std::set<int> exception_ids) {
  // Expands the natural elements of a material and returns a new material note
  // that this implementation relies on the fact that maps of ints are stored in
  // a sorted manner in C++.
  int n, nabund, znuc, zabund;
  comp_map newcomp;
  std::map<int, double>::iterator abund_itr, abund_end;
  if (pyne::natural_abund_map.empty())
    pyne::_load_atomic_mass_map();
  abund_itr = pyne::natural_abund_map.begin();
  abund_end = pyne::natural_abund_map.end();
  zabund = nucname::znum((*abund_itr).first);
  for (comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++) {
    // keep element as-is if in exception list
    if (0 < exception_ids.count(nuc->first)) {
      newcomp.insert(*nuc);
      continue;
    }

    if (abund_itr == abund_end)
      newcomp.insert(*nuc);
    else if (0 == nucname::anum((*nuc).first)) {
      n = (*nuc).first;
      znuc = nucname::znum(n);
      if (znuc < zabund) {
        newcomp.insert(*nuc);
        continue;
      }
      while (zabund <= znuc) {
        nabund = (*abund_itr).first;
        zabund = nucname::znum(nabund);

        if (zabund == znuc && 0 != nucname::anum(nabund) && 0.0 != (*abund_itr).second)
          newcomp[nabund] = (*abund_itr).second * (*nuc).second * \
                            atomic_mass(nabund) / atomic_mass(n);
        else if (n == nabund && 0.0 == (*abund_itr).second)
          newcomp.insert(*nuc);
        abund_itr++;
        if (abund_itr == abund_end) {
          zabund = INT_MAX;
          break;
        }
      }
    } else
      newcomp.insert(*nuc);
  }
  return Material(newcomp, mass, density, atoms_per_molecule, metadata);
}


// Wrapped version for calling from python
pyne::Material pyne::Material::expand_elements(int** int_ptr_arry) {
  std::set<int> nucvec;
  // Set first pointer to first int pointed to by arg
  if (int_ptr_arry != NULL) {
    int* int_ptr = *int_ptr_arry;
    while (int_ptr != NULL) {
      nucvec.insert(*int_ptr);
      int_ptr++;
    }
  }
  return expand_elements(nucvec);
}


pyne::Material pyne::Material::collapse_elements(std::set<int> exception_ids) {
  ////////////////////////////////////////////////////////////////////////
  // Assumptions
  //    - list passed in is of nucid's formed from the znum-anum of
  //      Fluka-named isotopes, since we want to preserve the full
  //      nucid of any such material in the problem
  // Algorithm
  // for each component listed in this material that has a nonzero frac or
  //    weight amount, look at its 'stripped' nucid, that is, the last four
  //    places replaced by zeros.
  //    if it's on the exception list, copy the component
  //    else it is to be collapsed
  //       => add its frac to the component of the znum
  //
  // * When from_hdf5 reads from a file the comp iterator will produce a
  //   hit for EVERY nucid in EVERY material in the file.  Only the nucids
  //   belonging to the CURRENT material have a nonzero fraction/mass amount
  /////////////////////////////////////////////////////////////////////////
  pyne::comp_map cm;

  for (pyne::comp_iter ptr = comp.begin(); ptr != comp.end(); ptr++) {
    if (0 < ptr->second) {
      // There is a nonzero amount of this nucid in the current material,
      // check if znum and anum are in the exception list,
      int cur_stripped_id = nucname::znum(ptr->first) * 10000000 +
                            nucname::anum(ptr->first) * 10000;
      if (0 < exception_ids.count(cur_stripped_id)) {
        // The znum/anum combination identify the current material as a
        // fluka-named exception list => copy, don't collapse
        cm[ptr->first] = (ptr->second) * mass;
      } else {
        // Not on exception list => add frac to id-component
        int znum_id = nucname::id(nucname::znum(ptr->first));
        cm[znum_id] += (ptr->second) * mass;
      }
    }
  }
  // Copy
  pyne::Material collapsed =
      pyne::Material(cm, mass, density, atoms_per_molecule, metadata);
  return collapsed;
}


// Wrapped version for calling from python
pyne::Material pyne::Material::collapse_elements(int** int_ptr_arry) {
  std::set<int> nucvec;
  // Set first pointer to first int pointed to by arg
  int* int_ptr = *int_ptr_arry;
  while (int_ptr != NULL) {
    nucvec.insert(*int_ptr);
    int_ptr++;
  }
  return collapse_elements(nucvec);
}


// Set up atom or mass frac map
std::map<int, double> pyne::Material::get_density_frac(std::string frac_type,
                                                       bool mult_den) {
  std::map<int, double> fracs;

  if ("atom" == frac_type) {
    if (density != -1.0 && mult_den) {
      fracs = to_atom_dens();
      for (comp_iter ci = fracs.begin(); ci != fracs.end(); ci++) {
        ci->second *= pyne::cm2_per_barn;  // unit requirememt is [10^24
        // atoms/cm3] = [atoms/b.cm]
      }
    } else {
      fracs = to_atom_frac();
    }
  } else {
    fracs = comp;
    if (density != -1.0 && mult_den) {
      for (comp_iter ci = fracs.begin(); ci != fracs.end(); ci++) {
        ci->second *= density;
      }
    }
  }
  return fracs;
}


double pyne::Material::mass_density(double num_dens, double apm) {
  if (0.0 <= num_dens) {
    double mw = molecular_mass(apm);
    density = num_dens * mw / pyne::N_A / atoms_per_molecule;
  }
  return density;
}


double pyne::Material::number_density(double mass_dens, double apm) {
  if (0 <= mass_dens)
    density = mass_dens;
  double mw = molecular_mass(apm);
  double num_dens = density * pyne::N_A * atoms_per_molecule / mw;
  return num_dens;
}


/*--- Stub-Stream Computation ---*/
pyne::Material pyne::Material::sub_mat(std::set<int> nucset) {
  // Grabs a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O,
  // 93 for Np, etc).

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if (0 < nucset.count(i->first))
      cm[i->first] = (i->second) * mass;
  }
  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::sub_mat(std::set<std::string> nucset) {
  // Grabs a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }
  return sub_mat(iset);
}


pyne::Material pyne::Material::set_mat(std::set<int> nucset, double value) {
  // Sets a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;

  // Add non-set components
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if (0 == nucset.count(i->first))
      cm[i->first] = (i->second) * mass;
  }

  // Add set component
  for (std::set<int>::iterator nuc = nucset.begin(); nuc != nucset.end(); nuc++)
    cm[*nuc] = value;

  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::set_mat(std::set<std::string> nucset, double value) {
  // Sets a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }
  return set_mat(iset, value);
}


pyne::Material pyne::Material::del_mat(std::set<int> nucset) {
  // Removes a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    // Only add to new comp if not in nucset
    if (0 == nucset.count(i->first))
      cm[i->first] = (i->second) * mass;
  }
  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::del_mat(std::set<std::string> nucset) {
  // Removes a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }
  return del_mat(iset);
}


pyne::Material pyne::Material::sub_range(int lower, int upper) {
  // Grabs a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = (i->second) * mass;
  }
  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::set_range(int lower, int upper, double value) {
  // Sets a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = value;
    else
      cm[i->first] = (i->second) * mass;
  }
  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::del_range(int lower, int upper) {
  // Removes a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((upper <= (i->first)) || ((i->first) < lower))
      cm[i->first] = (i->second) * mass;
  }
  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::sub_elem(int elem) {
  // Returns a material of the element that is a submaterial of this one.
  return sub_range(elem, elem + 10000000);
}


pyne::Material pyne::Material::sub_lan() {
  // Returns a material of Lanthanides that is a sub-material of this one.
  return sub_range(570000000, 720000000);
}


pyne::Material pyne::Material::sub_act() {
  //Returns a material of Actindes that is a sub-material of this one.
  return sub_range(890000000, 1040000000);
}


pyne::Material pyne::Material::sub_tru() {
  // Returns a material of Transuranics that is a sub-material of this one.
  return sub_range(930000000, INT_MAX);
}


pyne::Material pyne::Material::sub_ma() {
  // Returns a material of Minor Actinides that is a sub-material of this one.
  return sub_range(930000000, 1040000000).del_range(940000000, 950000000);
}


pyne::Material pyne::Material::sub_fp() {
  // Returns a material of Fission Products that is a sub-material of this one.
  return sub_range(0, 890000000);
}


/*--- Atom Frac Functions ---*/
std::map<int, double> pyne::Material::to_atom_frac() {
  // Returns an atom fraction map from this material's composition
  // the material's molecular mass
  double mat_mw = molecular_mass();

  std::map<int, double> atom_fracs = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_fracs[ci->first] = (ci->second) * mat_mw / pyne::atomic_mass(ci->first);

  return atom_fracs;
}


void pyne::Material::from_atom_frac(std::map<int, double> atom_fracs) {
  // atom frac must be of the form {nuc: af}, eg, water
  //  80160: 1.0
  //  10010: 2.0

  // clear existing components
  comp.clear();
  atoms_per_molecule = 0.0;

  for (std::map<int, double>::iterator afi = atom_fracs.begin(); afi != atom_fracs.end(); afi++) {
    comp[afi->first] = (afi->second) * pyne::atomic_mass(afi->first);
    atoms_per_molecule += (afi->second);
  }
  norm_comp();
}


std::map<int, double> pyne::Material::to_atom_dens() {
  // Returns an atom density map from this material's composition
  // the material's density

  std::map<int, double> atom_dens = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_dens[ci->first] = (ci->second) * density * pyne::N_A / pyne::atomic_mass(ci->first);

  return atom_dens;
}


std::vector<std::pair<double, double> > pyne::Material::gammas() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
      state_id = nucname::id_to_state_id(ci->first);
    else
      state_id = ci->first;

    std::vector<std::pair<double, double> > raw_gammas = pyne::gammas(state_id);
    for (int i = 0; i < raw_gammas.size(); ++i) {
      result.push_back(std::make_pair(raw_gammas[i].first,
                                      atom_fracs[ci->first]*raw_gammas[i].second));
    }
  }
  return result;
}


std::vector<std::pair<double, double> > pyne::Material::xrays() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
      state_id = nucname::id_to_state_id(ci->first);
    else
      state_id = ci->first;

    std::vector<std::pair<double, double> > raw_xrays = pyne::xrays(state_id);
    for (int i = 0; i < raw_xrays.size(); ++i) {
      result.push_back(std::make_pair(raw_xrays[i].first,
                                      atom_fracs[ci->first]*raw_xrays[i].second));
    }
  }
  return result;
}


std::vector<std::pair<double, double> > pyne::Material::photons(bool norm) {
  std::vector<std::pair<double, double> >  txray = this->xrays();
  std::vector<std::pair<double, double> >  tgammas = this->gammas();
  for (int i = 0; i < txray.size(); ++i)
    tgammas.push_back(txray[i]);
  if (norm)
    tgammas = normalize_radioactivity(tgammas);
  return tgammas;
}


std::vector<std::pair<double, double> > pyne::Material::normalize_radioactivity(
    std::vector<std::pair<double, double> > unnormed) {
  std::vector<std::pair<double, double> > normed;
  double sum = 0.0;
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second))
      sum = sum + unnormed[i].second;
  }
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second)) {
      normed.push_back(
          std::make_pair(unnormed[i].first, (unnormed[i].second) / sum));
    }
  }
  return normed;
}

void pyne::Material::from_activity(std::map<int, double> activities) {
  // activities must be of the form {nuc: act}, eg, tritium
  //  10030: 1.0

  // clear existing components
  comp.clear();

  for (std::map<int, double>::iterator acti = activities.begin();
       acti != activities.end(); acti++) {
    double dc = pyne::decay_const(acti->first);
    if (dc == 0.0 && (acti->second) == 0.0)
      continue;
    else if (dc == 0.0)
      throw std::invalid_argument("Activity keys must be radionuclides.");
    comp[acti->first] = (acti->second) * pyne::atomic_mass(acti->first) / \
                        pyne::N_A / dc;
  }
  norm_comp();
}


#ifdef PYNE_DECAY
pyne::Material pyne::Material::decay(double t) {
  throw pyne::ValueError("Material::decay is not supported in this amalgamated"
                         "version of PyNE.");
}
#endif //  PYNE_DECAY


pyne::Material pyne::Material::cram(std::vector<double> A,
                                    const int order) {
  throw pyne::ValueError("Material::cram is not supported in this amalgamated"
                         "version of PyNE.");
}


pyne::Material pyne::Material::operator+ (double y) {
  // Overloads x + y
  return pyne::Material(comp, mass + y, density);
}


pyne::Material pyne::Material::operator+ (Material y) {
  // Overloads x + y
  pyne::comp_map cm;
  pyne::comp_map xwgt = mult_by_mass();
  pyne::comp_map ywgt = y.mult_by_mass();

  for (pyne::comp_iter i = xwgt.begin(); i != xwgt.end(); i++) {
    if (0 < ywgt.count(i->first))
      cm[i->first] = xwgt[i->first] + ywgt[i->first];
    else
      cm[i->first] = xwgt[i->first];
  }

  for (pyne::comp_iter i = ywgt.begin(); i != ywgt.end(); i++) {
    if (0 == cm.count(i->first))
      cm[i->first] = ywgt[i->first];
  }

  return pyne::Material(cm, -1, -1);
}


pyne::Material pyne::Material::operator* (double y) {
  // Overloads x * y
  return pyne::Material(comp, mass * y, density);
}


pyne::Material pyne::Material::operator/ (double y) {
  // Overloads x / y
  return pyne::Material(comp, mass / y, density);
}


bool pyne::detect_nuclidelist(hid_t data_set, std::string& nucpath) {
  hid_t nuc_attr = H5Aopen(data_set, "nucpath", H5P_DEFAULT);

  // can't find the "nucpath" Attribute
  if (nuc_attr < 0)
    return false;

  H5A_info_t nuc_info;
  H5Aget_info(nuc_attr, &nuc_info);
  hsize_t nuc_attr_len = nuc_info.data_size;
  hid_t str_attr = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_attr, nuc_attr_len);
  char* nucpathbuf = new char[nuc_attr_len];
  H5Aread(nuc_attr, str_attr, nucpathbuf);
  nucpath = std::string(nucpathbuf, nuc_attr_len);
  delete[] nucpathbuf;
  H5Tclose(str_attr);
  H5Tclose(nuc_attr);
  return true;
}

//
// end of src/material.cpp
//


//
// start of src/material_library.cpp
//
#include <unistd.h>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

// Empty Constructor
pyne::MaterialLibrary::MaterialLibrary() {};

// Default constructor when loading from HDF5 file
pyne::MaterialLibrary::MaterialLibrary(const std::string& filename,
                                       const std::string& datapath) {
  // load materials

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename, datapath);
  else if (datapath != "") {
    std::string wrng =
        filename + " is not an hdf5 file but a datapath was prodided!";
    pyne::warning(wrng);
  }
  from_json(filename);
};

// Append Material to the library from hdf5 file
void pyne::MaterialLibrary::from_hdf5(const std::string& filename,
                                      const std::string& datapath) {
  if (!pyne::file_exists(filename)) {
    throw std::runtime_error("File " + filename +
                             " not found or no read permission");
  }
  std::string nucpath;
  std::string full_datapath = datapath;

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  // Get datapath status in the hdf5 file
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  hid_t matlib_group_id = db;
  H5O_info_t object_info;
  status = H5Oget_info_by_name(db, datapath.c_str(), &object_info, H5P_DEFAULT);

  // if datapath exist and is a dataset -> old hdf5 mat_lib format
  // else if datapath does not exist -> new hdf5 mat lib format
  // else don't know what to do with it: fail !
  if (status == 0 && object_info.type == H5O_TYPE_DATASET) {
    // Open dataset to retrieve nucclide list location form nucpath attribute
    hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    pyne::detect_nuclidelist(data_set, nucpath);
    H5Dclose(data_set);
  } else if (status != 0) {
    full_datapath = "/material_library" + datapath + "/composition";
    nucpath = "/material_library" + datapath + "/nuclidelist";
  } else {
    throw std::runtime_error(datapath + " exist but is not a dataset.");
  }

  int file_num_materials = get_length_of_table(db, full_datapath);

  for (int i = 0; i < file_num_materials; i++) {
    pyne::Material mat = pyne::Material();

    mat._load_comp_protocol1(db, full_datapath, nucpath, i);
    add_material(mat);
  }
  H5Fclose(db);
}

void pyne::MaterialLibrary::merge(const pyne::MaterialLibrary& mat_lib) {
  pyne::matname_set mats_to_add = mat_lib.get_keylist();
  for (auto it = mats_to_add.begin(); it != mats_to_add.end(); it++) {
    pyne::Material mat = Material(mat_lib.get_material(*it));
    (*this).add_material(*it, mat);
  }
}

void pyne::MaterialLibrary::merge(pyne::MaterialLibrary* mat_lib) {
  merge(*mat_lib);
}

void pyne::MaterialLibrary::load_json(Json::Value json) {
  Json::Value::Members keys = json.getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  for (; ikey != ikey_end; ++ikey) {
    pyne::Material mat = pyne::Material();
    mat.load_json(json[*ikey]);
    (*this).add_material(*ikey, Material(mat));
  }
}

Json::Value pyne::MaterialLibrary::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);

  for (auto name : keylist) {
    json[name] = material_library[name]->dump_json();
  }
  return json;
}

void pyne::MaterialLibrary::from_json(const std::string& filename) {
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
  f.seekg(0, std::ios::end);
  s.resize(f.tellg());
  f.seekg(0, std::ios::beg);
  f.read(&s[0], s.size());
  f.close();
  Json::Reader reader;
  Json::Value json;
  reader.parse(s, json);
  load_json(json);
}

void pyne::MaterialLibrary::write_json(const std::string& filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}

int pyne::MaterialLibrary::ensure_material_number(pyne::Material& mat) const {
  int mat_number = -1;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue) {
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
      mat.metadata["mat_number"] = mat_number;
    }
    std::set<int>::iterator mat_numb_it;
    mat_numb_it = mat_number_set.find(mat_number);
    if (mat_numb_it != mat_number_set.end()) {
      std::string msg = "Material number ";
      msg += mat_number;
      msg += " is already in the library.";
      warning(msg);
    }
  }
  return mat_number;
}

std::string pyne::MaterialLibrary::ensure_material_name_and_number(
    pyne::Material& mat) const {
  std::string mat_name = "";
  int mat_number = ensure_material_number(mat);

  if (mat.metadata.isMember("name")) {
    if (Json::intValue <= mat.metadata["name"].type() &&
        mat.metadata["name"].type() <= Json::realValue) {
      mat_name = std::to_string(mat.metadata["name"].asInt());
      mat.metadata["name"] = mat_name;
    } else {
      mat_name = mat.metadata["name"].asString();
    }
  } else {
    if (mat_number == -1) {
      if (!mat_number_set.empty())
        mat_number = *mat_number_set.rbegin() + 1;
      else {
        mat_number = 1;
      }
    }
    mat_name = "_" + std::to_string(mat_number);
    mat.metadata["name"] = mat_name;
  }
  return mat_name;
}

void pyne::MaterialLibrary::add_material(pyne::Material mat) {
  // if exists, get the material name from metadata make one instead
  std::string mat_name = ensure_material_name_and_number(mat);
  add_material(mat_name, mat);
}

void pyne::MaterialLibrary::add_material(const std::string& key,
                                         const pyne::Material& mat) {
  pyne::shr_mat_ptr new_mat = std::make_shared<pyne::Material>(mat);

  int mat_number = ensure_material_number(*new_mat);
  if (!new_mat->metadata.isMember("name")) {
    new_mat->metadata["name"] = key;
  }
  append_to_nuclist(*new_mat);
  if (mat_number > 0)
    mat_number_set.insert(mat_number);
  keylist.insert(key);
  material_library[key] = new_mat;
}

void pyne::MaterialLibrary::del_material(const std::string& key) {
  material_library.erase(key);
  keylist.erase(key);
}

pyne::Material pyne::MaterialLibrary::get_material(
    const std::string& mat_name) const {
  return *(this->get_material_ptr(mat_name));
}

pyne::shr_mat_ptr pyne::MaterialLibrary::get_material_ptr(
    const std::string& mat_name) const {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return it->second;
  } else {
    return std::make_shared<pyne::Material>(pyne::Material());
  }
}

void pyne::MaterialLibrary::write_hdf5(const std::string& filename,
                                       const std::string& datapath,
                                       bool h5_overwrite) const {
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  // Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  hid_t matlib_grp_id;
  hid_t data_id;
  hid_t db;
  if (!pyne::file_exists(filename)) {
    // Create the file
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  } else {
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  }
  // Check if root_path exist and what type it is
  std::string root_path = "/material_library";

  H5O_info_t object_info;
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check if root_path exist and what type it is
  status =
      H5Oget_info_by_name(db, root_path.c_str(), &object_info, H5P_DEFAULT);
  if (status == 0) {
    // "/material_library" not a group: fail!
    if (object_info.type != H5O_TYPE_GROUP) {
      throw std::runtime_error(
          "Non-group/non-dataset object /material_library already exists in "
          "the file. Can't write the Material");
    } else {  // "/material_library" is a group get its hid
      matlib_grp_id = H5Gopen2(db, root_path.c_str(), H5P_DEFAULT);
    }
  } else {  // "/material_library" does not exist -> create it !
    matlib_grp_id = H5Gcreate2(db, root_path.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT);
  }

  // Group "/root_path/datapath" does exist throw an error or remove it
  std::string full_path = root_path + datapath;
  if (h5wrap::path_exists(db, full_path)) {
    if (h5_overwrite) {
      H5Ldelete(db, full_path.c_str(), H5P_DEFAULT);
    } else {
      throw std::runtime_error(
          "datapath" + full_path +
          "already exist, use \"h5_overwrite=true\" option to "
          "overwrite existing datapath");
    }
  }

  // create "/root_path/datapath" Group
  data_id =
      H5Gcreate2(db, full_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::string compath = "composition";
  std::string nucpath = "nuclidelist";
  write_hdf5_nucpath(data_id, nucpath);
  std::vector<int> nuc_list;
  nuc_list.assign(nuclist.begin(), nuclist.end());
  for (auto mat : material_library) {
    mat.second->write_hdf5_datapath(data_id, compath, -0.0,
                                    DEFAULT_MAT_CHUNKSIZE, nuc_list);
  }

  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Gclose(data_id);
  H5Gclose(matlib_grp_id);
  H5Fclose(db);
}

void pyne::MaterialLibrary::write_hdf5_nucpath(hid_t db,
                                               std::string nucpath) const {
  int nuc_size;
  nuc_size = nuclist.size();
  int nuc_data[nuc_size];
  {
    int n = 0;
    for (auto nuc : nuclist) {
      nuc_data[n] = nuc;
      n++;
    }
  }
  hsize_t nuc_dims[1];
  nuc_dims[0] = nuc_size;
  hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
  hid_t nuc_set = H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuc_space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nuc_data);

  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Sclose(nuc_space);
  H5Dclose(nuc_set);
}

void pyne::MaterialLibrary::append_to_nuclist(const pyne::Material& mat) {
  pyne::comp_map mat_comp = mat.comp;
  for (auto nuclide : mat_comp) {
    nuclist.insert(nuclide.first);
  }
}

int pyne::MaterialLibrary::get_length_of_table(
    hid_t db, const std::string& datapath) const {
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);
  // failure to read the dapaspace return 0 for the size
  if (arr_space < 0)
    return 0;

  hsize_t arr_dims[1];
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  herr_t status = H5Eclear(H5E_DEFAULT);
  return arr_dims[0];
}
//
// end of src/material_library.cpp
//


//
// start of src/tally.cpp
//
// Tally.cpp
// Central Tally Class
// -- Andrew Davis

#include <iomanip>
#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
#include "particle.h"
#include "tally.h"
#endif

enum entity_type_enum {VOLUME, SURFACE, MESH}; // Enumeration for entity types
enum tally_type_enum {FLUX, CURRENT};  // Enumeration for tally types

const std::string tally_type_enum2string[] = {"Flux", "Current"};
const std::string entity_type_enum2string[] = {"Volume", "Surface", "Mesh"};
const std::string geometry_type_enum2string[] = {"Cartesian", "Cylinder"};

/***************************/
/*** Protected Functions ***/
/***************************/

// there are no protected functions currently
// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Tally::Tally() {
  // Empty Tally Constructor
  tally_type = "";
  entity_id = -1;
  entity_type = "";
  entity_name = "";
  tally_name = "";
  entity_size = -1.0;
  normalization = 1.0;
}

// Default constructor
pyne::Tally::Tally(std::string type, std::string part_name,
                   int ent, std::string ent_type,
                   std::string ent_name, std::string tal_name,
                   double size, double norm) {

  tally_type = type;
  particle_names.push_back(pyne::particle::name(part_name));
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string type, std::vector<std::string> part_names,
                   int ent, std::string ent_type,
                   std::string ent_name, std::string tal_name,
                   double size, double norm) {

  tally_type = type;
  particle_names = part_names;
  for (int i = 0; i < particle_names.size(); i++) {
    particle_names[i] = pyne::particle::name(particle_names[i]);
  }
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string part_name, std::string ent_geom,
                   std::vector<double> orgn,
                   std::vector<double> mesh_i, std::vector<double> mesh_j, std::vector<double> mesh_k,
                   std::vector<int> ints_i, std::vector<int> ints_j, std::vector<int> ints_k,
                   std::vector<double> e_bounds_, std::vector<int> e_ints_,
                   std::vector<double> axs_, std::vector<double> vec_,
                   std::string tal_name, double norm) {
  // Empty Tally Constructor
  entity_type = "Mesh";
  entity_name = "";
  particle_names.push_back(pyne::particle::name(part_name));
  entity_geometry = ent_geom;
  tally_name = tal_name;
  entity_size = -1;

  origin = orgn;
  vec = vec_;
  axs = axs_;
  meshes[0] = mesh_i;
  meshes[1] = mesh_j;
  meshes[2] = mesh_k;
  ints[0] = ints_i;
  ints[1] = ints_j;
  ints[2] = ints_k;
  e_bounds = e_bounds_;
  e_ints = e_ints_;
  normalization = norm;
}

// Destructor
pyne::Tally::~Tally() {}

/*--- Method definitions ---*/
//
void pyne::Tally::from_hdf5(char* filename, char* datapath, int row) {
  std::string fname(filename);
  std::string dpath(datapath);
  from_hdf5(fname, dpath, row);
}

//
void pyne::Tally::from_hdf5(std::string filename, std::string datapath,
                            int row) {
  // line of data to acces
  int data_row = row;

  // check for file existence
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());
  if (!is_h5)
    throw h5wrap::FileNotHDF5(filename);

  // Open file and dataset.
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);

  // Get dataspace and allocate memory for read buffer.
  hid_t space = H5Dget_space(dset);
  int rank  = H5Sget_simple_extent_ndims(space);
  hsize_t dims[1]; // for length of dataset

  // get the length of the dataset
  int ndims = H5Sget_simple_extent_dims(space, dims, NULL);

  // determine if chunked
  hid_t prop = H5Dget_create_plist(dset);

  hsize_t chunk_dimsr[1];
  int rank_chunk;

  if (H5D_CHUNKED == H5Pget_layout(prop))
    rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);

  // allocate memory for data from file
  tally_struct* read_data = new tally_struct[dims[0]];

  // if row number is larger than data set only give last element
  if (row >= dims[0])
    data_row = dims[0] - 1;

  // Create variable-length string datatype.
  hid_t strtype = H5Tcopy(H5T_C_S1);
  int status  = H5Tset_size(strtype, H5T_VARIABLE);

  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();

  // Create the compound datatype for the file
  hid_t filetype = create_filetype();

  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);

  // unpack the data and set values
  entity_id = read_data[data_row].entity_id;
  entity_type = entity_type_enum2string[read_data[data_row].entity_type];
  tally_type = tally_type_enum2string[read_data[data_row].tally_type];
  particle_names = pyne::split_string(read_data[data_row].particle_name, ",");
  tally_name = std::string(read_data[data_row].tally_name);
  entity_name = std::string(read_data[data_row].entity_name);
  entity_size = read_data[data_row].entity_size;
  normalization = read_data[data_row].normalization;

  // close the data sets
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(filetype);
  status = H5Fclose(file);

  // tidy up
  delete[] read_data;

}

// Dummy Wrapper around C Style Functions
void pyne::Tally::write_hdf5(char* filename, char* datapath) {
  std::string fname(filename);
  std::string groupname(datapath);
  write_hdf5(fname, groupname);
}

// create filetype
hid_t pyne::Tally::create_filetype() {
  herr_t status;  // iostatus

  // create string type
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);

  hid_t filetype = H5Tcreate(H5T_COMPOUND, 8 + 8 + 8 +
                             (3 * sizeof(hvl_t)) + 8 + 8);
  status = H5Tinsert(filetype, "entity_id", 0, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "entity_type", 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "tally_type", 8 + 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "particle_name", 8 + 8 + 8, strtype);
  status = H5Tinsert(filetype, "entity_name", 8 + 8 + 8 +
                     sizeof(hvl_t), strtype);
  status = H5Tinsert(filetype, "tally_name", 8 + 8 + 8 +
                     (2 * sizeof(hvl_t)), strtype);
  status = H5Tinsert(filetype, "entity_size", 8 + 8 + 8 +
                     (3 * sizeof(hvl_t)), H5T_IEEE_F64BE);
  status = H5Tinsert(filetype, "normalization", 8 + 8 + 8 +
                     (3 * sizeof(hvl_t)) + 8, H5T_IEEE_F64BE);
  return filetype;
}

// create memory type
hid_t pyne::Tally::create_memtype() {
  // iostatus
  herr_t status;

  //
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);
  // Create the compound datatype for memory.
  hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(tally_struct));
  status = H5Tinsert(memtype, "entity_id",
                     HOFFSET(tally_struct, entity_id), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "entity_type",
                     HOFFSET(tally_struct, entity_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "tally_type",
                     HOFFSET(tally_struct, tally_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "particle_name",
                     HOFFSET(tally_struct, particle_name), strtype);
  status = H5Tinsert(memtype, "entity_name", HOFFSET(tally_struct, entity_name),
                     strtype);
  status = H5Tinsert(memtype, "tally_name", HOFFSET(tally_struct, tally_name),
                     strtype);
  status = H5Tinsert(memtype, "entity_size",
                     HOFFSET(tally_struct, entity_size), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(memtype, "normalization",
                     HOFFSET(tally_struct, normalization), H5T_NATIVE_DOUBLE);
  return memtype;
}

hid_t pyne::Tally::create_dataspace(hid_t file, std::string datapath) {
  // enable chunking
  hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  // set chunk size
  hsize_t chunk_dimensions[1] = {1};
  herr_t status = H5Pset_chunk(prop, 1, chunk_dimensions);

  // allow varaible length strings
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);

  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();

  // Create the compound datatype for the file
  hid_t filetype = create_filetype();

  // max dims unlimted
  hsize_t max_dims[1] = {H5S_UNLIMITED};
  // only ever let 1 tally object be added
  hsize_t dims[1] = {1};
  // Create dataspace.  Setting maximum size to NULL sets the maximum
  hid_t space = H5Screate_simple(1, dims, max_dims);

  // Create the dataset and write the compound data to it.
  return H5Dcreate2(file, datapath.c_str(), filetype, space, H5P_DEFAULT, prop,
                    H5P_DEFAULT);
}

// Appends Tally object to dataset if file & datapath already exists
// if file exists & data path doesnt creates new datapath,
// otherwise creates new file
void pyne::Tally::write_hdf5(std::string filename, std::string datapath) {

  // turn of annoying hdf5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  tally_struct tally_data[1]; // storage for the tally to add

  // setup the data to write
  tally_data[0].entity_id = entity_id;
  // entity type
  if (entity_type.find("Volume") != std::string::npos)
    tally_data[0].entity_type = VOLUME;
  else if (entity_type.find("Surface") != std::string::npos)
    tally_data[0].entity_type = SURFACE;

  // tally kind
  if (tally_type.find("Flux") != std::string::npos)
    tally_data[0].tally_type = FLUX;
  else if (tally_type.find("Current") != std::string::npos)
    tally_data[0].tally_type = CURRENT;

  // unpack from class to struct array
  std::string particle_name = pyne::join_to_string(particle_names, ",");
  tally_data[0].entity_id = entity_id;
  tally_data[0].entity_name = entity_name.c_str();
  tally_data[0].particle_name = particle_name.c_str();
  tally_data[0].tally_name = tally_name.c_str();
  tally_data[0].entity_size = entity_size;
  tally_data[0].normalization = normalization;


  // check for file existence
  bool is_exist = pyne::file_exists(filename);
  // create new file

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());

  if (is_exist && !is_h5)
    throw h5wrap::FileNotHDF5(filename);

  if (!is_exist) {  // is a new file
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           H5P_DEFAULT);
    // create a dataspace
    hid_t dset = create_dataspace(file, datapath);

    hid_t memtype = create_memtype();

    herr_t status; // iostatus

    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

    // close the data sets
    status = H5Dclose(dset);
    //    status = H5Sclose(space);
    //    status = H5Tclose(filetype);
    status = H5Fclose(file);

  }  else if (is_exist && is_h5) {  // already exists and is an hdf file
    // then we append the data to the end
    herr_t data_status; // iostatus

    // Open file and dataset.
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // see if path exists
    data_status = H5Gget_objinfo(file, datapath.c_str(), 0, NULL);

    hid_t dset;
    // if fails neet to create dataset
    // still need to check that the datapath exists
    if (data_status != 0) { // doesnt exist
      dset = create_dataspace(file, datapath.c_str());
      hid_t memtype = create_memtype();
      herr_t status; // iostatus

      status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

      // close the data sets
      status = H5Dclose(dset);
      status = H5Fclose(file);
    } else {

      dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);

      // Get dataspace and allocate memory for read buffer.
      hid_t space = H5Dget_space(dset);
      int rank  = H5Sget_simple_extent_ndims(space);
      hsize_t dims[1]; // for length of dataset

      // get the length of the dataset
      int ndims = H5Sget_simple_extent_dims(space, dims, NULL);

      // determine if chunked
      hid_t prop = H5Dget_create_plist(dset);

      hsize_t chunk_dimsr[1];
      int rank_chunk;
      if (H5D_CHUNKED == H5Pget_layout(prop))
        rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);

      // allocate memory for data from file
      tally_struct* read_data = new tally_struct[dims[0]];

      // Create variable-length string datatype.
      hid_t strtype = H5Tcopy(H5T_C_S1);
      int status  = H5Tset_size(strtype, H5T_VARIABLE);

      // Create the compound datatype for memory.
      hid_t memtype = create_memtype();

      // Create the compound datatype for the file
      hid_t filetype = create_filetype();

      // Read the data.
      status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);

      // resize dims
      dims[0] += 1;

      // Extend the dataset
      status = H5Dextend(dset, dims);
      hid_t filespace = H5Dget_space(dset);
      // calculate the existing offset
      hsize_t offset[1] = {dims[0] - 1};

      // select hyerslab
      hsize_t new_length[1] = {1};
      status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                   new_length, NULL);

      // create dataspace for new data
      space = H5Screate_simple(1, new_length, NULL);

      // Write the dataset to memory
      status = H5Dwrite(dset, memtype, space, filespace, H5P_DEFAULT, tally_data);

      // tidy up
      status = H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, read_data);
      delete[] read_data;
      status = H5Dclose(dset);
      status = H5Sclose(space);
      status = H5Tclose(memtype);
      status = H5Tclose(strtype);
      status = H5Fclose(file);
    }
  }
}

std::ostream& operator<<(std::ostream& os, pyne::Tally tal) {
  //print the Tally to ostream
  os << "\t---------\n";
  os << "\t Tallying " << pyne::join_to_string(tal.particle_names, ", ") << " " << tal.tally_type << std::endl;
  os << "\t in/on " << tal.entity_type << " " << tal.entity_id << std::endl;
  return os;
}

// Sets string to valid mcnp formatted tally
// Takes mcnp version as arg, like 5 or 6
std::string pyne::Tally::mcnp(int tally_index, std::string mcnp_version,
                              std::string out) {
  std::stringstream output;  // output stream
  std::string particle_token = "";

  if (particle_names.size() == 0) {
    particle_token = "?";
  }

  // particle token
  for (int i = 0; i < particle_names.size(); i++) {
    if (i > 0)
      particle_token += ",";
    if (mcnp_version.find("mcnp5") != std::string::npos)
      particle_token += pyne::particle::mcnp(particle_names[i]);
    else if (mcnp_version.find("mcnp6") != std::string::npos)
      particle_token += pyne::particle::mcnp6(particle_names[i]);
    else
      particle_token += "?";
  }
  // print out comment line
  output << "C " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(6);

  if (normalization != 1.0)
    output << std::scientific;
  int tally_id = 0;

  // neednt check entity type
  if (entity_type.find("Surface") != std::string::npos) {
    if (tally_type.find("Current") != std::string::npos) {
      tally_id = 1;
    } else if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 2;
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Volume") != std::string::npos) {
    if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 4;
    } else if (tally_type.find("Current") != std::string::npos) {
      // makes no sense in mcnp
      return "";
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Mesh") != std::string::npos) {
    output << form_mcnp_meshtally(tally_index, particle_token, entity_geometry,
                                  axs, vec, origin, meshes, ints, e_bounds,
                                  e_ints, out);

  } else {
    std::cout << "tally/entity combination makes no sense for MCNP"
              << std::endl;
  }
  // print sd card if area/volume specified
  return output.str();
}

template <typename T>
bool pyne::Tally::is_zero(T vect) {
  int size = sizeof(vect) / sizeof(vect[0]);
  bool result = true;
  for (int i = 0; i < size; i++)
    result &= (vect[i] == 0);
  return result;
}


// Form the tally line as function of its properties
std::string pyne::Tally::form_mcnp_tally(int tally_index, int type,
                                         std::string particle_token,
                                         int entity_id, double entity_size,
                                         double normalization) {
  std::stringstream tally_stream;  // tally stream
  tally_stream << std::setiosflags(std::ios::fixed) << std::setprecision(6);
  if (normalization != 1.0)
    tally_stream << std::scientific;

  tally_stream << "F" << tally_index << type << ":" << particle_token << " "
               << entity_id << std::endl;

  if (entity_size > 0.0)
    tally_stream << "SD" << tally_index << type << " " << entity_size
                 << std::endl;

  if (normalization != 1.0)
    tally_stream << "FM" << tally_index << type << " " << normalization
                 << std::endl;

  return tally_stream.str();
}


// Form the mesh tally line as function of its properties
std::string pyne::Tally::form_mcnp_meshtally(
    int tally_index, std::string particle_token, std::string entity_geometry,
    std::vector<double> axs, std::vector<double> vec,
    std::vector<double> origin, std::vector<double> meshes[3],
    std::vector<int> ints[3], std::vector<double> e_bounds,
    std::vector<int> e_ints, std::string out) {
  std::stringstream mtally_stream;
  // indentation block
  std::string indent_block = "          ";

  mtally_stream << "FMESH" << tally_index << "4:" << particle_token << " ";
  mtally_stream << "GEOM=";

  if (entity_geometry.find("Cartesian") != std::string::npos) {
    mtally_stream << "XYZ" << std::endl;
    mtally_stream << indent_block;
  } else if (entity_geometry.find("Cylinder") != std::string::npos) {
    mtally_stream << "CYL" << std::endl;
    if (!is_zero(axs)) {
      mtally_stream << indent_block << "AXS=" << pyne::join_to_string(axs) << std::endl;
    }
    if (!is_zero(vec)) {
      mtally_stream << indent_block << "VEC=" << pyne::join_to_string(vec) << std::endl;
    }
    mtally_stream << indent_block;
  }

  mtally_stream << "ORIGIN=" << pyne::join_to_string(origin) << std::endl;
  std::string dir_name[3] = {"I", "J", "K"};

  for (int j = 0; j < 3; j++) {
    mtally_stream << indent_block;
    mtally_stream << dir_name[j] << "MESH=" << pyne::join_to_string(meshes[j]);
    if (ints[j].size() > 0) {
      mtally_stream << " " << dir_name[j]
                    << "INTS=" << pyne::join_to_string(ints[j]);
    }
    mtally_stream << std::endl;
  }
  if (e_bounds.size() > 0) {
    mtally_stream << indent_block << "EMESH=" << pyne::join_to_string(e_bounds);
  }
  mtally_stream << std::endl;
  if (e_ints.size() > 0) {
    mtally_stream << indent_block << "EINTS=" << pyne::join_to_string(e_ints);
  }
  if (out.size() > 0) {
    mtally_stream << std::endl << indent_block << "OUT=" << out;
  }
  return mtally_stream.str();
}


// Produces valid fluka tally
std::string pyne::Tally::fluka(std::string unit_number) {
  std::stringstream output; // output stream

  if (particle_names.size() != 1) {
    std::cout << "Fluka multi-particles tally have not been yet implemented!" << std:: endl;
    exit(1);
  }
  std::string particle_name = particle_names[0];

  // check entity type
  if (entity_type.find("Volume") != std::string::npos) {
    // ok
  }  else if (entity_type.find("Surface") != std::string::npos) {
    std::cout << "Surface tally not valid in FLUKA" << std::endl;
  } else {
    std::cout << "Unknown entity type" << std::endl;
  }

  output << "* " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(1);
  // check tally type
  if (tally_type.find("Flux") != std::string::npos) {
    output << std::setw(10) << std::left  << "USRTRACK";
    output << std::setw(10) << std::right << "     1.0";
    output << std::setw(10) << std::right
           << pyne::particle::fluka(particle_name);
    output << std::setw(10) << std::right << unit_number;
    output << std::setw(10) << std::right << entity_name;
    if (entity_size > 0.0) {
      output << std::scientific;
      output << std::setprecision(4);
      output << std::setw(10) << std::right << entity_size;
    } else
      output << std::setw(10) << std::right << 1.0;

    output << std::setw(10) << std::right << "   1000."; // number of eints
    tally_name.resize(8);
    output << std::setw(8) << std::left
           << tally_name; // may need to make sure less than 10 chars
    output << std::endl;
    output << std::setw(10) << std::left  << "USRTRACK";
    output << std::setw(10) << std::right << "   10.E1";
    output << std::setw(10) << std::right << "   1.E-3";
    output << std::setw(10) << std::right << "        ";
    output << std::setw(10) << std::right << "        ";
    output << std::setw(10) << std::right << "        ";
    output << std::setw(10) << std::right << "        ";
    output << std::setw(8) << std::left << "       &";
    // end of usrtrack
  } else if (tally_type.find("Current") != std::string::npos) {
    output << std::setw(10) << std::left  << "USRBDX  ";
    output << std::setw(10) << std::right << "   110.0";
    output << std::setw(10) << std::right
           << pyne::particle::fluka(particle_name);
    output << std::setw(10) << std::right << unit_number;
    output << std::setw(10) << std::right << entity_name; // upstream
    output << std::setw(10) << std::right << entity_name; // downstream
    if (entity_size > 0.0)
      output << std::setw(10) << std::right << entity_size; // area
    else
      output << std::setw(10) << std::right << 1.0;

    tally_name.resize(8);
    output << std::setw(8) << std::right
           << tally_name; // may need to make sure less than 10 chars
    output << std::endl;
    output << std::setw(10) << std::left  << "USRBDX  ";
    output << std::setw(10) << std::right << "  10.0E1";
    output << std::setw(10) << std::right << "     0.0";
    output << std::setw(10) << std::right << "  1000.0"; // number of ints
    output << std::setw(10) << std::right << "12.56637"; // 4pi
    output << std::setw(10) << std::right << "     0.0";
    output << std::setw(10) << std::right
           << "   240.0"; // number of angular ints
    output << std::setw(8) << std::left << "       &";
    // end of usrbdx
  } else {
    std::cout << "Unknown tally type" << std::endl;
  }
  return output.str();
}
//
// end of src/tally.cpp
//


//
// start of src/atomic_data.cpp
//
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This file has been autogenerated, modify atomicgen.py !!!!!!!!!!
!!!!!!          DO NOT MODIFY THIS FILE BY HAND!!            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
// Implements basic nuclear data functions.
#ifndef PYNE_IS_AMALGAMATED
#include "atomic_data.h"
#include "nucname.h"
#endif

void pyne::_load_atomic_mass_map_memory() {
  // header version of atomic weight table data
  //see if the data table is already loaded
  if (!atomic_mass_map.empty()) {
    return;
  } else {
    _insert_atomic_mass_map();
  }
  //see if the data table is already loaded
  if (!natural_abund_map.empty()) {
    return;
  } else {
    _insert_abund_map();
  }
  // calculate the atomic_masses of the elements
  std::map<int, double> :: iterator it;

  for (int z = 1; z <= 92 ; z++) {
    // loop through the natural abundance map
    double element_atomic_weight = 0.0;
    for (it = natural_abund_map.begin(); it != natural_abund_map.end() ; ++it) {
      // if the atomic number of the abudance matches the
      // that of index
      if (pyne::nucname::znum(it->first) == z) {
        // take atomic abundance and multiply by mass
        // to get the mass of that nuclide / 100 since abundance is in %
        element_atomic_weight += (it->second * atomic_mass_map[it->first] / 100.0);
      }
    }
    // insert the abundance of the element into the list
    atomic_mass_map[z * 10000000] = element_atomic_weight;
  }
}


void pyne::_insert_atomic_mass_map() {
  atomic_mass_map[10010000] = 1.00782503224;
  atomic_mass_map[10020000] = 2.01410177811;
  atomic_mass_map[10030000] = 3.01604928199;
  atomic_mass_map[20030000] = 3.01602932265;
  atomic_mass_map[30030000] = 3.030775;
  atomic_mass_map[10040000] = 4.026431868;
  atomic_mass_map[20040000] = 4.00260325413;
  atomic_mass_map[30040000] = 4.027185562;
  atomic_mass_map[10050000] = 5.035311493;
  atomic_mass_map[20050000] = 5.012057224;
  atomic_mass_map[30050000] = 5.0125378;
  atomic_mass_map[40050000] = 5.03987;
  atomic_mass_map[10060000] = 6.044955437;
  atomic_mass_map[20060000] = 6.018885891;
  atomic_mass_map[30060000] = 6.01512288742;
  atomic_mass_map[40060000] = 6.0197264090000004;
  atomic_mass_map[50060000] = 6.0508;
  atomic_mass_map[10070000] = 7.052749;
  atomic_mass_map[20070000] = 7.027990654;
  atomic_mass_map[30070000] = 7.01600343666;
  atomic_mass_map[40070000] = 7.016928717;
  atomic_mass_map[50070000] = 7.029712;
  atomic_mass_map[20080000] = 8.03393439;
  atomic_mass_map[30080000] = 8.022486246;
  atomic_mass_map[40080000] = 8.005305102;
  atomic_mass_map[50080000] = 8.024607316;
  atomic_mass_map[60080000] = 8.037643042;
  atomic_mass_map[20090000] = 9.043946419;
  atomic_mass_map[30090000] = 9.026790191;
  atomic_mass_map[40090000] = 9.012183066;
  atomic_mass_map[50090000] = 9.013329649;
  atomic_mass_map[60090000] = 9.031037207;
  atomic_mass_map[20100000] = 10.052815308;
  atomic_mass_map[30100000] = 10.035483453;
  atomic_mass_map[40100000] = 10.013534695;
  atomic_mass_map[50100000] = 10.012936862;
  atomic_mass_map[60100000] = 10.016853218;
  atomic_mass_map[70100000] = 10.041653543;
  atomic_mass_map[30110000] = 11.043723581;
  atomic_mass_map[40110000] = 11.021661081;
  atomic_mass_map[50110000] = 11.009305166;
  atomic_mass_map[60110000] = 11.011432597;
  atomic_mass_map[70110000] = 11.026090945;
  atomic_mass_map[30120000] = 12.052613941;
  atomic_mass_map[40120000] = 12.026922083;
  atomic_mass_map[50120000] = 12.014352638;
  atomic_mass_map[60120000] = 12.0;
  atomic_mass_map[70120000] = 12.018613182;
  atomic_mass_map[80120000] = 12.034261747;
  atomic_mass_map[30130000] = 13.061171503;
  atomic_mass_map[40130000] = 13.036134507;
  atomic_mass_map[50130000] = 13.017779981;
  atomic_mass_map[60130000] = 13.00335483521;
  atomic_mass_map[70130000] = 13.005738609;
  atomic_mass_map[80130000] = 13.024815437;
  atomic_mass_map[40140000] = 14.04289292;
  atomic_mass_map[50140000] = 14.025404012;
  atomic_mass_map[60140000] = 14.00324198843;
  atomic_mass_map[70140000] = 14.00307400446;
  atomic_mass_map[80140000] = 14.008596706;
  atomic_mass_map[90140000] = 14.034315199;
  atomic_mass_map[40150000] = 15.053490215;
  atomic_mass_map[50150000] = 15.031087953;
  atomic_mass_map[60150000] = 15.010599256;
  atomic_mass_map[70150000] = 15.00010889894;
  atomic_mass_map[80150000] = 15.003065618;
  atomic_mass_map[90150000] = 15.017785139;
  atomic_mass_map[100150000] = 15.04317298;
  atomic_mass_map[40160000] = 16.061672036;
  atomic_mass_map[50160000] = 16.03984192;
  atomic_mass_map[60160000] = 16.014701256;
  atomic_mass_map[70160000] = 16.006101925;
  atomic_mass_map[80160000] = 15.9949146196;
  atomic_mass_map[90160000] = 16.011465723;
  atomic_mass_map[100160000] = 16.025750864;
  atomic_mass_map[50170000] = 17.046931399;
  atomic_mass_map[60170000] = 17.022578672;
  atomic_mass_map[70170000] = 17.008448877;
  atomic_mass_map[80170000] = 16.99913175664;
  atomic_mass_map[90170000] = 17.002095238;
  atomic_mass_map[100170000] = 17.017713959;
  atomic_mass_map[110170000] = 17.03776;
  atomic_mass_map[50180000] = 18.055601682;
  atomic_mass_map[60180000] = 18.026751932;
  atomic_mass_map[70180000] = 18.014077565;
  atomic_mass_map[80180000] = 17.99915961284;
  atomic_mass_map[90180000] = 18.000937325;
  atomic_mass_map[100180000] = 18.005708693;
  atomic_mass_map[110180000] = 18.026879386;
  atomic_mass_map[50190000] = 19.064166;
  atomic_mass_map[60190000] = 19.034797596;
  atomic_mass_map[70190000] = 19.017022419;
  atomic_mass_map[80190000] = 19.00357797;
  atomic_mass_map[90190000] = 18.99840316288;
  atomic_mass_map[100190000] = 19.001880903;
  atomic_mass_map[110190000] = 19.013880272;
  atomic_mass_map[120190000] = 19.034169182;
  atomic_mass_map[50200000] = 20.073484;
  atomic_mass_map[60200000] = 20.040261732;
  atomic_mass_map[70200000] = 20.023367295;
  atomic_mass_map[80200000] = 20.004075358;
  atomic_mass_map[90200000] = 19.999981252;
  atomic_mass_map[100200000] = 19.99244017619;
  atomic_mass_map[110200000] = 20.007354426;
  atomic_mass_map[120200000] = 20.018763075;
  atomic_mass_map[50210000] = 21.083017;
  atomic_mass_map[60210000] = 21.049;
  atomic_mass_map[70210000] = 21.027087573;
  atomic_mass_map[80210000] = 21.00865495;
  atomic_mass_map[90210000] = 20.999948894;
  atomic_mass_map[100210000] = 20.993846685;
  atomic_mass_map[110210000] = 20.997654702;
  atomic_mass_map[120210000] = 21.011705764;
  atomic_mass_map[130210000] = 21.028975;
  atomic_mass_map[60220000] = 22.05755399;
  atomic_mass_map[70220000] = 22.034100918;
  atomic_mass_map[80220000] = 22.009965746;
  atomic_mass_map[90220000] = 22.002998809;
  atomic_mass_map[100220000] = 21.991385109;
  atomic_mass_map[110220000] = 21.994437418;
  atomic_mass_map[120220000] = 21.999570654;
  atomic_mass_map[130220000] = 22.01954;
  atomic_mass_map[140220000] = 22.03579;
  atomic_mass_map[60230000] = 23.06889;
  atomic_mass_map[70230000] = 23.039421;
  atomic_mass_map[80230000] = 23.015696686;
  atomic_mass_map[90230000] = 23.003526874;
  atomic_mass_map[100230000] = 22.9944669;
  atomic_mass_map[110230000] = 22.98976928199;
  atomic_mass_map[120230000] = 22.994123941;
  atomic_mass_map[130230000] = 23.007244351;
  atomic_mass_map[140230000] = 23.02544;
  atomic_mass_map[70240000] = 24.05039;
  atomic_mass_map[80240000] = 24.019861;
  atomic_mass_map[90240000] = 24.00809937;
  atomic_mass_map[100240000] = 23.993610645;
  atomic_mass_map[110240000] = 23.990963011;
  atomic_mass_map[120240000] = 23.985041697;
  atomic_mass_map[130240000] = 23.999947541;
  atomic_mass_map[140240000] = 24.011535441;
  atomic_mass_map[150240000] = 24.03577;
  atomic_mass_map[70250000] = 25.0601;
  atomic_mass_map[80250000] = 25.029338919;
  atomic_mass_map[90250000] = 25.012167727;
  atomic_mass_map[100250000] = 24.997814799;
  atomic_mass_map[110250000] = 24.989953973;
  atomic_mass_map[120250000] = 24.985836964;
  atomic_mass_map[130250000] = 24.990428306;
  atomic_mass_map[140250000] = 25.004108801;
  atomic_mass_map[150250000] = 25.02119;
  atomic_mass_map[80260000] = 26.037210155;
  atomic_mass_map[90260000] = 26.020020392;
  atomic_mass_map[100260000] = 26.000516496;
  atomic_mass_map[110260000] = 25.992634649;
  atomic_mass_map[120260000] = 25.982592971;
  atomic_mass_map[130260000] = 25.986891863;
  atomic_mass_map[140260000] = 25.992333804;
  atomic_mass_map[150260000] = 26.01178;
  atomic_mass_map[160260000] = 26.02907;
  atomic_mass_map[80270000] = 27.047955;
  atomic_mass_map[90270000] = 27.027322;
  atomic_mass_map[100270000] = 27.007569462;
  atomic_mass_map[110270000] = 26.994076408;
  atomic_mass_map[120270000] = 26.984340627999998;
  atomic_mass_map[130270000] = 26.981538408;
  atomic_mass_map[140270000] = 26.986704688;
  atomic_mass_map[150270000] = 26.999224409;
  atomic_mass_map[160270000] = 27.01828;
  atomic_mass_map[80280000] = 28.05591;
  atomic_mass_map[90280000] = 28.036223095;
  atomic_mass_map[100280000] = 28.012130767;
  atomic_mass_map[110280000] = 27.998939;
  atomic_mass_map[120280000] = 27.983876606;
  atomic_mass_map[130280000] = 27.981910087;
  atomic_mass_map[140280000] = 27.97692653499;
  atomic_mass_map[150280000] = 27.992326585;
  atomic_mass_map[160280000] = 28.004372766;
  atomic_mass_map[170280000] = 28.02954;
  atomic_mass_map[90290000] = 29.043103;
  atomic_mass_map[100290000] = 29.019753;
  atomic_mass_map[110290000] = 29.002877092;
  atomic_mass_map[120290000] = 28.988617393;
  atomic_mass_map[130290000] = 28.980453164;
  atomic_mass_map[140290000] = 28.97649466525;
  atomic_mass_map[150290000] = 28.981800368;
  atomic_mass_map[160290000] = 28.996611448;
  atomic_mass_map[170290000] = 29.014130178;
  atomic_mass_map[90300000] = 30.05165;
  atomic_mass_map[100300000] = 30.024992235;
  atomic_mass_map[110300000] = 30.009097932;
  atomic_mass_map[120300000] = 29.990462825999998;
  atomic_mass_map[130300000] = 29.982968388;
  atomic_mass_map[140300000] = 29.973770136;
  atomic_mass_map[150300000] = 29.978313489;
  atomic_mass_map[160300000] = 29.984906769;
  atomic_mass_map[170300000] = 30.00477;
  atomic_mass_map[180300000] = 30.022470511;
  atomic_mass_map[90310000] = 31.060272;
  atomic_mass_map[100310000] = 31.033474816;
  atomic_mass_map[110310000] = 31.013146656;
  atomic_mass_map[120310000] = 30.996648232;
  atomic_mass_map[130310000] = 30.983949756;
  atomic_mass_map[140310000] = 30.975363194;
  atomic_mass_map[150310000] = 30.97376199863;
  atomic_mass_map[160310000] = 30.979557007;
  atomic_mass_map[170310000] = 30.992448098;
  atomic_mass_map[180310000] = 31.012158;
  atomic_mass_map[100320000] = 32.03972;
  atomic_mass_map[110320000] = 32.020011026;
  atomic_mass_map[120320000] = 31.999110139;
  atomic_mass_map[130320000] = 31.988084339;
  atomic_mass_map[140320000] = 31.974151539;
  atomic_mass_map[150320000] = 31.973907643;
  atomic_mass_map[160320000] = 31.97207117443;
  atomic_mass_map[170320000] = 31.985684637;
  atomic_mass_map[180320000] = 31.997637826;
  atomic_mass_map[190320000] = 32.02265;
  atomic_mass_map[100330000] = 33.04938;
  atomic_mass_map[110330000] = 33.025529;
  atomic_mass_map[120330000] = 33.005327245;
  atomic_mass_map[130330000] = 32.990877687;
  atomic_mass_map[140330000] = 32.977976964;
  atomic_mass_map[150330000] = 32.971725694;
  atomic_mass_map[160330000] = 32.97145890985;
  atomic_mass_map[170330000] = 32.977451989;
  atomic_mass_map[180330000] = 32.989925547;
  atomic_mass_map[190330000] = 33.00756;
  atomic_mass_map[100340000] = 34.056728;
  atomic_mass_map[110340000] = 34.03401;
  atomic_mass_map[120340000] = 34.008935481;
  atomic_mass_map[130340000] = 33.996779057;
  atomic_mass_map[140340000] = 33.978575437;
  atomic_mass_map[150340000] = 33.973645887;
  atomic_mass_map[160340000] = 33.967867012;
  atomic_mass_map[170340000] = 33.973762491;
  atomic_mass_map[180340000] = 33.980270093;
  atomic_mass_map[190340000] = 33.99869;
  atomic_mass_map[200340000] = 34.01487;
  atomic_mass_map[110350000] = 35.041043;
  atomic_mass_map[120350000] = 35.01679;
  atomic_mass_map[130350000] = 34.999759817;
  atomic_mass_map[140350000] = 34.984550134;
  atomic_mass_map[150350000] = 34.973314053;
  atomic_mass_map[160350000] = 34.969032322;
  atomic_mass_map[170350000] = 34.968852694;
  atomic_mass_map[180350000] = 34.975257721;
  atomic_mass_map[190350000] = 34.988005407;
  atomic_mass_map[200350000] = 35.00514;
  atomic_mass_map[110360000] = 36.049708;
  atomic_mass_map[120360000] = 36.021879;
  atomic_mass_map[130360000] = 36.006388;
  atomic_mass_map[140360000] = 35.986649271;
  atomic_mass_map[150360000] = 35.978259619;
  atomic_mass_map[160360000] = 35.967080699;
  atomic_mass_map[170360000] = 35.968306822;
  atomic_mass_map[180360000] = 35.967545105;
  atomic_mass_map[190360000] = 35.98130201;
  atomic_mass_map[200360000] = 35.993074406;
  atomic_mass_map[210360000] = 36.01648;
  atomic_mass_map[110370000] = 37.057471;
  atomic_mass_map[120370000] = 37.030286265;
  atomic_mass_map[130370000] = 37.010531;
  atomic_mass_map[140370000] = 36.992945191;
  atomic_mass_map[150370000] = 36.979606956;
  atomic_mass_map[160370000] = 36.971125506999996;
  atomic_mass_map[170370000] = 36.965902584;
  atomic_mass_map[180370000] = 36.966776314;
  atomic_mass_map[190370000] = 36.973375889;
  atomic_mass_map[200370000] = 36.985897852;
  atomic_mass_map[210370000] = 37.003779;
  atomic_mass_map[120380000] = 38.03658;
  atomic_mass_map[130380000] = 38.017402;
  atomic_mass_map[140380000] = 37.995523;
  atomic_mass_map[150380000] = 37.984303105;
  atomic_mass_map[160380000] = 37.97116331;
  atomic_mass_map[170380000] = 37.968010418;
  atomic_mass_map[180380000] = 37.962732104;
  atomic_mass_map[190380000] = 37.969081116;
  atomic_mass_map[200380000] = 37.976319226;
  atomic_mass_map[210380000] = 37.995438;
  atomic_mass_map[220380000] = 38.011669;
  atomic_mass_map[120390000] = 39.045384;
  atomic_mass_map[130390000] = 39.022169;
  atomic_mass_map[140390000] = 39.002491;
  atomic_mass_map[150390000] = 38.986285865;
  atomic_mass_map[160390000] = 38.975133852;
  atomic_mass_map[170390000] = 38.968008162;
  atomic_mass_map[180390000] = 38.964313039;
  atomic_mass_map[190390000] = 38.96370648661;
  atomic_mass_map[200390000] = 38.970710813;
  atomic_mass_map[210390000] = 38.98478497;
  atomic_mass_map[220390000] = 39.002362;
  atomic_mass_map[120400000] = 40.051906;
  atomic_mass_map[130400000] = 40.029619;
  atomic_mass_map[140400000] = 40.005829;
  atomic_mass_map[150400000] = 39.991288865;
  atomic_mass_map[160400000] = 39.975482562;
  atomic_mass_map[170400000] = 39.970415469;
  atomic_mass_map[180400000] = 39.96238312378;
  atomic_mass_map[190400000] = 39.963998166;
  atomic_mass_map[200400000] = 39.962590865;
  atomic_mass_map[210400000] = 39.977967292;
  atomic_mass_map[220400000] = 39.990498721;
  atomic_mass_map[230400000] = 40.013065;
  atomic_mass_map[130410000] = 41.035878;
  atomic_mass_map[140410000] = 41.013011;
  atomic_mass_map[150410000] = 40.994654;
  atomic_mass_map[160410000] = 40.979593451;
  atomic_mass_map[170410000] = 40.970684525;
  atomic_mass_map[180410000] = 40.964500571;
  atomic_mass_map[190410000] = 40.96182525796;
  atomic_mass_map[200410000] = 40.962277921;
  atomic_mass_map[210410000] = 40.969251104;
  atomic_mass_map[220410000] = 40.983148;
  atomic_mass_map[230410000] = 41.000344;
  atomic_mass_map[130420000] = 42.043049;
  atomic_mass_map[140420000] = 42.017681;
  atomic_mass_map[150420000] = 42.001084;
  atomic_mass_map[160420000] = 41.9810651;
  atomic_mass_map[170420000] = 41.973342;
  atomic_mass_map[180420000] = 41.963045736;
  atomic_mass_map[190420000] = 41.962402306;
  atomic_mass_map[200420000] = 41.958617828;
  atomic_mass_map[210420000] = 41.965516522;
  atomic_mass_map[220420000] = 41.973049022;
  atomic_mass_map[230420000] = 41.99182;
  atomic_mass_map[240420000] = 42.007225;
  atomic_mass_map[130430000] = 43.050478;
  atomic_mass_map[140430000] = 43.0248;
  atomic_mass_map[150430000] = 43.005024;
  atomic_mass_map[160430000] = 42.986907635;
  atomic_mass_map[170430000] = 42.9740637;
  atomic_mass_map[180430000] = 42.965636055;
  atomic_mass_map[190430000] = 42.960734703;
  atomic_mass_map[200430000] = 42.95876643;
  atomic_mass_map[210430000] = 42.961150472;
  atomic_mass_map[220430000] = 42.968522521;
  atomic_mass_map[230430000] = 42.980766;
  atomic_mass_map[240430000] = 42.997885;
  atomic_mass_map[140440000] = 44.03061;
  atomic_mass_map[150440000] = 44.011219;
  atomic_mass_map[160440000] = 43.990118848;
  atomic_mass_map[170440000] = 43.978116312;
  atomic_mass_map[180440000] = 43.964923816;
  atomic_mass_map[190440000] = 43.961586986;
  atomic_mass_map[200440000] = 43.955481543;
  atomic_mass_map[210440000] = 43.959402867;
  atomic_mass_map[220440000] = 43.959689951;
  atomic_mass_map[230440000] = 43.97411;
  atomic_mass_map[240440000] = 43.985657;
  atomic_mass_map[250440000] = 44.007547;
  atomic_mass_map[140450000] = 45.040247;
  atomic_mass_map[150450000] = 45.016747;
  atomic_mass_map[160450000] = 44.995717;
  atomic_mass_map[170450000] = 44.980394353;
  atomic_mass_map[180450000] = 44.968039733;
  atomic_mass_map[190450000] = 44.960691493;
  atomic_mass_map[200450000] = 44.956186326;
  atomic_mass_map[210450000] = 44.955907503;
  atomic_mass_map[220450000] = 44.958121211;
  atomic_mass_map[230450000] = 44.965768951;
  atomic_mass_map[240450000] = 44.97905;
  atomic_mass_map[250450000] = 44.994364;
  atomic_mass_map[260450000] = 45.014774;
  atomic_mass_map[150460000] = 46.024659;
  atomic_mass_map[160460000] = 46.000365;
  atomic_mass_map[170460000] = 45.985121323;
  atomic_mass_map[180460000] = 45.968037446;
  atomic_mass_map[190460000] = 45.961981586;
  atomic_mass_map[200460000] = 45.953687988;
  atomic_mass_map[210460000] = 45.955167485;
  atomic_mass_map[220460000] = 45.952626856;
  atomic_mass_map[230460000] = 45.960197971;
  atomic_mass_map[240460000] = 45.96836097;
  atomic_mass_map[250460000] = 45.986506;
  atomic_mass_map[260460000] = 46.000977;
  atomic_mass_map[150470000] = 47.031895;
  atomic_mass_map[160470000] = 47.007912;
  atomic_mass_map[170470000] = 46.989501;
  atomic_mass_map[180470000] = 46.972768114;
  atomic_mass_map[190470000] = 46.961661614;
  atomic_mass_map[200470000] = 46.954541394;
  atomic_mass_map[210470000] = 46.952402704;
  atomic_mass_map[220470000] = 46.951757752;
  atomic_mass_map[230470000] = 46.954904038;
  atomic_mass_map[240470000] = 46.962895544;
  atomic_mass_map[250470000] = 46.975774;
  atomic_mass_map[260470000] = 46.992625;
  atomic_mass_map[270470000] = 47.011133;
  atomic_mass_map[160480000] = 48.0137;
  atomic_mass_map[170480000] = 47.995405;
  atomic_mass_map[180480000] = 47.97608;
  atomic_mass_map[190480000] = 47.965341186;
  atomic_mass_map[200480000] = 47.952522904;
  atomic_mass_map[210480000] = 47.952223157;
  atomic_mass_map[220480000] = 47.947940932;
  atomic_mass_map[230480000] = 47.952251229;
  atomic_mass_map[240480000] = 47.954028667;
  atomic_mass_map[250480000] = 47.968549085;
  atomic_mass_map[260480000] = 47.980676;
  atomic_mass_map[270480000] = 48.00161;
  atomic_mass_map[280480000] = 48.018028;
  atomic_mass_map[160490000] = 49.022644;
  atomic_mass_map[170490000] = 49.001009;
  atomic_mass_map[180490000] = 48.981546;
  atomic_mass_map[190490000] = 48.968210755;
  atomic_mass_map[200490000] = 48.955662875;
  atomic_mass_map[210490000] = 48.950014423;
  atomic_mass_map[220490000] = 48.947864627;
  atomic_mass_map[230490000] = 48.948510746;
  atomic_mass_map[240490000] = 48.951332955;
  atomic_mass_map[250490000] = 48.959612585;
  atomic_mass_map[260490000] = 48.973429;
  atomic_mass_map[270490000] = 48.989393;
  atomic_mass_map[280490000] = 49.008803;
  atomic_mass_map[170500000] = 50.008309;
  atomic_mass_map[180500000] = 49.98569;
  atomic_mass_map[190500000] = 49.972380017;
  atomic_mass_map[200500000] = 49.957499217;
  atomic_mass_map[210500000] = 49.952176415;
  atomic_mass_map[220500000] = 49.944785839;
  atomic_mass_map[230500000] = 49.947155845;
  atomic_mass_map[240500000] = 49.946041443;
  atomic_mass_map[250500000] = 49.954237391;
  atomic_mass_map[260500000] = 49.962988;
  atomic_mass_map[270500000] = 49.981073;
  atomic_mass_map[280500000] = 49.995577;
  atomic_mass_map[170510000] = 51.015341;
  atomic_mass_map[180510000] = 50.992818;
  atomic_mass_map[190510000] = 50.975827867;
  atomic_mass_map[200510000] = 50.960995665;
  atomic_mass_map[210510000] = 50.953592095;
  atomic_mass_map[220510000] = 50.9466096;
  atomic_mass_map[230510000] = 50.943956867;
  atomic_mass_map[240510000] = 50.944764652;
  atomic_mass_map[250510000] = 50.948208065;
  atomic_mass_map[260510000] = 50.956840779;
  atomic_mass_map[270510000] = 50.970647;
  atomic_mass_map[280510000] = 50.987225;
  atomic_mass_map[180520000] = 51.998626;
  atomic_mass_map[190520000] = 51.981602;
  atomic_mass_map[200520000] = 51.963213648;
  atomic_mass_map[210520000] = 51.956582351;
  atomic_mass_map[220520000] = 51.94689196;
  atomic_mass_map[230520000] = 51.944772839;
  atomic_mass_map[240520000] = 51.940504992;
  atomic_mass_map[250520000] = 51.945563488;
  atomic_mass_map[260520000] = 51.948115217;
  atomic_mass_map[270520000] = 51.963112;
  atomic_mass_map[280520000] = 51.976028;
  atomic_mass_map[290520000] = 51.997552;
  atomic_mass_map[180530000] = 53.00729;
  atomic_mass_map[190530000] = 52.9868;
  atomic_mass_map[200530000] = 52.968451;
  atomic_mass_map[210530000] = 52.958231821;
  atomic_mass_map[220530000] = 52.949724785;
  atomic_mass_map[230530000] = 52.944335593;
  atomic_mass_map[240530000] = 52.940646961;
  atomic_mass_map[250530000] = 52.941287742;
  atomic_mass_map[260530000] = 52.945305574;
  atomic_mass_map[270530000] = 52.954203217;
  atomic_mass_map[280530000] = 52.96819;
  atomic_mass_map[290530000] = 52.985754;
  atomic_mass_map[190540000] = 53.99463;
  atomic_mass_map[200540000] = 53.972989;
  atomic_mass_map[210540000] = 53.96361662;
  atomic_mass_map[220540000] = 53.951022786;
  atomic_mass_map[230540000] = 53.946437472;
  atomic_mass_map[240540000] = 53.938878012;
  atomic_mass_map[250540000] = 53.940356429;
  atomic_mass_map[260540000] = 53.939608306;
  atomic_mass_map[270540000] = 53.948459192;
  atomic_mass_map[280540000] = 53.957833;
  atomic_mass_map[290540000] = 53.977015;
  atomic_mass_map[300540000] = 53.993267;
  atomic_mass_map[190550000] = 55.00076;
  atomic_mass_map[200550000] = 54.9803;
  atomic_mass_map[210550000] = 54.967622601;
  atomic_mass_map[220550000] = 54.955267465;
  atomic_mass_map[230550000] = 54.947241114;
  atomic_mass_map[240550000] = 54.940837289;
  atomic_mass_map[250550000] = 54.938043172;
  atomic_mass_map[260550000] = 54.938291283;
  atomic_mass_map[270550000] = 54.941996531;
  atomic_mass_map[280550000] = 54.951329961;
  atomic_mass_map[290550000] = 54.966038;
  atomic_mass_map[300550000] = 54.984358;
  atomic_mass_map[190560000] = 56.00851;
  atomic_mass_map[200560000] = 55.985079999999996;
  atomic_mass_map[210560000] = 55.97332;
  atomic_mass_map[220560000] = 55.95778819;
  atomic_mass_map[230560000] = 55.950450694;
  atomic_mass_map[240560000] = 55.940649107;
  atomic_mass_map[250560000] = 55.938902947;
  atomic_mass_map[260560000] = 55.934935617;
  atomic_mass_map[270560000] = 55.93983815;
  atomic_mass_map[280560000] = 55.942127872;
  atomic_mass_map[290560000] = 55.958515;
  atomic_mass_map[300560000] = 55.972743;
  atomic_mass_map[310560000] = 55.996361;
  atomic_mass_map[200570000] = 56.99262;
  atomic_mass_map[210570000] = 56.97746;
  atomic_mass_map[220570000] = 56.963590068;
  atomic_mass_map[230570000] = 56.952320197;
  atomic_mass_map[240570000] = 56.943612408999996;
  atomic_mass_map[250570000] = 56.938285968;
  atomic_mass_map[260570000] = 56.935392134;
  atomic_mass_map[270570000] = 56.936289913;
  atomic_mass_map[280570000] = 56.939791525;
  atomic_mass_map[290570000] = 56.949211819;
  atomic_mass_map[300570000] = 56.965056;
  atomic_mass_map[310570000] = 56.983886;
  atomic_mass_map[200580000] = 57.99794;
  atomic_mass_map[210580000] = 57.98403;
  atomic_mass_map[220580000] = 57.966602;
  atomic_mass_map[230580000] = 57.956626932;
  atomic_mass_map[240580000] = 57.944184502;
  atomic_mass_map[250580000] = 57.940066646;
  atomic_mass_map[260580000] = 57.933273738;
  atomic_mass_map[270580000] = 57.935751429;
  atomic_mass_map[280580000] = 57.93534178;
  atomic_mass_map[290580000] = 57.944532413;
  atomic_mass_map[300580000] = 57.954590428;
  atomic_mass_map[310580000] = 57.974728999999996;
  atomic_mass_map[320580000] = 57.992399;
  atomic_mass_map[210590000] = 58.98894;
  atomic_mass_map[220590000] = 58.972614;
  atomic_mass_map[230590000] = 58.959385659;
  atomic_mass_map[240590000] = 58.94837781;
  atomic_mass_map[250590000] = 58.940391113;
  atomic_mass_map[260590000] = 58.934873649;
  atomic_mass_map[270590000] = 58.933193656;
  atomic_mass_map[280590000] = 58.934345571;
  atomic_mass_map[290590000] = 58.939496844;
  atomic_mass_map[300590000] = 58.949312017;
  atomic_mass_map[310590000] = 58.963757;
  atomic_mass_map[320590000] = 58.982963;
  atomic_mass_map[210600000] = 59.99565;
  atomic_mass_map[220600000] = 59.976028;
  atomic_mass_map[230600000] = 59.96431329;
  atomic_mass_map[240600000] = 59.949898146;
  atomic_mass_map[250600000] = 59.943136576;
  atomic_mass_map[260600000] = 59.934070411;
  atomic_mass_map[270600000] = 59.933815667;
  atomic_mass_map[280600000] = 59.930785256;
  atomic_mass_map[290600000] = 59.937363916;
  atomic_mass_map[300600000] = 59.94184145;
  atomic_mass_map[310600000] = 59.957498;
  atomic_mass_map[320600000] = 59.970918;
  atomic_mass_map[330600000] = 59.994128;
  atomic_mass_map[210610000] = 61.001;
  atomic_mass_map[220610000] = 60.982448;
  atomic_mass_map[230610000] = 60.96725;
  atomic_mass_map[240610000] = 60.954400963;
  atomic_mass_map[250610000] = 60.944452544;
  atomic_mass_map[260610000] = 60.936746244;
  atomic_mass_map[270610000] = 60.932476145;
  atomic_mass_map[280610000] = 60.931054945;
  atomic_mass_map[290610000] = 60.933457371;
  atomic_mass_map[300610000] = 60.93950696;
  atomic_mass_map[310610000] = 60.949398859;
  atomic_mass_map[320610000] = 60.964187;
  atomic_mass_map[330610000] = 60.981857;
  atomic_mass_map[220620000] = 61.986581;
  atomic_mass_map[230620000] = 61.97265;
  atomic_mass_map[240620000] = 61.956097451;
  atomic_mass_map[250620000] = 61.947907386;
  atomic_mass_map[260620000] = 61.936791812;
  atomic_mass_map[270620000] = 61.934058317;
  atomic_mass_map[280620000] = 61.928344871;
  atomic_mass_map[290620000] = 61.932594921;
  atomic_mass_map[300620000] = 61.934333477;
  atomic_mass_map[310620000] = 61.944189757;
  atomic_mass_map[320620000] = 61.95519;
  atomic_mass_map[330620000] = 61.973891;
  atomic_mass_map[220630000] = 62.993827;
  atomic_mass_map[230630000] = 62.9765;
  atomic_mass_map[240630000] = 62.961344384;
  atomic_mass_map[250630000] = 62.949664675;
  atomic_mass_map[260630000] = 62.9402727;
  atomic_mass_map[270630000] = 62.933599744;
  atomic_mass_map[280630000] = 62.929669139;
  atomic_mass_map[290630000] = 62.929597236;
  atomic_mass_map[300630000] = 62.933211167;
  atomic_mass_map[310630000] = 62.939294195;
  atomic_mass_map[320630000] = 62.949628;
  atomic_mass_map[330630000] = 62.964036;
  atomic_mass_map[220640000] = 63.9989;
  atomic_mass_map[230640000] = 63.98248;
  atomic_mass_map[240640000] = 63.964058;
  atomic_mass_map[250640000] = 63.95384937;
  atomic_mass_map[260640000] = 63.940987763;
  atomic_mass_map[270640000] = 63.935810291;
  atomic_mass_map[280640000] = 63.927966341;
  atomic_mass_map[290640000] = 63.929763857;
  atomic_mass_map[300640000] = 63.929141772;
  atomic_mass_map[310640000] = 63.936840365;
  atomic_mass_map[320640000] = 63.941689913;
  atomic_mass_map[330640000] = 63.95756;
  atomic_mass_map[340640000] = 63.971336;
  atomic_mass_map[230650000] = 64.987354;
  atomic_mass_map[240650000] = 64.969705;
  atomic_mass_map[250650000] = 64.95601975;
  atomic_mass_map[260650000] = 64.945015324;
  atomic_mass_map[270650000] = 64.936462073;
  atomic_mass_map[280650000] = 64.930084697;
  atomic_mass_map[290650000] = 64.927789487;
  atomic_mass_map[300650000] = 64.929240532;
  atomic_mass_map[310650000] = 64.932734395;
  atomic_mass_map[320650000] = 64.939368137;
  atomic_mass_map[330650000] = 64.949611;
  atomic_mass_map[340650000] = 64.964552;
  atomic_mass_map[230660000] = 65.993977;
  atomic_mass_map[240660000] = 65.973462;
  atomic_mass_map[250660000] = 65.960546834;
  atomic_mass_map[260660000] = 65.94624996;
  atomic_mass_map[270660000] = 65.939442945;
  atomic_mass_map[280660000] = 65.929139334;
  atomic_mass_map[290660000] = 65.928868814;
  atomic_mass_map[300660000] = 65.926033704;
  atomic_mass_map[310660000] = 65.931589832;
  atomic_mass_map[320660000] = 65.933862126;
  atomic_mass_map[330660000] = 65.944148779;
  atomic_mass_map[340660000] = 65.955276;
  atomic_mass_map[230670000] = 66.999302;
  atomic_mass_map[240670000] = 66.979946;
  atomic_mass_map[250670000] = 66.964079;
  atomic_mass_map[260670000] = 66.951035482;
  atomic_mass_map[270670000] = 66.940609628;
  atomic_mass_map[280670000] = 66.931569414;
  atomic_mass_map[290670000] = 66.927729526;
  atomic_mass_map[300670000] = 66.927127482;
  atomic_mass_map[310670000] = 66.928202384;
  atomic_mass_map[320670000] = 66.93273362;
  atomic_mass_map[330670000] = 66.939251111;
  atomic_mass_map[340670000] = 66.949994;
  atomic_mass_map[350670000] = 66.964798;
  atomic_mass_map[240680000] = 67.984112;
  atomic_mass_map[250680000] = 67.969533;
  atomic_mass_map[260680000] = 67.953314875;
  atomic_mass_map[270680000] = 67.944250135;
  atomic_mass_map[280680000] = 67.931868789;
  atomic_mass_map[290680000] = 67.929610889;
  atomic_mass_map[300680000] = 67.924844291;
  atomic_mass_map[310680000] = 67.927980221;
  atomic_mass_map[320680000] = 67.928095308;
  atomic_mass_map[330680000] = 67.93677413;
  atomic_mass_map[340680000] = 67.941825239;
  atomic_mass_map[350680000] = 67.958356;
  atomic_mass_map[240690000] = 68.990789;
  atomic_mass_map[250690000] = 68.973408;
  atomic_mass_map[260690000] = 68.9581;
  atomic_mass_map[270690000] = 68.946023102;
  atomic_mass_map[280690000] = 68.935610268;
  atomic_mass_map[290690000] = 68.929429268;
  atomic_mass_map[300690000] = 68.926550418;
  atomic_mass_map[310690000] = 68.925573531;
  atomic_mass_map[320690000] = 68.927964471;
  atomic_mass_map[330690000] = 68.932246294;
  atomic_mass_map[340690000] = 68.939414847;
  atomic_mass_map[350690000] = 68.950338413;
  atomic_mass_map[360690000] = 68.96518;
  atomic_mass_map[240700000] = 69.995191;
  atomic_mass_map[250700000] = 69.979066;
  atomic_mass_map[260700000] = 69.960805;
  atomic_mass_map[270700000] = 69.949941;
  atomic_mass_map[280700000] = 69.936431303;
  atomic_mass_map[290700000] = 69.932392079;
  atomic_mass_map[300700000] = 69.925319181;
  atomic_mass_map[310700000] = 69.926021917;
  atomic_mass_map[320700000] = 69.924248706;
  atomic_mass_map[330700000] = 69.930926151;
  atomic_mass_map[340700000] = 69.933515523;
  atomic_mass_map[350700000] = 69.944792323;
  atomic_mass_map[360700000] = 69.955877;
  atomic_mass_map[250710000] = 70.983285;
  atomic_mass_map[260710000] = 70.966259;
  atomic_mass_map[270710000] = 70.952366923;
  atomic_mass_map[280710000] = 70.940518964;
  atomic_mass_map[290710000] = 70.932676832;
  atomic_mass_map[300710000] = 70.92771958;
  atomic_mass_map[310710000] = 70.924702536;
  atomic_mass_map[320710000] = 70.924952284;
  atomic_mass_map[330710000] = 70.927113758;
  atomic_mass_map[340710000] = 70.932209432;
  atomic_mass_map[350710000] = 70.939342156;
  atomic_mass_map[360710000] = 70.950265696;
  atomic_mass_map[370710000] = 70.965582;
  atomic_mass_map[250720000] = 71.989372;
  atomic_mass_map[260720000] = 71.969479;
  atomic_mass_map[270720000] = 71.956844;
  atomic_mass_map[280720000] = 71.941785926;
  atomic_mass_map[290720000] = 71.935820307;
  atomic_mass_map[300720000] = 71.926842807;
  atomic_mass_map[310720000] = 71.926367434;
  atomic_mass_map[320720000] = 71.922075826;
  atomic_mass_map[330720000] = 71.926752295;
  atomic_mass_map[340720000] = 71.927140507;
  atomic_mass_map[350720000] = 71.936594607;
  atomic_mass_map[360720000] = 71.942092407;
  atomic_mass_map[370720000] = 71.958851;
  atomic_mass_map[260730000] = 72.975416;
  atomic_mass_map[270730000] = 72.95983;
  atomic_mass_map[280730000] = 72.946206683;
  atomic_mass_map[290730000] = 72.936674378;
  atomic_mass_map[300730000] = 72.929582582;
  atomic_mass_map[310730000] = 72.925174682;
  atomic_mass_map[320730000] = 72.923458956;
  atomic_mass_map[330730000] = 72.923829089;
  atomic_mass_map[340730000] = 72.926754883;
  atomic_mass_map[350730000] = 72.931671621;
  atomic_mass_map[360730000] = 72.939289195;
  atomic_mass_map[370730000] = 72.950529;
  atomic_mass_map[380730000] = 72.9657;
  atomic_mass_map[260740000] = 73.978969;
  atomic_mass_map[270740000] = 73.964766;
  atomic_mass_map[280740000] = 73.94798;
  atomic_mass_map[290740000] = 73.939874862;
  atomic_mass_map[300740000] = 73.929407262;
  atomic_mass_map[310740000] = 73.926945726;
  atomic_mass_map[320740000] = 73.921177762;
  atomic_mass_map[330740000] = 73.923928598;
  atomic_mass_map[340740000] = 73.922475935;
  atomic_mass_map[350740000] = 73.929910281;
  atomic_mass_map[360740000] = 73.933084017;
  atomic_mass_map[370740000] = 73.944265868;
  atomic_mass_map[380740000] = 73.95617;
  atomic_mass_map[260750000] = 74.985357;
  atomic_mass_map[270750000] = 74.96817;
  atomic_mass_map[280750000] = 74.952732;
  atomic_mass_map[290750000] = 74.941522606;
  atomic_mass_map[300750000] = 74.932840246;
  atomic_mass_map[310750000] = 74.926500246;
  atomic_mass_map[320750000] = 74.922858371;
  atomic_mass_map[330750000] = 74.921594562;
  atomic_mass_map[340750000] = 74.922522871;
  atomic_mass_map[350750000] = 74.92581057;
  atomic_mass_map[360750000] = 74.930945746;
  atomic_mass_map[370750000] = 74.938573201;
  atomic_mass_map[380750000] = 74.94995277;
  atomic_mass_map[390750000] = 74.96584;
  atomic_mass_map[270760000] = 75.973687;
  atomic_mass_map[280760000] = 75.955308;
  atomic_mass_map[290760000] = 75.945275025;
  atomic_mass_map[300760000] = 75.933114957;
  atomic_mass_map[310760000] = 75.928827625;
  atomic_mass_map[320760000] = 75.921402726;
  atomic_mass_map[330760000] = 75.92239201;
  atomic_mass_map[340760000] = 75.919213704;
  atomic_mass_map[350760000] = 75.924541577;
  atomic_mass_map[360760000] = 75.925910726;
  atomic_mass_map[370760000] = 75.935073032;
  atomic_mass_map[380760000] = 75.941762761;
  atomic_mass_map[390760000] = 75.95869;
  atomic_mass_map[270770000] = 76.97744;
  atomic_mass_map[280770000] = 76.960494;
  atomic_mass_map[290770000] = 76.9478;
  atomic_mass_map[300770000] = 76.936887199;
  atomic_mass_map[310770000] = 76.9291543;
  atomic_mass_map[320770000] = 76.923549844;
  atomic_mass_map[330770000] = 76.920647564;
  atomic_mass_map[340770000] = 76.91991415;
  atomic_mass_map[350770000] = 76.921379194;
  atomic_mass_map[360770000] = 76.92467;
  atomic_mass_map[370770000] = 76.9304016;
  atomic_mass_map[380770000] = 76.937945455;
  atomic_mass_map[390770000] = 76.950146;
  atomic_mass_map[400770000] = 76.965604;
  atomic_mass_map[280780000] = 77.963618;
  atomic_mass_map[290780000] = 77.95223;
  atomic_mass_map[300780000] = 77.938289205;
  atomic_mass_map[310780000] = 77.931608845;
  atomic_mass_map[320780000] = 77.922852912;
  atomic_mass_map[330780000] = 77.921827795;
  atomic_mass_map[340780000] = 77.917309243;
  atomic_mass_map[350780000] = 77.921145859;
  atomic_mass_map[360780000] = 77.920366341;
  atomic_mass_map[370780000] = 77.928141868;
  atomic_mass_map[380780000] = 77.93217998;
  atomic_mass_map[390780000] = 77.94399;
  atomic_mass_map[400780000] = 77.956146;
  atomic_mass_map[280790000] = 78.970402;
  atomic_mass_map[290790000] = 78.95519;
  atomic_mass_map[300790000] = 78.942638068;
  atomic_mass_map[310790000] = 78.932852301;
  atomic_mass_map[320790000] = 78.925360129;
  atomic_mass_map[330790000] = 78.920948445;
  atomic_mass_map[340790000] = 78.918499251;
  atomic_mass_map[350790000] = 78.918337601;
  atomic_mass_map[360790000] = 78.920082945;
  atomic_mass_map[370790000] = 78.923989864;
  atomic_mass_map[380790000] = 78.929707664;
  atomic_mass_map[390790000] = 78.93793;
  atomic_mass_map[400790000] = 78.94979;
  atomic_mass_map[410790000] = 78.966022;
  atomic_mass_map[280800000] = 79.975706;
  atomic_mass_map[290800000] = 79.961138;
  atomic_mass_map[300800000] = 79.94455293;
  atomic_mass_map[310800000] = 79.936420774;
  atomic_mass_map[320800000] = 79.925350774;
  atomic_mass_map[330800000] = 79.922474548;
  atomic_mass_map[340800000] = 79.916521785;
  atomic_mass_map[350800000] = 79.91852981;
  atomic_mass_map[360800000] = 79.916378048;
  atomic_mass_map[370800000] = 79.922516444;
  atomic_mass_map[380800000] = 79.92451754;
  atomic_mass_map[390800000] = 79.934354755;
  atomic_mass_map[400800000] = 79.941642;
  atomic_mass_map[410800000] = 79.958754;
  atomic_mass_map[290810000] = 80.966269;
  atomic_mass_map[300810000] = 80.950402619;
  atomic_mass_map[310810000] = 80.938133842;
  atomic_mass_map[320810000] = 80.928832942;
  atomic_mass_map[330810000] = 80.92213229;
  atomic_mass_map[340810000] = 80.917993044;
  atomic_mass_map[350810000] = 80.916288206;
  atomic_mass_map[360810000] = 80.916589714;
  atomic_mass_map[370810000] = 80.918993927;
  atomic_mass_map[380810000] = 80.923211394;
  atomic_mass_map[390810000] = 80.929454283;
  atomic_mass_map[400810000] = 80.938314;
  atomic_mass_map[410810000] = 80.95023;
  atomic_mass_map[420810000] = 80.965915;
  atomic_mass_map[290820000] = 81.972818;
  atomic_mass_map[300820000] = 81.954574099;
  atomic_mass_map[310820000] = 81.943176533;
  atomic_mass_map[320820000] = 81.929774033;
  atomic_mass_map[330820000] = 81.924738733;
  atomic_mass_map[340820000] = 81.916699537;
  atomic_mass_map[350820000] = 81.91680176;
  atomic_mass_map[360820000] = 81.9134811552;
  atomic_mass_map[370820000] = 81.918209024;
  atomic_mass_map[380820000] = 81.918399847;
  atomic_mass_map[390820000] = 81.926930188;
  atomic_mass_map[400820000] = 81.931689;
  atomic_mass_map[410820000] = 81.944079;
  atomic_mass_map[420820000] = 81.956661;
  atomic_mass_map[300830000] = 82.961041;
  atomic_mass_map[310830000] = 82.947120301;
  atomic_mass_map[320830000] = 82.934539101;
  atomic_mass_map[330830000] = 82.925206901;
  atomic_mass_map[340830000] = 82.919118609;
  atomic_mass_map[350830000] = 82.915175289;
  atomic_mass_map[360830000] = 82.914126518;
  atomic_mass_map[370830000] = 82.915114182;
  atomic_mass_map[380830000] = 82.917554374;
  atomic_mass_map[390830000] = 82.922484025;
  atomic_mass_map[400830000] = 82.929240925;
  atomic_mass_map[410830000] = 82.938211;
  atomic_mass_map[420830000] = 82.950252;
  atomic_mass_map[430830000] = 82.966377;
  atomic_mass_map[300840000] = 83.965722;
  atomic_mass_map[310840000] = 83.95267;
  atomic_mass_map[320840000] = 83.937575091;
  atomic_mass_map[330840000] = 83.929303291;
  atomic_mass_map[340840000] = 83.918466762;
  atomic_mass_map[350840000] = 83.916496419;
  atomic_mass_map[360840000] = 83.91149772863;
  atomic_mass_map[370840000] = 83.914375225;
  atomic_mass_map[380840000] = 83.91341912;
  atomic_mass_map[390840000] = 83.920671061;
  atomic_mass_map[400840000] = 83.923325662;
  atomic_mass_map[410840000] = 83.934279;
  atomic_mass_map[420840000] = 83.941846;
  atomic_mass_map[430840000] = 83.959527;
  atomic_mass_map[300850000] = 84.972914;
  atomic_mass_map[310850000] = 84.95722;
  atomic_mass_map[320850000] = 84.942969659;
  atomic_mass_map[330850000] = 84.932163659;
  atomic_mass_map[340850000] = 84.922260759;
  atomic_mass_map[350850000] = 84.915645759;
  atomic_mass_map[360850000] = 84.912527262;
  atomic_mass_map[370850000] = 84.9117897376;
  atomic_mass_map[380850000] = 84.912932043;
  atomic_mass_map[390850000] = 84.916433039;
  atomic_mass_map[400850000] = 84.921443198;
  atomic_mass_map[410850000] = 84.928845837;
  atomic_mass_map[420850000] = 84.938260737;
  atomic_mass_map[430850000] = 84.950778;
  atomic_mass_map[440850000] = 84.966774;
  atomic_mass_map[310860000] = 85.963414;
  atomic_mass_map[320860000] = 85.946967;
  atomic_mass_map[330860000] = 85.936701533;
  atomic_mass_map[340860000] = 85.924311733;
  atomic_mass_map[350860000] = 85.918805433;
  atomic_mass_map[360860000] = 85.91061062627;
  atomic_mass_map[370860000] = 85.911167443;
  atomic_mass_map[380860000] = 85.90926072631;
  atomic_mass_map[390860000] = 85.914886098;
  atomic_mass_map[400860000] = 85.916296815;
  atomic_mass_map[410860000] = 85.925781535;
  atomic_mass_map[420860000] = 85.931174817;
  atomic_mass_map[430860000] = 85.944637;
  atomic_mass_map[440860000] = 85.957305;
  atomic_mass_map[310870000] = 86.968599;
  atomic_mass_map[320870000] = 86.95268;
  atomic_mass_map[330870000] = 86.940291718;
  atomic_mass_map[340870000] = 86.928688618;
  atomic_mass_map[350870000] = 86.920674018;
  atomic_mass_map[360870000] = 86.913354759;
  atomic_mass_map[370870000] = 86.909180531;
  atomic_mass_map[380870000] = 86.90887749615;
  atomic_mass_map[390870000] = 86.910876102;
  atomic_mass_map[400870000] = 86.914817339;
  atomic_mass_map[410870000] = 86.920692472;
  atomic_mass_map[420870000] = 86.928196201;
  atomic_mass_map[430870000] = 86.938067187;
  atomic_mass_map[440870000] = 86.951132;
  atomic_mass_map[320880000] = 87.95691;
  atomic_mass_map[330880000] = 87.94555;
  atomic_mass_map[340880000] = 87.931417491;
  atomic_mass_map[350880000] = 87.924083291;
  atomic_mass_map[360880000] = 87.914447881;
  atomic_mass_map[370880000] = 87.911315591;
  atomic_mass_map[380880000] = 87.90561225561;
  atomic_mass_map[390880000] = 87.909501276;
  atomic_mass_map[400880000] = 87.910220709;
  atomic_mass_map[410880000] = 87.918224287;
  atomic_mass_map[420880000] = 87.921967781;
  atomic_mass_map[430880000] = 87.933782381;
  atomic_mass_map[440880000] = 87.941664;
  atomic_mass_map[450880000] = 87.960429;
  atomic_mass_map[320890000] = 88.96379;
  atomic_mass_map[330890000] = 88.94976;
  atomic_mass_map[340890000] = 88.936669059;
  atomic_mass_map[350890000] = 88.926704559;
  atomic_mass_map[360890000] = 88.91783545;
  atomic_mass_map[370890000] = 88.912278137;
  atomic_mass_map[380890000] = 88.907450808;
  atomic_mass_map[390890000] = 88.905841205;
  atomic_mass_map[400890000] = 88.908882332;
  atomic_mass_map[410890000] = 88.913445272;
  atomic_mass_map[420890000] = 88.91946815;
  atomic_mass_map[430890000] = 88.92764865;
  atomic_mass_map[440890000] = 88.937455;
  atomic_mass_map[450890000] = 88.950767;
  atomic_mass_map[320900000] = 89.96863;
  atomic_mass_map[330900000] = 89.95563;
  atomic_mass_map[340900000] = 89.940096;
  atomic_mass_map[350900000] = 89.93129285;
  atomic_mass_map[360900000] = 89.91952793;
  atomic_mass_map[370900000] = 89.914798803;
  atomic_mass_map[380900000] = 89.907730885;
  atomic_mass_map[390900000] = 89.9071448;
  atomic_mass_map[400900000] = 89.904698758;
  atomic_mass_map[410900000] = 89.911259204;
  atomic_mass_map[420900000] = 89.913931272;
  atomic_mass_map[430900000] = 89.924073921;
  atomic_mass_map[440900000] = 89.930344379;
  atomic_mass_map[450900000] = 89.944498;
  atomic_mass_map[460900000] = 89.95737;
  atomic_mass_map[330910000] = 90.96039;
  atomic_mass_map[340910000] = 90.9457;
  atomic_mass_map[350910000] = 90.934398618;
  atomic_mass_map[360910000] = 90.92380631;
  atomic_mass_map[370910000] = 90.916537265;
  atomic_mass_map[380910000] = 90.910195958;
  atomic_mass_map[390910000] = 90.907298066;
  atomic_mass_map[400910000] = 90.905640223;
  atomic_mass_map[410910000] = 90.906990274;
  atomic_mass_map[420910000] = 90.911745195;
  atomic_mass_map[430910000] = 90.918424975;
  atomic_mass_map[440910000] = 90.926741532;
  atomic_mass_map[450910000] = 90.937123;
  atomic_mass_map[460910000] = 90.950692;
  atomic_mass_map[330920000] = 91.96674;
  atomic_mass_map[340920000] = 91.94984;
  atomic_mass_map[350920000] = 91.939631597;
  atomic_mass_map[360920000] = 91.926173094;
  atomic_mass_map[370920000] = 91.919728481;
  atomic_mass_map[380920000] = 91.911038224;
  atomic_mass_map[390920000] = 91.908945745;
  atomic_mass_map[400920000] = 91.905035322;
  atomic_mass_map[410920000] = 91.907188568;
  atomic_mass_map[420920000] = 91.906807155;
  atomic_mass_map[430920000] = 91.915269779;
  atomic_mass_map[440920000] = 91.920234375;
  atomic_mass_map[450920000] = 91.932367694;
  atomic_mass_map[460920000] = 91.941406;
  atomic_mass_map[470920000] = 91.960139;
  atomic_mass_map[340930000] = 92.95629;
  atomic_mass_map[350930000] = 92.94322;
  atomic_mass_map[360930000] = 92.931147174;
  atomic_mass_map[370930000] = 92.922039325;
  atomic_mass_map[380930000] = 92.914024311;
  atomic_mass_map[390930000] = 92.909578422;
  atomic_mass_map[400930000] = 92.906470646;
  atomic_mass_map[410930000] = 92.906373161;
  atomic_mass_map[420930000] = 92.906808773;
  atomic_mass_map[430930000] = 92.910245149;
  atomic_mass_map[440930000] = 92.917104444;
  atomic_mass_map[450930000] = 92.925912781;
  atomic_mass_map[460930000] = 92.93666;
  atomic_mass_map[470930000] = 92.95033;
  atomic_mass_map[340940000] = 93.96049;
  atomic_mass_map[350940000] = 93.949114;
  atomic_mass_map[360940000] = 93.934140454;
  atomic_mass_map[370940000] = 93.926394818;
  atomic_mass_map[380940000] = 93.915355643;
  atomic_mass_map[390940000] = 93.911592063;
  atomic_mass_map[400940000] = 93.906312524;
  atomic_mass_map[410940000] = 93.907278992;
  atomic_mass_map[420940000] = 93.905083592;
  atomic_mass_map[430940000] = 93.909652325;
  atomic_mass_map[440940000] = 93.911342863;
  atomic_mass_map[450940000] = 93.921730453;
  atomic_mass_map[460940000] = 93.929036292;
  atomic_mass_map[470940000] = 93.943736;
  atomic_mass_map[480940000] = 93.956908;
  atomic_mass_map[340950000] = 94.9673;
  atomic_mass_map[350950000] = 94.95301;
  atomic_mass_map[360950000] = 94.939710923;
  atomic_mass_map[370950000] = 94.929262568;
  atomic_mass_map[380950000] = 94.91935584;
  atomic_mass_map[390950000] = 94.912818711;
  atomic_mass_map[400950000] = 94.908040267;
  atomic_mass_map[410950000] = 94.906831115;
  atomic_mass_map[420950000] = 94.905837442;
  atomic_mass_map[430950000] = 94.907652287;
  atomic_mass_map[440950000] = 94.91040442;
  atomic_mass_map[450950000] = 94.915897895;
  atomic_mass_map[460950000] = 94.924888512;
  atomic_mass_map[470950000] = 94.93602;
  atomic_mass_map[480950000] = 94.94994;
  atomic_mass_map[350960000] = 95.95903;
  atomic_mass_map[360960000] = 95.943016618;
  atomic_mass_map[370960000] = 95.934133393;
  atomic_mass_map[380960000] = 95.921712692;
  atomic_mass_map[390960000] = 95.915902953;
  atomic_mass_map[400960000] = 95.908277621;
  atomic_mass_map[410960000] = 95.908101591;
  atomic_mass_map[420960000] = 95.904674774;
  atomic_mass_map[430960000] = 95.907866681;
  atomic_mass_map[440960000] = 95.907588914;
  atomic_mass_map[450960000] = 95.91445171;
  atomic_mass_map[460960000] = 95.918213744;
  atomic_mass_map[470960000] = 95.930743906;
  atomic_mass_map[480960000] = 95.94034;
  atomic_mass_map[490960000] = 95.959323;
  atomic_mass_map[350970000] = 96.96344;
  atomic_mass_map[360970000] = 96.949088784;
  atomic_mass_map[370970000] = 96.937177118;
  atomic_mass_map[380970000] = 96.926374776;
  atomic_mass_map[390970000] = 96.918280286;
  atomic_mass_map[400970000] = 96.910957386;
  atomic_mass_map[410970000] = 96.908098414;
  atomic_mass_map[420970000] = 96.906016903;
  atomic_mass_map[430970000] = 96.906360723;
  atomic_mass_map[440970000] = 96.907545779;
  atomic_mass_map[450970000] = 96.911327876;
  atomic_mass_map[460970000] = 96.916471987;
  atomic_mass_map[470970000] = 96.923965326;
  atomic_mass_map[480970000] = 96.9351;
  atomic_mass_map[490970000] = 96.94934;
  atomic_mass_map[350980000] = 97.969672;
  atomic_mass_map[360980000] = 97.95243;
  atomic_mass_map[370980000] = 97.941632317;
  atomic_mass_map[380980000] = 97.92869186;
  atomic_mass_map[390980000] = 97.92238836;
  atomic_mass_map[400980000] = 97.912735124;
  atomic_mass_map[410980000] = 97.91033265;
  atomic_mass_map[420980000] = 97.905403608;
  atomic_mass_map[430980000] = 97.907211205;
  atomic_mass_map[440980000] = 97.905286713;
  atomic_mass_map[450980000] = 97.91070774;
  atomic_mass_map[460980000] = 97.912698337;
  atomic_mass_map[470980000] = 97.921559972;
  atomic_mass_map[480980000] = 97.927389317;
  atomic_mass_map[490980000] = 97.94214;
  atomic_mass_map[360990000] = 98.95839;
  atomic_mass_map[370990000] = 98.94511919199999;
  atomic_mass_map[380990000] = 98.932880511;
  atomic_mass_map[390990000] = 98.924154288;
  atomic_mass_map[400990000] = 98.916670835;
  atomic_mass_map[410990000] = 98.911609371;
  atomic_mass_map[420990000] = 98.907707298;
  atomic_mass_map[430990000] = 98.906249678;
  atomic_mass_map[440990000] = 98.905930278;
  atomic_mass_map[450990000] = 98.90812469;
  atomic_mass_map[460990000] = 98.91177329;
  atomic_mass_map[470990000] = 98.917645768;
  atomic_mass_map[480990000] = 98.924925847;
  atomic_mass_map[490990000] = 98.93411;
  atomic_mass_map[500990000] = 98.94853;
  atomic_mass_map[361000000] = 99.96236999999999;
  atomic_mass_map[371000000] = 99.950351731;
  atomic_mass_map[381000000] = 99.935779615;
  atomic_mass_map[391000000] = 99.927721063;
  atomic_mass_map[401000000] = 99.918005444;
  atomic_mass_map[411000000] = 99.914333963;
  atomic_mass_map[421000000] = 99.907467976;
  atomic_mass_map[431000000] = 99.907652711;
  atomic_mass_map[441000000] = 99.904210452;
  atomic_mass_map[451000000] = 99.908114141;
  atomic_mass_map[461000000] = 99.908520315;
  atomic_mass_map[471000000] = 99.916115445;
  atomic_mass_map[481000000] = 99.92034882;
  atomic_mass_map[491000000] = 99.93095718000001;
  atomic_mass_map[501000000] = 99.938504196;
  atomic_mass_map[361010000] = 100.96873;
  atomic_mass_map[371010000] = 100.954004;
  atomic_mass_map[381010000] = 100.940606266;
  atomic_mass_map[391010000] = 100.930154138;
  atomic_mass_map[401010000] = 100.92145311;
  atomic_mass_map[411010000] = 100.915306496;
  atomic_mass_map[421010000] = 100.910337641;
  atomic_mass_map[431010000] = 100.90730526;
  atomic_mass_map[441010000] = 100.905573075;
  atomic_mass_map[451010000] = 100.906158905;
  atomic_mass_map[461010000] = 100.908284828;
  atomic_mass_map[471010000] = 100.912683953;
  atomic_mass_map[481010000] = 100.918586211;
  atomic_mass_map[491010000] = 100.92634;
  atomic_mass_map[501010000] = 100.935259244;
  atomic_mass_map[371020000] = 101.95952;
  atomic_mass_map[381020000] = 101.94400468;
  atomic_mass_map[391020000] = 101.934327889;
  atomic_mass_map[401020000] = 101.923147431;
  atomic_mass_map[411020000] = 101.918083697;
  atomic_mass_map[421020000] = 101.910288138;
  atomic_mass_map[431020000] = 101.909207275;
  atomic_mass_map[441020000] = 101.9043403;
  atomic_mass_map[451020000] = 101.90683427;
  atomic_mass_map[461020000] = 101.905632058;
  atomic_mass_map[471020000] = 101.91170454;
  atomic_mass_map[481020000] = 101.914481799;
  atomic_mass_map[491020000] = 101.924105916;
  atomic_mass_map[501020000] = 101.93028953;
  atomic_mass_map[371030000] = 102.96392;
  atomic_mass_map[381030000] = 102.94909;
  atomic_mass_map[391030000] = 102.937243208;
  atomic_mass_map[401030000] = 102.92719724;
  atomic_mass_map[411030000] = 102.919453403;
  atomic_mass_map[421030000] = 102.91308514;
  atomic_mass_map[431030000] = 102.909174008;
  atomic_mass_map[441030000] = 102.906314833;
  atomic_mass_map[451030000] = 102.905494068;
  atomic_mass_map[461030000] = 102.90611084;
  atomic_mass_map[471030000] = 102.90896056;
  atomic_mass_map[481030000] = 102.913416923;
  atomic_mass_map[491030000] = 102.919878613;
  atomic_mass_map[501030000] = 102.928101962;
  atomic_mass_map[511030000] = 102.93969;
  atomic_mass_map[381040000] = 103.95265;
  atomic_mass_map[391040000] = 103.94196;
  atomic_mass_map[401040000] = 103.929442315;
  atomic_mass_map[411040000] = 103.922899115;
  atomic_mass_map[421040000] = 103.913740756;
  atomic_mass_map[431040000] = 103.911428905;
  atomic_mass_map[441040000] = 103.90542536;
  atomic_mass_map[451040000] = 103.906645295;
  atomic_mass_map[461040000] = 103.904030401;
  atomic_mass_map[471040000] = 103.908623725;
  atomic_mass_map[481040000] = 103.90985623;
  atomic_mass_map[491040000] = 103.91821454;
  atomic_mass_map[501040000] = 103.923105197;
  atomic_mass_map[511040000] = 103.936474502;
  atomic_mass_map[381050000] = 104.95855;
  atomic_mass_map[391050000] = 104.944959;
  atomic_mass_map[401050000] = 104.93401489;
  atomic_mass_map[411050000] = 104.924942564;
  atomic_mass_map[421050000] = 104.916975159;
  atomic_mass_map[431050000] = 104.911657952;
  atomic_mass_map[441050000] = 104.907745525;
  atomic_mass_map[451050000] = 104.905687806;
  atomic_mass_map[461050000] = 104.905079487;
  atomic_mass_map[471050000] = 104.906525607;
  atomic_mass_map[481050000] = 104.909463895;
  atomic_mass_map[491050000] = 104.914502324;
  atomic_mass_map[501050000] = 104.921268423;
  atomic_mass_map[511050000] = 104.931276549;
  atomic_mass_map[521050000] = 104.943304508;
  atomic_mass_map[381060000] = 105.962651;
  atomic_mass_map[391060000] = 105.95056;
  atomic_mass_map[401060000] = 105.937144;
  atomic_mass_map[411060000] = 105.928927768;
  atomic_mass_map[421060000] = 105.918266218;
  atomic_mass_map[431060000] = 105.914356697;
  atomic_mass_map[441060000] = 105.907328203;
  atomic_mass_map[451060000] = 105.907285901;
  atomic_mass_map[461060000] = 105.903480293;
  atomic_mass_map[471060000] = 105.906663507;
  atomic_mass_map[481060000] = 105.906459797;
  atomic_mass_map[491060000] = 105.913463603;
  atomic_mass_map[501060000] = 105.916957396;
  atomic_mass_map[511060000] = 105.928637982;
  atomic_mass_map[521060000] = 105.937498526;
  atomic_mass_map[381070000] = 106.968975;
  atomic_mass_map[391070000] = 106.95452;
  atomic_mass_map[401070000] = 106.941621;
  atomic_mass_map[411070000] = 106.931589672;
  atomic_mass_map[421070000] = 106.922112692;
  atomic_mass_map[431070000] = 106.915458485;
  atomic_mass_map[441070000] = 106.909969885;
  atomic_mass_map[451070000] = 106.906747974;
  atomic_mass_map[461070000] = 106.905128064;
  atomic_mass_map[471070000] = 106.905091531;
  atomic_mass_map[481070000] = 106.906612108;
  atomic_mass_map[491070000] = 106.910290071;
  atomic_mass_map[501070000] = 106.915713651;
  atomic_mass_map[511070000] = 106.924150624;
  atomic_mass_map[521070000] = 106.935008356;
  atomic_mass_map[531070000] = 106.946935;
  atomic_mass_map[391080000] = 107.95996;
  atomic_mass_map[401080000] = 107.94487;
  atomic_mass_map[411080000] = 107.936074988;
  atomic_mass_map[421080000] = 107.924040367;
  atomic_mass_map[431080000] = 107.918493541;
  atomic_mass_map[441080000] = 107.910185841;
  atomic_mass_map[451080000] = 107.908714688;
  atomic_mass_map[461080000] = 107.903891805;
  atomic_mass_map[471080000] = 107.905950266;
  atomic_mass_map[481080000] = 107.904183587;
  atomic_mass_map[491080000] = 107.909693655;
  atomic_mass_map[501080000] = 107.911894292;
  atomic_mass_map[511080000] = 107.922226734;
  atomic_mass_map[521080000] = 107.929380471;
  atomic_mass_map[531080000] = 107.943478321;
  atomic_mass_map[391090000] = 108.964358;
  atomic_mass_map[401090000] = 108.95041;
  atomic_mass_map[411090000] = 108.939141;
  atomic_mass_map[421090000] = 108.928431106;
  atomic_mass_map[431090000] = 108.920254156;
  atomic_mass_map[441090000] = 108.913323756;
  atomic_mass_map[451090000] = 108.908749326;
  atomic_mass_map[461090000] = 108.905950574;
  atomic_mass_map[471090000] = 108.904755773;
  atomic_mass_map[481090000] = 108.904986698;
  atomic_mass_map[491090000] = 108.907149685;
  atomic_mass_map[501090000] = 108.911292843;
  atomic_mass_map[511090000] = 108.918141204;
  atomic_mass_map[521090000] = 108.927304534;
  atomic_mass_map[531090000] = 108.938086025;
  atomic_mass_map[541090000] = 108.950434948;
  atomic_mass_map[401100000] = 109.95396;
  atomic_mass_map[411100000] = 109.943843;
  atomic_mass_map[421100000] = 109.93071068;
  atomic_mass_map[431100000] = 109.923741312;
  atomic_mass_map[441100000] = 109.914038548;
  atomic_mass_map[451100000] = 109.911079742;
  atomic_mass_map[461100000] = 109.905172868;
  atomic_mass_map[471100000] = 109.906110719;
  atomic_mass_map[481100000] = 109.90300746;
  atomic_mass_map[491100000] = 109.907170665;
  atomic_mass_map[501100000] = 109.907844835;
  atomic_mass_map[511100000] = 109.916854286;
  atomic_mass_map[521100000] = 109.922458104;
  atomic_mass_map[531100000] = 109.935089033;
  atomic_mass_map[541100000] = 109.944258765;
  atomic_mass_map[401110000] = 110.959678;
  atomic_mass_map[411110000] = 110.94753;
  atomic_mass_map[421110000] = 110.935652016;
  atomic_mass_map[431110000] = 110.925899016;
  atomic_mass_map[441110000] = 110.917567616;
  atomic_mass_map[451110000] = 110.911642531;
  atomic_mass_map[461110000] = 110.907690347;
  atomic_mass_map[471110000] = 110.905296816;
  atomic_mass_map[481110000] = 110.904183766;
  atomic_mass_map[491110000] = 110.905107233;
  atomic_mass_map[501110000] = 110.907741126;
  atomic_mass_map[511110000] = 110.913218189;
  atomic_mass_map[521110000] = 110.921000589;
  atomic_mass_map[531110000] = 110.930269239;
  atomic_mass_map[541110000] = 110.941603989;
  atomic_mass_map[551110000] = 110.95403;
  atomic_mass_map[401120000] = 111.963703;
  atomic_mass_map[411120000] = 111.95247;
  atomic_mass_map[421120000] = 111.93831;
  atomic_mass_map[431120000] = 111.929941644;
  atomic_mass_map[441120000] = 111.918806972;
  atomic_mass_map[451120000] = 111.914404705;
  atomic_mass_map[461120000] = 111.907329986;
  atomic_mass_map[471120000] = 111.90704855;
  atomic_mass_map[481120000] = 111.902763883;
  atomic_mass_map[491120000] = 111.905538704;
  atomic_mass_map[501120000] = 111.904824877;
  atomic_mass_map[511120000] = 111.912399903;
  atomic_mass_map[521120000] = 111.91672785;
  atomic_mass_map[531120000] = 111.92800455;
  atomic_mass_map[541120000] = 111.935559071;
  atomic_mass_map[551120000] = 111.950305341;
  atomic_mass_map[411130000] = 112.95651;
  atomic_mass_map[421130000] = 112.94365;
  atomic_mass_map[431130000] = 112.932569033;
  atomic_mass_map[441130000] = 112.922846396;
  atomic_mass_map[451130000] = 112.91543956699999;
  atomic_mass_map[461130000] = 112.910261267;
  atomic_mass_map[471130000] = 112.906572858;
  atomic_mass_map[481130000] = 112.904408097;
  atomic_mass_map[491130000] = 112.904060448;
  atomic_mass_map[501130000] = 112.905175845;
  atomic_mass_map[511130000] = 112.909374652;
  atomic_mass_map[521130000] = 112.915891;
  atomic_mass_map[531130000] = 112.923650064;
  atomic_mass_map[541130000] = 112.933221666;
  atomic_mass_map[551130000] = 112.944428488;
  atomic_mass_map[561130000] = 112.95729;
  atomic_mass_map[411140000] = 113.96201;
  atomic_mass_map[421140000] = 113.94653;
  atomic_mass_map[431140000] = 113.93709;
  atomic_mass_map[441140000] = 113.92461378;
  atomic_mass_map[451140000] = 113.918721296;
  atomic_mass_map[461140000] = 113.91036878;
  atomic_mass_map[471140000] = 113.908823031;
  atomic_mass_map[481140000] = 113.90336499;
  atomic_mass_map[491140000] = 113.904916402;
  atomic_mass_map[501140000] = 113.902780132;
  atomic_mass_map[511140000] = 113.909289191;
  atomic_mass_map[521140000] = 113.912089;
  atomic_mass_map[531140000] = 113.92185;
  atomic_mass_map[541140000] = 113.927980331;
  atomic_mass_map[551140000] = 113.941296175;
  atomic_mass_map[561140000] = 113.950718495;
  atomic_mass_map[411150000] = 114.96634;
  atomic_mass_map[421150000] = 114.95196;
  atomic_mass_map[431150000] = 114.939538;
  atomic_mass_map[441150000] = 114.928942393;
  atomic_mass_map[451150000] = 114.920310993;
  atomic_mass_map[461150000] = 114.913658718;
  atomic_mass_map[471150000] = 114.908767363;
  atomic_mass_map[481150000] = 114.905437417;
  atomic_mass_map[491150000] = 114.903878773;
  atomic_mass_map[501150000] = 114.903344697;
  atomic_mass_map[511150000] = 114.906598;
  atomic_mass_map[521150000] = 114.911902;
  atomic_mass_map[531150000] = 114.918048;
  atomic_mass_map[541150000] = 114.926293945;
  atomic_mass_map[551150000] = 114.93591;
  atomic_mass_map[561150000] = 114.947375;
  atomic_mass_map[421160000] = 115.955448;
  atomic_mass_map[431160000] = 115.94476;
  atomic_mass_map[441160000] = 115.931219193;
  atomic_mass_map[451160000] = 115.924061645;
  atomic_mass_map[461160000] = 115.91429721;
  atomic_mass_map[471160000] = 115.911386812;
  atomic_mass_map[481160000] = 115.90476323;
  atomic_mass_map[491160000] = 115.905259992;
  atomic_mass_map[501160000] = 115.901742824;
  atomic_mass_map[511160000] = 115.906792583;
  atomic_mass_map[521160000] = 115.90846;
  atomic_mass_map[531160000] = 115.916808658;
  atomic_mass_map[541160000] = 115.921581112;
  atomic_mass_map[551160000] = 115.933395;
  atomic_mass_map[561160000] = 115.941406;
  atomic_mass_map[571160000] = 115.956365;
  atomic_mass_map[421170000] = 116.96117;
  atomic_mass_map[431170000] = 116.94806;
  atomic_mass_map[441170000] = 116.936135;
  atomic_mass_map[451170000] = 116.926035623;
  atomic_mass_map[461170000] = 116.917954944;
  atomic_mass_map[471170000] = 116.911773974;
  atomic_mass_map[481170000] = 116.907226038;
  atomic_mass_map[491170000] = 116.904515712;
  atomic_mass_map[501170000] = 116.902954017;
  atomic_mass_map[511170000] = 116.904841535;
  atomic_mass_map[521170000] = 116.908646313;
  atomic_mass_map[531170000] = 116.913648314;
  atomic_mass_map[541170000] = 116.92035876;
  atomic_mass_map[551170000] = 116.928616726;
  atomic_mass_map[561170000] = 116.938316561;
  atomic_mass_map[571170000] = 116.950111;
  atomic_mass_map[421180000] = 117.96497;
  atomic_mass_map[431180000] = 117.95299;
  atomic_mass_map[441180000] = 117.938529;
  atomic_mass_map[451180000] = 117.930340443;
  atomic_mass_map[461180000] = 117.919066847;
  atomic_mass_map[471180000] = 117.914595487;
  atomic_mass_map[481180000] = 117.906921955;
  atomic_mass_map[491180000] = 117.906356659;
  atomic_mass_map[501180000] = 117.901606609;
  atomic_mass_map[511180000] = 117.905532174;
  atomic_mass_map[521180000] = 117.905853839;
  atomic_mass_map[531180000] = 117.913074;
  atomic_mass_map[541180000] = 117.91617868;
  atomic_mass_map[551180000] = 117.926559519;
  atomic_mass_map[561180000] = 117.93306;
  atomic_mass_map[571180000] = 117.946795;
  atomic_mass_map[431190000] = 118.95666;
  atomic_mass_map[441190000] = 118.94357;
  atomic_mass_map[451190000] = 118.932556952;
  atomic_mass_map[461190000] = 118.923340459;
  atomic_mass_map[471190000] = 118.915570293;
  atomic_mass_map[481190000] = 118.909846903;
  atomic_mass_map[491190000] = 118.905850944;
  atomic_mass_map[501190000] = 118.903311216;
  atomic_mass_map[511190000] = 118.903945512;
  atomic_mass_map[521190000] = 118.906407148;
  atomic_mass_map[531190000] = 118.910074;
  atomic_mass_map[541190000] = 118.915410713;
  atomic_mass_map[551190000] = 118.92237733;
  atomic_mass_map[561190000] = 118.930659686;
  atomic_mass_map[571190000] = 118.941181;
  atomic_mass_map[581190000] = 118.952828;
  atomic_mass_map[431200000] = 119.96187;
  atomic_mass_map[441200000] = 119.94631;
  atomic_mass_map[451200000] = 119.93686;
  atomic_mass_map[461200000] = 119.924551258;
  atomic_mass_map[471200000] = 119.918784767;
  atomic_mass_map[481200000] = 119.909868067;
  atomic_mass_map[491200000] = 119.907966805;
  atomic_mass_map[501200000] = 119.902201873;
  atomic_mass_map[511200000] = 119.905079624;
  atomic_mass_map[521200000] = 119.904059514;
  atomic_mass_map[531200000] = 119.910087465;
  atomic_mass_map[541200000] = 119.91178427;
  atomic_mass_map[551200000] = 119.920677279;
  atomic_mass_map[561200000] = 119.926045;
  atomic_mass_map[571200000] = 119.938196;
  atomic_mass_map[581200000] = 119.946752;
  atomic_mass_map[431210000] = 120.965883;
  atomic_mass_map[441210000] = 120.95164;
  atomic_mass_map[451210000] = 120.939613;
  atomic_mass_map[461210000] = 120.928950343;
  atomic_mass_map[471210000] = 120.920125282;
  atomic_mass_map[481210000] = 120.912963663;
  atomic_mass_map[491210000] = 120.907851286;
  atomic_mass_map[501210000] = 120.904242792;
  atomic_mass_map[511210000] = 120.903810093;
  atomic_mass_map[521210000] = 120.904942488;
  atomic_mass_map[531210000] = 120.907405255;
  atomic_mass_map[541210000] = 120.911453014;
  atomic_mass_map[551210000] = 120.917227238;
  atomic_mass_map[561210000] = 120.924052289;
  atomic_mass_map[571210000] = 120.933236;
  atomic_mass_map[581210000] = 120.943435;
  atomic_mass_map[591210000] = 120.955532;
  atomic_mass_map[441220000] = 121.95475;
  atomic_mass_map[451220000] = 121.94409;
  atomic_mass_map[461220000] = 121.930631694;
  atomic_mass_map[471220000] = 121.923664448;
  atomic_mass_map[481220000] = 121.913459052;
  atomic_mass_map[491220000] = 121.910280966;
  atomic_mass_map[501220000] = 121.903444001;
  atomic_mass_map[511220000] = 121.905168074;
  atomic_mass_map[521220000] = 121.903043434;
  atomic_mass_map[531220000] = 121.90758882;
  atomic_mass_map[541220000] = 121.908367658;
  atomic_mass_map[551220000] = 121.916108145;
  atomic_mass_map[561220000] = 121.919904;
  atomic_mass_map[571220000] = 121.93071;
  atomic_mass_map[581220000] = 121.93787;
  atomic_mass_map[591220000] = 121.951927;
  atomic_mass_map[441230000] = 122.960193;
  atomic_mass_map[451230000] = 122.94701;
  atomic_mass_map[461230000] = 122.935126;
  atomic_mass_map[471230000] = 122.925337062;
  atomic_mass_map[481230000] = 122.916892453;
  atomic_mass_map[491230000] = 122.910433826;
  atomic_mass_map[501230000] = 122.905725446;
  atomic_mass_map[511230000] = 122.904214016;
  atomic_mass_map[521230000] = 122.904269747;
  atomic_mass_map[531230000] = 122.90558852;
  atomic_mass_map[541230000] = 122.90848174999999;
  atomic_mass_map[551230000] = 122.912996062;
  atomic_mass_map[561230000] = 122.918781062;
  atomic_mass_map[571230000] = 122.9263;
  atomic_mass_map[581230000] = 122.93528;
  atomic_mass_map[591230000] = 122.946076;
  atomic_mass_map[441240000] = 123.963542;
  atomic_mass_map[451240000] = 123.951809;
  atomic_mass_map[461240000] = 123.937316;
  atomic_mass_map[471240000] = 123.928931229;
  atomic_mass_map[481240000] = 123.917657363;
  atomic_mass_map[491240000] = 123.913182263;
  atomic_mass_map[501240000] = 123.905276692;
  atomic_mass_map[511240000] = 123.905935789;
  atomic_mass_map[521240000] = 123.902817064;
  atomic_mass_map[531240000] = 123.906209021;
  atomic_mass_map[541240000] = 123.905891588;
  atomic_mass_map[551240000] = 123.912257798;
  atomic_mass_map[561240000] = 123.915093629;
  atomic_mass_map[571240000] = 123.924574275;
  atomic_mass_map[581240000] = 123.93031;
  atomic_mass_map[591240000] = 123.94294;
  atomic_mass_map[601240000] = 123.9522;
  atomic_mass_map[451250000] = 124.954911;
  atomic_mass_map[461250000] = 124.9419;
  atomic_mass_map[471250000] = 124.930735;
  atomic_mass_map[481250000] = 124.921257577;
  atomic_mass_map[491250000] = 124.913604591;
  atomic_mass_map[501250000] = 124.907786442;
  atomic_mass_map[511250000] = 124.905252987;
  atomic_mass_map[521250000] = 124.9044299;
  atomic_mass_map[531250000] = 124.904629333;
  atomic_mass_map[541250000] = 124.90639405;
  atomic_mass_map[551250000] = 124.909727867;
  atomic_mass_map[561250000] = 124.914471843;
  atomic_mass_map[571250000] = 124.920815932;
  atomic_mass_map[581250000] = 124.92844;
  atomic_mass_map[591250000] = 124.937799;
  atomic_mass_map[601250000] = 124.9489;
  atomic_mass_map[451260000] = 125.959957;
  atomic_mass_map[461260000] = 125.944326;
  atomic_mass_map[471260000] = 125.934857;
  atomic_mass_map[481260000] = 125.922429127;
  atomic_mass_map[491260000] = 125.916507344;
  atomic_mass_map[501260000] = 125.907658836;
  atomic_mass_map[511260000] = 125.907253036;
  atomic_mass_map[521260000] = 125.903310866;
  atomic_mass_map[531260000] = 125.905623313;
  atomic_mass_map[541260000] = 125.904296794;
  atomic_mass_map[551260000] = 125.909445655;
  atomic_mass_map[561260000] = 125.911250204;
  atomic_mass_map[571260000] = 125.919512667;
  atomic_mass_map[581260000] = 125.923971;
  atomic_mass_map[591260000] = 125.93524;
  atomic_mass_map[601260000] = 125.94311;
  atomic_mass_map[611260000] = 125.957756;
  atomic_mass_map[451270000] = 126.963467;
  atomic_mass_map[461270000] = 126.94935;
  atomic_mass_map[471270000] = 126.937262;
  atomic_mass_map[481270000] = 126.926196624;
  atomic_mass_map[491270000] = 126.917448546;
  atomic_mass_map[501270000] = 126.910390401;
  atomic_mass_map[511270000] = 126.906924277;
  atomic_mass_map[521270000] = 126.905225714;
  atomic_mass_map[531270000] = 126.904471838;
  atomic_mass_map[541270000] = 126.905182899;
  atomic_mass_map[551270000] = 126.907417381;
  atomic_mass_map[561270000] = 126.911091275;
  atomic_mass_map[571270000] = 126.916375084;
  atomic_mass_map[581270000] = 126.922727;
  atomic_mass_map[591270000] = 126.93071;
  atomic_mass_map[601270000] = 126.94038;
  atomic_mass_map[611270000] = 126.95192;
  atomic_mass_map[461280000] = 127.952238;
  atomic_mass_map[471280000] = 127.941363;
  atomic_mass_map[481280000] = 127.927812857;
  atomic_mass_map[491280000] = 127.920401053;
  atomic_mass_map[501280000] = 127.910507197;
  atomic_mass_map[511280000] = 127.909145645;
  atomic_mass_map[521280000] = 127.904461311;
  atomic_mass_map[531280000] = 127.9058086;
  atomic_mass_map[541280000] = 127.903530996;
  atomic_mass_map[551280000] = 127.907748648;
  atomic_mass_map[561280000] = 127.908342408;
  atomic_mass_map[571280000] = 127.915592123;
  atomic_mass_map[581280000] = 127.918911;
  atomic_mass_map[591280000] = 127.928791;
  atomic_mass_map[601280000] = 127.93525;
  atomic_mass_map[611280000] = 127.9487;
  atomic_mass_map[621280000] = 127.958486;
  atomic_mass_map[461290000] = 128.959624;
  atomic_mass_map[471290000] = 128.944197;
  atomic_mass_map[481290000] = 128.932304399;
  atomic_mass_map[491290000] = 128.921805486;
  atomic_mass_map[501290000] = 128.913482102;
  atomic_mass_map[511290000] = 128.909146696;
  atomic_mass_map[521290000] = 128.906596492;
  atomic_mass_map[531290000] = 128.904983687;
  atomic_mass_map[541290000] = 128.90478085892;
  atomic_mass_map[551290000] = 128.90606569;
  atomic_mass_map[561290000] = 128.908680896;
  atomic_mass_map[571290000] = 128.912694475;
  atomic_mass_map[581290000] = 128.918102;
  atomic_mass_map[591290000] = 128.925095;
  atomic_mass_map[601290000] = 128.933102;
  atomic_mass_map[611290000] = 128.94323;
  atomic_mass_map[621290000] = 128.954911;
  atomic_mass_map[471300000] = 129.950942;
  atomic_mass_map[481300000] = 129.934387566;
  atomic_mass_map[491300000] = 129.924977288;
  atomic_mass_map[501300000] = 129.913974533;
  atomic_mass_map[511300000] = 129.911662688;
  atomic_mass_map[521300000] = 129.906222747;
  atomic_mass_map[531300000] = 129.906670211;
  atomic_mass_map[541300000] = 129.903509349;
  atomic_mass_map[551300000] = 129.906709283;
  atomic_mass_map[561300000] = 129.906320874;
  atomic_mass_map[571300000] = 129.912369413;
  atomic_mass_map[581300000] = 129.914736;
  atomic_mass_map[591300000] = 129.92359;
  atomic_mass_map[601300000] = 129.928506;
  atomic_mass_map[611300000] = 129.94053;
  atomic_mass_map[621300000] = 129.949;
  atomic_mass_map[631300000] = 129.96384;
  atomic_mass_map[471310000] = 130.95665;
  atomic_mass_map[481310000] = 130.94072;
  atomic_mass_map[491310000] = 130.926972122;
  atomic_mass_map[501310000] = 130.917053066;
  atomic_mass_map[511310000] = 130.911989341;
  atomic_mass_map[521310000] = 130.908522211;
  atomic_mass_map[531310000] = 130.906126384;
  atomic_mass_map[541310000] = 130.905084136;
  atomic_mass_map[551310000] = 130.905464999;
  atomic_mass_map[561310000] = 130.906941181;
  atomic_mass_map[571310000] = 130.91007;
  atomic_mass_map[581310000] = 130.914429465;
  atomic_mass_map[591310000] = 130.92023496;
  atomic_mass_map[601310000] = 130.92724802;
  atomic_mass_map[611310000] = 130.935952;
  atomic_mass_map[621310000] = 130.94618;
  atomic_mass_map[631310000] = 130.957842;
  atomic_mass_map[471320000] = 131.963725;
  atomic_mass_map[481320000] = 131.94604;
  atomic_mass_map[491320000] = 131.932998449;
  atomic_mass_map[501320000] = 131.917823902;
  atomic_mass_map[511320000] = 131.914508015;
  atomic_mass_map[521320000] = 131.908546716;
  atomic_mass_map[531320000] = 131.907993514;
  atomic_mass_map[541320000] = 131.90415508697;
  atomic_mass_map[551320000] = 131.906437743;
  atomic_mass_map[561320000] = 131.905061098;
  atomic_mass_map[571320000] = 131.910118959;
  atomic_mass_map[581320000] = 131.911463846;
  atomic_mass_map[591320000] = 131.91924;
  atomic_mass_map[601320000] = 131.923321237;
  atomic_mass_map[611320000] = 131.93384;
  atomic_mass_map[621320000] = 131.94087;
  atomic_mass_map[631320000] = 131.954696;
  atomic_mass_map[481330000] = 132.95285;
  atomic_mass_map[491330000] = 132.93831;
  atomic_mass_map[501330000] = 132.923913756;
  atomic_mass_map[511330000] = 132.91527213;
  atomic_mass_map[521330000] = 132.910963332;
  atomic_mass_map[531330000] = 132.907827361;
  atomic_mass_map[541330000] = 132.90591075;
  atomic_mass_map[551330000] = 132.905451961;
  atomic_mass_map[561330000] = 132.906007325;
  atomic_mass_map[571330000] = 132.908218;
  atomic_mass_map[581330000] = 132.911520402;
  atomic_mass_map[591330000] = 132.916330561;
  atomic_mass_map[601330000] = 132.922348;
  atomic_mass_map[611330000] = 132.929782;
  atomic_mass_map[621330000] = 132.93856;
  atomic_mass_map[631330000] = 132.94929;
  atomic_mass_map[641330000] = 132.961503;
  atomic_mass_map[481340000] = 133.958218;
  atomic_mass_map[491340000] = 133.94454;
  atomic_mass_map[501340000] = 133.928680433;
  atomic_mass_map[511340000] = 133.920535675;
  atomic_mass_map[521340000] = 133.911396379;
  atomic_mass_map[531340000] = 133.909775663;
  atomic_mass_map[541340000] = 133.905393033;
  atomic_mass_map[551340000] = 133.906718503;
  atomic_mass_map[561340000] = 133.904508399;
  atomic_mass_map[571340000] = 133.908514011;
  atomic_mass_map[581340000] = 133.908928142;
  atomic_mass_map[591340000] = 133.915696729;
  atomic_mass_map[601340000] = 133.91879021;
  atomic_mass_map[611340000] = 133.928353;
  atomic_mass_map[621340000] = 133.93411;
  atomic_mass_map[631340000] = 133.9464;
  atomic_mass_map[641340000] = 133.95566;
  atomic_mass_map[491350000] = 134.95005;
  atomic_mass_map[501350000] = 134.934908605;
  atomic_mass_map[511350000] = 134.925184357;
  atomic_mass_map[521350000] = 134.916554718;
  atomic_mass_map[531350000] = 134.910059382;
  atomic_mass_map[541350000] = 134.907231661;
  atomic_mass_map[551350000] = 134.905977234;
  atomic_mass_map[561350000] = 134.905688606;
  atomic_mass_map[571350000] = 134.906984568;
  atomic_mass_map[581350000] = 134.909160799;
  atomic_mass_map[591350000] = 134.913111774;
  atomic_mass_map[601350000] = 134.91818132;
  atomic_mass_map[611350000] = 134.924796;
  atomic_mass_map[621350000] = 134.93252;
  atomic_mass_map[631350000] = 134.94187;
  atomic_mass_map[641350000] = 134.952345;
  atomic_mass_map[651350000] = 134.96476;
  atomic_mass_map[491360000] = 135.956511;
  atomic_mass_map[501360000] = 135.93999;
  atomic_mass_map[511360000] = 135.930749011;
  atomic_mass_map[521360000] = 135.920101182;
  atomic_mass_map[531360000] = 135.914604695;
  atomic_mass_map[541360000] = 135.907214476;
  atomic_mass_map[551360000] = 135.90731159;
  atomic_mass_map[561360000] = 135.904575959;
  atomic_mass_map[571360000] = 135.907634962;
  atomic_mass_map[581360000] = 135.907129438;
  atomic_mass_map[591360000] = 135.912677532;
  atomic_mass_map[601360000] = 135.914976064;
  atomic_mass_map[611360000] = 135.923595949;
  atomic_mass_map[621360000] = 135.928275555;
  atomic_mass_map[631360000] = 135.93962;
  atomic_mass_map[641360000] = 135.9473;
  atomic_mass_map[651360000] = 135.961213;
  atomic_mass_map[491370000] = 136.962383;
  atomic_mass_map[501370000] = 136.94655;
  atomic_mass_map[511370000] = 136.935522522;
  atomic_mass_map[521370000] = 136.925599357;
  atomic_mass_map[531370000] = 136.91802818;
  atomic_mass_map[541370000] = 136.911557773;
  atomic_mass_map[551370000] = 136.907089464;
  atomic_mass_map[561370000] = 136.905827375;
  atomic_mass_map[571370000] = 136.906450618;
  atomic_mass_map[581370000] = 136.907762596;
  atomic_mass_map[591370000] = 136.910679304;
  atomic_mass_map[601370000] = 136.914562448;
  atomic_mass_map[611370000] = 136.920479522;
  atomic_mass_map[621370000] = 136.926970517;
  atomic_mass_map[631370000] = 136.935430722;
  atomic_mass_map[641370000] = 136.94502;
  atomic_mass_map[651370000] = 136.95602;
  atomic_mass_map[501380000] = 137.95184;
  atomic_mass_map[511380000] = 137.941792;
  atomic_mass_map[521380000] = 137.929472454;
  atomic_mass_map[531380000] = 137.922726394;
  atomic_mass_map[541380000] = 137.914146271;
  atomic_mass_map[551380000] = 137.911017207;
  atomic_mass_map[561380000] = 137.905247229;
  atomic_mass_map[571380000] = 137.907117834;
  atomic_mass_map[581380000] = 137.905988743;
  atomic_mass_map[591380000] = 137.910752059;
  atomic_mass_map[601380000] = 137.911949717;
  atomic_mass_map[611380000] = 137.919548077;
  atomic_mass_map[621380000] = 137.92324399;
  atomic_mass_map[631380000] = 137.933709;
  atomic_mass_map[641380000] = 137.940096;
  atomic_mass_map[651380000] = 137.95312;
  atomic_mass_map[661380000] = 137.9625;
  atomic_mass_map[501390000] = 138.958733;
  atomic_mass_map[511390000] = 138.94655;
  atomic_mass_map[521390000] = 138.935367193;
  atomic_mass_map[531390000] = 138.926493403;
  atomic_mass_map[541390000] = 138.918792203;
  atomic_mass_map[551390000] = 138.913363992;
  atomic_mass_map[561390000] = 138.908841334;
  atomic_mass_map[571390000] = 138.906358804;
  atomic_mass_map[581390000] = 138.906657625;
  atomic_mass_map[591390000] = 138.90894327;
  atomic_mass_map[601390000] = 138.911954407;
  atomic_mass_map[611390000] = 138.916799806;
  atomic_mass_map[621390000] = 138.922296634;
  atomic_mass_map[631390000] = 138.92979231;
  atomic_mass_map[641390000] = 138.93813;
  atomic_mass_map[651390000] = 138.94833;
  atomic_mass_map[661390000] = 138.95959;
  atomic_mass_map[511400000] = 139.95283;
  atomic_mass_map[521400000] = 139.939262917;
  atomic_mass_map[531400000] = 139.931715917;
  atomic_mass_map[541400000] = 139.921645817;
  atomic_mass_map[551400000] = 139.917283305;
  atomic_mass_map[561400000] = 139.910606666;
  atomic_mass_map[571400000] = 139.909483184;
  atomic_mass_map[581400000] = 139.905446424;
  atomic_mass_map[591400000] = 139.909083592;
  atomic_mass_map[601400000] = 139.909544332;
  atomic_mass_map[611400000] = 139.916034122;
  atomic_mass_map[621400000] = 139.918994717;
  atomic_mass_map[631400000] = 139.928087637;
  atomic_mass_map[641400000] = 139.933674;
  atomic_mass_map[651400000] = 139.945805049;
  atomic_mass_map[661400000] = 139.95402;
  atomic_mass_map[671400000] = 139.968589;
  atomic_mass_map[511410000] = 140.958014;
  atomic_mass_map[521410000] = 140.9458;
  atomic_mass_map[531410000] = 140.935666084;
  atomic_mass_map[541410000] = 140.926787184;
  atomic_mass_map[551410000] = 140.920045086;
  atomic_mass_map[561410000] = 140.9144035;
  atomic_mass_map[571410000] = 140.910969222;
  atomic_mass_map[581410000] = 140.908283987;
  atomic_mass_map[591410000] = 140.907658403;
  atomic_mass_map[601410000] = 140.909615488;
  atomic_mass_map[611410000] = 140.913555084;
  atomic_mass_map[621410000] = 140.918481591;
  atomic_mass_map[631410000] = 140.924931745;
  atomic_mass_map[641410000] = 140.932126;
  atomic_mass_map[651410000] = 140.941448;
  atomic_mass_map[661410000] = 140.95128;
  atomic_mass_map[671410000] = 140.963108;
  atomic_mass_map[521420000] = 141.95022;
  atomic_mass_map[531420000] = 141.941202;
  atomic_mass_map[541420000] = 141.929973098;
  atomic_mass_map[551420000] = 141.924299512;
  atomic_mass_map[561420000] = 141.916432888;
  atomic_mass_map[571420000] = 141.914090454;
  atomic_mass_map[581420000] = 141.909249884;
  atomic_mass_map[591420000] = 141.91005044;
  atomic_mass_map[601420000] = 141.907728895;
  atomic_mass_map[611420000] = 141.912890428;
  atomic_mass_map[621420000] = 141.915204532;
  atomic_mass_map[631420000] = 141.923441836;
  atomic_mass_map[641420000] = 141.928116;
  atomic_mass_map[651420000] = 141.939280859;
  atomic_mass_map[661420000] = 141.946194;
  atomic_mass_map[671420000] = 141.96001;
  atomic_mass_map[681420000] = 141.969909;
  atomic_mass_map[521430000] = 142.95676;
  atomic_mass_map[531430000] = 142.945646;
  atomic_mass_map[541430000] = 142.935369553;
  atomic_mass_map[551430000] = 142.927347348;
  atomic_mass_map[561430000] = 142.92062515;
  atomic_mass_map[571430000] = 142.916079422;
  atomic_mass_map[581430000] = 142.91239163;
  atomic_mass_map[591430000] = 142.910822564;
  atomic_mass_map[601430000] = 142.909819887;
  atomic_mass_map[611430000] = 142.910938073;
  atomic_mass_map[621430000] = 142.914634821;
  atomic_mass_map[631430000] = 142.920298681;
  atomic_mass_map[641430000] = 142.926750682;
  atomic_mass_map[651430000] = 142.935137335;
  atomic_mass_map[661430000] = 142.943994335;
  atomic_mass_map[671430000] = 142.95486;
  atomic_mass_map[681430000] = 142.966441;
  atomic_mass_map[531440000] = 143.95139;
  atomic_mass_map[541440000] = 143.938945079;
  atomic_mass_map[551440000] = 143.932075404;
  atomic_mass_map[561440000] = 143.922954821;
  atomic_mass_map[571440000] = 143.919645589;
  atomic_mass_map[581440000] = 143.91365283;
  atomic_mass_map[591440000] = 143.91331075;
  atomic_mass_map[601440000] = 143.910092865;
  atomic_mass_map[611440000] = 143.912596224;
  atomic_mass_map[621440000] = 143.912006373;
  atomic_mass_map[631440000] = 143.918819517;
  atomic_mass_map[641440000] = 143.922963;
  atomic_mass_map[651440000] = 143.933045;
  atomic_mass_map[661440000] = 143.939269514;
  atomic_mass_map[671440000] = 143.952109714;
  atomic_mass_map[681440000] = 143.9607;
  atomic_mass_map[691440000] = 143.976104;
  atomic_mass_map[531450000] = 144.95605;
  atomic_mass_map[541450000] = 144.944719634;
  atomic_mass_map[551450000] = 144.93552893;
  atomic_mass_map[561450000] = 144.9275184;
  atomic_mass_map[571450000] = 144.921808066;
  atomic_mass_map[581450000] = 144.917265144;
  atomic_mass_map[591450000] = 144.914518033;
  atomic_mass_map[601450000] = 144.912579199;
  atomic_mass_map[611450000] = 144.912755773;
  atomic_mass_map[621450000] = 144.913417244;
  atomic_mass_map[631450000] = 144.916272668;
  atomic_mass_map[641450000] = 144.92171037;
  atomic_mass_map[651450000] = 144.928729105;
  atomic_mass_map[661450000] = 144.937473994;
  atomic_mass_map[671450000] = 144.947267394;
  atomic_mass_map[681450000] = 144.957874;
  atomic_mass_map[691450000] = 144.970389;
  atomic_mass_map[541460000] = 145.948518248;
  atomic_mass_map[551460000] = 145.94062187;
  atomic_mass_map[561460000] = 145.930276431;
  atomic_mass_map[571460000] = 145.925871468;
  atomic_mass_map[581460000] = 145.918802065;
  atomic_mass_map[591460000] = 145.917679549;
  atomic_mass_map[601460000] = 145.913122503;
  atomic_mass_map[611460000] = 145.914702286;
  atomic_mass_map[621460000] = 145.913046881;
  atomic_mass_map[631460000] = 145.917210909;
  atomic_mass_map[641460000] = 145.918318548;
  atomic_mass_map[651460000] = 145.927252768;
  atomic_mass_map[661460000] = 145.932844529;
  atomic_mass_map[671460000] = 145.944993506;
  atomic_mass_map[681460000] = 145.952418359;
  atomic_mass_map[691460000] = 145.966661;
  atomic_mass_map[541470000] = 146.954525;
  atomic_mass_map[551470000] = 146.944261515;
  atomic_mass_map[561470000] = 146.9353039;
  atomic_mass_map[571470000] = 146.9284178;
  atomic_mass_map[581470000] = 146.922689903;
  atomic_mass_map[591470000] = 146.919007458;
  atomic_mass_map[601470000] = 146.91610601;
  atomic_mass_map[611470000] = 146.915144638;
  atomic_mass_map[621470000] = 146.914904064;
  atomic_mass_map[631470000] = 146.916752276;
  atomic_mass_map[641470000] = 146.919100987;
  atomic_mass_map[651470000] = 146.92405462;
  atomic_mass_map[661470000] = 146.931082715;
  atomic_mass_map[671470000] = 146.940142295;
  atomic_mass_map[681470000] = 146.949964458;
  atomic_mass_map[691470000] = 146.96137989;
  atomic_mass_map[541480000] = 147.958561;
  atomic_mass_map[551480000] = 147.949639029;
  atomic_mass_map[561480000] = 147.938170578;
  atomic_mass_map[571480000] = 147.9326794;
  atomic_mass_map[581480000] = 147.924424196;
  atomic_mass_map[591480000] = 147.922130015;
  atomic_mass_map[601480000] = 147.916899093;
  atomic_mass_map[611480000] = 147.917481255;
  atomic_mass_map[621480000] = 147.914829012;
  atomic_mass_map[631480000] = 147.918089294;
  atomic_mass_map[641480000] = 147.918121503;
  atomic_mass_map[651480000] = 147.924275323;
  atomic_mass_map[661480000] = 147.927149772;
  atomic_mass_map[671480000] = 147.937743928;
  atomic_mass_map[681480000] = 147.944735029;
  atomic_mass_map[691480000] = 147.958384029;
  atomic_mass_map[701480000] = 147.967439;
  atomic_mass_map[551490000] = 148.953569;
  atomic_mass_map[561490000] = 148.942973;
  atomic_mass_map[571490000] = 148.93535126;
  atomic_mass_map[581490000] = 148.9284269;
  atomic_mass_map[591490000] = 148.9237361;
  atomic_mass_map[601490000] = 148.920154648;
  atomic_mass_map[611490000] = 148.918341658;
  atomic_mass_map[621490000] = 148.917191375;
  atomic_mass_map[631490000] = 148.917937086;
  atomic_mass_map[641490000] = 148.919347831;
  atomic_mass_map[651490000] = 148.923253753;
  atomic_mass_map[661490000] = 148.927325448;
  atomic_mass_map[671490000] = 148.933819672;
  atomic_mass_map[681490000] = 148.942306;
  atomic_mass_map[691490000] = 148.95289;
  atomic_mass_map[701490000] = 148.96436;
  atomic_mass_map[551500000] = 149.959023;
  atomic_mass_map[561500000] = 149.94643;
  atomic_mass_map[571500000] = 149.939742;
  atomic_mass_map[581500000] = 149.930384035;
  atomic_mass_map[591500000] = 149.926676415;
  atomic_mass_map[601500000] = 149.920901525;
  atomic_mass_map[611500000] = 149.920990217;
  atomic_mass_map[621500000] = 149.917282195;
  atomic_mass_map[631500000] = 149.919707229;
  atomic_mass_map[641500000] = 149.918664066;
  atomic_mass_map[651500000] = 149.923664864;
  atomic_mass_map[661500000] = 149.92559308;
  atomic_mass_map[671500000] = 149.933498358;
  atomic_mass_map[681500000] = 149.937915528;
  atomic_mass_map[691500000] = 149.95009;
  atomic_mass_map[701500000] = 149.95852;
  atomic_mass_map[711500000] = 149.973548;
  atomic_mass_map[551510000] = 150.963253;
  atomic_mass_map[561510000] = 150.951755;
  atomic_mass_map[571510000] = 150.942769;
  atomic_mass_map[581510000] = 150.9342722;
  atomic_mass_map[591510000] = 150.928309114;
  atomic_mass_map[601510000] = 150.923839565;
  atomic_mass_map[611510000] = 150.921216817;
  atomic_mass_map[621510000] = 150.919939066;
  atomic_mass_map[631510000] = 150.91985686;
  atomic_mass_map[641510000] = 150.920355109;
  atomic_mass_map[651510000] = 150.923109001;
  atomic_mass_map[661510000] = 150.926191253;
  atomic_mass_map[671510000] = 150.931698177;
  atomic_mass_map[681510000] = 150.937448567;
  atomic_mass_map[691510000] = 150.945493201;
  atomic_mass_map[701510000] = 150.955402458;
  atomic_mass_map[711510000] = 150.967677;
  atomic_mass_map[551520000] = 151.968942;
  atomic_mass_map[561520000] = 151.955222;
  atomic_mass_map[571520000] = 151.947085;
  atomic_mass_map[581520000] = 151.936682;
  atomic_mass_map[591520000] = 151.9315529;
  atomic_mass_map[601520000] = 151.924691509;
  atomic_mass_map[611520000] = 151.923505481;
  atomic_mass_map[621520000] = 151.91973904;
  atomic_mass_map[631520000] = 151.921751235;
  atomic_mass_map[641520000] = 151.919798822;
  atomic_mass_map[651520000] = 151.924082263;
  atomic_mass_map[661520000] = 151.924725363;
  atomic_mass_map[671520000] = 151.931717465;
  atomic_mass_map[681520000] = 151.935050169;
  atomic_mass_map[691520000] = 151.944476;
  atomic_mass_map[701520000] = 151.9503267;
  atomic_mass_map[711520000] = 151.96412;
  atomic_mass_map[561530000] = 152.960848;
  atomic_mass_map[571530000] = 152.950553;
  atomic_mass_map[581530000] = 152.941052;
  atomic_mass_map[591530000] = 152.933903532;
  atomic_mass_map[601530000] = 152.927717949;
  atomic_mass_map[611530000] = 152.924156436;
  atomic_mass_map[621530000] = 152.922103969;
  atomic_mass_map[631530000] = 152.921237043;
  atomic_mass_map[641530000] = 152.921757359;
  atomic_mass_map[651530000] = 152.923441978;
  atomic_mass_map[661530000] = 152.925771992;
  atomic_mass_map[671530000] = 152.930206632;
  atomic_mass_map[681530000] = 152.935084279;
  atomic_mass_map[691530000] = 152.942057244;
  atomic_mass_map[701530000] = 152.94932;
  atomic_mass_map[711530000] = 152.958805054;
  atomic_mass_map[721530000] = 152.970692;
  atomic_mass_map[561540000] = 153.964766;
  atomic_mass_map[571540000] = 153.955416;
  atomic_mass_map[581540000] = 153.94394;
  atomic_mass_map[591540000] = 153.937621738;
  atomic_mass_map[601540000] = 153.929333977;
  atomic_mass_map[611540000] = 153.926449364;
  atomic_mass_map[621540000] = 153.922216164;
  atomic_mass_map[631540000] = 153.922985955;
  atomic_mass_map[641540000] = 153.920873398;
  atomic_mass_map[651540000] = 153.924684106;
  atomic_mass_map[661540000] = 153.924429028;
  atomic_mass_map[671540000] = 153.930606841;
  atomic_mass_map[681540000] = 153.932790743;
  atomic_mass_map[691540000] = 153.941570067;
  atomic_mass_map[701540000] = 153.946395701;
  atomic_mass_map[711540000] = 153.957364;
  atomic_mass_map[721540000] = 153.964927;
  atomic_mass_map[571550000] = 154.95928;
  atomic_mass_map[581550000] = 154.948706;
  atomic_mass_map[591550000] = 154.940509259;
  atomic_mass_map[601550000] = 154.933135668;
  atomic_mass_map[611550000] = 154.928137024;
  atomic_mass_map[621550000] = 154.924647051;
  atomic_mass_map[631550000] = 154.922900102;
  atomic_mass_map[641550000] = 154.922629796;
  atomic_mass_map[651550000] = 154.923509921;
  atomic_mass_map[661550000] = 154.925758459;
  atomic_mass_map[671550000] = 154.929103634;
  atomic_mass_map[681550000] = 154.933215684;
  atomic_mass_map[691550000] = 154.939209578;
  atomic_mass_map[701550000] = 154.945783217;
  atomic_mass_map[711550000] = 154.954326011;
  atomic_mass_map[721550000] = 154.963317;
  atomic_mass_map[731550000] = 154.974312;
  atomic_mass_map[571560000] = 155.964519;
  atomic_mass_map[581560000] = 155.951884;
  atomic_mass_map[591560000] = 155.94464;
  atomic_mass_map[601560000] = 155.935078868;
  atomic_mass_map[611560000] = 155.93111749;
  atomic_mass_map[621560000] = 155.925538511;
  atomic_mass_map[631560000] = 155.924763285;
  atomic_mass_map[641560000] = 155.922130562;
  atomic_mass_map[651560000] = 155.92475443;
  atomic_mass_map[661560000] = 155.924284038;
  atomic_mass_map[671560000] = 155.929705436;
  atomic_mass_map[681560000] = 155.93106589;
  atomic_mass_map[691560000] = 155.938985597;
  atomic_mass_map[701560000] = 155.942816893;
  atomic_mass_map[711560000] = 155.953086606;
  atomic_mass_map[721560000] = 155.959401889;
  atomic_mass_map[731560000] = 155.972237;
  atomic_mass_map[581570000] = 156.957133;
  atomic_mass_map[591570000] = 156.94789;
  atomic_mass_map[601570000] = 156.939386037;
  atomic_mass_map[611570000] = 156.93312137;
  atomic_mass_map[621570000] = 156.928418673;
  atomic_mass_map[631570000] = 156.925432791;
  atomic_mass_map[641570000] = 156.92396787;
  atomic_mass_map[651570000] = 156.924032328;
  atomic_mass_map[661570000] = 156.925469667;
  atomic_mass_map[671570000] = 156.928251999;
  atomic_mass_map[681570000] = 156.931922655;
  atomic_mass_map[691570000] = 156.936973;
  atomic_mass_map[701570000] = 156.94264923;
  atomic_mass_map[711570000] = 156.950144045;
  atomic_mass_map[721570000] = 156.958236;
  atomic_mass_map[731570000] = 156.968230251;
  atomic_mass_map[741570000] = 156.979098;
  atomic_mass_map[581580000] = 157.960644;
  atomic_mass_map[591580000] = 157.95241;
  atomic_mass_map[601580000] = 157.94197;
  atomic_mass_map[611580000] = 157.936565121;
  atomic_mass_map[621580000] = 157.929950979;
  atomic_mass_map[631580000] = 157.927798581;
  atomic_mass_map[641580000] = 157.924111646;
  atomic_mass_map[651580000] = 157.925420166;
  atomic_mass_map[661580000] = 157.924414597;
  atomic_mass_map[671580000] = 157.928944692;
  atomic_mass_map[681580000] = 157.929893474;
  atomic_mass_map[691580000] = 157.936979525;
  atomic_mass_map[701580000] = 157.939870534;
  atomic_mass_map[711580000] = 157.949315626;
  atomic_mass_map[721580000] = 157.954801222;
  atomic_mass_map[731580000] = 157.966541;
  atomic_mass_map[741580000] = 157.974629;
  atomic_mass_map[591590000] = 158.95589;
  atomic_mass_map[601590000] = 158.94653;
  atomic_mass_map[611590000] = 158.939286479;
  atomic_mass_map[621590000] = 158.933217202;
  atomic_mass_map[631590000] = 158.929099612;
  atomic_mass_map[641590000] = 158.926396267;
  atomic_mass_map[651590000] = 158.925353933;
  atomic_mass_map[661590000] = 158.925746023;
  atomic_mass_map[671590000] = 158.927718768;
  atomic_mass_map[681590000] = 158.930690875;
  atomic_mass_map[691590000] = 158.934975;
  atomic_mass_map[701590000] = 158.940054787;
  atomic_mass_map[711590000] = 158.946635615;
  atomic_mass_map[721590000] = 158.953995838;
  atomic_mass_map[731590000] = 158.963028052;
  atomic_mass_map[741590000] = 158.972845;
  atomic_mass_map[751590000] = 158.984171;
  atomic_mass_map[591600000] = 159.960794;
  atomic_mass_map[601600000] = 159.9494;
  atomic_mass_map[611600000] = 159.9431;
  atomic_mass_map[621600000] = 159.935335286;
  atomic_mass_map[631600000] = 159.931850916;
  atomic_mass_map[641600000] = 159.927061537;
  atomic_mass_map[651600000] = 159.927174778;
  atomic_mass_map[661600000] = 159.925203244;
  atomic_mass_map[671600000] = 159.928735204;
  atomic_mass_map[681600000] = 159.92907713;
  atomic_mass_map[691600000] = 159.935263106;
  atomic_mass_map[701600000] = 159.937559763;
  atomic_mass_map[711600000] = 159.946033;
  atomic_mass_map[721600000] = 159.950682513;
  atomic_mass_map[731600000] = 159.961541679;
  atomic_mass_map[741600000] = 159.968516753;
  atomic_mass_map[751600000] = 159.98203;
  atomic_mass_map[601610000] = 160.95428;
  atomic_mass_map[611610000] = 160.94607;
  atomic_mass_map[621610000] = 160.939160143;
  atomic_mass_map[631610000] = 160.933664066;
  atomic_mass_map[641610000] = 160.929676602;
  atomic_mass_map[651610000] = 160.927577001;
  atomic_mass_map[661610000] = 160.926939088;
  atomic_mass_map[671610000] = 160.927860759;
  atomic_mass_map[681610000] = 160.930003191;
  atomic_mass_map[691610000] = 160.933549;
  atomic_mass_map[701610000] = 160.937906846;
  atomic_mass_map[711610000] = 160.943572;
  atomic_mass_map[721610000] = 160.950279151;
  atomic_mass_map[731610000] = 160.958369031;
  atomic_mass_map[741610000] = 160.967197;
  atomic_mass_map[751610000] = 160.977627121;
  atomic_mass_map[761610000] = 160.989287;
  atomic_mass_map[601620000] = 161.957541;
  atomic_mass_map[611620000] = 161.95022;
  atomic_mass_map[621620000] = 161.94146;
  atomic_mass_map[631620000] = 161.936979303;
  atomic_mass_map[641620000] = 161.930992146;
  atomic_mass_map[651620000] = 161.929493955;
  atomic_mass_map[661620000] = 161.926804168;
  atomic_mass_map[671620000] = 161.929101485;
  atomic_mass_map[681620000] = 161.92878696;
  atomic_mass_map[691620000] = 161.934000872;
  atomic_mass_map[701620000] = 161.935773771;
  atomic_mass_map[711620000] = 161.943282776;
  atomic_mass_map[721620000] = 161.947214896;
  atomic_mass_map[731620000] = 161.957294202;
  atomic_mass_map[741620000] = 161.963500347;
  atomic_mass_map[751620000] = 161.975844;
  atomic_mass_map[761620000] = 161.984498;
  atomic_mass_map[611630000] = 162.95357;
  atomic_mass_map[621630000] = 162.94555;
  atomic_mass_map[631630000] = 162.939360977;
  atomic_mass_map[641630000] = 162.934176832;
  atomic_mass_map[651630000] = 162.930653261;
  atomic_mass_map[661630000] = 162.928736879;
  atomic_mass_map[671630000] = 162.928739921;
  atomic_mass_map[681630000] = 162.930039567;
  atomic_mass_map[691630000] = 162.932657941;
  atomic_mass_map[701630000] = 162.9363398;
  atomic_mass_map[711630000] = 162.941179;
  atomic_mass_map[721630000] = 162.947113258;
  atomic_mass_map[731630000] = 162.954337195;
  atomic_mass_map[741630000] = 162.962524511;
  atomic_mass_map[751630000] = 162.972085441;
  atomic_mass_map[761630000] = 162.982617;
  atomic_mass_map[611640000] = 163.958271;
  atomic_mass_map[621640000] = 163.94836;
  atomic_mass_map[631640000] = 163.942693;
  atomic_mass_map[641640000] = 163.93583;
  atomic_mass_map[651640000] = 163.933356559;
  atomic_mass_map[661640000] = 163.929180472;
  atomic_mass_map[671640000] = 163.930239483;
  atomic_mass_map[681640000] = 163.929207392;
  atomic_mass_map[691640000] = 163.933543281;
  atomic_mass_map[701640000] = 163.934495103;
  atomic_mass_map[711640000] = 163.941339;
  atomic_mass_map[721640000] = 163.944370544;
  atomic_mass_map[731640000] = 163.953534;
  atomic_mass_map[741640000] = 163.958952222;
  atomic_mass_map[751640000] = 163.970507124;
  atomic_mass_map[761640000] = 163.978075966;
  atomic_mass_map[771640000] = 163.992116;
  atomic_mass_map[621650000] = 164.95297;
  atomic_mass_map[631650000] = 164.945546;
  atomic_mass_map[641650000] = 164.939395;
  atomic_mass_map[651650000] = 164.93498;
  atomic_mass_map[661650000] = 164.931709054;
  atomic_mass_map[671650000] = 164.930328047;
  atomic_mass_map[681650000] = 164.930733198;
  atomic_mass_map[691650000] = 164.932442269;
  atomic_mass_map[701650000] = 164.935270241;
  atomic_mass_map[711650000] = 164.939406758;
  atomic_mass_map[721650000] = 164.944567;
  atomic_mass_map[731650000] = 164.950780303;
  atomic_mass_map[741650000] = 164.958280974;
  atomic_mass_map[751650000] = 164.967085375;
  atomic_mass_map[761650000] = 164.976602;
  atomic_mass_map[771650000] = 164.987555;
  atomic_mass_map[621660000] = 165.956275;
  atomic_mass_map[631660000] = 165.94932;
  atomic_mass_map[641660000] = 165.94146;
  atomic_mass_map[651660000] = 165.937858119;
  atomic_mass_map[661660000] = 165.932812461;
  atomic_mass_map[671660000] = 165.932290139;
  atomic_mass_map[681660000] = 165.930299023;
  atomic_mass_map[691660000] = 165.933560092;
  atomic_mass_map[701660000] = 165.933874249;
  atomic_mass_map[711660000] = 165.939859;
  atomic_mass_map[721660000] = 165.94218;
  atomic_mass_map[731660000] = 165.950512;
  atomic_mass_map[741660000] = 165.955031346;
  atomic_mass_map[751660000] = 165.96576094;
  atomic_mass_map[761660000] = 165.972698141;
  atomic_mass_map[771660000] = 165.985664;
  atomic_mass_map[781660000] = 165.994923;
  atomic_mass_map[631670000] = 166.952753;
  atomic_mass_map[641670000] = 166.94545;
  atomic_mass_map[651670000] = 166.93996;
  atomic_mass_map[661670000] = 166.935661823;
  atomic_mass_map[671670000] = 166.933138994;
  atomic_mass_map[681670000] = 166.932054119;
  atomic_mass_map[691670000] = 166.932856635;
  atomic_mass_map[701670000] = 166.934953337;
  atomic_mass_map[711670000] = 166.93827;
  atomic_mass_map[721670000] = 166.9426;
  atomic_mass_map[731670000] = 166.948093;
  atomic_mass_map[741670000] = 166.954805873;
  atomic_mass_map[751670000] = 166.962607;
  atomic_mass_map[761670000] = 166.971548938;
  atomic_mass_map[771670000] = 166.981671981;
  atomic_mass_map[781670000] = 166.992901;
  atomic_mass_map[631680000] = 167.957337;
  atomic_mass_map[641680000] = 167.94808;
  atomic_mass_map[651680000] = 167.9434;
  atomic_mass_map[661680000] = 167.937133716;
  atomic_mass_map[671680000] = 167.935521676;
  atomic_mass_map[681680000] = 167.932376192;
  atomic_mass_map[691680000] = 167.934177868;
  atomic_mass_map[701680000] = 167.933889106;
  atomic_mass_map[711680000] = 167.938735139;
  atomic_mass_map[721680000] = 167.940568;
  atomic_mass_map[731680000] = 167.948047;
  atomic_mass_map[741680000] = 167.951805262;
  atomic_mass_map[751680000] = 167.961572608;
  atomic_mass_map[761680000] = 167.967798812;
  atomic_mass_map[771680000] = 167.979960981;
  atomic_mass_map[781680000] = 167.988183004;
  atomic_mass_map[641690000] = 168.9526;
  atomic_mass_map[651690000] = 168.94597;
  atomic_mass_map[661690000] = 168.940313971;
  atomic_mass_map[671690000] = 168.93687863;
  atomic_mass_map[681690000] = 168.934596353;
  atomic_mass_map[691690000] = 168.93421835;
  atomic_mass_map[701690000] = 168.935182016;
  atomic_mass_map[711690000] = 168.937643653;
  atomic_mass_map[721690000] = 168.941259;
  atomic_mass_map[731690000] = 168.946011;
  atomic_mass_map[741690000] = 168.951778677;
  atomic_mass_map[751690000] = 168.958765991;
  atomic_mass_map[761690000] = 168.967017833;
  atomic_mass_map[771690000] = 168.976281287;
  atomic_mass_map[781690000] = 168.986567;
  atomic_mass_map[791690000] = 168.99808;
  atomic_mass_map[641700000] = 169.955577;
  atomic_mass_map[651700000] = 169.94984;
  atomic_mass_map[661700000] = 169.94239;
  atomic_mass_map[671700000] = 169.939625289;
  atomic_mass_map[681700000] = 169.935470673;
  atomic_mass_map[691700000] = 169.935806507;
  atomic_mass_map[701700000] = 169.934767245;
  atomic_mass_map[711700000] = 169.938479234;
  atomic_mass_map[721700000] = 169.939609;
  atomic_mass_map[731700000] = 169.946175;
  atomic_mass_map[741700000] = 169.9492312;
  atomic_mass_map[751700000] = 169.958224966;
  atomic_mass_map[761700000] = 169.963578673;
  atomic_mass_map[771700000] = 169.974922;
  atomic_mass_map[781700000] = 169.982502095;
  atomic_mass_map[791700000] = 169.995972;
  atomic_mass_map[651710000] = 170.95273;
  atomic_mass_map[661710000] = 170.94612;
  atomic_mass_map[671710000] = 170.94147149;
  atomic_mass_map[681710000] = 170.938036148;
  atomic_mass_map[691710000] = 170.936435126;
  atomic_mass_map[701710000] = 170.936331517;
  atomic_mass_map[711710000] = 170.93791866;
  atomic_mass_map[721710000] = 170.940492;
  atomic_mass_map[731710000] = 170.944476;
  atomic_mass_map[741710000] = 170.949451;
  atomic_mass_map[751710000] = 170.955716;
  atomic_mass_map[761710000] = 170.963175348;
  atomic_mass_map[771710000] = 170.971645522;
  atomic_mass_map[781710000] = 170.981245502;
  atomic_mass_map[791710000] = 170.991881542;
  atomic_mass_map[801710000] = 171.003736;
  atomic_mass_map[651720000] = 171.957219;
  atomic_mass_map[661720000] = 171.94846;
  atomic_mass_map[671720000] = 171.94473;
  atomic_mass_map[681720000] = 171.939362344;
  atomic_mass_map[691720000] = 171.938406067;
  atomic_mass_map[701720000] = 171.936386658;
  atomic_mass_map[711720000] = 171.939091417;
  atomic_mass_map[721720000] = 171.939449716;
  atomic_mass_map[731720000] = 171.944895;
  atomic_mass_map[741720000] = 171.947292;
  atomic_mass_map[751720000] = 171.955408079;
  atomic_mass_map[761720000] = 171.960017088;
  atomic_mass_map[771720000] = 171.970607036;
  atomic_mass_map[781720000] = 171.977340788;
  atomic_mass_map[791720000] = 171.989996708;
  atomic_mass_map[801720000] = 171.998863391;
  atomic_mass_map[661730000] = 172.95283;
  atomic_mass_map[671730000] = 172.94702;
  atomic_mass_map[681730000] = 172.9424;
  atomic_mass_map[691730000] = 172.939606632;
  atomic_mass_map[701730000] = 172.938216215;
  atomic_mass_map[711730000] = 172.938935822;
  atomic_mass_map[721730000] = 172.940513;
  atomic_mass_map[731730000] = 172.94375;
  atomic_mass_map[741730000] = 172.947689;
  atomic_mass_map[751730000] = 172.953243;
  atomic_mass_map[761730000] = 172.959808375;
  atomic_mass_map[771730000] = 172.967505496;
  atomic_mass_map[781730000] = 172.976443315;
  atomic_mass_map[791730000] = 172.986223808;
  atomic_mass_map[801730000] = 172.997091;
  atomic_mass_map[661740000] = 173.955587;
  atomic_mass_map[671740000] = 173.95095;
  atomic_mass_map[681740000] = 173.94423;
  atomic_mass_map[691740000] = 173.942174064;
  atomic_mass_map[701740000] = 173.938867548;
  atomic_mass_map[711740000] = 173.940342938;
  atomic_mass_map[721740000] = 173.94004848;
  atomic_mass_map[731740000] = 173.944454;
  atomic_mass_map[741740000] = 173.946079;
  atomic_mass_map[751740000] = 173.953115;
  atomic_mass_map[761740000] = 173.957063152;
  atomic_mass_map[771740000] = 173.966866676;
  atomic_mass_map[781740000] = 173.972819832;
  atomic_mass_map[791740000] = 173.984718;
  atomic_mass_map[801740000] = 173.992870583;
  atomic_mass_map[671750000] = 174.95362;
  atomic_mass_map[681750000] = 174.94777;
  atomic_mass_map[691750000] = 174.943842313;
  atomic_mass_map[701750000] = 174.94128191;
  atomic_mass_map[711750000] = 174.940777308;
  atomic_mass_map[721750000] = 174.941511527;
  atomic_mass_map[731750000] = 174.943737;
  atomic_mass_map[741750000] = 174.946717;
  atomic_mass_map[751750000] = 174.951381;
  atomic_mass_map[761750000] = 174.956945105;
  atomic_mass_map[771750000] = 174.964149521;
  atomic_mass_map[781750000] = 174.972395457;
  atomic_mass_map[791750000] = 174.981316085;
  atomic_mass_map[801750000] = 174.991441086;
  atomic_mass_map[671760000] = 175.95782;
  atomic_mass_map[681760000] = 175.94994;
  atomic_mass_map[691760000] = 175.946997711;
  atomic_mass_map[701760000] = 175.942574708;
  atomic_mass_map[711760000] = 175.942691809;
  atomic_mass_map[721760000] = 175.941409905;
  atomic_mass_map[731760000] = 175.944857;
  atomic_mass_map[741760000] = 175.945634;
  atomic_mass_map[751760000] = 175.951623;
  atomic_mass_map[761760000] = 175.954806;
  atomic_mass_map[771760000] = 175.963630119;
  atomic_mass_map[781760000] = 175.968938214;
  atomic_mass_map[791760000] = 175.980116927;
  atomic_mass_map[801760000] = 175.987348335;
  atomic_mass_map[811760000] = 176.000624367;
  atomic_mass_map[681770000] = 176.95399;
  atomic_mass_map[691770000] = 176.94904;
  atomic_mass_map[701770000] = 176.945263848;
  atomic_mass_map[711770000] = 176.943763668;
  atomic_mass_map[721770000] = 176.94323032;
  atomic_mass_map[731770000] = 176.944482073;
  atomic_mass_map[741770000] = 176.946643;
  atomic_mass_map[751770000] = 176.950328;
  atomic_mass_map[761770000] = 176.954957882;
  atomic_mass_map[771770000] = 176.9613015;
  atomic_mass_map[781770000] = 176.968469529;
  atomic_mass_map[791770000] = 176.976870379;
  atomic_mass_map[801770000] = 176.986277376;
  atomic_mass_map[811770000] = 176.996413797;
  atomic_mass_map[681780000] = 177.956779;
  atomic_mass_map[691780000] = 177.95264;
  atomic_mass_map[701780000] = 177.94664971;
  atomic_mass_map[711780000] = 177.945960162;
  atomic_mass_map[721780000] = 177.943708456;
  atomic_mass_map[731780000] = 177.945681;
  atomic_mass_map[741780000] = 177.945885925;
  atomic_mass_map[751780000] = 177.950989;
  atomic_mass_map[761780000] = 177.9532533;
  atomic_mass_map[771780000] = 177.961082;
  atomic_mass_map[781780000] = 177.965649248;
  atomic_mass_map[791780000] = 177.976055945;
  atomic_mass_map[801780000] = 177.982484158;
  atomic_mass_map[811780000] = 177.994857;
  atomic_mass_map[821780000] = 178.003837163;
  atomic_mass_map[691790000] = 178.95534;
  atomic_mass_map[701790000] = 178.95004;
  atomic_mass_map[711790000] = 178.947333082;
  atomic_mass_map[721790000] = 178.945825838;
  atomic_mass_map[731790000] = 178.945939187;
  atomic_mass_map[741790000] = 178.947079501;
  atomic_mass_map[751790000] = 178.949989715;
  atomic_mass_map[761790000] = 178.953816669;
  atomic_mass_map[771790000] = 178.959117596;
  atomic_mass_map[781790000] = 178.965358719;
  atomic_mass_map[791790000] = 178.973173668;
  atomic_mass_map[801790000] = 178.981826899;
  atomic_mass_map[811790000] = 178.991123405;
  atomic_mass_map[821790000] = 179.002201452;
  atomic_mass_map[691800000] = 179.959291;
  atomic_mass_map[701800000] = 179.95212;
  atomic_mass_map[711800000] = 179.949890876;
  atomic_mass_map[721800000] = 179.946559669;
  atomic_mass_map[731800000] = 179.947468392;
  atomic_mass_map[741800000] = 179.946713435;
  atomic_mass_map[751800000] = 179.950791568;
  atomic_mass_map[761800000] = 179.95237993;
  atomic_mass_map[771800000] = 179.959229446;
  atomic_mass_map[781800000] = 179.963031563;
  atomic_mass_map[791800000] = 179.972489883;
  atomic_mass_map[801800000] = 179.978260249;
  atomic_mass_map[811800000] = 179.989923019;
  atomic_mass_map[821800000] = 179.997915842;
  atomic_mass_map[691810000] = 180.962243;
  atomic_mass_map[701810000] = 180.95589;
  atomic_mass_map[711810000] = 180.951908;
  atomic_mass_map[721810000] = 180.949110965;
  atomic_mass_map[731810000] = 180.947999331;
  atomic_mass_map[741810000] = 180.948218863;
  atomic_mass_map[751810000] = 180.950061523;
  atomic_mass_map[761810000] = 180.953247188;
  atomic_mass_map[771810000] = 180.957634694;
  atomic_mass_map[781810000] = 180.963089927;
  atomic_mass_map[791810000] = 180.970079103;
  atomic_mass_map[801810000] = 180.977819357;
  atomic_mass_map[811810000] = 180.986259992;
  atomic_mass_map[821810000] = 180.996653386;
  atomic_mass_map[701820000] = 181.958325;
  atomic_mass_map[711820000] = 181.95504;
  atomic_mass_map[721820000] = 181.950563816;
  atomic_mass_map[731820000] = 181.950155413;
  atomic_mass_map[741820000] = 181.948205721;
  atomic_mass_map[751820000] = 181.951211645;
  atomic_mass_map[761820000] = 181.952110153;
  atomic_mass_map[771820000] = 181.958076296;
  atomic_mass_map[781820000] = 181.961171571;
  atomic_mass_map[791820000] = 181.969617874;
  atomic_mass_map[801820000] = 181.974689132;
  atomic_mass_map[811820000] = 181.98569188;
  atomic_mass_map[821820000] = 181.99267294;
  atomic_mass_map[701830000] = 182.962319;
  atomic_mass_map[711830000] = 182.957363;
  atomic_mass_map[721830000] = 182.953534004;
  atomic_mass_map[731830000] = 182.95137618;
  atomic_mass_map[741830000] = 182.9502245;
  atomic_mass_map[751830000] = 182.95082139;
  atomic_mass_map[761830000] = 182.953124719;
  atomic_mass_map[771830000] = 182.956839968;
  atomic_mass_map[781830000] = 182.961596653;
  atomic_mass_map[791830000] = 182.967588108;
  atomic_mass_map[801830000] = 182.974444629;
  atomic_mass_map[811830000] = 182.982192846;
  atomic_mass_map[821830000] = 182.991867668;
  atomic_mass_map[701840000] = 183.965067;
  atomic_mass_map[711840000] = 183.96091;
  atomic_mass_map[721840000] = 183.955448587;
  atomic_mass_map[731840000] = 183.954010038;
  atomic_mass_map[741840000] = 183.95093326;
  atomic_mass_map[751840000] = 183.952528267;
  atomic_mass_map[761840000] = 183.952492949;
  atomic_mass_map[771840000] = 183.957476;
  atomic_mass_map[781840000] = 183.959920039;
  atomic_mass_map[791840000] = 183.967451524;
  atomic_mass_map[801840000] = 183.971713221;
  atomic_mass_map[811840000] = 183.981875093;
  atomic_mass_map[821840000] = 183.988135702;
  atomic_mass_map[831840000] = 184.00114125;
  atomic_mass_map[701850000] = 184.969404;
  atomic_mass_map[711850000] = 184.96362;
  atomic_mass_map[721850000] = 184.958862;
  atomic_mass_map[731850000] = 184.955561396;
  atomic_mass_map[741850000] = 184.953421286;
  atomic_mass_map[751850000] = 184.952958337;
  atomic_mass_map[761850000] = 184.954045995;
  atomic_mass_map[771850000] = 184.956698;
  atomic_mass_map[781850000] = 184.960613659;
  atomic_mass_map[791850000] = 184.965798874;
  atomic_mass_map[801850000] = 184.971890676;
  atomic_mass_map[811850000] = 184.978789191;
  atomic_mass_map[821850000] = 184.987609989;
  atomic_mass_map[831850000] = 184.9976;
  atomic_mass_map[711860000] = 185.967568;
  atomic_mass_map[721860000] = 185.960897;
  atomic_mass_map[731860000] = 185.958553111;
  atomic_mass_map[741860000] = 185.954365215;
  atomic_mass_map[751860000] = 185.954989419;
  atomic_mass_map[761860000] = 185.95383766;
  atomic_mass_map[771860000] = 185.957946754;
  atomic_mass_map[781860000] = 185.959350846;
  atomic_mass_map[791860000] = 185.965952703;
  atomic_mass_map[801860000] = 185.969362017;
  atomic_mass_map[811860000] = 185.978650841;
  atomic_mass_map[821860000] = 185.984238196;
  atomic_mass_map[831860000] = 185.996622402;
  atomic_mass_map[841860000] = 186.004402577;
  atomic_mass_map[711870000] = 186.970392;
  atomic_mass_map[721870000] = 186.96477;
  atomic_mass_map[731870000] = 186.960391;
  atomic_mass_map[741870000] = 186.957161323;
  atomic_mass_map[751870000] = 186.955752288;
  atomic_mass_map[761870000] = 186.95574964;
  atomic_mass_map[771870000] = 186.957542;
  atomic_mass_map[781870000] = 186.960616976;
  atomic_mass_map[791870000] = 186.964543155;
  atomic_mass_map[801870000] = 186.969814158;
  atomic_mass_map[811870000] = 186.975904743;
  atomic_mass_map[821870000] = 186.983910836;
  atomic_mass_map[831870000] = 186.993147276;
  atomic_mass_map[841870000] = 187.003036624;
  atomic_mass_map[711880000] = 187.97446;
  atomic_mass_map[721880000] = 187.96685;
  atomic_mass_map[731880000] = 187.963916;
  atomic_mass_map[741880000] = 187.958488395;
  atomic_mass_map[751880000] = 187.958113728;
  atomic_mass_map[761880000] = 187.955837361;
  atomic_mass_map[771880000] = 187.958835046;
  atomic_mass_map[781880000] = 187.95939756;
  atomic_mass_map[791880000] = 187.965247969;
  atomic_mass_map[801880000] = 187.96757691;
  atomic_mass_map[811880000] = 187.976020886;
  atomic_mass_map[821880000] = 187.980874592;
  atomic_mass_map[831880000] = 187.992276184;
  atomic_mass_map[841880000] = 187.999415655;
  atomic_mass_map[721890000] = 188.97084;
  atomic_mass_map[731890000] = 188.96583;
  atomic_mass_map[741890000] = 188.961763;
  atomic_mass_map[751890000] = 188.959227817;
  atomic_mass_map[761890000] = 188.958146005;
  atomic_mass_map[771890000] = 188.958722669;
  atomic_mass_map[781890000] = 188.960848542;
  atomic_mass_map[791890000] = 188.963948286;
  atomic_mass_map[801890000] = 188.968194748;
  atomic_mass_map[811890000] = 188.973573527;
  atomic_mass_map[821890000] = 188.980843639;
  atomic_mass_map[831890000] = 188.989195141;
  atomic_mass_map[841890000] = 188.998473415;
  atomic_mass_map[721900000] = 189.973129;
  atomic_mass_map[731900000] = 189.96939;
  atomic_mass_map[741900000] = 189.963089066;
  atomic_mass_map[751900000] = 189.96174336;
  atomic_mass_map[761900000] = 189.958445496;
  atomic_mass_map[771900000] = 189.960543445;
  atomic_mass_map[781900000] = 189.959949876;
  atomic_mass_map[791900000] = 189.96475175;
  atomic_mass_map[801900000] = 189.966322169;
  atomic_mass_map[811900000] = 189.973835551;
  atomic_mass_map[821900000] = 189.978081828;
  atomic_mass_map[831900000] = 189.988620883;
  atomic_mass_map[841900000] = 189.995100519;
  atomic_mass_map[731910000] = 190.97156;
  atomic_mass_map[741910000] = 190.966531;
  atomic_mass_map[751910000] = 190.963123437;
  atomic_mass_map[761910000] = 190.960928159;
  atomic_mass_map[771910000] = 190.960591527;
  atomic_mass_map[781910000] = 190.961676363;
  atomic_mass_map[791910000] = 190.963716455;
  atomic_mass_map[801910000] = 190.967158247;
  atomic_mass_map[811910000] = 190.971784096;
  atomic_mass_map[821910000] = 190.978281;
  atomic_mass_map[831910000] = 190.985786975;
  atomic_mass_map[841910000] = 190.994558488;
  atomic_mass_map[851910000] = 191.004148086;
  atomic_mass_map[731920000] = 191.97524;
  atomic_mass_map[741920000] = 191.96817;
  atomic_mass_map[751920000] = 191.966088;
  atomic_mass_map[761920000] = 191.961478881;
  atomic_mass_map[771920000] = 191.962602485;
  atomic_mass_map[781920000] = 191.961042736;
  atomic_mass_map[791920000] = 191.964817684;
  atomic_mass_map[801920000] = 191.965634182;
  atomic_mass_map[811920000] = 191.972225;
  atomic_mass_map[821920000] = 191.975785115;
  atomic_mass_map[831920000] = 191.985470078;
  atomic_mass_map[841920000] = 191.991335788;
  atomic_mass_map[851920000] = 192.003141034;
  atomic_mass_map[731930000] = 192.977595;
  atomic_mass_map[741930000] = 192.97178;
  atomic_mass_map[751930000] = 192.967545;
  atomic_mass_map[761930000] = 192.964149753;
  atomic_mass_map[771930000] = 192.962923824;
  atomic_mass_map[781930000] = 192.962984616;
  atomic_mass_map[791930000] = 192.964138447;
  atomic_mass_map[801930000] = 192.966653377;
  atomic_mass_map[811930000] = 192.970501997;
  atomic_mass_map[821930000] = 192.976173234;
  atomic_mass_map[831930000] = 192.982947223;
  atomic_mass_map[841930000] = 192.991062403;
  atomic_mass_map[851930000] = 192.999927728;
  atomic_mass_map[861930000] = 193.009707964;
  atomic_mass_map[731940000] = 193.981428;
  atomic_mass_map[741940000] = 193.97367;
  atomic_mass_map[751940000] = 193.97076;
  atomic_mass_map[761940000] = 193.965179477;
  atomic_mass_map[771940000] = 193.965075773;
  atomic_mass_map[781940000] = 193.962683527;
  atomic_mass_map[791940000] = 193.965419062;
  atomic_mass_map[801940000] = 193.965449111;
  atomic_mass_map[811940000] = 193.971081411;
  atomic_mass_map[821940000] = 193.974011706;
  atomic_mass_map[831940000] = 193.982792362;
  atomic_mass_map[841940000] = 193.988186015;
  atomic_mass_map[851940000] = 193.999226872;
  atomic_mass_map[861940000] = 194.006144424;
  atomic_mass_map[741950000] = 194.977445;
  atomic_mass_map[751950000] = 194.97254;
  atomic_mass_map[761950000] = 194.968318;
  atomic_mass_map[771950000] = 194.965976967;
  atomic_mass_map[781950000] = 194.964794353;
  atomic_mass_map[791950000] = 194.965037851;
  atomic_mass_map[801950000] = 194.966705751;
  atomic_mass_map[811950000] = 194.969774096;
  atomic_mass_map[821950000] = 194.974548743;
  atomic_mass_map[831950000] = 194.980648762;
  atomic_mass_map[841950000] = 194.988130617;
  atomic_mass_map[851950000] = 194.996274485;
  atomic_mass_map[861950000] = 195.005421699;
  atomic_mass_map[741960000] = 195.979731;
  atomic_mass_map[751960000] = 195.9758;
  atomic_mass_map[761960000] = 195.969643277;
  atomic_mass_map[771960000] = 195.968399696;
  atomic_mass_map[781960000] = 195.964954675;
  atomic_mass_map[791960000] = 195.966571221;
  atomic_mass_map[801960000] = 195.965833444;
  atomic_mass_map[811960000] = 195.970481192;
  atomic_mass_map[821960000] = 195.972787466;
  atomic_mass_map[831960000] = 195.980666509;
  atomic_mass_map[841960000] = 195.985536094;
  atomic_mass_map[851960000] = 195.995797421;
  atomic_mass_map[861960000] = 196.002115945;
  atomic_mass_map[741970000] = 196.983747;
  atomic_mass_map[751970000] = 196.97799;
  atomic_mass_map[761970000] = 196.97283;
  atomic_mass_map[771970000] = 196.969657233;
  atomic_mass_map[781970000] = 196.967343053;
  atomic_mass_map[791970000] = 196.966570114;
  atomic_mass_map[801970000] = 196.967213713;
  atomic_mass_map[811970000] = 196.969573986;
  atomic_mass_map[821970000] = 196.973434717;
  atomic_mass_map[831970000] = 196.978864929;
  atomic_mass_map[841970000] = 196.985659607;
  atomic_mass_map[851970000] = 196.993177357;
  atomic_mass_map[861970000] = 197.00162143;
  atomic_mass_map[871970000] = 197.01100809;
  atomic_mass_map[751980000] = 197.9816;
  atomic_mass_map[761980000] = 197.97441;
  atomic_mass_map[771980000] = 197.97228;
  atomic_mass_map[781980000] = 197.967896734;
  atomic_mass_map[791980000] = 197.968243724;
  atomic_mass_map[801980000] = 197.966769179;
  atomic_mass_map[811980000] = 197.970446673;
  atomic_mass_map[821980000] = 197.972015397;
  atomic_mass_map[831980000] = 197.979206;
  atomic_mass_map[841980000] = 197.983388672;
  atomic_mass_map[851980000] = 197.992791673;
  atomic_mass_map[861980000] = 197.998679156;
  atomic_mass_map[871980000] = 198.010278138;
  atomic_mass_map[751990000] = 198.984047;
  atomic_mass_map[761990000] = 198.97801;
  atomic_mass_map[771990000] = 198.973807115;
  atomic_mass_map[781990000] = 198.970597038;
  atomic_mass_map[791990000] = 198.968766582;
  atomic_mass_map[801990000] = 198.968280989;
  atomic_mass_map[811990000] = 198.969877;
  atomic_mass_map[821990000] = 198.972912542;
  atomic_mass_map[831990000] = 198.977672893;
  atomic_mass_map[841990000] = 198.983673021;
  atomic_mass_map[851990000] = 198.990527719;
  atomic_mass_map[861990000] = 198.998390273;
  atomic_mass_map[871990000] = 199.007269389;
  atomic_mass_map[762000000] = 199.97984;
  atomic_mass_map[772000000] = 199.9768;
  atomic_mass_map[782000000] = 199.971444625;
  atomic_mass_map[792000000] = 199.970756556;
  atomic_mass_map[802000000] = 199.968326934;
  atomic_mass_map[812000000] = 199.970963602;
  atomic_mass_map[822000000] = 199.971818332;
  atomic_mass_map[832000000] = 199.978131093;
  atomic_mass_map[842000000] = 199.98181227;
  atomic_mass_map[852000000] = 199.9903511;
  atomic_mass_map[862000000] = 199.995700707;
  atomic_mass_map[872000000] = 200.006583507;
  atomic_mass_map[762010000] = 200.98364;
  atomic_mass_map[772010000] = 200.97864;
  atomic_mass_map[782010000] = 200.974513293;
  atomic_mass_map[792010000] = 200.971657665;
  atomic_mass_map[802010000] = 200.970303038;
  atomic_mass_map[812010000] = 200.970820168;
  atomic_mass_map[822010000] = 200.972870425;
  atomic_mass_map[832010000] = 200.977008512;
  atomic_mass_map[842010000] = 200.98226377700001;
  atomic_mass_map[852010000] = 200.988417061;
  atomic_mass_map[862010000] = 200.995628179;
  atomic_mass_map[872010000] = 201.003852496;
  atomic_mass_map[882010000] = 201.012814683;
  atomic_mass_map[762020000] = 201.98595;
  atomic_mass_map[772020000] = 201.98199;
  atomic_mass_map[782020000] = 201.975639;
  atomic_mass_map[792020000] = 201.973856;
  atomic_mass_map[802020000] = 201.970643585;
  atomic_mass_map[812020000] = 201.972109089;
  atomic_mass_map[822020000] = 201.972151604;
  atomic_mass_map[832020000] = 201.9777331;
  atomic_mass_map[842020000] = 201.980738881;
  atomic_mass_map[852020000] = 201.98863038;
  atomic_mass_map[862020000] = 201.993263902;
  atomic_mass_map[872020000] = 202.003323946;
  atomic_mass_map[882020000] = 202.009742264;
  atomic_mass_map[762030000] = 202.991798;
  atomic_mass_map[772030000] = 202.98423;
  atomic_mass_map[782030000] = 202.97893;
  atomic_mass_map[792030000] = 202.975154498;
  atomic_mass_map[802030000] = 202.972872326;
  atomic_mass_map[812030000] = 202.972344022;
  atomic_mass_map[822030000] = 202.973390535;
  atomic_mass_map[832030000] = 202.976892145;
  atomic_mass_map[842030000] = 202.981415995;
  atomic_mass_map[852030000] = 202.986942957;
  atomic_mass_map[862030000] = 202.993393732;
  atomic_mass_map[872030000] = 203.000940872;
  atomic_mass_map[882030000] = 203.009298745;
  atomic_mass_map[772040000] = 203.9896;
  atomic_mass_map[782040000] = 203.98076;
  atomic_mass_map[792040000] = 203.977831;
  atomic_mass_map[802040000] = 203.973494037;
  atomic_mass_map[812040000] = 203.973863337;
  atomic_mass_map[822040000] = 203.97304342;
  atomic_mass_map[832040000] = 203.977835717;
  atomic_mass_map[842040000] = 203.980309863;
  atomic_mass_map[852040000] = 203.987251197;
  atomic_mass_map[862040000] = 203.991443644;
  atomic_mass_map[872040000] = 204.000651974;
  atomic_mass_map[882040000] = 204.006502228;
  atomic_mass_map[772050000] = 204.993602;
  atomic_mass_map[782050000] = 204.98608;
  atomic_mass_map[792050000] = 204.97985;
  atomic_mass_map[802050000] = 204.976073125;
  atomic_mass_map[812050000] = 204.974427237;
  atomic_mass_map[822050000] = 204.974481597;
  atomic_mass_map[832050000] = 204.977386323;
  atomic_mass_map[842050000] = 204.981190004;
  atomic_mass_map[852050000] = 204.986074041;
  atomic_mass_map[862050000] = 204.991723204;
  atomic_mass_map[872050000] = 204.998593858;
  atomic_mass_map[882050000] = 205.006268415;
  atomic_mass_map[892050000] = 205.015144158;
  atomic_mass_map[782060000] = 205.98966;
  atomic_mass_map[792060000] = 205.98474;
  atomic_mass_map[802060000] = 205.977513756;
  atomic_mass_map[812060000] = 205.976110026;
  atomic_mass_map[822060000] = 205.974465124;
  atomic_mass_map[832060000] = 205.978498757;
  atomic_mass_map[842060000] = 205.980473654;
  atomic_mass_map[852060000] = 205.986656148;
  atomic_mass_map[862060000] = 205.990195358;
  atomic_mass_map[872060000] = 205.998666211;
  atomic_mass_map[882060000] = 206.003827763;
  atomic_mass_map[892060000] = 206.014470787;
  atomic_mass_map[782070000] = 206.995126;
  atomic_mass_map[792070000] = 206.988395;
  atomic_mass_map[802070000] = 206.9823;
  atomic_mass_map[812070000] = 206.977418586;
  atomic_mass_map[822070000] = 206.975896735;
  atomic_mass_map[832070000] = 206.978470471;
  atomic_mass_map[842070000] = 206.981593252;
  atomic_mass_map[852070000] = 206.985799783;
  atomic_mass_map[862070000] = 206.9907302;
  atomic_mass_map[872070000] = 206.996946474;
  atomic_mass_map[882070000] = 207.003805161;
  atomic_mass_map[892070000] = 207.011965973;
  atomic_mass_map[782080000] = 207.998937;
  atomic_mass_map[792080000] = 207.99345;
  atomic_mass_map[802080000] = 207.985759;
  atomic_mass_map[812080000] = 207.982017992;
  atomic_mass_map[822080000] = 207.976651918;
  atomic_mass_map[832080000] = 207.979741981;
  atomic_mass_map[842080000] = 207.981245616;
  atomic_mass_map[852080000] = 207.986613042;
  atomic_mass_map[862080000] = 207.989634295;
  atomic_mass_map[872080000] = 207.997138018;
  atomic_mass_map[882080000] = 208.001854929;
  atomic_mass_map[892080000] = 208.011544073;
  atomic_mass_map[902080000] = 208.017910722;
  atomic_mass_map[792090000] = 208.997273;
  atomic_mass_map[802090000] = 208.99072;
  atomic_mass_map[812090000] = 208.98535175;
  atomic_mass_map[822090000] = 208.981089898;
  atomic_mass_map[832090000] = 208.980398519;
  atomic_mass_map[842090000] = 208.982430276;
  atomic_mass_map[852090000] = 208.986169944;
  atomic_mass_map[862090000] = 208.990401388;
  atomic_mass_map[872090000] = 208.995953197;
  atomic_mass_map[882090000] = 209.001994879;
  atomic_mass_map[892090000] = 209.00949422;
  atomic_mass_map[902090000] = 209.017571;
  atomic_mass_map[792100000] = 210.0025;
  atomic_mass_map[802100000] = 209.99424;
  atomic_mass_map[812100000] = 209.99007297;
  atomic_mass_map[822100000] = 209.984188301;
  atomic_mass_map[832100000] = 209.984120156;
  atomic_mass_map[842100000] = 209.982873601;
  atomic_mass_map[852100000] = 209.987147338;
  atomic_mass_map[862100000] = 209.989688854;
  atomic_mass_map[872100000] = 209.996421657;
  atomic_mass_map[882100000] = 210.000475356;
  atomic_mass_map[892100000] = 210.00943613;
  atomic_mass_map[902100000] = 210.015093437;
  atomic_mass_map[802110000] = 210.99933;
  atomic_mass_map[812110000] = 210.993475;
  atomic_mass_map[822110000] = 210.988735356;
  atomic_mass_map[832110000] = 210.987268698;
  atomic_mass_map[842110000] = 210.986653085;
  atomic_mass_map[852110000] = 210.987496147;
  atomic_mass_map[862110000] = 210.990600686;
  atomic_mass_map[872110000] = 210.995555259;
  atomic_mass_map[882110000] = 211.000893213;
  atomic_mass_map[892110000] = 211.007731894;
  atomic_mass_map[902110000] = 211.014933183;
  atomic_mass_map[912110000] = 211.023704;
  atomic_mass_map[802120000] = 212.00296;
  atomic_mass_map[812120000] = 211.998335;
  atomic_mass_map[822120000] = 211.991895975;
  atomic_mass_map[832120000] = 211.991285016;
  atomic_mass_map[842120000] = 211.988867896;
  atomic_mass_map[852120000] = 211.990737223;
  atomic_mass_map[862120000] = 211.990703528;
  atomic_mass_map[872120000] = 211.996225453;
  atomic_mass_map[882120000] = 211.999786399;
  atomic_mass_map[892120000] = 212.007812501;
  atomic_mass_map[902120000] = 212.013001487;
  atomic_mass_map[912120000] = 212.023181425;
  atomic_mass_map[802130000] = 213.00823;
  atomic_mass_map[812130000] = 213.001915;
  atomic_mass_map[822130000] = 212.996560867;
  atomic_mass_map[832130000] = 212.994383608;
  atomic_mass_map[842130000] = 212.992857083;
  atomic_mass_map[852130000] = 212.992936514;
  atomic_mass_map[862130000] = 212.993885064;
  atomic_mass_map[872130000] = 212.996185861;
  atomic_mass_map[882130000] = 213.00037097;
  atomic_mass_map[892130000] = 213.006607333;
  atomic_mass_map[902130000] = 213.013011447;
  atomic_mass_map[912130000] = 213.021108697;
  atomic_mass_map[802140000] = 214.012;
  atomic_mass_map[812140000] = 214.00694;
  atomic_mass_map[822140000] = 213.999803788;
  atomic_mass_map[832140000] = 213.998710938;
  atomic_mass_map[842140000] = 213.995201208;
  atomic_mass_map[852140000] = 213.996371601;
  atomic_mass_map[862140000] = 213.995362566;
  atomic_mass_map[872140000] = 213.998970785;
  atomic_mass_map[882140000] = 214.000099554;
  atomic_mass_map[892140000] = 214.006917762;
  atomic_mass_map[902140000] = 214.011481431;
  atomic_mass_map[912140000] = 214.020918561;
  atomic_mass_map[802150000] = 215.0174;
  atomic_mass_map[812150000] = 215.01064;
  atomic_mass_map[822150000] = 215.00466159;
  atomic_mass_map[832150000] = 215.001749149;
  atomic_mass_map[842150000] = 214.999418454;
  atomic_mass_map[852150000] = 214.99865189;
  atomic_mass_map[862150000] = 214.998745498;
  atomic_mass_map[872150000] = 215.000341456;
  atomic_mass_map[882150000] = 215.00272008;
  atomic_mass_map[892150000] = 215.006474132;
  atomic_mass_map[902150000] = 215.011724805;
  atomic_mass_map[912150000] = 215.019177728;
  atomic_mass_map[922150000] = 215.026756035;
  atomic_mass_map[802160000] = 216.02132;
  atomic_mass_map[812160000] = 216.0158;
  atomic_mass_map[822160000] = 216.00803;
  atomic_mass_map[832160000] = 216.006305989;
  atomic_mass_map[842160000] = 216.001913506;
  atomic_mass_map[852160000] = 216.002422631;
  atomic_mass_map[862160000] = 216.000271464;
  atomic_mass_map[872160000] = 216.003189445;
  atomic_mass_map[882160000] = 216.003533117;
  atomic_mass_map[892160000] = 216.008743367;
  atomic_mass_map[902160000] = 216.011055714;
  atomic_mass_map[912160000] = 216.019108242;
  atomic_mass_map[922160000] = 216.024762747;
  atomic_mass_map[812170000] = 217.01966;
  atomic_mass_map[822170000] = 217.01314;
  atomic_mass_map[832170000] = 217.009372;
  atomic_mass_map[842170000] = 217.006316216;
  atomic_mass_map[852170000] = 217.004717835;
  atomic_mass_map[862170000] = 217.003927562;
  atomic_mass_map[872170000] = 217.004631902;
  atomic_mass_map[882170000] = 217.006322806;
  atomic_mass_map[892170000] = 217.009343777;
  atomic_mass_map[902170000] = 217.013103444;
  atomic_mass_map[912170000] = 217.018323692;
  atomic_mass_map[922170000] = 217.024663;
  atomic_mass_map[812180000] = 218.024885;
  atomic_mass_map[822180000] = 218.01659;
  atomic_mass_map[832180000] = 218.014188;
  atomic_mass_map[842180000] = 218.008971502;
  atomic_mass_map[852180000] = 218.008693735;
  atomic_mass_map[862180000] = 218.005601052;
  atomic_mass_map[872180000] = 218.007578274;
  atomic_mass_map[882180000] = 218.007140325;
  atomic_mass_map[892180000] = 218.011641093;
  atomic_mass_map[902180000] = 218.013276242;
  atomic_mass_map[912180000] = 218.020057853;
  atomic_mass_map[922180000] = 218.023504829;
  atomic_mass_map[822190000] = 219.02177;
  atomic_mass_map[832190000] = 219.01748;
  atomic_mass_map[842190000] = 219.013614;
  atomic_mass_map[852190000] = 219.011160647;
  atomic_mass_map[862190000] = 219.009478753;
  atomic_mass_map[872190000] = 219.009251553;
  atomic_mass_map[882190000] = 219.010085176;
  atomic_mass_map[892190000] = 219.012420348;
  atomic_mass_map[902190000] = 219.015535677;
  atomic_mass_map[912190000] = 219.01990365;
  atomic_mass_map[922190000] = 219.024999161;
  atomic_mass_map[932190000] = 219.031623021;
  atomic_mass_map[822200000] = 220.02541;
  atomic_mass_map[832200000] = 220.02235;
  atomic_mass_map[842200000] = 220.016386;
  atomic_mass_map[852200000] = 220.015433;
  atomic_mass_map[862200000] = 220.011392534;
  atomic_mass_map[872200000] = 220.012326778;
  atomic_mass_map[882200000] = 220.011025562;
  atomic_mass_map[892200000] = 220.01475445;
  atomic_mass_map[902200000] = 220.015747926;
  atomic_mass_map[912200000] = 220.021705;
  atomic_mass_map[922200000] = 220.02462;
  atomic_mass_map[932200000] = 220.03254;
  atomic_mass_map[832210000] = 221.02587;
  atomic_mass_map[842210000] = 221.021228;
  atomic_mass_map[852210000] = 221.018017;
  atomic_mass_map[862210000] = 221.015535709;
  atomic_mass_map[872210000] = 221.014253757;
  atomic_mass_map[882210000] = 221.013917224;
  atomic_mass_map[892210000] = 221.015591199;
  atomic_mass_map[902210000] = 221.018186236;
  atomic_mass_map[912210000] = 221.021874846;
  atomic_mass_map[922210000] = 221.026323299;
  atomic_mass_map[932210000] = 221.032045;
  atomic_mass_map[832220000] = 222.030842;
  atomic_mass_map[842220000] = 222.02414;
  atomic_mass_map[852220000] = 222.022494;
  atomic_mass_map[862220000] = 222.017576286;
  atomic_mass_map[872220000] = 222.01758262;
  atomic_mass_map[882220000] = 222.015373355;
  atomic_mass_map[892220000] = 222.017843887;
  atomic_mass_map[902220000] = 222.0184683;
  atomic_mass_map[912220000] = 222.023784;
  atomic_mass_map[922220000] = 222.026057953;
  atomic_mass_map[932220000] = 222.0333;
  atomic_mass_map[832230000] = 223.0345;
  atomic_mass_map[842230000] = 223.02907;
  atomic_mass_map[852230000] = 223.025151;
  atomic_mass_map[862230000] = 223.021889285;
  atomic_mass_map[872230000] = 223.019734313;
  atomic_mass_map[882230000] = 223.018500719;
  atomic_mass_map[892230000] = 223.019136872;
  atomic_mass_map[902230000] = 223.020811546;
  atomic_mass_map[912230000] = 223.023962232;
  atomic_mass_map[922230000] = 223.027737168;
  atomic_mass_map[932230000] = 223.03285;
  atomic_mass_map[832240000] = 224.039539;
  atomic_mass_map[842240000] = 224.03211;
  atomic_mass_map[852240000] = 224.029749;
  atomic_mass_map[862240000] = 224.024095804;
  atomic_mass_map[872240000] = 224.0233481;
  atomic_mass_map[882240000] = 224.020210453;
  atomic_mass_map[892240000] = 224.021722239;
  atomic_mass_map[902240000] = 224.021464157;
  atomic_mass_map[912240000] = 224.02561721;
  atomic_mass_map[922240000] = 224.027613974;
  atomic_mass_map[932240000] = 224.03422;
  atomic_mass_map[842250000] = 225.03707;
  atomic_mass_map[852250000] = 225.03263;
  atomic_mass_map[862250000] = 225.028485574;
  atomic_mass_map[872250000] = 225.025572478;
  atomic_mass_map[882250000] = 225.023610574;
  atomic_mass_map[892250000] = 225.023228647;
  atomic_mass_map[902250000] = 225.023950907;
  atomic_mass_map[912250000] = 225.026130844;
  atomic_mass_map[922250000] = 225.029393555;
  atomic_mass_map[932250000] = 225.033910797;
  atomic_mass_map[842260000] = 226.04031;
  atomic_mass_map[852260000] = 226.03716;
  atomic_mass_map[862260000] = 226.030861382;
  atomic_mass_map[872260000] = 226.029544515;
  atomic_mass_map[882260000] = 226.025408455;
  atomic_mass_map[892260000] = 226.026097069;
  atomic_mass_map[902260000] = 226.024903686;
  atomic_mass_map[912260000] = 226.027947872;
  atomic_mass_map[922260000] = 226.029338749;
  atomic_mass_map[932260000] = 226.035188;
  atomic_mass_map[842270000] = 227.04539;
  atomic_mass_map[852270000] = 227.04024;
  atomic_mass_map[862270000] = 227.035304396;
  atomic_mass_map[872270000] = 227.031865417;
  atomic_mass_map[882270000] = 227.029176474;
  atomic_mass_map[892270000] = 227.027750666;
  atomic_mass_map[902270000] = 227.027702618;
  atomic_mass_map[912270000] = 227.028804477;
  atomic_mass_map[922270000] = 227.031181587;
  atomic_mass_map[932270000] = 227.034956832;
  atomic_mass_map[942270000] = 227.039474;
  atomic_mass_map[852280000] = 228.04475;
  atomic_mass_map[862280000] = 228.037835418;
  atomic_mass_map[872280000] = 228.035839437;
  atomic_mass_map[882280000] = 228.031068657;
  atomic_mass_map[892280000] = 228.031019767;
  atomic_mass_map[902280000] = 228.028739835;
  atomic_mass_map[912280000] = 228.031050748;
  atomic_mass_map[922280000] = 228.031371351;
  atomic_mass_map[932280000] = 228.036066462;
  atomic_mass_map[942280000] = 228.038741387;
  atomic_mass_map[852290000] = 229.04812;
  atomic_mass_map[862290000] = 229.042257276;
  atomic_mass_map[872290000] = 229.038291455;
  atomic_mass_map[882290000] = 229.034956707;
  atomic_mass_map[892290000] = 229.032947;
  atomic_mass_map[902290000] = 229.031761431;
  atomic_mass_map[912290000] = 229.032095652;
  atomic_mass_map[922290000] = 229.033505909;
  atomic_mass_map[932290000] = 229.036263974;
  atomic_mass_map[942290000] = 229.040145819;
  atomic_mass_map[952290000] = 229.045249909;
  atomic_mass_map[862300000] = 230.04514;
  atomic_mass_map[872300000] = 230.042390791;
  atomic_mass_map[882300000] = 230.03705478;
  atomic_mass_map[892300000] = 230.036327;
  atomic_mass_map[902300000] = 230.033132358;
  atomic_mass_map[912300000] = 230.034539789;
  atomic_mass_map[922300000] = 230.033940102;
  atomic_mass_map[932300000] = 230.037827716;
  atomic_mass_map[942300000] = 230.039650703;
  atomic_mass_map[952300000] = 230.046089;
  atomic_mass_map[862310000] = 231.04987;
  atomic_mass_map[872310000] = 231.045175357;
  atomic_mass_map[882310000] = 231.041027086;
  atomic_mass_map[892310000] = 231.038393;
  atomic_mass_map[902310000] = 231.036302853;
  atomic_mass_map[912310000] = 231.035882575;
  atomic_mass_map[922310000] = 231.036292252;
  atomic_mass_map[932310000] = 231.03824449;
  atomic_mass_map[942310000] = 231.04112641;
  atomic_mass_map[952310000] = 231.045529;
  atomic_mass_map[962310000] = 231.050746;
  atomic_mass_map[872320000] = 232.049461224;
  atomic_mass_map[882320000] = 232.04347527;
  atomic_mass_map[892320000] = 232.042034;
  atomic_mass_map[902320000] = 232.038053689;
  atomic_mass_map[912320000] = 232.0385903;
  atomic_mass_map[922320000] = 232.03715486;
  atomic_mass_map[932320000] = 232.040107;
  atomic_mass_map[942320000] = 232.041184526;
  atomic_mass_map[952320000] = 232.046527;
  atomic_mass_map[962320000] = 232.049718;
  atomic_mass_map[872330000] = 233.052517838;
  atomic_mass_map[882330000] = 233.047594573;
  atomic_mass_map[892330000] = 233.044346;
  atomic_mass_map[902330000] = 233.041580208;
  atomic_mass_map[912330000] = 233.040246605;
  atomic_mass_map[922330000] = 233.039634367;
  atomic_mass_map[932330000] = 233.040739489;
  atomic_mass_map[942330000] = 233.042997345;
  atomic_mass_map[952330000] = 233.046445;
  atomic_mass_map[962330000] = 233.050772206;
  atomic_mass_map[972330000] = 233.056748;
  atomic_mass_map[882340000] = 234.050382104;
  atomic_mass_map[892340000] = 234.048139;
  atomic_mass_map[902340000] = 234.04359986;
  atomic_mass_map[912340000] = 234.043305615;
  atomic_mass_map[922340000] = 234.04095037;
  atomic_mass_map[932340000] = 234.04289332;
  atomic_mass_map[942340000] = 234.043317478;
  atomic_mass_map[952340000] = 234.047731;
  atomic_mass_map[962340000] = 234.050160959;
  atomic_mass_map[972340000] = 234.057387;
  atomic_mass_map[882350000] = 235.05489;
  atomic_mass_map[892350000] = 235.05084;
  atomic_mass_map[902350000] = 235.047255;
  atomic_mass_map[912350000] = 235.045399;
  atomic_mass_map[922350000] = 235.04392819;
  atomic_mass_map[932350000] = 235.044061591;
  atomic_mass_map[942350000] = 235.045284682;
  atomic_mass_map[952350000] = 235.047907371;
  atomic_mass_map[962350000] = 235.051567;
  atomic_mass_map[972350000] = 235.05658;
  atomic_mass_map[892360000] = 236.054988;
  atomic_mass_map[902360000] = 236.049657;
  atomic_mass_map[912360000] = 236.048668;
  atomic_mass_map[922360000] = 236.045566201;
  atomic_mass_map[932360000] = 236.046568392;
  atomic_mass_map[942360000] = 236.046056756;
  atomic_mass_map[952360000] = 236.049427;
  atomic_mass_map[962360000] = 236.051374506;
  atomic_mass_map[972360000] = 236.05748;
  atomic_mass_map[892370000] = 237.057993;
  atomic_mass_map[902370000] = 237.053629;
  atomic_mass_map[912370000] = 237.051023;
  atomic_mass_map[922370000] = 237.04872838;
  atomic_mass_map[932370000] = 237.04817171;
  atomic_mass_map[942370000] = 237.048407957;
  atomic_mass_map[952370000] = 237.049995;
  atomic_mass_map[962370000] = 237.052868923;
  atomic_mass_map[972370000] = 237.0571;
  atomic_mass_map[982370000] = 237.062199993;
  atomic_mass_map[902380000] = 238.056388;
  atomic_mass_map[912380000] = 238.054637;
  atomic_mass_map[922380000] = 238.050786996;
  atomic_mass_map[932380000] = 238.050944671;
  atomic_mass_map[942380000] = 238.04955825;
  atomic_mass_map[952380000] = 238.051982607;
  atomic_mass_map[962380000] = 238.053081595;
  atomic_mass_map[972380000] = 238.058203;
  atomic_mass_map[982380000] = 238.06149;
  atomic_mass_map[902390000] = 239.060602;
  atomic_mass_map[912390000] = 239.05726;
  atomic_mass_map[922390000] = 239.054292048;
  atomic_mass_map[932390000] = 239.052937599;
  atomic_mass_map[942390000] = 239.052161669;
  atomic_mass_map[952390000] = 239.053022803;
  atomic_mass_map[962390000] = 239.054908593;
  atomic_mass_map[972390000] = 239.05824;
  atomic_mass_map[982390000] = 239.062554;
  atomic_mass_map[992390000] = 239.06823;
  atomic_mass_map[912400000] = 240.061095;
  atomic_mass_map[922400000] = 240.056592425;
  atomic_mass_map[932400000] = 240.05616383;
  atomic_mass_map[942400000] = 240.053811812;
  atomic_mass_map[952400000] = 240.055298444;
  atomic_mass_map[962400000] = 240.055528329;
  atomic_mass_map[972400000] = 240.059758;
  atomic_mass_map[982400000] = 240.062255842;
  atomic_mass_map[992400000] = 240.06892;
  atomic_mass_map[912410000] = 241.064026;
  atomic_mass_map[922410000] = 241.06033;
  atomic_mass_map[932410000] = 241.058250697;
  atomic_mass_map[942410000] = 241.056849722;
  atomic_mass_map[952410000] = 241.056827413;
  atomic_mass_map[962410000] = 241.057651288;
  atomic_mass_map[972410000] = 241.060153;
  atomic_mass_map[982410000] = 241.06369;
  atomic_mass_map[992410000] = 241.06856;
  atomic_mass_map[1002410000] = 241.07421;
  atomic_mass_map[922420000] = 242.062931;
  atomic_mass_map[932420000] = 242.061639615;
  atomic_mass_map[942420000] = 242.058741045;
  atomic_mass_map[952420000] = 242.059547428;
  atomic_mass_map[962420000] = 242.058834263;
  atomic_mass_map[972420000] = 242.06198;
  atomic_mass_map[982420000] = 242.063754533;
  atomic_mass_map[992420000] = 242.069567;
  atomic_mass_map[1002420000] = 242.07343;
  atomic_mass_map[922430000] = 243.066946;
  atomic_mass_map[932430000] = 243.064279;
  atomic_mass_map[942430000] = 243.062002119;
  atomic_mass_map[952430000] = 243.06137994;
  atomic_mass_map[962430000] = 243.061387403;
  atomic_mass_map[972430000] = 243.06300598;
  atomic_mass_map[982430000] = 243.065475;
  atomic_mass_map[992430000] = 243.069509;
  atomic_mass_map[1002430000] = 243.07449;
  atomic_mass_map[932440000] = 244.06785;
  atomic_mass_map[942440000] = 244.064204415;
  atomic_mass_map[952440000] = 244.064282964;
  atomic_mass_map[962440000] = 244.062750694;
  atomic_mass_map[972440000] = 244.065179039;
  atomic_mass_map[982440000] = 244.065999543;
  atomic_mass_map[992440000] = 244.070881;
  atomic_mass_map[1002440000] = 244.074038;
  atomic_mass_map[932450000] = 245.070736;
  atomic_mass_map[942450000] = 245.067824568;
  atomic_mass_map[952450000] = 245.06645289;
  atomic_mass_map[962450000] = 245.065491113;
  atomic_mass_map[972450000] = 245.066359885;
  atomic_mass_map[982450000] = 245.068046825;
  atomic_mass_map[992450000] = 245.071247;
  atomic_mass_map[1002450000] = 245.075349;
  atomic_mass_map[1012450000] = 245.080808;
  atomic_mass_map[942460000] = 246.070204209;
  atomic_mass_map[952460000] = 246.069774;
  atomic_mass_map[962460000] = 246.067222082;
  atomic_mass_map[972460000] = 246.068671367;
  atomic_mass_map[982460000] = 246.068803762;
  atomic_mass_map[992460000] = 246.072894;
  atomic_mass_map[1002460000] = 246.075350815;
  atomic_mass_map[1012460000] = 246.081713;
  atomic_mass_map[942470000] = 247.07419;
  atomic_mass_map[952470000] = 247.072092;
  atomic_mass_map[962470000] = 247.070352726;
  atomic_mass_map[972470000] = 247.07030594;
  atomic_mass_map[982470000] = 247.070965462;
  atomic_mass_map[992470000] = 247.073621932;
  atomic_mass_map[1002470000] = 247.076944;
  atomic_mass_map[1012470000] = 247.081521;
  atomic_mass_map[952480000] = 248.075752;
  atomic_mass_map[962480000] = 248.072349101;
  atomic_mass_map[972480000] = 248.073087;
  atomic_mass_map[982480000] = 248.072182978;
  atomic_mass_map[992480000] = 248.075469;
  atomic_mass_map[1002480000] = 248.077185528;
  atomic_mass_map[1012480000] = 248.082822;
  atomic_mass_map[1022480000] = 248.08655;
  atomic_mass_map[952490000] = 249.07848;
  atomic_mass_map[962490000] = 249.075954006;
  atomic_mass_map[972490000] = 249.074983182;
  atomic_mass_map[982490000] = 249.074850491;
  atomic_mass_map[992490000] = 249.076409;
  atomic_mass_map[1002490000] = 249.078926098;
  atomic_mass_map[1012490000] = 249.082912;
  atomic_mass_map[1022490000] = 249.087797;
  atomic_mass_map[962500000] = 250.078357556;
  atomic_mass_map[972500000] = 250.078315027;
  atomic_mass_map[982500000] = 250.076404561;
  atomic_mass_map[992500000] = 250.078611;
  atomic_mass_map[1002500000] = 250.079519828;
  atomic_mass_map[1012500000] = 250.084413;
  atomic_mass_map[1022500000] = 250.087562;
  atomic_mass_map[962510000] = 251.082285036;
  atomic_mass_map[972510000] = 251.080760603;
  atomic_mass_map[982510000] = 251.079587219;
  atomic_mass_map[992510000] = 251.079992224;
  atomic_mass_map[1002510000] = 251.081539889;
  atomic_mass_map[1012510000] = 251.084774291;
  atomic_mass_map[1022510000] = 251.088942;
  atomic_mass_map[1032510000] = 251.09418;
  atomic_mass_map[962520000] = 252.08487;
  atomic_mass_map[972520000] = 252.08431;
  atomic_mass_map[982520000] = 252.081626523;
  atomic_mass_map[992520000] = 252.082979189;
  atomic_mass_map[1002520000] = 252.082464972;
  atomic_mass_map[1012520000] = 252.086432;
  atomic_mass_map[1022520000] = 252.088966141;
  atomic_mass_map[1032520000] = 252.095263;
  atomic_mass_map[972530000] = 253.08688;
  atomic_mass_map[982530000] = 253.085133738;
  atomic_mass_map[992530000] = 253.084821305;
  atomic_mass_map[1002530000] = 253.08518116;
  atomic_mass_map[1012530000] = 253.087143;
  atomic_mass_map[1022530000] = 253.090562831;
  atomic_mass_map[1032530000] = 253.095089;
  atomic_mass_map[1042530000] = 253.100438;
  atomic_mass_map[972540000] = 254.0906;
  atomic_mass_map[982540000] = 254.08732359;
  atomic_mass_map[992540000] = 254.088020527;
  atomic_mass_map[1002540000] = 254.086852726;
  atomic_mass_map[1012540000] = 254.08959;
  atomic_mass_map[1022540000] = 254.090954259;
  atomic_mass_map[1032540000] = 254.096481;
  atomic_mass_map[1042540000] = 254.100053;
  atomic_mass_map[982550000] = 255.091047;
  atomic_mass_map[992550000] = 255.090273553;
  atomic_mass_map[1002550000] = 255.089962633;
  atomic_mass_map[1012550000] = 255.091082787;
  atomic_mass_map[1022550000] = 255.093191404;
  atomic_mass_map[1032550000] = 255.096562404;
  atomic_mass_map[1042550000] = 255.101267;
  atomic_mass_map[1052550000] = 255.106918;
  atomic_mass_map[982560000] = 256.093442;
  atomic_mass_map[992560000] = 256.093599;
  atomic_mass_map[1002560000] = 256.091773878;
  atomic_mass_map[1012560000] = 256.093888;
  atomic_mass_map[1022560000] = 256.094280866;
  atomic_mass_map[1032560000] = 256.098494029;
  atomic_mass_map[1042560000] = 256.101151535;
  atomic_mass_map[1052560000] = 256.107889;
  atomic_mass_map[992570000] = 257.095979;
  atomic_mass_map[1002570000] = 257.095105317;
  atomic_mass_map[1012570000] = 257.095537977;
  atomic_mass_map[1022570000] = 257.096884419;
  atomic_mass_map[1032570000] = 257.09948;
  atomic_mass_map[1042570000] = 257.102916848;
  atomic_mass_map[1052570000] = 257.107576;
  atomic_mass_map[992580000] = 258.09952;
  atomic_mass_map[1002580000] = 258.097077;
  atomic_mass_map[1012580000] = 258.098429825;
  atomic_mass_map[1022580000] = 258.098205;
  atomic_mass_map[1032580000] = 258.101753;
  atomic_mass_map[1042580000] = 258.103426362;
  atomic_mass_map[1052580000] = 258.109284;
  atomic_mass_map[1062580000] = 258.112984;
  atomic_mass_map[1002590000] = 259.100596;
  atomic_mass_map[1012590000] = 259.10051;
  atomic_mass_map[1022590000] = 259.100997503;
  atomic_mass_map[1032590000] = 259.102901;
  atomic_mass_map[1042590000] = 259.105596;
  atomic_mass_map[1052590000] = 259.109491865;
  atomic_mass_map[1062590000] = 259.114353;
  atomic_mass_map[1002600000] = 260.102809;
  atomic_mass_map[1012600000] = 260.103653;
  atomic_mass_map[1022600000] = 260.102643;
  atomic_mass_map[1032600000] = 260.105504;
  atomic_mass_map[1042600000] = 260.106439;
  atomic_mass_map[1052600000] = 260.111297;
  atomic_mass_map[1062600000] = 260.114383508;
  atomic_mass_map[1072600000] = 260.121658;
  atomic_mass_map[1012610000] = 261.105828;
  atomic_mass_map[1022610000] = 261.105696;
  atomic_mass_map[1032610000] = 261.10688;
  atomic_mass_map[1042610000] = 261.10876999;
  atomic_mass_map[1052610000] = 261.11198;
  atomic_mass_map[1062610000] = 261.115948188;
  atomic_mass_map[1072610000] = 261.121454;
  atomic_mass_map[1012620000] = 262.109101;
  atomic_mass_map[1022620000] = 262.107463;
  atomic_mass_map[1032620000] = 262.109611;
  atomic_mass_map[1042620000] = 262.109923;
  atomic_mass_map[1052620000] = 262.114068;
  atomic_mass_map[1062620000] = 262.116335446;
  atomic_mass_map[1072620000] = 262.122965;
  atomic_mass_map[1022630000] = 263.110714;
  atomic_mass_map[1032630000] = 263.111358;
  atomic_mass_map[1042630000] = 263.11246;
  atomic_mass_map[1052630000] = 263.114988;
  atomic_mass_map[1062630000] = 263.118294;
  atomic_mass_map[1072630000] = 263.122916;
  atomic_mass_map[1082630000] = 263.12848;
  atomic_mass_map[1022640000] = 264.112734;
  atomic_mass_map[1032640000] = 264.1142;
  atomic_mass_map[1042640000] = 264.113878;
  atomic_mass_map[1052640000] = 264.117405;
  atomic_mass_map[1062640000] = 264.118929;
  atomic_mass_map[1072640000] = 264.124593;
  atomic_mass_map[1082640000] = 264.128356405;
  atomic_mass_map[1032650000] = 265.116193;
  atomic_mass_map[1042650000] = 265.116683;
  atomic_mass_map[1052650000] = 265.118608;
  atomic_mass_map[1062650000] = 265.12109;
  atomic_mass_map[1072650000] = 265.124977;
  atomic_mass_map[1082650000] = 265.129791799;
  atomic_mass_map[1092650000] = 265.135995;
  atomic_mass_map[1032660000] = 266.119831;
  atomic_mass_map[1042660000] = 266.118172;
  atomic_mass_map[1052660000] = 266.121028;
  atomic_mass_map[1062660000] = 266.121973;
  atomic_mass_map[1072660000] = 266.12679;
  atomic_mass_map[1082660000] = 266.130045252;
  atomic_mass_map[1092660000] = 266.137373;
  atomic_mass_map[1042670000] = 267.121787;
  atomic_mass_map[1052670000] = 267.122464;
  atomic_mass_map[1062670000] = 267.124322;
  atomic_mass_map[1072670000] = 267.1275;
  atomic_mass_map[1082670000] = 267.131673;
  atomic_mass_map[1092670000] = 267.137189;
  atomic_mass_map[1102670000] = 267.143726;
  atomic_mass_map[1042680000] = 268.123968;
  atomic_mass_map[1052680000] = 268.125671;
  atomic_mass_map[1062680000] = 268.125392;
  atomic_mass_map[1072680000] = 268.129691;
  atomic_mass_map[1082680000] = 268.131863;
  atomic_mass_map[1092680000] = 268.138649;
  atomic_mass_map[1102680000] = 268.143477;
  atomic_mass_map[1052690000] = 269.127911;
  atomic_mass_map[1062690000] = 269.12857;
  atomic_mass_map[1072690000] = 269.130412;
  atomic_mass_map[1082690000] = 269.133725;
  atomic_mass_map[1092690000] = 269.138884;
  atomic_mass_map[1102690000] = 269.144751021;
  atomic_mass_map[1052700000] = 270.131302;
  atomic_mass_map[1062700000] = 270.130426;
  atomic_mass_map[1072700000] = 270.133362;
  atomic_mass_map[1082700000] = 270.134314;
  atomic_mass_map[1092700000] = 270.140323;
  atomic_mass_map[1102700000] = 270.14458309;
  atomic_mass_map[1062710000] = 271.133932;
  atomic_mass_map[1072710000] = 271.135182;
  atomic_mass_map[1082710000] = 271.137135;
  atomic_mass_map[1092710000] = 271.140742;
  atomic_mass_map[1102710000] = 271.145946;
  atomic_mass_map[1062720000] = 272.13589;
  atomic_mass_map[1072720000] = 272.138261;
  atomic_mass_map[1082720000] = 272.138494;
  atomic_mass_map[1092720000] = 272.143406;
  atomic_mass_map[1102720000] = 272.146018;
  atomic_mass_map[1112720000] = 272.153273;
  atomic_mass_map[1062730000] = 273.13958;
  atomic_mass_map[1072730000] = 273.14024;
  atomic_mass_map[1082730000] = 273.14159;
  atomic_mass_map[1092730000] = 273.14462;
  atomic_mass_map[1102730000] = 273.148531;
  atomic_mass_map[1112730000] = 273.153189;
  atomic_mass_map[1072740000] = 274.143513;
  atomic_mass_map[1082740000] = 274.143303;
  atomic_mass_map[1092740000] = 274.147339;
  atomic_mass_map[1102740000] = 274.149434;
  atomic_mass_map[1112740000] = 274.155249;
  atomic_mass_map[1072750000] = 275.14567;
  atomic_mass_map[1082750000] = 275.146667;
  atomic_mass_map[1092750000] = 275.149039;
  atomic_mass_map[1102750000] = 275.151976;
  atomic_mass_map[1112750000] = 275.155981;
  atomic_mass_map[1082760000] = 276.148455;
  atomic_mass_map[1092760000] = 276.151708;
  atomic_mass_map[1102760000] = 276.153024;
  atomic_mass_map[1112760000] = 276.158333;
  atomic_mass_map[1122760000] = 276.16141;
  atomic_mass_map[1082770000] = 277.151899;
  atomic_mass_map[1092770000] = 277.153483;
  atomic_mass_map[1102770000] = 277.155815;
  atomic_mass_map[1112770000] = 277.159247;
  atomic_mass_map[1122770000] = 277.163611;
  atomic_mass_map[1092780000] = 278.156454;
  atomic_mass_map[1102780000] = 278.157146;
  atomic_mass_map[1112780000] = 278.161587;
  atomic_mass_map[1122780000] = 278.164179;
  atomic_mass_map[1132780000] = 278.170574;
  atomic_mass_map[1092790000] = 279.158343;
  atomic_mass_map[1102790000] = 279.160093;
  atomic_mass_map[1112790000] = 279.162937;
  atomic_mass_map[1122790000] = 279.166432;
  atomic_mass_map[1132790000] = 279.17095;
  atomic_mass_map[1102800000] = 280.16159;
  atomic_mass_map[1112800000] = 280.165203;
  atomic_mass_map[1122800000] = 280.167147;
  atomic_mass_map[1132800000] = 280.172991;
  atomic_mass_map[1102810000] = 281.164715;
  atomic_mass_map[1112810000] = 281.166718;
  atomic_mass_map[1122810000] = 281.169641;
  atomic_mass_map[1132810000] = 281.17371;
  atomic_mass_map[1112820000] = 282.169405;
  atomic_mass_map[1122820000] = 282.170668;
  atomic_mass_map[1132820000] = 282.175766;
  atomic_mass_map[1112830000] = 283.170995;
  atomic_mass_map[1122830000] = 283.173362;
  atomic_mass_map[1132830000] = 283.17682;
  atomic_mass_map[1122840000] = 284.174499;
  atomic_mass_map[1132840000] = 284.178843;
  atomic_mass_map[1142840000] = 284.181344;
  atomic_mass_map[1122850000] = 285.177321;
  atomic_mass_map[1132850000] = 285.180066;
  atomic_mass_map[1142850000] = 285.183579;
  atomic_mass_map[1132860000] = 286.182518;
  atomic_mass_map[1142860000] = 286.184406;
  atomic_mass_map[1132870000] = 287.18384;
  atomic_mass_map[1142870000] = 287.186875;
  atomic_mass_map[1152870000] = 287.190978;
  atomic_mass_map[1142880000] = 288.187916;
  atomic_mass_map[1152880000] = 288.192992;
  atomic_mass_map[1142890000] = 289.190623;
  atomic_mass_map[1152890000] = 289.193953;
  atomic_mass_map[1162890000] = 289.198099;
  atomic_mass_map[1152900000] = 290.196345;
  atomic_mass_map[1162900000] = 290.198818;
  atomic_mass_map[1152910000] = 291.197522;
  atomic_mass_map[1162910000] = 291.201169;
  atomic_mass_map[1172910000] = 291.205906;
  atomic_mass_map[1162920000] = 292.202086;
  atomic_mass_map[1172920000] = 292.207812;
  atomic_mass_map[1162930000] = 293.204691;
  atomic_mass_map[1172930000] = 293.20868;
  atomic_mass_map[1182930000] = 293.213498;
  atomic_mass_map[1172940000] = 294.210974;
  atomic_mass_map[1182940000] = 294.214132;
  atomic_mass_map[1182950000] = 295.216332;
}


void pyne::_insert_abund_map() {
  natural_abund_map[10010000] = 99.9885;
  natural_abund_map[10020000] = 0.0115;
  natural_abund_map[20030000] = 0.000134;
  natural_abund_map[20040000] = 99.999866;
  natural_abund_map[30060000] = 7.59;
  natural_abund_map[30070000] = 92.41;
  natural_abund_map[40090000] = 100.0;
  natural_abund_map[50100000] = 19.9;
  natural_abund_map[50110000] = 80.1;
  natural_abund_map[60120000] = 98.93;
  natural_abund_map[60130000] = 1.07;
  natural_abund_map[70140000] = 99.636;
  natural_abund_map[70150000] = 0.364;
  natural_abund_map[80160000] = 99.757;
  natural_abund_map[80170000] = 0.038;
  natural_abund_map[80180000] = 0.205;
  natural_abund_map[90190000] = 100.0;
  natural_abund_map[100200000] = 90.48;
  natural_abund_map[100210000] = 0.27;
  natural_abund_map[100220000] = 9.25;
  natural_abund_map[110230000] = 100.0;
  natural_abund_map[120240000] = 78.99;
  natural_abund_map[120250000] = 10.00;
  natural_abund_map[120260000] = 11.01;
  natural_abund_map[130270000] = 100.0;
  natural_abund_map[140280000] = 92.223;
  natural_abund_map[140290000] = 4.685;
  natural_abund_map[140300000] = 3.092;
  natural_abund_map[150310000] = 100.0;
  natural_abund_map[160320000] = 94.99;
  natural_abund_map[160330000] = 0.75;
  natural_abund_map[160340000] = 4.25;
  natural_abund_map[160360000] = 0.01;
  natural_abund_map[170350000] = 75.76;
  natural_abund_map[170370000] = 24.24;
  natural_abund_map[180360000] = 0.3336;
  natural_abund_map[180380000] = 0.0629;
  natural_abund_map[180400000] = 99.6035;
  natural_abund_map[190390000] = 93.2581;
  natural_abund_map[190400000] = 0.0117;
  natural_abund_map[190410000] = 6.7302;
  natural_abund_map[200400000] = 96.941;
  natural_abund_map[200420000] = 0.647;
  natural_abund_map[200430000] = 0.135;
  natural_abund_map[200440000] = 2.086;
  natural_abund_map[200460000] = 0.004;
  natural_abund_map[200480000] = 0.187;
  natural_abund_map[210450000] = 100.0;
  natural_abund_map[220460000] = 8.25;
  natural_abund_map[220470000] = 7.44;
  natural_abund_map[220480000] = 73.72;
  natural_abund_map[220490000] = 5.41;
  natural_abund_map[220500000] = 5.18;
  natural_abund_map[230500000] = 0.250;
  natural_abund_map[230510000] = 99.750;
  natural_abund_map[240500000] = 4.345;
  natural_abund_map[240520000] = 83.789;
  natural_abund_map[240530000] = 9.501;
  natural_abund_map[240540000] = 2.365;
  natural_abund_map[250550000] = 100.0;
  natural_abund_map[260540000] = 5.845;
  natural_abund_map[260560000] = 91.754;
  natural_abund_map[260570000] = 2.119;
  natural_abund_map[260580000] = 0.282;
  natural_abund_map[270590000] = 100.0;
  natural_abund_map[280580000] = 68.077;
  natural_abund_map[280600000] = 26.223;
  natural_abund_map[280610000] = 1.1399;
  natural_abund_map[280620000] = 3.6346;
  natural_abund_map[280640000] = 0.9255;
  natural_abund_map[290630000] = 69.15;
  natural_abund_map[290650000] = 30.85;
  natural_abund_map[300640000] = 49.17;
  natural_abund_map[300660000] = 27.73;
  natural_abund_map[300670000] = 4.04;
  natural_abund_map[300680000] = 18.45;
  natural_abund_map[300700000] = 0.61;
  natural_abund_map[310690000] = 60.108;
  natural_abund_map[310710000] = 39.892;
  natural_abund_map[320700000] = 20.57;
  natural_abund_map[320720000] = 27.45;
  natural_abund_map[320730000] = 7.75;
  natural_abund_map[320740000] = 36.50;
  natural_abund_map[320760000] = 7.73;
  natural_abund_map[330750000] = 100.0;
  natural_abund_map[340740000] = 0.89;
  natural_abund_map[340760000] = 9.37;
  natural_abund_map[340770000] = 7.63;
  natural_abund_map[340780000] = 23.77;
  natural_abund_map[340800000] = 49.61;
  natural_abund_map[340820000] = 8.73;
  natural_abund_map[350790000] = 50.69;
  natural_abund_map[350810000] = 49.31;
  natural_abund_map[360780000] = 0.355;
  natural_abund_map[360800000] = 2.286;
  natural_abund_map[360820000] = 11.593;
  natural_abund_map[360830000] = 11.500;
  natural_abund_map[360840000] = 56.987;
  natural_abund_map[360860000] = 17.279;
  natural_abund_map[370850000] = 72.17;
  natural_abund_map[370870000] = 27.83;
  natural_abund_map[380840000] = 0.56;
  natural_abund_map[380860000] = 9.86;
  natural_abund_map[380870000] = 7.00;
  natural_abund_map[380880000] = 82.58;
  natural_abund_map[390890000] = 100.0;
  natural_abund_map[400900000] = 51.45;
  natural_abund_map[400910000] = 11.22;
  natural_abund_map[400920000] = 17.15;
  natural_abund_map[400940000] = 17.38;
  natural_abund_map[400960000] = 2.80;
  natural_abund_map[410930000] = 100.0;
  natural_abund_map[420920000] = 14.53;
  natural_abund_map[420940000] = 9.15;
  natural_abund_map[420950000] = 15.84;
  natural_abund_map[420960000] = 16.67;
  natural_abund_map[420970000] = 9.60;
  natural_abund_map[420980000] = 24.39;
  natural_abund_map[421000000] = 9.82;
  natural_abund_map[440960000] = 5.54;
  natural_abund_map[440980000] = 1.87;
  natural_abund_map[440990000] = 12.76;
  natural_abund_map[441000000] = 12.60;
  natural_abund_map[441010000] = 17.06;
  natural_abund_map[441020000] = 31.55;
  natural_abund_map[441040000] = 18.62;
  natural_abund_map[451030000] = 100.0;
  natural_abund_map[461020000] = 1.02;
  natural_abund_map[461040000] = 11.14;
  natural_abund_map[461050000] = 22.33;
  natural_abund_map[461060000] = 27.33;
  natural_abund_map[461080000] = 26.46;
  natural_abund_map[461100000] = 11.72;
  natural_abund_map[471070000] = 51.839;
  natural_abund_map[471090000] = 48.161;
  natural_abund_map[481060000] = 1.25;
  natural_abund_map[481080000] = 0.89;
  natural_abund_map[481100000] = 12.49;
  natural_abund_map[481110000] = 12.80;
  natural_abund_map[481120000] = 24.13;
  natural_abund_map[481130000] = 12.22;
  natural_abund_map[481140000] = 28.73;
  natural_abund_map[481160000] = 7.49;
  natural_abund_map[491130000] = 4.29;
  natural_abund_map[491150000] = 95.71;
  natural_abund_map[501120000] = 0.97;
  natural_abund_map[501140000] = 0.66;
  natural_abund_map[501150000] = 0.34;
  natural_abund_map[501160000] = 14.54;
  natural_abund_map[501170000] = 7.68;
  natural_abund_map[501180000] = 24.22;
  natural_abund_map[501190000] = 8.59;
  natural_abund_map[501200000] = 32.58;
  natural_abund_map[501220000] = 4.63;
  natural_abund_map[501240000] = 5.79;
  natural_abund_map[511210000] = 57.21;
  natural_abund_map[511230000] = 42.79;
  natural_abund_map[521200000] = 0.09;
  natural_abund_map[521220000] = 2.55;
  natural_abund_map[521230000] = 0.89;
  natural_abund_map[521240000] = 4.74;
  natural_abund_map[521250000] = 7.07;
  natural_abund_map[521260000] = 18.84;
  natural_abund_map[521280000] = 31.74;
  natural_abund_map[521300000] = 34.08;
  natural_abund_map[531270000] = 100.0;
  natural_abund_map[541240000] = 0.0952;
  natural_abund_map[541260000] = 0.0890;
  natural_abund_map[541280000] = 1.9102;
  natural_abund_map[541290000] = 26.4006;
  natural_abund_map[541300000] = 4.0710;
  natural_abund_map[541310000] = 21.2324;
  natural_abund_map[541320000] = 26.9086;
  natural_abund_map[541340000] = 10.4357;
  natural_abund_map[541360000] = 8.8573;
  natural_abund_map[551330000] = 100.0;
  natural_abund_map[561300000] = 0.106;
  natural_abund_map[561320000] = 0.101;
  natural_abund_map[561340000] = 2.417;
  natural_abund_map[561350000] = 6.592;
  natural_abund_map[561360000] = 7.854;
  natural_abund_map[561370000] = 11.232;
  natural_abund_map[561380000] = 71.698;
  natural_abund_map[571380000] = 0.08881;
  natural_abund_map[571390000] = 99.91119;
  natural_abund_map[581360000] = 0.185;
  natural_abund_map[581380000] = 0.251;
  natural_abund_map[581400000] = 88.450;
  natural_abund_map[581420000] = 11.114;
  natural_abund_map[591410000] = 100.0;
  natural_abund_map[601420000] = 27.152;
  natural_abund_map[601430000] = 12.174;
  natural_abund_map[601440000] = 23.798;
  natural_abund_map[601450000] = 8.293;
  natural_abund_map[601460000] = 17.189;
  natural_abund_map[601480000] = 5.756;
  natural_abund_map[601500000] = 5.638;
  natural_abund_map[621440000] = 3.07;
  natural_abund_map[621470000] = 14.99;
  natural_abund_map[621480000] = 11.24;
  natural_abund_map[621490000] = 13.82;
  natural_abund_map[621500000] = 7.38;
  natural_abund_map[621520000] = 26.75;
  natural_abund_map[621540000] = 22.75;
  natural_abund_map[631510000] = 47.81;
  natural_abund_map[631530000] = 52.19;
  natural_abund_map[641520000] = 0.20;
  natural_abund_map[641540000] = 2.18;
  natural_abund_map[641550000] = 14.80;
  natural_abund_map[641560000] = 20.47;
  natural_abund_map[641570000] = 15.65;
  natural_abund_map[641580000] = 24.84;
  natural_abund_map[641600000] = 21.86;
  natural_abund_map[651590000] = 100.0;
  natural_abund_map[661560000] = 0.056;
  natural_abund_map[661580000] = 0.095;
  natural_abund_map[661600000] = 2.329;
  natural_abund_map[661610000] = 18.889;
  natural_abund_map[661620000] = 25.475;
  natural_abund_map[661630000] = 24.896;
  natural_abund_map[661640000] = 28.260;
  natural_abund_map[671650000] = 100.0;
  natural_abund_map[681620000] = 0.139;
  natural_abund_map[681640000] = 1.601;
  natural_abund_map[681660000] = 33.503;
  natural_abund_map[681670000] = 22.869;
  natural_abund_map[681680000] = 26.978;
  natural_abund_map[681700000] = 14.910;
  natural_abund_map[691690000] = 100.0;
  natural_abund_map[701680000] = 0.123;
  natural_abund_map[701700000] = 2.982;
  natural_abund_map[701710000] = 14.09;
  natural_abund_map[701720000] = 21.68;
  natural_abund_map[701730000] = 16.103;
  natural_abund_map[701740000] = 32.026;
  natural_abund_map[701760000] = 12.996;
  natural_abund_map[711750000] = 97.401;
  natural_abund_map[711760000] = 2.599;
  natural_abund_map[721740000] = 0.16;
  natural_abund_map[721760000] = 5.26;
  natural_abund_map[721770000] = 18.60;
  natural_abund_map[721780000] = 27.28;
  natural_abund_map[721790000] = 13.62;
  natural_abund_map[721800000] = 35.08;
  natural_abund_map[731800000] = 0.01201;
  natural_abund_map[731810000] = 99.98799;
  natural_abund_map[741800000] = 0.12;
  natural_abund_map[741820000] = 26.50;
  natural_abund_map[741830000] = 14.31;
  natural_abund_map[741840000] = 30.64;
  natural_abund_map[741860000] = 28.43;
  natural_abund_map[751850000] = 37.40;
  natural_abund_map[751870000] = 62.60;
  natural_abund_map[761840000] = 0.02;
  natural_abund_map[761860000] = 1.59;
  natural_abund_map[761870000] = 1.96;
  natural_abund_map[761880000] = 13.24;
  natural_abund_map[761890000] = 16.15;
  natural_abund_map[761900000] = 26.26;
  natural_abund_map[761920000] = 40.78;
  natural_abund_map[771910000] = 37.3;
  natural_abund_map[771930000] = 62.7;
  natural_abund_map[781900000] = 0.012;
  natural_abund_map[781920000] = 0.782;
  natural_abund_map[781940000] = 32.86;
  natural_abund_map[781950000] = 33.78;
  natural_abund_map[781960000] = 25.21;
  natural_abund_map[781980000] = 7.356;
  natural_abund_map[791970000] = 100.0;
  natural_abund_map[801960000] = 0.15;
  natural_abund_map[801980000] = 9.97;
  natural_abund_map[801990000] = 16.87;
  natural_abund_map[802000000] = 23.10;
  natural_abund_map[802010000] = 13.18;
  natural_abund_map[802020000] = 29.86;
  natural_abund_map[802040000] = 6.87;
  natural_abund_map[812030000] = 29.52;
  natural_abund_map[812050000] = 70.48;
  natural_abund_map[822040000] = 1.4;
  natural_abund_map[822060000] = 24.1;
  natural_abund_map[822070000] = 22.1;
  natural_abund_map[822080000] = 52.4;
  natural_abund_map[832090000] = 100.0;
  natural_abund_map[902320000] = 100.0;
  natural_abund_map[912310000] = 100.0;
  natural_abund_map[922340000] = 0.0054;
  natural_abund_map[922350000] = 0.7204;
  natural_abund_map[922380000] = 99.2742;
}
//
// end of src/atomic_data.cpp
//


//
// start of src/measure.cpp
//
/*
 * Copyright (c) 2005 Lawrence Livermore National Laboratory under
 * contract number B545069 with the University of Wisconsin - Madison.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>

#ifndef PYNE_IS_AMALGAMATED
#include "measure.h"
#endif

class CartVect {
 private:
  double coords[3];

 public:

  inline CartVect() {}

  inline CartVect(double tx, double ty, double tz) { set(tx, ty, tz); }

  inline CartVect(const CartVect& other) { set(other.coords); }

  inline void set(double tx, double ty, double tz)
  { coords[0] = tx; coords[1] = ty; coords[2] = tz; }

  inline void set(const double* c)
  { coords[0] = c[0]; coords[1] = c[1]; coords[2] = c[2]; }

  inline double x() const { return coords[0]; }
  inline double y() const { return coords[1]; }
  inline double z() const { return coords[2]; }

  inline CartVect& operator+=(const CartVect& other) {
    coords[0] += other.coords[0];
    coords[1] += other.coords[1];
    coords[2] += other.coords[2];
    return *this;
  }

  inline CartVect& operator-=(const CartVect& other) {
    coords[0] -= other.coords[0];
    coords[1] -= other.coords[1];
    coords[2] -= other.coords[2];
    return *this;
  }

  inline CartVect& operator*=(const CartVect& other);

  inline double lensqr() const;

  inline double len() const;

  inline CartVect operator~() const;


  inline CartVect& operator*=(double a) {
    coords[0] *= a;
    coords[1] *= a;
    coords[2] *= a;
    return *this;
  }

  inline CartVect& operator/=(double a) {
    coords[0] /= a;
    coords[1] /= a;
    coords[2] /= a;
    return *this;
  }


};

inline CartVect operator+(const CartVect& v1, const CartVect& v2) {
  CartVect rval(v1);
  rval += v2;
  return rval;
}

inline CartVect operator-(const CartVect& v1, const CartVect& v2) {
  CartVect rval(v1);
  rval -= v2;
  return rval;
}

inline double operator%(const CartVect& v1, const CartVect& v2) {
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline CartVect operator*(const CartVect& v1, const CartVect& v2) {
  return CartVect(v1.y() * v2.z() - v1.z() * v2.y(),
                  v1.z() * v2.x() - v1.x() * v2.z(),
                  v1.x() * v2.y() - v1.y() * v2.x());
}

inline CartVect CartVect::operator~() const {
  double invlen = 1.0 / len();
  return CartVect(invlen * x(), invlen * y(), invlen * z());
}

inline CartVect& CartVect::operator*=(const CartVect& other)
{ return *this = *this * other; }

inline double CartVect::lensqr() const
{ return *this % *this; }

inline double CartVect::len() const
{ return sqrt(lensqr()); }

inline static double tet_volume(const CartVect& v0,
                                const CartVect& v1,
                                const CartVect& v2,
                                const CartVect& v3) {
  return 1. / 6. * (((v1 - v0) * (v2 - v0)) % (v3 - v0));
}

double edge_length(const double* start_vtx_coords,
                   const double*   end_vtx_coords) {
  const CartVect* start = reinterpret_cast<const CartVect*>(start_vtx_coords);
  const CartVect*   end = reinterpret_cast<const CartVect*>(end_vtx_coords);
  return (*start - *end).len();
}

double measure(moab::EntityType type,
               int num_vertices,
               const double* vertex_coordinates) {
  const CartVect* coords = reinterpret_cast<const CartVect*>(vertex_coordinates);
  switch (type) {
    case moab::MBEDGE:
      return (coords[0] - coords[1]).len();
    case moab::MBTRI:
      return 0.5 * ((coords[1] - coords[0]) * (coords[2] - coords[0])).len();
    case moab::MBQUAD:
      num_vertices = 4;
    case moab::MBPOLYGON: {
      CartVect mid(0, 0, 0);
      for (int i = 0; i < num_vertices; ++i)
        mid += coords[i];
      mid /= num_vertices;

      double sum = 0.0;
      for (int i = 0; i < num_vertices; ++i) {
        int j = (i + 1) % num_vertices;
        sum += ((mid - coords[i]) * (mid - coords[j])).len();
      }
      return 0.5 * sum;
    }
    case moab::MBTET:
      return tet_volume(coords[0], coords[1], coords[2], coords[3]) ;
    case moab::MBPYRAMID:
      return tet_volume(coords[0], coords[1], coords[2], coords[4]) +
             tet_volume(coords[0], coords[2], coords[3], coords[4]) ;
    case moab::MBPRISM:
      return tet_volume(coords[0], coords[1], coords[2], coords[5]) +
             tet_volume(coords[3], coords[5], coords[4], coords[0]) +
             tet_volume(coords[1], coords[4], coords[5], coords[0]) ;
    case moab::MBHEX:
      return tet_volume(coords[0], coords[1], coords[3], coords[4]) +
             tet_volume(coords[7], coords[3], coords[6], coords[4]) +
             tet_volume(coords[4], coords[5], coords[1], coords[6]) +
             tet_volume(coords[1], coords[6], coords[3], coords[4]) +
             tet_volume(coords[2], coords[6], coords[3], coords[1]) ;
    default:
      return 0.0;
  }
}

//
// end of src/measure.cpp
//


//
// start of src/source_sampling.cpp
//
#ifndef PYNE_IS_AMALGAMATED
#include "source_sampling.h"
#endif

// Global sampler instance
static pyne::Sampler* sampler = NULL;
// Global variable for mode range
const int SUBVOXEL_START = 3;

// Fortran API
void pyne::sampling_setup_(int* mode, int* cell_list_size) {
  if (sampler == NULL) {
    std::string filename("source.h5m");
    std::string src_tag_name("source_density");
    std::string e_bounds_file("e_bounds");
    std::vector<double> e_bounds = read_e_bounds(e_bounds_file);
    std::map<std::string, std::string> tag_names;
    tag_names.insert(std::pair<std::string, std::string> ("src_tag_name",
                                                          "source_density"));
    tag_names.insert(std::pair<std::string, std::string> ("bias_tag_name",
                                                          "biased_source_density"));
    tag_names.insert(std::pair<std::string, std::string> ("cell_number_tag_name",
                                                          "cell_number"));
    tag_names.insert(std::pair<std::string, std::string> ("cell_fracs_tag_name",
                                                          "cell_fracs"));
    sampler = new pyne::Sampler(filename, tag_names, e_bounds, *mode);
    *cell_list_size = sampler->get_cell_list_size();
  }
}

void pyne::particle_birth_(double* rands,
                           double* x,
                           double* y,
                           double* z,
                           double* e,
                           double* w,
                           int* cell_list) {
  std::vector<double> rands2(rands, rands + 6);
  pyne::SourceParticle src = sampler->particle_birth(rands2);
  *x = src.get_x();
  *y = src.get_y();
  *z = src.get_z();
  *e = src.get_e();
  *w = src.get_w();
  std::vector<int> c_list = src.get_cell_list();
  int cell_list_size = sampler->get_cell_list_size();
  std::copy(c_list.begin(), c_list.end(), cell_list);
}

std::vector<double> pyne::read_e_bounds(std::string e_bounds_file) {
  std::vector<double> e_bounds;
  std::ifstream inputFile(e_bounds_file.c_str());
  double value;
  if (inputFile) {
    while (inputFile >> value)
      e_bounds.push_back(value);
  } else {
    throw std::runtime_error("File " + e_bounds_file + " not found or no read permission");
  }
  return e_bounds;
}


// C++ API
pyne::Sampler::Sampler(std::string filename,
                       std::string src_tag_name,
                       std::vector<double> e_bounds,
                       bool uniform)
  : filename(filename), src_tag_name(src_tag_name), e_bounds(e_bounds) {
  bias_mode = (uniform) ? UNIFORM : ANALOG;
  setup();
}

pyne::Sampler::Sampler(std::string filename,
                       std::string src_tag_name,
                       std::vector<double> e_bounds,
                       std::string bias_tag_name)
  : filename(filename),
    src_tag_name(src_tag_name),
    e_bounds(e_bounds),
    bias_tag_name(bias_tag_name) {
  bias_mode = USER;
  setup();
}

pyne::Sampler::Sampler(std::string filename,
                       std::map<std::string, std::string> tag_names,
                       std::vector<double> e_bounds,
                       int mode)
  : filename(filename),
    tag_names(tag_names),
    e_bounds(e_bounds),
    mode(mode) {
  // determine the bias_mode
  if (mode == 0) {
    bias_mode = ANALOG;
  } else if (mode == 1) {
    bias_mode = UNIFORM;
  } else if (mode == 2) {
    bias_mode = USER;
  } else if (mode == 3) {
    bias_mode = ANALOG;
  } else if (mode == 4) {
    bias_mode = UNIFORM;
  } else if (mode == 5) {
    bias_mode = USER;
  }

  // find out the src_tag_name and bias_tag_name
  if (tag_names.find("src_tag_name") == tag_names.end()) {
    // src_tag_name not found
    throw std::invalid_argument("src_tag_name not found");
  } else {
    // found src_tag_name
    src_tag_name = tag_names["src_tag_name"];
  }
  if (bias_mode == USER) {
    // bias_tag_name required
    if (tag_names.find("bias_tag_name") == tag_names.end()) {
      // bias_tag_name not found
      throw std::invalid_argument("bias_tag_name not found");
    } else {
      // found bias_tag_name
      bias_tag_name = tag_names["bias_tag_name"];
    }
  }
  setup();
}

pyne::SourceParticle pyne::Sampler::particle_birth(std::vector<double> rands) {
  // select mesh volume and energy group
  // For Unstructured mesh, p_src_num_cells and max_num_cells are set to 1
  // For Cartisian mesh, max_num_cells is obtained via cell_fracs tag
  int pdf_idx = at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx / p_src_num_cells / num_e_groups;
  int c_idx = (pdf_idx / num_e_groups) % p_src_num_cells;
  int e_idx = pdf_idx % num_e_groups;

  // Sample uniformly within the selected mesh volume element and energy
  // group.
  std::vector<double> xyz_rands;
  xyz_rands.push_back(rands[2]);
  xyz_rands.push_back(rands[3]);
  xyz_rands.push_back(rands[4]);
  moab::CartVect pos = sample_xyz(ve_idx, xyz_rands);
  cell_list.resize(0);
  if (ve_type == moab::MBHEX) {
    if (mesh_mode == SUBVOXEL) {
      cell_list.emplace_back(cell_number[ve_idx * max_num_cells + c_idx]);
    } else { // Voxel
      for (int c = 0; c < max_num_cells; c++) {
        cell_list.emplace_back(cell_number[ve_idx * max_num_cells + c]);
      }
    }
  }
  pyne::SourceParticle src = SourceParticle(pos[0], pos[1], pos[2],
                                            sample_e(e_idx, rands[5]), sample_w(pdf_idx), cell_list);
  return src;
}


void pyne::Sampler::setup() {
  moab::ErrorCode rval;
  moab::EntityHandle loaded_file_set;
  // Create MOAB instance
  mesh = new moab::Core();
  rval = mesh->create_meshset(moab::MESHSET_SET, loaded_file_set);
  rval = mesh->load_file(filename.c_str(), &loaded_file_set);
  if (rval != moab::MB_SUCCESS)
    throw std::invalid_argument("Could not load mesh file.");

  // Get mesh volume elements
  moab::Range ves;
  rval = mesh->get_entities_by_dimension(loaded_file_set, 3, ves);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem entities of dimension 3");
  num_ves = ves.size();
  int num_hex, num_tet;
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBHEX, num_hex);
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBTET, num_tet);
  if (num_hex == num_ves) {
    ve_type = moab::MBHEX;
    verts_per_ve = 8;
  } else if (num_tet == num_ves) {
    ve_type = moab::MBTET;
    verts_per_ve = 4;
  } else
    throw std::invalid_argument("Mesh file must contain only tets or hexes.");

  // Assign MeshMode: VOXEL, SUBVOXEL, TET
  // Accept mode: 0, 1, 2, 3, 4, 5
  if (ve_type == moab::MBHEX) {
    if (mode < SUBVOXEL_START) {
      mesh_mode = VOXEL;
    } else {
      mesh_mode = SUBVOXEL;
    }
  } else {
    mesh_mode = TET;
  }

  if (ve_type == moab::MBHEX) {
    // cell_number_tag
    if (tag_names.find("cell_number_tag_name") == tag_names.end()) {
      // cell_number_tag_name not found
      throw std::invalid_argument("cell_number_tag_name not found");
    } else {
      // found cell_number_tag_name
      cell_number_tag_name = tag_names["cell_number_tag_name"];
    }
    // cell_fracs_tag
    if (tag_names.find("cell_fracs_tag_name") == tag_names.end()) {
      // cell_fracs_tag_name not found
      throw std::invalid_argument("cell_fracs_tag_name not found");
    } else {
      // found cell_fracs_tag_name
      cell_fracs_tag_name = tag_names["cell_fracs_tag_name"];
    }
  }

  // Process all the spatial and tag data and create an alias table.
  std::vector<double> volumes(num_ves);
  mesh_geom_data(ves, volumes);
  mesh_tag_data(ves, volumes);
}

void pyne::Sampler::mesh_geom_data(moab::Range ves, std::vector<double>& volumes) {
  // Get connectivity.
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> connect;
  rval = mesh->get_connectivity_by_type(ve_type, connect);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting mesh connectivity.");

  // Grab the coordinates that define 4 connected points within a mesh volume
  // element and setup a data structure to allow uniform sampling with each
  // mesh volume element.
  double coords[verts_per_ve * 3];
  int v;
  for (v = 0; v < num_ves; ++v) {
    rval = mesh->get_coords(&connect[verts_per_ve * v], verts_per_ve, &coords[0]);
    if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem vertex coordinates.");
    volumes[v] = measure(ve_type, verts_per_ve, &coords[0]);
    if (ve_type == moab::MBHEX) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[9], coords[10], coords[11]);
      moab::CartVect z(coords[12], coords[13], coords[14]);
      edge_points ep = {o, x - o, y - o, z - o};
      all_edge_points.push_back(ep);
    } else if (ve_type == moab::MBTET) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[6], coords[7], coords[8]);
      moab::CartVect z(coords[9], coords[10], coords[11]);
      edge_points ep = {o, x - o, y - o, z - o};
      all_edge_points.push_back(ep);
    }
  }
}

void pyne::Sampler::mesh_tag_data(moab::Range ves,
                                  const std::vector<double> volumes) {
  moab::ErrorCode rval;
  moab::Tag src_tag;
  moab::Tag cell_number_tag;
  moab::Tag cell_fracs_tag;
  rval = mesh->tag_get_handle(src_tag_name.c_str(),
                              moab::MB_TAG_VARLEN,
                              moab::MB_TYPE_DOUBLE,
                              src_tag);
  // THIS rval FAILS because we do not know number of energy groups a priori.
  // That's okay. That's what the next line is all about:
  num_e_groups = num_groups(src_tag);

  // Set the default value of max_num_cells to 1, so that the structured mesh
  // and unstructured mesh r2s can use the same form of pdf size description.
  max_num_cells = 1;
  p_src_num_cells = 1;
  int tag_size;
  // set the default value of cell_number to -1, cell_fracs to 1.0.
  cell_number.resize(num_ves, -1);
  cell_fracs.resize(num_ves, 1.0);
  if (ve_type == moab::MBHEX) {
    // Read the cell_number tag and cell_fracs tag
    rval = mesh->tag_get_handle(cell_number_tag_name.c_str(),
                                cell_number_tag);
    rval = mesh->tag_get_handle(cell_fracs_tag_name.c_str(),
                                cell_fracs_tag);
    has_cell_fracs = check_cell_fracs(cell_fracs_tag);
    max_num_cells = get_max_num_cells(cell_fracs_tag);
    if (mesh_mode == SUBVOXEL) {
      p_src_num_cells = max_num_cells;
      num_e_groups /= p_src_num_cells;
      // cell_fracs must exist in SUBVOXEL mode
      if (has_cell_fracs == false) {
        throw std::runtime_error("No cell_fracs tag found in sub-voxel R2S. Wrong source file used.");
      }
      cell_fracs.resize(num_ves * p_src_num_cells);
      rval = mesh->tag_get_data(cell_fracs_tag, ves, &cell_fracs[0]);
    }
    // check and set cell_number
    if (has_cell_fracs) {
      cell_number.resize(num_ves * max_num_cells);
      rval = mesh->tag_get_data(cell_number_tag, ves, &cell_number[0]);
    }
  }
  std::cout << " comment. max_num_cells=" << max_num_cells << std::endl;
  std::vector<double> pdf(num_ves * num_e_groups * p_src_num_cells);
  rval = mesh->tag_get_data(src_tag, ves, &pdf[0]);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting source tag data.");

  // Multiply the source densities by the VE volumes
  int v, c, e;
  for (v = 0; v < num_ves; ++v) {
    for (c = 0; c < p_src_num_cells; ++c) {
      for (e = 0; e < num_e_groups; ++e) {
        pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] *=
            volumes[v] * cell_fracs[v * p_src_num_cells + c];
      }
    }
  }
  normalize_pdf(pdf);

  // Setup alias table based off PDF or biased PDF
  if (bias_mode == ANALOG) {
    at = new AliasTable(pdf);
  } else {
    std::vector<double> bias_pdf = read_bias_pdf(ves, volumes, pdf);
    normalize_pdf(bias_pdf);
    //  Create alias table based off biased pdf and calculate birth weights.
    biased_weights.resize(num_ves * p_src_num_cells * num_e_groups);
    for (int i = 0; i < biased_weights.size(); ++i) {
      biased_weights[i] = pdf[i] / bias_pdf[i];
    }
    at = new AliasTable(bias_pdf);
  }
}

std::vector<double> pyne::Sampler::read_bias_pdf(moab::Range ves,
                                                 std::vector<double> volumes,
                                                 std::vector<double> pdf) {
  std::vector<double> bias_pdf(num_ves * p_src_num_cells * num_e_groups);
  int v, c, e;
  moab::ErrorCode rval;
  if (bias_mode == UNIFORM) {
    // Sub-voxel Uniform sampling: uniform in space, analog in energy. Biased PDF is
    // found by normalizing the total photon emission density to 1 in each
    // mesh volume element and multiplying by the volume of the element.
    double q_in_group;
    for (v = 0; v < num_ves; ++v) {
      for (c = 0; c < p_src_num_cells; ++c) {
        q_in_group = 0.0;
        for (e = 0; e < num_e_groups; ++e) {
          q_in_group += pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e];
        }

        if (q_in_group > 0) {
          for (e = 0; e < num_e_groups; ++e) {
            bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] =
                volumes[v] * cell_fracs[v * p_src_num_cells + c] *
                pdf[v * p_src_num_cells * num_e_groups +
                      c * num_e_groups + e] / q_in_group;
          }
        } else {
          for (e = 0; e < num_e_groups; ++e) {
            bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] = 0.0;
          }
        }

      }
    }
  } else if (bias_mode == USER) {
    // Get the biased PDF from the mesh
    moab::Tag bias_tag;
    rval = mesh->tag_get_handle(bias_tag_name.c_str(),
                                moab::MB_TAG_VARLEN,
                                moab::MB_TYPE_DOUBLE,
                                bias_tag);
    num_bias_groups = num_groups(bias_tag);
    if (num_bias_groups == num_e_groups * p_src_num_cells) {
      // Spatial, cell and energy biasing. The supplied bias PDF values are
      // applied to each specific energy group and sub-voxels in a mesh
      // volume element.
      rval = mesh->tag_get_data(bias_tag, ves, &bias_pdf[0]);
      if (rval != moab::MB_SUCCESS)
        throw std::runtime_error("Problem getting bias tag data.");
      for (v = 0; v < num_ves; ++v) {
        for (c = 0; c < p_src_num_cells; c++) {
          for (e = 0; e < num_e_groups; ++e)
            bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] *=
                volumes[v] * cell_fracs[v * p_src_num_cells + c];
        }
      }
    } else if (num_bias_groups == 1) {
      // Spatial biasing only: the supplied bias PDF values are applied
      // to all energy groups within a mesh volume element, which are
      // sampled in analog.
      std::vector<double> spatial_pdf(num_ves);
      rval = mesh->tag_get_data(bias_tag, ves, &spatial_pdf[0]);
      if (rval != moab::MB_SUCCESS)
        throw std::runtime_error("Problem getting bias tag data.");
      double q_in_group;
      for (v = 0; v < num_ves; ++v) {
        q_in_group = 0;
        for (c = 0; c < p_src_num_cells; ++c) {
          for (e = 0; e < num_e_groups; ++e) {
            q_in_group += pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e];
          }
        }
        if (q_in_group > 0) {
          for (c = 0; c < p_src_num_cells; ++c) {
            for (e = 0; e < num_e_groups; ++e) {
              bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] =
                  spatial_pdf[v] * volumes[v] * cell_fracs[v * p_src_num_cells + c] *
                  pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] /
                  q_in_group;
            }
          }
        } else {
          for (c = 0; c < p_src_num_cells; ++c)
            for (e = 0; e < num_e_groups; ++e) {
              bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] =  0;
            }
        }
      }
    } else if (num_bias_groups == num_e_groups) {
      // Voxel and energy biasing. Apply the energy bias to all the sub-voxel in the voxel
      std::vector<double> spa_erg_pdf(num_ves * num_e_groups);
      rval = mesh->tag_get_data(bias_tag, ves, &spa_erg_pdf[0]);
      if (rval != moab::MB_SUCCESS)
        throw std::runtime_error("Problem getting bias tag data.");
      double q_in_group;
      for (v = 0; v < num_ves; ++v) {
        for (e = 0; e < num_e_groups; ++e) {
          q_in_group = 0.0;
          for (c = 0; c < p_src_num_cells; ++c) {
            q_in_group += pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e];
          }
          if (q_in_group > 0) {
            for (c = 0; c < p_src_num_cells; ++c) {
              bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] =
                  spa_erg_pdf[v * num_e_groups + e] * volumes[v] * cell_fracs[v * p_src_num_cells + c] *
                  pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] / q_in_group;
            }
          } else {
            for (c = 0; c < p_src_num_cells; ++c) {
              bias_pdf[v * p_src_num_cells * num_e_groups + c * num_e_groups + e] = 0.0;
            }
          }
        }
      }
    } else {
      throw std::length_error("Length of bias tag must equal length of the"
                              "  p_src_num_cells*num_e_group, num_e_groups, or 1.");
    }
  }
  double q_in_all = 0.0;
  for (int i = 0; i < bias_pdf.size(); i++)
    q_in_all += bias_pdf[i];
  if (q_in_all <= 0.0) {
    throw std::runtime_error("Bias data are ALL ZERO!");
  }
  return bias_pdf;
}


moab::CartVect pyne::Sampler::sample_xyz(int ve_idx, std::vector<double> rands) {
  double s = rands[0];
  double t = rands[1];
  double u = rands[2];

  // Transform s, t, u to uniformly sample a tetrahedron. See:
  // C. Rocchini and P. Cignoni, “Generating Random Points in a Tetrahedron,”
  //  Journal of Graphics Tools, 5, 200–202 (2001).
  if (ve_type == moab::MBTET) {
    if (s + t > 1) {
      s = 1.0 - s;
      t = 1.0 - t;
    }
    if (s + t + u > 1) {
      if (t + u > 1) {
        double old_t = t;
        t = 1.0 - u;
        u = 1.0 - s - old_t;
      } else if (t + u <= 1) {
        double old_s = s;
        s = 1.0 - t - u;
        u = old_s + t + u - 1;
      }
    }
  }

  return s * all_edge_points[ve_idx].x_vec + \
         t * all_edge_points[ve_idx].y_vec + \
         u * all_edge_points[ve_idx].z_vec + \
         all_edge_points[ve_idx].o_point;
}

double pyne::Sampler::sample_e(int e_idx, double rand) {
  double e_min = e_bounds[e_idx];
  double e_max = e_bounds[e_idx + 1];
  return rand * (e_max - e_min) + e_min;
}

double pyne::Sampler::sample_w(int pdf_idx) {
  return (bias_mode == ANALOG) ? 1.0 : biased_weights[pdf_idx];
}

void pyne::Sampler::normalize_pdf(std::vector<double>& pdf) {
  double sum = 0;
  for (int i = 0; i < pdf.size(); ++i)
    sum += pdf[i];
  for (int i = 0; i < pdf.size(); ++i)
    pdf[i] /= sum;
}

int pyne::Sampler::num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mesh->tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting tag size.");
  return tag_size / sizeof(double);
}

bool pyne::Sampler::check_cell_fracs(moab::Tag cell_fracs_tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mesh->tag_get_bytes(cell_fracs_tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS) {
    std::cout << "Warning: Old version source file used. No cell_number and cell_fracs tag found!" << std::endl;
    std::cout << "Default cell_number [-1] and cell_fracs [1.0] will be used." << std::endl;
    has_cell_fracs = false;
  } else {
    has_cell_fracs = true;
  }
  return has_cell_fracs;
}

int pyne::Sampler::get_max_num_cells(moab::Tag cell_fracs_tag) {
  moab::ErrorCode rval;
  int tag_size;
  if (has_cell_fracs) {
    rval = mesh->tag_get_bytes(cell_fracs_tag, *(&tag_size));
    return tag_size / sizeof(double);
  } else {
    return 1;
  }
}

int pyne::Sampler::get_cell_list_size() {
  // cell_list_size should be:
  // 0: for unstructured mesh
  // 1: for sub-voxel R2S
  // max_num_cells: for voxel R2S
  if (ve_type == moab::MBTET or has_cell_fracs == false) {
    return 0;
  } else {
    if (mesh_mode == VOXEL) {
      return max_num_cells;
    } else {
      return 1;
    }
  }
}

// Random-number sampling using the Walker-Vose alias method,
// Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
// M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
// A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)

pyne::AliasTable::AliasTable(std::vector<double> p) {
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  std::vector<double> small(n);
  std::vector<double> large(n);
  int i, a, g;

  for (i = 0; i < n; ++i)
    p[i] *= n;

  // Set separate index lists for small and large probabilities:
  int n_s = 0;
  int n_l = 0;
  for (i = n - 1; i >= 0; --i) {
    // at variance from Schwarz, we revert the index order
    if (p[i] < 1)
      small[n_s++] = i;
    else
      large[n_l++] = i;
  }

  // Work through index lists
  while (n_s && n_l) {
    a = small[--n_s]; // Schwarz's l
    g = large[--n_l]; // Schwarz's g
    prob[a] = p[a];
    alias[a] = g;
    p[g] = p[g] + p[a] - 1;
    if (p[g] < 1)
      small[n_s++] = g;
    else
      large[n_l++] = g;
  }

  while (n_l)
    prob[large[--n_l]] = 1;

  while (n_s)
    // can only happen through numeric instability
    prob[small[--n_s] ] = 1;
}

int pyne::AliasTable::sample_pdf(double rand1, double rand2) {
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}

pyne::SourceParticle::SourceParticle() {
  x = -1.0;
  y = -1.0;
  z = -1.0;
  e = -1.0;
  w = -1.0;
}

pyne::SourceParticle::SourceParticle(double _x, double _y, double _z,
                                     double _e, double _w, std::vector<int> _cell_list)
  : x(_x), y(_y), z(_z),
    e(_e), w(_w), cell_list(_cell_list) {}

pyne::SourceParticle::~SourceParticle() {};
//
// end of src/source_sampling.cpp
//


