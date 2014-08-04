#ifndef DAGMC_FLUKA_FLUKA_FUNCS_H
#define DAGMC_FLUKA_FLUKA_FUNCS_H

#include <iostream>    // cout, endl
#include <stdlib.h>
#include <string>      // string, cout
#include <vector>
#include <set>
#include <utility>     // makepair

#include "moab/Types.hpp"
#include "MBInterface.hpp"
#include "MBCartVect.hpp"
#include "DagMC.hpp"
#include "pyne/pyne.h"


/// Called from mainFludag when only one argument is given to the program.
//  This function writes out a simple numerical material assignment to the named argument file
//  Example usage:  mainFludag dagmc.html
//  Outputs
//           mat.inp  contains MATERIAL and ASSIGNMAt records for the input geometry.
//                    The MATERIAL is gotten by parsing the Cubit volume name on underscores.  
//                    The string after "M_" is considered to be the material for that volume.
//                    There are no MATERIAL cards for the materials in the FLUKA_mat_set list
//                    For the remaining materials, there is one MATERIAL card apiece (no dups)
//                    User-named (not predefined) materials are TRUNCATED to 8 chars.
//                    User-named material id's start at 25 and increment by 1 for each MATERIAL card
//           index-id.txt  Map of FluDAG volume index vs Cubit volume ids, for info only.
//  Note that a preprocessing step to this call sets up the the DAG object that contains 
//  all the geometry information contained in dagmc.html.  
//  the name of the (currently hardcoded) output file is "mat.inp"
//  The graveyard is assumed to be the last region.
//  Overall function call; calls other fludag functions
void fludag_write(std::string matfile, std::string lfname);

// Gets the number of volumes via MOAB entities and DAGMC calls;
// parses properties
int fludag_setup(std::map<int, std::string>& map_name);
// Make a string from the groupname
std::string makeMaterialName (int index);

// Load the PyNE material objects in the named file
void pyne_get_materials(std::string mat_file, std::list<pyne::Material>& pyne_list);

// Extract PyNE nucids from a known list of elements
std::set<int> make_exception_set();
/*
 * Write the material assignment for each volume to an output stream
 */
void fludagwrite_assignma(std::ostringstream& ostr, int num_vols, 
                                           std::list<pyne::Material> pyne_list, 
					   std::map<int, std::string>& map_name);
/*
 * Write material cards
 */
// std::list<pyne::Material> fludag_write_material(std::ostringstream& ostr, 
void fludag_write_material(std::ostringstream& ostr, 
                                                  int& last_id,
                                                  std::set<int> exception_set,
                                                  std::list<pyne::Material> pyne_list);
/*
 * Write compound cards
 */
void fludag_write_compound(std::ostringstream& cstr, int& last_id, 
                            pyne::Material& material);
/*
 *  Global mapping from z-a nucid to fluka name.
 *  WAS filled by readElements()
 */
std::map<int, std::string> fluka_name_map;
/*
 * Convenience function to make a different fluka_name; may not be needed
 */
std::string modify_fluka_name(std::string& origName);
/*
 * Creates a Fluka MATERIAL card with a fake density
 */
void define_material(std::ostringstream &cstr, int &last_id, std::string dname);
/* 
 * Return a vector of structs, one for each line in the file
 */
std::map<int, std::string> readElements(std::ifstream& fin);
/*
 * Convenience function
 */
void print_material( pyne::Material test_mat, std::string xtraTitle=" ");

// Defined but not called
void slow_check(double pos[3], const double dir[3], int &oldReg);

/*
 * Prepare a descriptive string that creates the properties of the volume whose index is index
 */
  std::string mat_property_string (int index, std::vector<std::string> &properties);

// This function is defined but never called
void dagmc_version_(double* dagmcVersion);


///////////////////////////////////////////////////////////////////
// Data for Material and Compound Cards
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//
// Interface for Dag Wrappers
//
///////////////////////////////////////////////////////////////////

#define f_idnr idnrwr_
#define g_step g1wr_
#define f_g1rt g1rtwr_
#define inihwr inihwr_
#define jomiwr jomiwr_
#define f_lookdb lkdbwr_
#define lkfxwr lkfxwr_
#define lkmgwr lkmgwr_
#define f_look lkwr_
#define fldwr fldwr_
#define flgfwr flgfwr_
#define f_normal nrmlwr_
#define rgrpwr rgrpwr_
#define isvhwr isvhwr_
#define rg2nwr rg2nwr_

////////////////////////////////////////////////////////////////////
// Start of functions that could be called from Fortran 
//  as such they cannot be overloaded and should be declared 
//  external so they will be compiled without mangling
////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif

// WrapFlgfwr.cc stubs this
void flgfwr(int& flkflg);

// The FLUKA internal function is used.
int f_idnr(const int & nreg, const int & mlat);

// The function is defined in fluka_funcs.cpp.  It calls g_fire.
void  g_step(double& pSx, double& pSy, double& pSz, double* pV,
                      int& oldReg, const int& oldLttc, double& propStep,
                      int& nascFlag, double& retStep, int& newReg,
	              double& saf, int& newLttc, int& LttcFlag,
                      double* sLt, int* jrLt);

// Stub function
void f_g1rt(void);

// WrapInit.cc - has been deleted, the function is now
// defined in fluka_funcs.cpp
void jomiwr(int & nge, const int& lin, const int& lou,
                       int& flukaReg);

// WrapLookDB.cc has been deleted, the function is now
// defined in fluka_funcs.cpp.  It sets some of its return values
void f_lookdb(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
               	       int& newReg, int& flagErr, int& newLttc);

// WrapLookFX
// Stubbed in WrapLookFX.cc and linked in.
void lkfxwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
                       int& newReg, int& flagErr, int& newLttc);
	    
// WrapLookMG.cc stubs this function and is linked in.
void lkmgwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
		       int& flagErr, int& newReg, int& newLttc);
	    
// WrapLookZ has been deleted.  This function is defined in
// fluka_funcs.cc.  It is called by look(..).
void f_look(double& pSx, double& pSy, double& pSz,
                     double* pV, const int& oldReg, const int& oldLttc,
	             int& newReg, int& flagErr, int& newLttc);

// WrapMag.cc stubs this function and is linked in
void fldwr(const double& pX, const double& pY, const double& pZ,
                       double& cosBx, double& cosBy, double& cosBz, 
                       double& Bmag, int& reg, int& idiscflag);
	    
// WrapNorml.cc has been removed.  This function is defined in fluka_funcs.
//  It is called by normal(..).
void f_normal(double& pSx, double& pSy, double& pSz,
                       double& pVx, double& pVy, double& pVz,
	               double* norml, const int& oldReg, 
	               const int& newReg, int& flagErr);

// WrapReg
void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
                       int* indMother, int* repMother, int& depthFluka);

// WrapSavHist.cc has been removed.  This function has been commented
// out of an implementation in fluka_funcs.cpp, so the internal FLUKA 
// version, if any, is called.
int isvhwr(const int& fCheck, const int& intHist);

// WrapReg2name.cc defines this and is linked in.  It calls
// region2name, which is defined in fluka_funcs
// ToDo:  Remove WrapReg2name.cc and implement rg2nwr(..) in fluka_funcs.cpp
void rg2nwr(const int& mreg, char* Vname);

#ifdef __cplusplus
} // extern "C"
#endif
/////////////////////////  End of Fortran-called Functions /////////////////////

////////////////////////////////////////////////////////////////////////////////
// Helper functions for the wrapper functions declared above
////////////////////////////////////////////////////////////////////////////////

// Wrapper for f_look clarifying which arguments are used.
int look( double& posx, double& posy, double& posz, double* dir, int& region);

// Wrapper for f_normal clarifying which arguments are used
int  normal (double& posx, double& posy, double& posz, double *norml, int& regionForSense);

// Protoype for boundary test funct
int boundary_test(MBEntityHandle vol, double xyz[3], double uvw[3]);

// protoype for mat props
std::string mat_property_string (int index, std::vector<std::string> &properties);

// Defined in fluka_funcs, called by rg2nwr
void region2name(int volindex, char * vname );


/**
 * g_fire is called by fludag's implementation of g_step.  It calls DAG->ray_fire(...).
 * oldRegion region of start point
 * point     start point
 * dir       direction vector
 * propStep
 * retStep   set to distance to next surface
 * newRegion region ofter step
 */
void g_fire(int& oldRegion, double point[], double dir[], 
	      double &propStep, double& retStep, double &safety,
	      int& newRegion);

#endif /* DAGMC_FLUKA_FLUKA_FUNCS_H */
