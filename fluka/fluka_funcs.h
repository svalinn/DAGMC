#ifndef DAGMC_FLUKA_FLUKA_FUNCS_H
#define DAGMC_FLUKA_FLUKA_FUNCS_H

#include <iostream>    // cout, endl
#include <stdlib.h>
#include <string>      // string, cout
#include <vector>
#include <set>
#include <utility>     // makepair

#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "DagMC.hpp"
#include "../pyne/pyne.h"
#include "../uwuw/uwuw.hpp"

// defines for the flkstk common block structure
#define STACK_SIZE 40001 // because the fortran array goes from [0:40000]
#define MKBMX1 11
#define MKBMX2 11

//
#define MULBOU_SIZE 2001 // because the fortran array goes from [0:2000]

// flkstk common block
extern "C" {
  extern struct {
    // all the doubles
    double xflk[STACK_SIZE];
    double yflk[STACK_SIZE];
    double zflk[STACK_SIZE];
    double txflk[STACK_SIZE];
    double tyflk[STACK_SIZE];
    double tzflk[STACK_SIZE];
    double txpol[STACK_SIZE];
    double typol[STACK_SIZE];
    double tzpol[STACK_SIZE];
    double txnor[STACK_SIZE];
    double tynor[STACK_SIZE];
    double tznor[STACK_SIZE];
    double wtflk[STACK_SIZE];
    double pmoflk[STACK_SIZE];
    double tkeflk[STACK_SIZE];
    double dfnear[STACK_SIZE];
    double agestk[STACK_SIZE];
    double aknshr[STACK_SIZE];
    double raddly[STACK_SIZE];
    double cmpath[STACK_SIZE];
    double frcphn[STACK_SIZE];
    double dchflk[STACK_SIZE];
    // spared doubles
    double sparek[STACK_SIZE][MKBMX1]; // fortran arrays other way round
    // ints
    int ispark[STACK_SIZE][MKBMX2];
    int iloflk[STACK_SIZE];
    int igroup[STACK_SIZE];
    int loflk[STACK_SIZE];
    int louse[STACK_SIZE];
    int nrgflk[STACK_SIZE];
    int nlattc[STACK_SIZE];
    int nhspnt[STACK_SIZE];
    int nevent[STACK_SIZE];
    int numpar[STACK_SIZE];
    int lraddc[STACK_SIZE];
    int lfrphn[STACK_SIZE];
    int lchflk[STACK_SIZE];
    int nparma;
    int nstmax;
    int npflka;
    int nstaol;
    int igroun;
  } flkstk_;
}

// needed for the lsense flag
// mulbou common block
extern "C" {
  extern struct {
    // all the doubles
    double xold;
    double yold;
    double zold;
    double xmiddl;
    double ymiddl;
    double zmiddl;
    double umiddl;
    double vmiddl;
    double wmiddl;
    double pstep1;
    double pstep2;
    double uold;
    double vold;
    double wold;
    double umag;
    double vmag;
    double wmag;
    double unorml;
    double vnorml;
    double wnorml;
    double usense;
    double vsense;
    double wsense;
    double xnorml;
    double ynorml;
    double znorml;
    double tsense;
    double ddsens;
    double dsmall;
    double tslttc[MULBOU_SIZE];
    // integer vars
    int multtc[MULBOU_SIZE];
    int nssens;
    int nulttc;
    int iplgnl;
    int nrgbef;
    int nrgaft;
    // logical vars
    int llda;
    int lagain;
    int lstnew;
    int lartef;
    int lnorml;
    int lsense;
    int lmgnor;
    int lsnsct;
    int lplgnl;
    int lnwghs;
    int lmagea;
    int lmgnmv;
    int lbndrx;
  } mulbou_;
}

// struct to hold particle state
struct particle_state {
  moab::DagMC::RayHistory history;
  bool on_boundary;
  double old_direction[3];
  double old_position[3];
  moab::EntityHandle next_surface; // the next suface the ray will hit
  moab::EntityHandle prev_surface; // the last value of next surface
  moab::EntityHandle PrevRegion; // the integer region that the particle was in previously
  int stack_count;
};

// Interface for Dag Wrappers
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
#define flabrt flabrt_

/**** Start of directly called Fluka functions ****/

#ifdef __cplusplus
extern "C" {
#endif

  void flabrt(const char* function, const char* message, int function_strlen, int message_strlen);

  /**
   * \brief does nothing
   *
   * \param[in/out] fklflg indicates something is on/off
   */
  void flgfwr(int& flkflg);

  /**
   * \brief Sets if we would like to use dnear, internally the
   *        safety variable, doesnt use the inputs to set dnear
   *
   * \param[in] nreg the number of regions in the problem
   * \param[in] mlat the number oflattices in the problem
   */
  int f_idnr(const int & nreg, const int & mlat);

  /**
   * \brief This is a the main Fluka tracking call made to FluDAG,
   *        this gives access to the variables passed by argument
   *        this function is expected to return the approved particle step
   *        distance (either limited by proposed step or the geometry, and
   *        the geometry ID number after the step
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] *pV the particle direction vector (double[3])
   * \param[in] *oldReg the current region number of the particle
   * \param[in] oldLttc the integer ID of the lattice of the particle
   * \param[in] propStep, the proposed step length as dictated by the
   *             physics of the current collision
   * \param[in] nascFlag (unknown)
   * \param[out] retStep, the approved step length either equal to propStep
   *             or limited by geometry
   * \param[out] newReg the region number after the particle step
   * \param[out] saf, the new approved safety (distance to the nearest surface)
   * \param[out] newLttc, the new lattice ID (ignored)
   * \param[out] LttcFlag unknown (ignored)
   * \param[out] sLt, unknown (ignored)
   * \param[out] jrLt, unknown (ignored)
   */
  void  g_step(double& pSx, double& pSy, double& pSz, double* pV,
               int& oldReg, const int& oldLttc, double& propStep,
               int& nascFlag, double& retStep, int& newReg,
               double& saf, int& newLttc, int& LttcFlag,
               double* sLt, int* jrLt);

  /**
   * \brief Called by Fluka directly, used to indicate when a history has been
   *         terminated. We use this to reset the accumulated particle state
   */
  void f_g1rt(void);

  /**
   * \brief Sets the number of volumes in the problem, this should include the
   *         implicit compliment, the other arguments are unused
   *
   * \param[in] nge unknown (unused)
   * \param[in] lin logical input unit number for fortran
   * \param[in] lou logical output unit number for fortran
   * \param[out] flukaReg the number of regions in the problem
   */
  void jomiwr(int & nge, const int& lin, const int& lou,
              int& flukaReg);

  /**
   * \brief In our case this functions sets newReg and newLttc to 0, and flagErr to -1
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] *pV the particle direction vector (double[3])
   * \param[in] oldReg the current region number of the particle
   * \param[in] oldLttc the integer ID of the lattice of the particle
   * \param[out] newReg the new region number
   * \param[out] flagErr used to indicate an error
   * \param[out] newLttc used to indicate the new lattice id number
   */
  void f_lookdb(double& pSx, double& pSy, double& pSz,
                double* pV, const int& oldReg, const int& oldLttc,
                int& newReg, int& flagErr, int& newLttc);
  /**
   * \brief Does nothing.
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] *pV the particle direction vector (double[3])
   * \param[in] oldReg the current region number of the particle
   * \param[in] oldLttc the integer ID of the lattice of the particle
   * \param[out] newReg the new region number
   * \param[out] flagErr used to indicate an error
   * \param[out] newLttc used to indicate the new lattice id number
   */
  void lkfxwr(double& pSx, double& pSy, double& pSz,
              double* pV, const int& oldReg, int& oldLttc,
              int& newReg, int& flagErr, int& newLttc);

  /**
   * \brief Determines what region number the current particle is
   *        in, bears a strong similarity to lkwr, but this routine
   *        is only called when magnetic field tracking is on
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] *pV the particle direction vector (double[3])
   * \param[in] oldReg the current region number of the particle
   * \param[in] oldLttc the integer ID of the lattice of the particle
   * \param[out] newReg the new region number
   * \param[out] flagErr used to indicate an error
   * \param[out] newLttc used to indicate the new lattice id number
   */
  void lkmgwr(double& pSx, double& pSy, double& pSz,
              double* pV, const int& oldReg, const int& oldLttc,
              int& flagErr, int& newReg, int& newLttc);

  /**
   * \brief Determines what region number the current particle is
   *        in. Uses the global variable state to assist. Sets newReg
   *        to the integer ID of the region the particle is in
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] *pV the particle direction vector (double[3])
   * \param[in] oldReg the current region number of the particle
   * \param[in] oldLttc the integer ID of the lattice of the particle
   * \param[out] newReg the new region number
   * \param[out] flagErr used to indicate an error
   * \param[out] newLttc used to indicate the new lattice id number
   */
  void f_look(double& pSx, double& pSy, double& pSz,
              double* pV, const int& oldReg, const int& oldLttc,
              int& newReg, int& flagErr, int& newLttc);

  /**
   * \brief For a given position, determines the magnetic field strength
   *        and direction.
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[out] cosBx magnetic field cosine direction vector in x direction
   * \param[out] cosBy magnetic field cosine direction vector in y direction
   * \param[out] cosBz magnetic field cosine direction vector in z direction
   * \param[out] Bmag the magnetic field strength (tesla?)
   * \param[in] reg the region id number
   * \param[?] idiscflag (unknown)
   */
  void fldwr(const double& pX, const double& pY, const double& pZ,
             double& cosBx, double& cosBy, double& cosBz,
             double& Bmag, int& reg, int& idiscflag);

  /**
   * \brief For a given position and particle direction, determines the normal
   *        of the surface between newReg and oldReg
   *
   * \param[in] pSx the particle x coordinate (cm)
   * \param[in] pSy the particle y coordinate (cm)
   * \param[in] pSz the particle z coordinate (cm)
   * \param[in] pVx the particle direction vector in the x direction (cm)
   * \param[in] pVy the particle direction vector in the y direction (cm)
   * \param[in] pVz the particle direction vector in the x direction (cm)
   * \param[out] *norml the normal vector to the surface, normal should always point towards
                oldReg. i.e. the inward pointing normal

   * \param[in] oldReg the current region number
   * \param[in] newReg the region number the particle will step to next
   * \param[out] flagErr flag some error condition
   */
  void f_normal(double& pSx, double& pSy, double& pSz,
                double& pVx, double& pVy, double& pVz,
                double* norml, const int& oldReg,
                const int& newReg, int& flagErr);

  /**
   * \brief Does nothing in FluDAG.
   *
   */
  void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
              int* indMother, int* repMother, int& depthFluka);

  /**
   * \brief Used to translate fluka region number ino names, in FluDAG this function
   *        calls region2name, and takes the Fluka ID number and turns it into a string
   *
   * \param[in] mreg the Fluka region id number
   * \param[out] Vname the name of the Fluka region with id mreg
   */
  void rg2nwr(const int& mreg, char* Vname);

#ifdef __cplusplus
} // extern "C"
#endif

/**** END OF Fluka called functions ****/

/*** start of testable wrappers ***/

// The testable interfaces to the direct Fluka interface calls

/**
 * \brief The testable interface to allow fludag to interface with the fluka
 *        abort function. When called, passes the messages to flabrt_ which
 *        raises an error call, dumps simulation data and writes messages to file
 *
 * \param[in] function_name, char* denoting the name of the function that raised the error
 * \param[in] message, char* denoting the error message you want to pass
 * \param[in] error_code, int, the error code you are raising from DAGMC
 */
void fludag_abort(const char* function_name, const char* message, int error_code);


double dot_product(moab::EntityHandle surface, double point[3], double direction[3]);

/**
 * \brief The testable interface to f_look, the function which given a position
 *       , direction determines the region number the position is in. Indicates
 *        an error when the particle is nowhere.
 *
 * \param[in] posx the x position of the ray
 * \param[in] posy the y position of the ray
 * \param[in] posz the z position of the ray
 * \param[in] *dir the direction vector of the ray (double[3])
 * \param[out] the region index that the position belong to
 */
int look( double& posx, double& posy, double& posz, double* dir, int& region);

/**
 * \brief g_fire is called by g_step, which is the external interface to FLuka.
 *
 * \param[in] oldRegion region of start point
 * \param[in] point[3] start position
 * \param[in] dir[3] direction vector
 * \param[in] propStep physics proposed step length
 * \param[out] retStep actual returned distance, governed by geometry or physics
 * \param[out] newRegion region after step
 **/
void g_fire(int& oldRegion, double point[], double dir[],
            double &propStep, double& retStep, double &safety,
            int& newRegion);

/**
 * \brief normal is the testable wrapper to f_normal, the external interface to FLuka.
 *        it expects to be only called in certain circumstances, it expects to be only called
 *        when the position is on a boundary, and that boundary belongs to curRegion.
 *
 * \param[in] posx the x coordinate of the ray
 * \param[in] posy the y coordinate of the ray
 * \param[in] posz the z coordinate of the ray
 * \param[out] *norml the returned normal vector
 * \param[int] curRegion the region contained by posx,posy, posz
 **/
int normal (double& posx, double& posy, double& posz, double *norml, int& curRegion);

/**
 * \brief Tests if a particle is on the boundary of the volume vol, returns 1 true or 0 false,
 *        expects to use the global particle state variable state.next_surf, so must be called
 *        after a succesful ray fire, like g_fire
 *
 * \param[in] vol the volume to test
 * \param[in] xyz[3] the position to test
 * \param[in] uvw[3] the direction vector
 * \param[returns] either yes (1) or no (0)
 **/
int boundary_test(moab::EntityHandle vol, double xyz[3], double uvw[3]);

/**
 * \brief turns an integer volume index into its string form with a period appended
 *         i.e turns int 1 into "1."
 * \param[in] volindex any integer
 * \param[out] vname string of the integer form
 */
void region2name(int volindex, std::string &vname );

void print_state(particle_state &state);

/**
 * \brief resets the particle state associated with the particle.
 *
 * \param[in/out] state the particle state structure
 */
void reset_state(particle_state &state);

/**
 * \brief Called from mainFludag when only one argument is given to the program.
 *  This function writes out a simple numerical material assignment to the named argument file
 *  Example usage:  mainFludag dagmc.html
 *  Outputs
 *           mat.inp  contains MATERIAL and ASSIGNMAt records for the input geometry.
 *                    The MATERIAL is gotten by parsing the Cubit volume name on underscores.
 *                    The string after "M_" is considered to be the material for that volume.
 *                    There are no MATERIAL cards for the materials in the FLUKA_mat_set list
 *                    For the remaining materials, there is one MATERIAL card apiece (no dups)
 *                    User-named (not predefined) materials are TRUNCATED to 8 chars.
 *                    User-named material id's start at 25 and increment by 1 for each MATERIAL card
 *           index-id.txt  Map of FluDAG volume index vs Cubit volume ids, for info only.
 *  Note that a preprocessing step to this call sets up the the DAG object that contains
 *  all the geometry information contained in dagmc.html.
 *  the name of the (currently hardcoded) output file is "mat.inp"
 *  The graveyard is assumed to be the last region.
 *  Overall function call; calls other fludag functions
 */
void fludag_write(std::string matfile, std::string lfname);

/**
 * \brief makes a standard set of nuclides which we would like to retain in nucid form
 */
std::set<int> make_exception_set();

/*** end of testable wrappers ***/

/*** start of uwuw interface functions ***/

void fludag_write_ididx(std::string ididx);

/**
 * \brief Writes the importance assignment for each volume to an output stream
 *
 * \param[out] ostr the output stream of which the text is printed
 */
void fludagwrite_importances(std::ostringstream& ostr);

/**
 * \brief Writes the material assignment for each volume to an output stream
 *
 * \param[out] ostr the output stream of which the text is printed
 * \param[in] pyne_map the map of group name vs PyNE material object
 * \param[in] map_name ???
 */
void fludagwrite_assignma(std::ostringstream& ostr,
                          std::map<std::string, pyne::Material> pyne_map);
/**
 * \brief writes all the material compisitions from the map of materials to output stream
 *
 * \param[out] ostr the output stream to which the text is printed
 * \param[in] pyne_map the map of all material objects
 *
 */
void fludag_all_materials(std::ostringstream& mstr,
                          std::map<std::string, pyne::Material> pyne_list);

/**
 * \brief Writes all the PyNE tally objects out that are contained in the tally map
 *
 * \param[out] mstr the output stream to which the tallies are to be printed
 * \param[in] tally_map the map of all tallies in the problem indexed by the tally name
 */
void fludag_all_tallies(std::ostringstream& mstr, std::map<std::string,pyne::Tally> tally_map);


/**
 * \brief get all properties (metadata) stored on all entities of a given dimension.
 *
 * \param[in] property the property you are looking for
 * \param[in] dimension of the entities to be queried
 * \param[in] delimiters the possible characters used as delimiters
 * \return map of vector of property values in an entity handlewise map
 */
std::map<moab::EntityHandle,std::vector<std::string> > get_property_assignments(std::string property,
    int dimension, std::string delimiters);

/*** end of uwuw functions ***/

#endif /* DAGMC_FLUKA_FLUKA_FUNCS_H */
