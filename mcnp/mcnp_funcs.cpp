#include "mcnp_funcs.h"

#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
using moab::DagMC;

#include <limits>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef CUBIT_LIBS_PRESENT
#include <fenv.h>
#endif

// globals

moab::DagMC *DAG;
dagmcMetaData *DMD;
UWUW *workflow_data;

#define DGFM_SEQ   0
#define DGFM_READ  1
#define DGFM_BCAST 2

#ifdef ENABLE_RAYSTAT_DUMPS

#include <fstream>
#include <numeric>

static std::ostream* raystat_dump = NULL;
#endif


/* Static values used by dagmctrack_ */

static DagMC::RayHistory history;
static int last_nps = 0;
static double last_uvw[3] = {0,0,0};
static std::vector< DagMC::RayHistory > history_bank;
static std::vector< DagMC::RayHistory > pblcm_history_stack;
static bool visited_surface = false;

static bool use_dist_limit = false;
static double dist_limit; // needs to be thread-local


void dagmcinit_(char *cfile, int *clen,  // geom
                char *ftol,  int *ftlen, // faceting tolerance
                int *parallel_file_mode, // parallel read mode
                double* dagmc_version, int* moab_version, int* max_pbl )
{

  moab::ErrorCode rval;

  // make new DagMC
  DAG = new moab::DagMC();

#ifdef ENABLE_RAYSTAT_DUMPS
  // the file to which ray statistics dumps will be written
  raystat_dump = new std::ofstream("dagmc_raystat_dump.csv");
#endif

  *dagmc_version = DAG->version();
  *moab_version = DAG->interface_revision();

  // terminate all filenames with null char
  cfile[*clen] = ftol[*ftlen] = '\0';

  // read geometry
  rval = DAG->load_file(cfile);
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to read input file: " << cfile << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef CUBIT_LIBS_PRESENT
  // The Cubit 10.2 libraries enable floating point exceptions.
  // This is bad because MOAB may divide by zero and expect to continue executing.
  // See MOAB mailing list discussion on April 28 2010.
  // As a workaround, put a hold exceptions when Cubit is present.

  fenv_t old_fenv;
  if ( feholdexcept( &old_fenv ) ) {
    std::cerr << "Warning: could not hold floating-point exceptions!" << std::endl;
  }
#endif


  // initialize geometry
  rval = DAG->init_OBBTree();
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // intialize the metadata
  DMD = new dagmcMetaData(DAG);
  DMD->load_property_data();
  // all metadata now loaded

  pblcm_history_stack.resize( *max_pbl+1 ); // fortran will index from 1

}

void dagmcwritefacets_(char *ffile, int *flen)  // facet file
{
  // terminate all filenames with null char
  ffile[*flen]  = '\0';

  moab::ErrorCode rval = DAG->write_mesh(ffile,*flen);
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to write mesh file: " << ffile <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // if there a uwuw material library write it
  // out to the fcad file
  if(workflow_data->material_library.size() != 0) {
    // get the mat lib
    std::map<std::string, pyne::Material> mat_lib = workflow_data->material_library;
    std::map<std::string, pyne::Material> ::iterator it;

    // iterate over the map
    for ( it = mat_lib.begin() ; it != mat_lib.end() ; ++it ) {
      // check to make sure we find "mat in the key"
      if(it->first.find("mat:") != std::string::npos) {
        // write the hdf5 file data
        it->second.write_hdf5(ffile,"/materials");
      } else {
        std::cout << "Warning: Spurious material " << it->first << std::endl;
        std::cout << "not written to file " << std::endl;
      }
    }
  }

  return;

}

void dagmcwritemcnp_(char* dagfile, char *lfile, int *llen, char *mcnp_version_major)  // file with cell/surface cards
{
  workflow_data = new UWUW(dagfile);
  std::string full_dagfilename = workflow_data->full_filepath;

  lfile[*llen]  = '\0';

  std::string lfname(lfile, *llen);

  std::cout << "Going to write an lcad file = " << lfname << std::endl;
  // Before opening file for writing, check for an existing file
  if( lfname != "lcad" ) {
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( lfname.c_str(), R_OK ) == 0 ) {
      std::cout << "DagMC: reading from existing lcad file " << lfname << std::endl;
      return;
    }
  }

  // by default overwrites the existing file at lfname.c_str()
  std::ofstream lcadfile( lfname.c_str(), std::ios::out );

  // string stream for output
  std::ostringstream lcadfile_str;

  write_cell_cards(lcadfile_str,mcnp_version_major);
  lcadfile_str << std::endl;
  write_surface_cards(lcadfile_str);
  lcadfile_str << std::endl;

  if(workflow_data->material_library.size() > 0)
    write_material_data(lcadfile_str);

  if(workflow_data->tally_library.size() > 0)
    write_tally_data(lcadfile_str);

  // all done
  lcadfile << lcadfile_str.str();
  lcadfile.close();
  return;
}

// write all cell related data
void write_cell_cards(std::ostringstream &lcadfile, char* mcnp_version_major)
{
  int num_cells = DAG->num_entities( 3 );

  std::string mat_num, density;

  // loop over all cells
  for( int i = 1; i <= num_cells; ++i ) {
    int cellid = DAG->id_by_index( 3, i );
    moab::EntityHandle entity = DAG->entity_by_index( 3, i );

    // deal with material number & density
    if(workflow_data->material_library.size() == 0) {
      // assuming simplified naming scheme check to make sure
      // that material numbers are assigned
      mat_num = DMD->volume_material_data_eh[entity];
      // if we cant make an int from the mat_num
      if(mat_num.find("Graveyard") == std::string::npos &&
         mat_num.find("Vacuum") == std::string::npos ) {
        if(!DMD->try_to_make_int(mat_num)) {
          std::cout << "Failed to cast material number to an integer"  << std::endl;
          std::cout << "volume with ID " << cellid << " has material assignment" << std::endl;
          std::cout << mat_num << " which appears to be a name rather than an integer" << std::endl;
          std::cout << "Did you forget to run uwuw_preproc?" << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      density = DMD->volume_density_data_eh[entity];
      // if we have a vacuum problem
      if(mat_num == "Graveyard" || mat_num == "Vacuum") {
        mat_num = "0";
        density = "";
      }
    } else {
      std::string mat_name = DMD->volume_material_property_data_eh[entity];
      // if we not vacuum or graveyard
      if(mat_name.find("Vacuum") == std::string::npos && mat_name.find("Graveyard") == std::string::npos) {
        if(workflow_data->material_library.count(mat_name) == 0) {
          std::cout << "Material with name " << mat_name <<  " not found " << std::endl;
          std::cout << "In the material library" << std::endl;
          exit(EXIT_FAILURE);
        }

        pyne::Material material = workflow_data->material_library[mat_name];
        int matnumber = material.metadata["mat_number"].asInt();
        mat_num = _to_string(matnumber);
        density = "-"+_to_string(material.density);
      } else {
        mat_num = "0";
        density = "";
      }
    }
    // write out cell, matnum and density
    lcadfile  << cellid << " " << mat_num << " " << density << " ";

    // string to collect importance data
    std::string importances = "";
    std::set<std::string>::iterator it;
    std::set<std::string> set = DMD->imp_particles;
    // deal with importances;
    std::string mat_name = DMD->volume_material_property_data_eh[entity];
    for ( it = set.begin() ; it != set.end() ; ++it) {
      std::string particle_name = *it;
      std::string mcnp_name;
      if(mcnp_version_major[0] == '5') {
        mcnp_name = pyne::particle::mcnp(particle_name);
      } else if(mcnp_version_major[0] == '6') {
        mcnp_name = pyne::particle::mcnp6(particle_name);
      } else {
        std::cout << "Unknown MCNP verison: " << mcnp_version_major << std::endl;
        exit(EXIT_FAILURE);
      }
      double imp = 1.0;
      // if we find graveyard always have importance 0.0
      if(mat_name.find("Graveyard") != std::string::npos) {
        imp = 0.0;
        // no splitting can happenin vacuum set to 1
      } else if (mat_name.find("Vacuum") != std::string::npos) {
        imp = 1.0;
        // otherwise as the map says
      } else {
        imp = DMD->importance_map[entity][particle_name];
      }
      importances += "imp:"+mcnp_name+"="+_to_string(imp)+" ";
    }
    // its possible no importances were assigned
    if(set.size() == 0) {
      if (mat_name.find("Graveyard") == std::string::npos ) {
        importances = "imp:n=1";
      } else {
        importances = "imp:n=0";
      }
    }

    // add descriptive comments for special volumes
    if (mat_name.find("Graveyard") != std::string::npos) {
      importances += "  $ graveyard";
    } else if (DAG->is_implicit_complement(entity)) {
      importances += "  $ implicit complement";
    }

    // write out importances to lcadfile
    lcadfile << importances << std::endl;

  }
  // all done
  return;
}

// write the surface data as appropriate
void write_surface_cards(std::ostringstream &lcadfile)
{
  int num_surfaces = DAG->num_entities( 2 );

  std::string surface_property = "";
  // loop over all cells
  for( int i = 1; i <= num_surfaces; ++i ) {
    int surfid = DAG->id_by_index( 2, i );
    moab::EntityHandle entity = DAG->entity_by_index(2,i);
    std::string boundary_prop = DMD->surface_boundary_data_eh[entity];
    if(boundary_prop.find("Reflecting") != std::string::npos)
      surface_property = "*";
    if(boundary_prop.find("White") != std::string::npos)
      surface_property = "+";
    lcadfile  << surface_property << _to_string(surfid) << std::endl;
  }
  return;
}

// write out all the tally data from the uwuw file
void write_material_data(std::ostringstream &lcadfile)
{
  std::map<std::string,pyne::Material> material_library = workflow_data->material_library;
  // loop over all tallies
  std::cout << "Writing Materials ..." << std::endl;

  lcadfile << "C materials from library" << std::endl;
  // loop over the material and print them out
  for(std::map<std::string,pyne::Material>::const_iterator it = material_library.begin() ;
      it != material_library.end() ; ++it ) {
    pyne::Material new_material = (it->second);
    std::string material_card = new_material.mcnp();
    lcadfile << material_card;
  }
  return;
}

// write out all the tally data from the uwuw file
void write_tally_data(std::ostringstream &lcadfile)
{
  std::map<std::string,pyne::Tally> tally_library = workflow_data->tally_library;
  // loop over all tallies
  std::cout << "Writing Tallies ..." << std::endl;
  int count = 1;
  for( std::map<std::string,pyne::Tally>::iterator it = tally_library.begin() ; it != tally_library.end() ; ++it ) {
    std::string tally_card = (it->second).mcnp(count,"mcnp5");
    lcadfile << tally_card;
    count++;
  }
}

void dagmcangl_(int *jsu, double *xxx, double *yyy, double *zzz, double *ang)
{
  moab::EntityHandle surf = DAG->entity_by_index( 2, *jsu );
  double xyz[3] = {*xxx, *yyy, *zzz};
  moab::ErrorCode rval = DAG->get_angle(surf, xyz, ang, &history );
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC: failed in calling get_angle" <<  std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "angl: " << *xxx << ", " << *yyy << ", " << *zzz << " --> "
            << ang[0] <<", " << ang[1] << ", " << ang[2] << std::endl;
  CartVect uvw(last_uvw);
  CartVect norm(ang);
  double aa = angle(uvw,norm) * (180.0/M_PI);
  std::cout << "    : " << aa << " deg to uvw" << (aa>90.0? " (!)":"")  << std::endl;
#endif

}

void dagmcchkcel_by_angle_( double *uuu, double *vvv, double *www,
                            double *xxx, double *yyy, double *zzz,
                            int *jsu, int *i1, int *j)
{


#ifdef TRACE_DAGMC_CALLS
  std::cout<< " " << std::endl;
  std::cout<< "chkcel_by_angle: vol=" << DAG->id_by_index(3,*i1) << " surf=" << DAG->id_by_index(2,*jsu)
           << " xyz=" << *xxx  << " " << *yyy << " " << *zzz << std::endl;
  std::cout<< "               : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};

  moab::EntityHandle surf = DAG->entity_by_index( 2, *jsu );
  moab::EntityHandle vol  = DAG->entity_by_index( 3, *i1 );

  int result;
  moab::ErrorCode rval = DAG->test_volume_boundary( vol, surf, xyz, uvw, result, &history );
  if( moab::MB_SUCCESS != rval ) {
    std::cerr << "DAGMC: failed calling test_volume_boundary" << std::endl;
    exit(EXIT_FAILURE);
  }

  switch (result) {
  case 1:
    *j = 0; // inside==  1 -> inside volume -> j=0
    break;
  case 0:
    *j = 1; // outside== 0  -> outside volume -> j=1
    break;
  default:
    std::cerr << "Impossible result in dagmcchkcel_by_angle" << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout<< "chkcel_by_angle: j=" << *j << std::endl;
#endif

}

void dagmcchkcel_(double *uuu,double *vvv,double *www,double *xxx,
                  double *yyy,double *zzz, int *i1, int *j)
{


#ifdef TRACE_DAGMC_CALLS
  std::cout<< " " << std::endl;
  std::cout<< "chkcel: vol=" << DAG->id_by_index(3,*i1) << " xyz=" << *xxx
           << " " << *yyy << " " << *zzz << std::endl;
  std::cout<< "      : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  int inside;
  moab::EntityHandle vol = DAG->entity_by_index( 3, *i1 );
  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};
  moab::ErrorCode rval = DAG->point_in_volume( vol, xyz, inside, uvw );

  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC: failed in point_in_volume" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  if (moab::MB_SUCCESS != rval) *j = -2;
  else
    switch (inside) {
    case 1:
      *j = 0; // inside==  1 -> inside volume -> j=0
      break;
    case 0:
      *j = 1; // outside== 0  -> outside volume -> j=1
      break;
    case -1:
      *j = 1; // onboundary== -1 -> on boundary -> j=1 (assume leaving volume)
      break;
    default:
      std::cerr << "Impossible result in dagmcchkcel" << std::endl;
      exit(EXIT_FAILURE);
    }

#ifdef TRACE_DAGMC_CALLS
  std::cout<< "chkcel: j=" << *j << std::endl;
#endif

}


void dagmcdbmin_( int *ih, double *xxx, double *yyy, double *zzz, double *huge, double* dbmin)
{
  double point[3] = {*xxx, *yyy, *zzz};

  // get handle for this volume (*ih)
  moab::EntityHandle vol  = DAG->entity_by_index( 3, *ih );

  // get distance to closest surface
  moab::ErrorCode rval = DAG->closest_to_location(vol,point,*dbmin);

  // if failed, return 'huge'
  if (moab::MB_SUCCESS != rval) {
    *dbmin = *huge;
    std::cerr << "DAGMC: error in closest_to_location, returning huge value from dbmin_" <<  std::endl;
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "dbmin " << DAG->id_by_index( 3, *ih ) << " dist = " << *dbmin << std::endl;
#endif

}

void dagmcnewcel_( int *jsu, int *icl, int *iap )
{

  moab::EntityHandle surf = DAG->entity_by_index( 2, *jsu );
  moab::EntityHandle vol  = DAG->entity_by_index( 3, *icl );
  moab::EntityHandle newvol = 0;

  moab::ErrorCode rval = DAG->next_vol( surf, vol, newvol );
  if( moab::MB_SUCCESS != rval ) {
    *iap = -1;
    std::cerr << "DAGMC: error calling next_vol, newcel_ returning -1" << std::endl;
  }

  *iap = DAG->index_by_handle( newvol );

  visited_surface = true;

#ifdef TRACE_DAGMC_CALLS
  std::cout<< "newcel: prev_vol=" << DAG->id_by_index(3,*icl) << " surf= "
           << DAG->id_by_index(2,*jsu) << " next_vol= " << DAG->id_by_index(3,*iap) <<std::endl;

#endif
}

void dagmc_surf_reflection_( double *uuu, double *vvv, double *www, int* verify_dir_change )
{


#ifdef TRACE_DAGMC_CALLS
  // compute and report the angle between old and new
  CartVect oldv(last_uvw);
  CartVect newv( *uuu, *vvv, *www );

  std::cout << "surf_reflection: " << angle(oldv,newv)*(180.0/M_PI) << std::endl;;
#endif

  // a surface was visited
  visited_surface = true;

  bool update = true;
  if( *verify_dir_change ) {
    if( last_uvw[0] == *uuu && last_uvw[1] == *vvv && last_uvw[2] == *www  )
      update = false;
  }

  if( update ) {
    last_uvw[0] = *uuu;
    last_uvw[1] = *vvv;
    last_uvw[2] = *www;
    history.reset_to_last_intersection();
  }

#ifdef TRACE_DAGMC_CALLS
  else {
    // mark it in the log if nothing happened
    std::cout << "(noop)";
  }

  std::cout << std::endl;
#endif

}

void dagmc_particle_terminate_( )
{
  history.reset();

#ifdef TRACE_DAGMC_CALLS
  std::cout << "particle_terminate:" << std::endl;
#endif
}

// *ih              - volue index
// *uuu, *vvv, *www - ray direction
// *xxx, *yyy, *zzz - ray point
// *huge            - passed to ray_fire as 'huge'
// *dls             - output from ray_fire as 'dist_traveled'
// *jap             - intersected surface index, or zero if none
// *jsu             - previous surface index
void dagmctrack_(int *ih, double *uuu,double *vvv,double *www,double *xxx,
                 double *yyy,double *zzz,double *huge,double *dls,int *jap,int *jsu,
                 int *nps )
{
  // Get data from IDs
  moab::EntityHandle vol = DAG->entity_by_index( 3, *ih );
  moab::EntityHandle prev = DAG->entity_by_index( 2, *jsu );
  moab::EntityHandle next_surf = 0;
  double next_surf_dist;

#ifdef ENABLE_RAYSTAT_DUMPS
  moab::OrientedBoxTreeTool::TrvStats trv;
#endif

  double point[3] = {*xxx,*yyy,*zzz};
  double dir[3]   = {*uuu,*vvv,*www};

  /* detect streaming or reflecting situations */
  if( last_nps != *nps || prev == 0 ) {
    // not streaming or reflecting: reset history
    history.reset();
#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: new history" << std::endl;
#endif

  } else if( last_uvw[0] == *uuu && last_uvw[1] == *vvv && last_uvw[2] == *www ) {
    // streaming -- use history without change
    // unless a surface was not visited
    if( !visited_surface ) {
      history.rollback_last_intersection();
#ifdef TRACE_DAGMC_CALLS
      std::cout << "     : (rbl)" << std::endl;
#endif
    }
#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: streaming " << history.size() << std::endl;
#endif
  } else {
    // not streaming or reflecting
    history.reset();

#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: reset" << std::endl;
#endif

  }

  moab::ErrorCode result = DAG->ray_fire(vol, point, dir,
                                         next_surf, next_surf_dist, &history,
                                         (use_dist_limit ? dist_limit : 0 )
#ifdef ENABLE_RAYSTAT_DUMPS
                                         , raystat_dump ? &trv : NULL
#endif
                                        );


  if(moab::MB_SUCCESS != result) {
    std::cerr << "DAGMC: failed in ray_fire" << std::endl;
    exit( EXIT_FAILURE );
  }


  for( int i = 0; i < 3; ++i ) {
    last_uvw[i] = dir[i];
  }
  last_nps = *nps;

  // Return results: if next_surf exists, then next_surf_dist will be nearer than dist_limit (if any)
  if( next_surf != 0 ) {
    *jap = DAG->index_by_handle( next_surf );
    *dls = next_surf_dist;
  } else {
    // no next surface
    *jap = 0;
    if( use_dist_limit ) {
      // Dist limit on: return a number bigger than dist_limit
      *dls = dist_limit * 2.0;
    } else {
      // Dist limit off: return huge value, triggering lost particle
      *dls = *huge;
    }
  }

  visited_surface = false;

#ifdef ENABLE_RAYSTAT_DUMPS
  if( raystat_dump ) {

    *raystat_dump << *ih << ",";
    *raystat_dump << trv.ray_tri_tests() << ",";
    *raystat_dump << std::accumulate( trv.nodes_visited().begin(), trv.nodes_visited().end(), 0 ) << ",";
    *raystat_dump << std::accumulate( trv.leaves_visited().begin(), trv.leaves_visited().end(), 0 ) << std::endl;

  }
#endif

#ifdef TRACE_DAGMC_CALLS

  std::cout<< "track: vol=" << DAG->id_by_index(3,*ih) << " prev_surf=" << DAG->id_by_index(2,*jsu)
           << " next_surf=" << DAG->id_by_index(2,*jap) << " nps=" << *nps <<std::endl;
  std::cout<< "     : xyz=" << *xxx << " " << *yyy << " "<< *zzz << " dist = " << *dls << std::flush;
  if( use_dist_limit && *jap == 0 ) std::cout << " > distlimit" << std::flush;
  std::cout << std::endl;
  std::cout<< "     : uvw=" << *uuu << " " << *vvv << " "<< *www << std::endl;
#endif

}

void dagmc_bank_push_( int* nbnk )
{
  if( ((unsigned)*nbnk) != history_bank.size() ) {
    std::cerr << "bank push size mismatch: F" << *nbnk << " C" << history_bank.size() << std::endl;
  }
  history_bank.push_back( history );

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_push (" << *nbnk+1 << ")" << std::endl;
#endif
}

void dagmc_bank_usetop_( )
{

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_usetop" << std::endl;
#endif

  if( history_bank.size() ) {
    history = history_bank.back();
  } else {
    std::cerr << "dagmc_bank_usetop_() called without bank history!" << std::endl;
  }
}

void dagmc_bank_pop_( int* nbnk )
{

  if( ((unsigned)*nbnk) != history_bank.size() ) {
    std::cerr << "bank pop size mismatch: F" << *nbnk << " C" << history_bank.size() << std::endl;
  }

  if( history_bank.size() ) {
    history_bank.pop_back( );
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_pop (" << *nbnk-1 << ")" << std::endl;
#endif

}

void dagmc_bank_clear_( )
{
  history_bank.clear();
#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_clear" << std::endl;
#endif
}

void dagmc_savpar_( int* n )
{
#ifdef TRACE_DAGMC_CALLS
  std::cout << "savpar: " << *n << " ("<< history.size() << ")" << std::endl;
#endif
  pblcm_history_stack[*n] = history;
}

void dagmc_getpar_( int* n )
{
#ifdef TRACE_DAGMC_CALLS
  std::cout << "getpar: " << *n << " (" << pblcm_history_stack[*n].size() << ")" << std::endl;
#endif
  history = pblcm_history_stack[*n];
}


void dagmcvolume_(int* mxa, double* vols, int* mxj, double* aras)
{
  moab::ErrorCode rval;

  // get size of each volume
  int num_vols = DAG->num_entities(3);
  for (int i = 0; i < num_vols; ++i) {
    rval = DAG->measure_volume( DAG->entity_by_index(3, i+1), vols[i*2] );
    if( moab::MB_SUCCESS != rval ) {
      std::cerr << "DAGMC: could not measure volume " << i+1 << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  // get size of each surface
  int num_surfs = DAG->num_entities(2);
  for (int i = 0; i < num_surfs; ++i) {
    rval = DAG->measure_area( DAG->entity_by_index(2, i+1), aras[i*2] );
    if( moab::MB_SUCCESS != rval ) {
      std::cerr << "DAGMC: could not measure surface " << i+1 << std::endl;
      exit( EXIT_FAILURE );
    }
  }

}

void dagmc_setdis_(double *d)
{
  dist_limit = *d;
#ifdef TRACE_DAGMC_CALLS
  std::cout << "setdis: " << *d << std::endl;
#endif
}

void dagmc_set_settings_(int* fort_use_dist_limit, int* use_cad, double* overlap_thickness, int* srccell_mode )
{

  if( *fort_use_dist_limit ) {
    std::cout << "DAGMC distance limit optimization is ENABLED" << std::endl;
    use_dist_limit = true;
  }

  if( *srccell_mode ) {
    std::cout << "DAGMC source cell optimization is ENABLED (warning: experimental!)" << std::endl;
  }

  DAG->set_overlap_thickness( *overlap_thickness );

}

void dagmc_init_settings_(int* fort_use_dist_limit, int* use_cad,
                          double* overlap_thickness, double* facet_tol, int* srccell_mode )
{

  *fort_use_dist_limit = use_dist_limit ? 1 : 0;

  *overlap_thickness = DAG->overlap_thickness();

  *facet_tol = DAG->faceting_tolerance();


  if( *srccell_mode ) {
    std::cout << "DAGMC source cell optimization is ENABLED (warning: experimental!)" << std::endl;
  }
}

// delete the stored data
void dagmc_teardown_()
{
  delete DMD;
  delete DAG;
}

// these functions should be replaced when we adopt C++11
// int to string
std::string _to_string(int var)
{
  std::ostringstream outstr;
  outstr << var;
  std::string ret_string = outstr.str();
  return ret_string;
}

// double to string
std::string _to_string(double var)
{
  std::ostringstream outstr;
  outstr << var;
  std::string ret_string = outstr.str();
  return ret_string;
}
