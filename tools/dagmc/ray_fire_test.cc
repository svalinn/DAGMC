#include "MBInterface.hpp"
#include "MBCore.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "MBCartVect.hpp"

#include <vector>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <limits>
#include <fcntl.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cfloat>
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

// define following macro for verbose debugging of random ray generation
//#define DEBUG

void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem);

static const double PI = acos(-1.0);
static const double denom = 1.0 / ((double) RAND_MAX);
static const double denomPI = PI * denom;
  
inline void RNDVEC(double &u, double &v, double &w, double &az) 
{
  
  double theta = denom * az * rand();
  double phi = denomPI * rand();
  u = cos(theta)*sin(phi);
  v = sin(theta)*sin(phi);
  w = cos(phi);
}

/* program global data, including settings with their defaults*/
typedef struct{ MBCartVect p; MBCartVect v; } ray_t;
std::vector< ray_t > rays; // list of user-specified rays (given with -f flag)


static double facet_tol = 1e-4;
static double source_rad = 0;
static int vol_index = 1;
static int num_random_rays = 1000;
static int randseed = 12345;
static bool do_stat_report = false;
static bool do_trv_stats   = false;
static double location_az = 2.0 * PI;
static double direction_az = location_az;

/* Most of the argument handling code was stolen/adapted from MOAB/test/obb/obb_test.cpp */
static void usage( const char* error, const char* opt, const char* name = "ray_fire_test" )
{
  const char* default_message = "Invalid option";
  if (opt && !error)
    error = default_message;

  std::ostream& str = error ? std::cerr : std::cout;
  if (error) {
    str << error;
    if (opt)
      str << ": " << opt;
    str << std::endl;
  }

  str << "Usage: " << name << " [options] input_file" << std::endl;
  str << "       " << name << " -h" << std::endl; 

  if( !error ){
    str << "-h  print this help" << std::endl;
    str << "-s  print OBB tree structural statistics" << std::endl;
    str << "-S  track and print OBB tree traversal statistics (to be implemented soon)" << std::endl;
    str << "-i <int>   specify volume to upon which to test ray intersections (default 1)" << std::endl;
    str << "-t <real>  specify faceting tolerance (default 1e-4)" << std::endl;
    str << "-n <int>   specify number of random rays to fire (default 1000)" << std::endl;
    str << "-r <real>  random ray radius.  Random rays begin at this distance from the origin." << std::endl;
    str << "           if < 0, fire rays inward through the origin" << std::endl;
    str << "           if >= 0, fire random rays outward.  (default 0)" << std::endl;
    str << "-f <x> <y> <z> <u> <v> <w>  Fire one given ray and report result." << std::endl;
    str << "           (May be given multiple times.  -s implies -n 0)" << std::endl;
    str << "-z <int>   seed the random number generator (default 12345)" << std::endl;
    str << "-L <real>  if present, limit random ray Location to between +-<value> degrees" << std::endl;
    str << "-D <real>  if present, limit random ray Direction to between +-<value> degrees" << std::endl;
  }

  exit( error ? 1 : 0 );
}

static const char* get_option( int& i, int argc, char* argv[] ) {
  ++i;
  if (i == argc)
    usage( "Expected argument following option", argv[i-1] );
  return argv[i];
}

static int get_int_option( int& i, int argc, char* argv[] ) {
  const char* str = get_option( i, argc, argv );
  char* end_ptr;
  long val = strtol( str, &end_ptr, 0 );
  if (!*str || *end_ptr) 
    usage( "Expected integer following option", argv[i-1] );
  return val;
}

static double get_double_option( int& i, int argc, char* argv[] ) {
  const char* str = get_option( i, argc, argv );
  char* end_ptr;
  double val = strtod( str, &end_ptr );
  if (!*str || *end_ptr) 
    usage( "Expected real number following option", argv[i-1] );
  return val;
}

static void parse_ray( int& i, int argc, char* argv[] )
{
  double params[6]; bool err = false;
  for( int j = 0; j<6 && !err; ++j ){
    if( ++i == argc ){
      err = true; break;
    }    
    char* end_ptr;
    params[j] = strtod( argv[i], &end_ptr );
    if( !argv[i][0] || *end_ptr ){
      err = true; break;
    }
  }
  if( err ){
      usage( "Expected ray specified as <x> <y> <z> <u> <v> <w>", 0, argv[0] );
  }

  MBCartVect point(params), direction(params+3); 
  direction.normalize(); 
  ray_t ray; ray.p = point; ray.v = direction;
  rays.push_back( ray );
}


int main( int argc, char* argv[] )
{

  char* filename = NULL;
  bool flags = true;
  for (int i = 1; i < argc; ++i) {
    if (flags && argv[i][0] =='-') {
      if (!argv[i][1] || argv[i][2])
        usage(0,argv[i],argv[0]);
      switch (argv[i][1]) {
        default:  usage( 0, argv[i], argv[0] );   break;
        case '-': flags = false;         break;
        case 'h': usage(0,0,argv[0]);    break;
        case 's': do_stat_report = true; break;
        case 'S': do_trv_stats   = true; break;
        case 'i': 
          vol_index = get_int_option( i, argc, argv );
          break;
        case 't': 
	  facet_tol = get_double_option( i, argc, argv );
	  break;
        case 'n':
	  num_random_rays = get_int_option( i, argc, argv );
	  break;
        case 'r':
	  source_rad = get_double_option( i, argc, argv );
	  break;
        case 'f':
	  parse_ray( i, argc, argv );
          num_random_rays = 0;
	  break;
        case 'z':
          randseed = get_int_option( i, argc, argv );
          break;
        case 'L':
          location_az = get_double_option( i, argc, argv ) * (PI / 180.0);
          break;
        case 'D':
          direction_az = get_double_option( i, argc, argv ) * (PI / 180.0);
          break;        
      }
    }
    else {
      if( !filename ){ filename =  argv[i]; }
      else{ usage( "Unexpected parameter", 0, argv[0] ); }
    }
  }

  if( !filename ){
    usage("No filename specified", 0, argv[0] );
  }
     
  MBErrorCode rval;
  MBEntityHandle surf = 0, vol = 0;
  double x, y, z, u, v, w, dist;


  /* Initialize DAGMC and find the appropriate volume */
  std::cout << "Initializing DagMC, facet_tol = " << facet_tol << std::endl;
  DagMC& dagmc = *DagMC::instance();
  rval = dagmc.load_file( filename, facet_tol );
  if(MB_SUCCESS != rval) {
    std::cerr << "Failed to load file '" << filename << "'" << std::endl;
    return 2;
  }
  
  rval = dagmc.init_OBBTree( );
  if(MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC." << std::endl;
    return 2;
  }
  
  vol = dagmc.entity_by_id(3, vol_index);
  if(0 == vol) {
    std::cerr << "Problem getting volume " << vol_index << std::endl;
    return 2;
  }

  /* Fire any rays specified with -f flag */
  if( rays.size() > 0 ){
    
    std:: cout << "Firing user-specified rays at volume " << vol_index << std::endl;
    for( unsigned i = 0; i < rays.size(); ++i ){

      ray_t ray = rays[i];
      std::cout << " Ray: point = " << ray.p << " dir = " << ray.v << std::endl;

      rval = dagmc.ray_fire( vol, 0, 1, 
                             ray.v[0], ray.v[1], ray.v[2], 
                             ray.p[0], ray.p[1], ray.p[2], 
                             DBL_MAX, dist, surf ); 

      if(MB_SUCCESS != rval) {
        std::cerr << "ERROR: ray_fire() failed!" << std::endl;
        return 2;
      }      
      if(0 == surf) {
        std::cerr << "ERROR: Ray finds no surface.  Particle is lost." << std::endl;
        // might as well keep going here, in case other user specified rays were given
        continue;
      }
      
      int surf_id = dagmc.id_by_index( 2, dagmc.index_by_handle( surf ) );
      std::cout << "       hits surf_id " << surf_id << " dist=" << dist << " new_xyz="
                << ray.p+(ray.v*dist) << std::endl;
    }
  }


  /* Fire and time random rays */
  if( num_random_rays > 0 ){
    std::cout << "Firing " << num_random_rays 
              << " random rays at volume " << vol_index << "..." << std::flush;
  }

  double ttime1, utime1, stime1, tmem1, ttime2, utime2, stime2, tmem2;
  get_time_mem(ttime1, utime1, stime1, tmem1);

  srand( randseed );

#ifdef DEBUG
  double uavg = 0.0, vavg = 0.0, wavg = 0.0;
#endif
  
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(u, v, w, location_az);

    x = -u*source_rad;
    y = -v*source_rad;
    z = -w*source_rad;

    if (source_rad >= 0.0) {
      RNDVEC(u, v, w, direction_az);
    }

#ifdef DEBUG
      std::cout << "x, y, z, u, v, w, u^2 + v^2 + w^2 = "
                << u << " " << v << " " << w << " " << (u*u + v*v + w*w) << std::endl;
      uavg += u; vavg += v; wavg += w;
#endif
    
    dagmc.ray_fire(vol, 0, 1, u, v, w, x, y, z, DBL_MAX,
                   dist, surf);

  }
  get_time_mem(ttime2, utime2, stime2, tmem1);
  double timewith = ttime2 - ttime1;

    // now without ray fire call, to subtract out overhead
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(u, v, w, location_az);

    x = -u*source_rad;
    y = -v*source_rad;
    z = -w*source_rad;

    if (source_rad >= 0.0) {
      RNDVEC(u, v, w, direction_az);
    }
  }
  
  get_time_mem(ttime1, utime1, stime1, tmem2);
  double timewithout = ttime1 - ttime2;

  std::cout << " done." << std::endl;
  
   
  std::cout << "Total time per ray fire: " << timewith/num_random_rays 
            << " sec" << std::endl;
  std::cout << "Estimated time per call (excluding ray generation): " 
            << (timewith - timewithout) / num_random_rays << " sec" << std::endl;
  std::cout << "Program memory used: " 
            << tmem2 << " bytes (" << tmem2/(1024*1024) << " MB)" << std::endl;

  /* Gather OBB tree stats and make final reports */
  MBEntityHandle root;
  MBErrorCode result = dagmc.get_root(vol, root);
  if (MB_SUCCESS != result) {
    std::cerr << "Trouble getting tree stats." << std::endl;
    return 2;
  }

  if (do_stat_report) {
    std::cout << "Tree statistics: " << std::endl;
    dagmc.obb_tree()->stats(root, std::cout);
  }
  else{
    unsigned int entities_in_tree, tree_height, node_count, num_leaves;
    double root_volume, tot_node_volume, tot_to_root_volume;
    dagmc.obb_tree()->stats(root, entities_in_tree, root_volume, tot_node_volume,
                            tot_to_root_volume, tree_height, node_count, num_leaves);

    std::cout << "Tree dimensions:" << std::endl;
    std::cout << "   facets: " << entities_in_tree << ", height: " << tree_height << std::endl;
    std::cout << "   num leaves: " << num_leaves << ", num nodes: " << node_count << std::endl;
  }

#ifdef DEBUG
  std::cout << "uavg, vavg, wavg = " << uavg << " " << vavg << " " << wavg << std::endl;
#endif

  return 0;

}

void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  user_time = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  sys_time = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = user_time + sys_time;
  tot_mem = 0;
  if (0 != r_usage.ru_maxrss) {
    tot_mem = r_usage.ru_idrss; 
  }
  else {
      // this machine doesn't return rss - try going to /proc
      // print the file name to open
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24)
      tot_mem = ((double)vm_size);
  }
}
