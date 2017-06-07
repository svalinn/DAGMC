#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"

#include <vector>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <limits>
#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
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

using namespace moab;

// define following macro for verbose debugging of random ray generation
//#define DEBUG

void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem);

void dump_pyfile( char* filename, double timewith, double timewithout, double tmem, DagMC& dagmc,
		  OrientedBoxTreeTool::TrvStats* trv_stats, EntityHandle tree_root );

static const double PI = acos(-1.0);
static const double denom = 1.0 / ((double) RAND_MAX);
static const double denomPI = PI * denom;
  
inline void RNDVEC(CartVect& uvw, double &az) 
{
  // denom normalizes rand values (see global defines)
  double theta = az * denom * rand(); // randomly samples from 0 to az. (Default az is 2PI)
  double u = 2 * denom * rand() - 1; // randomly samples from -1 to 1.
  uvw[0] = sqrt(1-u*u)*cos(theta);
  uvw[1] = sqrt(1-u*u)*sin(theta);
  uvw[2] = u;

}

/* program global data, including settings with their defaults*/
typedef struct{ CartVect p; CartVect v; } ray_t;
std::vector< ray_t > rays; // list of user-specified rays (given with -f flag)

static CartVect ray_source(0,0,0);

static double source_rad = 0;
static int vol_index = 1;
static int num_random_rays = 1000;
static int randseed = 12345;
static bool do_stat_report = false;
static bool do_trv_stats   = false;
static double location_az = 2.0 * PI;
static double direction_az = location_az;
static const char* pyfile = NULL;

static int random_rays_missed = 0; // count of random rays that did not hit a surface

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
    str << "-S  track and print OBB tree traversal statistics" << std::endl;
    str << "-i <int>   specify volume to upon which to test ray intersections (default 1)" << std::endl;
    str << "-n <int>   specify number of random rays to fire (default 1000)" << std::endl;
    str << "-c <x> <y> <z>  Specify center of of random ray generation (default origin)." << std::endl;
    str << "-r <real>  random ray radius.  Random rays begin at this distance from the center." << std::endl;
    str << "           if < 0, fire rays inward through the center" << std::endl;
    str << "           if >= 0, fire randomly directed rays.  (default 0)" << std::endl;
    str << "-f <x> <y> <z> <u> <v> <w>  Fire one given ray and report result." << std::endl;
    str << "           (May be given multiple times.  -f implies -n 0)" << std::endl;
    str << "-z <int>   seed the random number generator (default 12345)" << std::endl;
    str << "-L <real>  if present, limit random ray Location to between +-<value> degrees" << std::endl;
    str << "-D <real>  if present, limit random ray Direction to between +-<value> degrees" << std::endl;
    str << "           (unused if random ray radius < 0)" << std::endl;
    str << "-p <filename>  if present, save parameters and results to a python dictionary" << std::endl;
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

static CartVect parse_vector( int& i, int argc, char* argv[] ){
  double params[3]; bool err = false;
  for( int j = 0; j<3 && !err; ++j ){
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
      usage( "Expected vector specified as <x> <y> <z>", 0, argv[0] );
  }

  return CartVect(params);
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

  CartVect point(params), direction(params+3); 
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
        case 'c':
          ray_source = parse_vector( i, argc, argv );
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
        case 'p':
	  pyfile = get_option( i, argc, argv );
	  break;
      }
    }
    else {
      if( !filename ){ filename =  argv[i]; }
      else{ usage( "Unexpected parameter", 0, argv[0] ); }
    }
  }

  std::cout << "DagMC Version: " << DagMC::version()
	    << ", svn revision: " << DagMC::interface_revision() << std::endl;

  if( !filename ){
    usage("No filename specified", 0, argv[0] );
  }
     
  ErrorCode rval;
  EntityHandle surf = 0, vol = 0;
  double dist;

  OrientedBoxTreeTool::TrvStats* trv_stats = NULL;
  if( do_trv_stats ){ trv_stats = new OrientedBoxTreeTool::TrvStats; }

  DagMC dagmc = DagMC();
  rval = dagmc.load_file( filename );
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

      // added ray orientation
      rval = dagmc.ray_fire( vol, ray.p.array(), ray.v.array(), surf, dist, NULL, 0, 1, trv_stats );

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

  CartVect xyz, uvw;

  double ttime1, utime1, stime1, tmem1, ttime2, utime2, stime2, tmem2;
  get_time_mem(ttime1, utime1, stime1, tmem1);

  srand( randseed );

#ifdef DEBUG
  double uavg = 0.0, vavg = 0.0, wavg = 0.0;
#endif
  
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(uvw, location_az);

    xyz = uvw * source_rad + ray_source;
    if (source_rad >= 0.0) {
      RNDVEC(uvw, direction_az);
    }

#ifdef DEBUG
    std::cout << "x,y,z,u,v,w,u^2 + v^2 + w^2 = " << xyz 
              << " " << uvw << " " << uvw%uvw << std::endl;
    uavg += uvw[0]; vavg += uvw[1]; wavg += uvw[2];
#endif
    // added ray orientation
    dagmc.ray_fire(vol, xyz.array(), uvw.array(), surf, dist, NULL, 0, 1, trv_stats );

    if( surf == 0){ random_rays_missed++; }

  }
  get_time_mem(ttime2, utime2, stime2, tmem1);
  double timewith = ttime2 - ttime1;

  srand(randseed); // reseed to generate the same values as before

    // now without ray fire call, to subtract out overhead
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(uvw, location_az);

    xyz = uvw * source_rad + ray_source;
    if (source_rad >= 0.0) {
      RNDVEC(uvw, direction_az);
    }
  }
  
  get_time_mem(ttime1, utime1, stime1, tmem2);
  double timewithout = ttime1 - ttime2;

  std::cout << " done." << std::endl;

  if( random_rays_missed ){
    std::cout << "Warning: " << random_rays_missed << " random rays did not hit the target volume" << std::endl;
  }
  
  if( num_random_rays > 0 ){
    std::cout << "Total time per ray fire: " << timewith/num_random_rays 
	      << " sec" << std::endl;
    std::cout << "Estimated time per call (excluding ray generation): " 
	      << (timewith - timewithout) / num_random_rays << " sec" << std::endl;
  }
  std::cout << "Program memory used: " 
            << tmem2 << " bytes (" << tmem2/(1024*1024) << " MB)" << std::endl;

  /* Gather OBB tree stats and make final reports */
  EntityHandle root;
  ErrorCode result = dagmc.get_root(vol, root);
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

  if( do_trv_stats ){
    std::cout << "Traversal statistics:" << std::endl;
    trv_stats->print( std::cout );
  }

  if( pyfile ){
    dump_pyfile( filename, timewith, timewithout, tmem2, dagmc, trv_stats, root );
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

  // try going to /proc to estimate total memory
    char file_str[4096], dum_str[4096];
    int file_ptr = open("/proc/self/stat", O_RDONLY);
    int file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) {
      close(file_ptr);
      return;
    }

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


class HistogramBuilder : public OrientedBoxTreeTool::Op {

protected:
  int last_depth;

  void ensure_depth( unsigned depth ){
    while( node_counts.size() < depth+1 ){
      node_counts.push_back(0);
      leaf_counts.push_back(0);
    }
  }

public:
  std::vector<int> node_counts;
  std::vector<int> leaf_counts;

  virtual ErrorCode visit( EntityHandle /*node*/, int depth, bool& descend ){
    ensure_depth( depth );
    last_depth = depth;

    node_counts[last_depth] += 1;
    
    descend = true;
    return MB_SUCCESS;
  }
       
  virtual ErrorCode leaf( EntityHandle /*node*/ ){
    leaf_counts[last_depth] += 1;
    return MB_SUCCESS;
  }
};

void write_obbtree_histogram( EntityHandle root, OrientedBoxTreeTool& tree, std::ostream& out ){

  HistogramBuilder hb;
  tree.preorder_traverse( root, hb ); 

  out << "[";
  for( unsigned i = 0; i < hb.node_counts.size(); ++i ){
    out << "(" << hb.node_counts[i] << ", " << hb.leaf_counts.at(i) << "),";
  }
  out << "]";

}

void moab_memory_estimates( Interface* mbi, unsigned long long& moab_data_bytes, unsigned long long& moab_alldata_est_bytes ){
  
  unsigned long long storage, amortized_storage;
  mbi->estimated_memory_use( NULL, 0, &storage, &amortized_storage );

  moab_data_bytes = storage;
  moab_alldata_est_bytes = amortized_storage;

}

#define DICT_VAL(X) out << "'" #X "':" << X << "," << std::endl;
#define DICT_VAL_STR(X) out << "'" #X  "':'" << X << "'," << std::endl;
void dump_pyfile( char* filename, double timewith, double timewithout, double tmem, DagMC& dagmc,
		  OrientedBoxTreeTool::TrvStats* trv_stats, EntityHandle tree_root ){
  std::ofstream out(pyfile);
  out.precision( 14 );

  out << "# This file was automatically generated by ray_fire_test, and contains" << std::endl;
  out << "# the results of a single program run.  It is formatted as a python dictionary," << std::endl;
  out << "# and its contents may be loaded into python with a command such as: " << std::endl;
  out << "#  data = eval( compile( open('filename').read(), 'filename','eval'))" << std::endl;
  out << "# (with the name of this file substituted for 'filename')" << std::endl;
  out << "{" << std::endl;

  time_t now = time(NULL);
  std::string finish_time = ctime(&now);
  // remove trailing \n
  finish_time.resize( finish_time.length()-1 );
  DICT_VAL_STR(finish_time);
  DICT_VAL_STR(filename);
  out << "'ray_source':[" << ray_source[0] << "," << ray_source[1] << "," << ray_source[2] << "]," << std::endl;
  DICT_VAL(source_rad);
  unsigned num_user_rays = rays.size();
  DICT_VAL(num_user_rays);
  out << "'user_rays':[";
  for( unsigned i = 0; i < num_user_rays; ++i ){
    out << "(("   << rays[i].p[0] << "," << rays[i].p[1] << "," << rays[i].p[2] 
	<< "),(" << rays[i].v[0] << "," << rays[i].v[1] << "," << rays[i].v[2] << ")),";
  }
  out << "]," << std::endl;

  DICT_VAL(num_random_rays);
  DICT_VAL(random_rays_missed);
  if( num_random_rays > 0 ){
    DICT_VAL(randseed);
    DICT_VAL(timewith);
    DICT_VAL(timewith-timewithout);
  }
  DICT_VAL(tmem);
  unsigned long long moab_data_bytes, moab_alldata_est_bytes;
  moab_memory_estimates( dagmc.moab_instance(), moab_data_bytes, moab_alldata_est_bytes );
  DICT_VAL( moab_data_bytes );
  DICT_VAL( moab_alldata_est_bytes );

  unsigned int entities_in_tree, tree_height, node_count, num_leaves;
  double root_volume, tot_node_volume, tot_to_root_volume;
  dagmc.obb_tree()->stats(tree_root, entities_in_tree, root_volume, tot_node_volume,
			  tot_to_root_volume, tree_height, node_count, num_leaves);
  
  DICT_VAL( entities_in_tree );
  DICT_VAL( tree_height );
  DICT_VAL( node_count );
  DICT_VAL( num_leaves );

  out << "'tree_structure':";
  write_obbtree_histogram( tree_root, *dagmc.obb_tree(), out );
  out << "," <<  std::endl;

  if(trv_stats){
    unsigned stat_depth = trv_stats->nodes_visited().size();
    out << "'nodes_visited':[";
    for( unsigned i = 0; i < stat_depth; ++i){
      out << trv_stats->nodes_visited()[i] << ",";
    }
    out << "]," << std::endl;

    out << "'leaves_visited':[";
    for( unsigned i = 0; i < stat_depth; ++i){
      out << trv_stats->leaves_visited()[i] << ",";
    }
    out << "]," << std::endl;

    out << "'traversals_ended':[";
    for( unsigned i = 0; i < stat_depth; ++i){
      out << trv_stats->traversals_ended()[i] << ",";
    }
    out << "]," << std::endl;

    unsigned int tri_test_count = trv_stats->ray_tri_tests(); 
    DICT_VAL( tri_test_count );
  }

  out << "'stat_string':\"\"\"";
  dagmc.obb_tree()->stats(tree_root, out);
  out << "\"\"\"" << std::endl;

  out << "}" << std::endl;
  
  
}
