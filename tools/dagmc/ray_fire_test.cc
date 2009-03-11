#include "MBInterface.hpp"
#include "MBCore.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"

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

#define CHKERR if (MB_SUCCESS != rval) return rval

void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem);

#define RNDVEC(u, v, w) \
    {u = 2.0 * (denom * rand() - 0.5);          \
    v = 2.0 * (denom * rand() - 0.5); \
    w = 2.0 * (denom * rand() - 0.5); \
    double normal = 1.0 / (u*u + v*v + w*w); \
    u *= normal; \
    v *= normal; \
    w *= normal;}

int main( int argc, char* argv[] )
{
  MBErrorCode rval;

  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " [-f] <filename> "
              << " <facet_tol> <source_rad> <vol_index> <#calls> " << std::endl;
    std::cerr << "-f: report full tree statistics" << std::endl;
    std::cerr << "filename: mesh or geometry filename" << std::endl;
    std::cerr << "facet_tol: facet tolerance" << std::endl;
    std::cerr << "source_rad: if < 0, ray at source_rad and random angle pointed at origin;" << std::endl;
    std::cerr << "            otherwise, ray at radius and random angle, pointed in random direction;" << std::endl;
    std::cerr << "vol_index: index of the volume at which to fire rays" << std::endl;
    std::cerr << "#calls: # iterations of ray-tracing loop" << std::endl;
    return 1;
  }
  
  double facet_tol;
  int ncalls;
  bool full = false;
  int i = 1;
  if (!strcmp(argv[i], "-f")) {
    full = true;
    i++;
  }
  
  char* filename = argv[i++];
  facet_tol = atof(argv[i++]);
  ncalls = atoi(argv[i++]);
  double rad = atof(argv[i++]);
  int vol_idx = atoi(argv[i++]);
  ncalls = atoi(argv[i++]);
  
  DagMC& dagmc = *DagMC::instance();
  rval = dagmc.load_file_and_init( filename, strlen(filename), 0, 0, facet_tol);
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC." << std::endl;
    return 2;
  }

  MBEntityHandle vol = dagmc.entity_by_index(3, vol_idx);
  if (0 == vol) {
    std::cerr << "Problem getting first volume." << std::endl;
    return 2;
  }
  
  double ttime1, utime1, stime1, tmem1, ttime2, utime2, stime2, tmem2;
  get_time_mem(ttime1, utime1, stime1, tmem1);

    // initialize random number generator using ttime1
  srand((unsigned int) ttime1);
  double denom = 1.0 / ((double) RAND_MAX);
  double x, y, z, u, v, w, dist;
  MBEntityHandle nsurf;
  
  for (int j = 0; j < ncalls; j++) {
    RNDVEC(u, v, w);

    x = -u*rad;
    y = -v*rad;
    z = -w*rad;

    if (rad >= 0.0) {
      RNDVEC(u, v, w);
    }

    dagmc.ray_fire(vol, 0, 1, u, v, w, x, y, z, DBL_MAX,
                   dist, nsurf);
  }
  get_time_mem(ttime2, utime2, stime2, tmem1);
  double timewith = ttime2 - ttime1;

    // now without ray fire call, to subtract out overhead
  for (int j = 0; j < ncalls; j++) {
    RNDVEC(u, v, w);

    x = -u*rad;
    y = -v*rad;
    z = -w*rad;

    if (rad >= 0.0) {
      RNDVEC(u, v, w);
    }
  }
  
  get_time_mem(ttime1, utime1, stime1, tmem2);
  double timewithout = ttime1 - ttime2;
  
  MBEntityHandle root;
  MBErrorCode result = dagmc.get_root(vol, root);
  if (MB_SUCCESS != result) {
    std::cerr << "Trouble getting tree stats." << std::endl;
    return 2;
  }
    
  unsigned int entities_in_tree, tree_height, node_count, num_leaves;
  double root_volume, tot_node_volume, tot_to_root_volume;
  dagmc.obb_tree()->stats(root, entities_in_tree, root_volume, tot_node_volume,
                          tot_to_root_volume, tree_height, node_count, num_leaves);

  std::cout << "Gross time per call, net time per call, memory = " 
            << timewith/ncalls << " " << (timewith - timewithout)/ncalls
            << " " << tmem2 << std::endl;
  std::cout << "Tree stats: facets, tree_height, num_leaves = " 
            << entities_in_tree << " " << tree_height << " " << num_leaves << std::endl;

  if (full) {
    std::cout << "Tree data: " << std::endl;
    dagmc.obb_tree()->stats(root, std::cout);
  }
  
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
