#include "DagMC.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"

#include <vector>
#include <deque>
#include <iostream>
#include <math.h>
#include <limits>

#define DEFAULT_NUMRAYS 1000

void get_rand_dir(moab::CartVect &dir ) {

  static const double denom = 1.0 / ((double) RAND_MAX);
  static const double pi = acos(-1.0);

  double theta;

  dir[0] = dir[1] = dir[2] = 0;

  dir[2] = rand()*denom;
  dir[2] = 2*dir[2] -1;
  theta = 2 * pi * rand() * denom;

  dir[0] = dir[1] = sqrt(1-dir[2]*dir[2]);
  
  dir[0] *= cos(theta);
  dir[1] *= sin(theta);

}

int main(int argc, char* argv[]) {

  moab::DagMC* dag = moab::DagMC::instance();
  moab::ErrorCode rval;

  std::string geom_file;
  int num_rays;
  moab::CartVect ref_point;
  moab::EntityHandle start_vol = 0;
  moab::EntityHandle graveyard = 0;
  double huge = std::numeric_limits<double>::max();

  std::deque<double> slab_length;
  std::deque<double> slab_density;
  std::deque<std::string> slab_mat_name;

  // parse command-line for inputs
  // - REQUIRED: geometry file: no default
  // - reference point: default = origin
  // - number of rays: default = 1000
  ProgOptions po("hzetrn_paths: a tool to generate rays for HZETRN");
  po.addRequiredArg<std::string>("geometry_file","Path to the DAGMC-compliant geometry file", &geom_file);
  po.addOpt<double>("ref_point_x,x","X location of point to which rays should be fired from boundary");
  po.addOpt<double>("ref_point_y,y","Y location of point to which rays should be fired from boundary");
  po.addOpt<double>("ref_point_z,z","Z location of point to which rays should be fired from boundary");
  po.addOpt<int>("num_rays,n","Number of rays to fire (default=1000)");

  po.parseCommandLine( argc, argv );
  
  // load geometry
  rval = dag->load_file(geom_file.c_str());
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "Could not load file " << geom_file << std::endl;
    exit(rval);
  }

  // initialize OBB tree
  rval = dag->init_OBBTree();

  // find volume ID of reference point
  if (!po.getOpt("ref_point_x",&ref_point[0])) {
    ref_point[0] = 0;
  }
  if (!po.getOpt("ref_point_y",&ref_point[1])) {
    ref_point[1] = 0;
  }
  if (!po.getOpt("ref_point_z",&ref_point[2])) {
    ref_point[2] = 0;
  }

  // for volume in list of volumes check point_in_volume for reference point
  //   & find graveyard
  for (int vol_num=dag->num_entities(3); vol_num > 0; vol_num--) {
    int inside;
    moab::EntityHandle test_vol = dag->entity_by_index(3,vol_num);
    rval = dag->point_in_volume(test_vol,&ref_point[0],inside);
    
    // find graveyard
    if (dag->has_prop(test_vol,"graveyard")){
      graveyard = test_vol;
    }

    // find start volume
    if (inside) {
      start_vol = test_vol;
      std::cout << "Found start volume: idx=" << dag->index_by_handle(start_vol) 
                << " id=" << dag->get_entity_id(start_vol) << std::endl;
    }
    
  }
  
  if (0 == start_vol) {
    std::cerr << "Error: Can't find starting volume" << std::endl;
    exit(0);
  }
  
  if (0 == graveyard) {
    std::cerr << "Warning: No graveryard found.  Using lost ray for boundary." << std::endl;
  }
  
  // read the number of rays from command line
  if (!po.getOpt("num_rays",&num_rays)) {
    num_rays = DEFAULT_NUMRAYS;
  }

  // for each ray requested in the input
  for (int ray_num = 0; ray_num < num_rays; ray_num++) {

    std::cout << "Ray " << ray_num << std::endl;

    // initialize the new list
    moab::EntityHandle vol = start_vol;
    moab::EntityHandle surf = 0;
    moab::CartVect current_pt = ref_point;
    moab::CartVect dir;
    get_rand_dir(dir);
    double dist = 0;
    slab_length.clear();
    slab_density.clear();
    slab_mat_name.clear();

    // while not at the graveyard
    while (vol != graveyard) {
      // std::cout << current_pt << "\t" << dir << std::endl;
      rval = dag->ray_fire(vol,current_pt.array(),dir.array(),surf,dist);
      if (dist < huge && surf != 0) {
        // std::cout << dist << "\t" << surf << std::endl;
        slab_length.push_front(dist);
        std::string mat_name;
        rval = dag->prop_value(vol,"mat",mat_name);
        std::cout << mat_name << std::endl;
        slab_mat_name.push_front(mat_name);
        slab_density.push_front(-1);
        moab::EntityHandle new_vol;
        rval = dag->next_vol(surf,vol,new_vol);
        vol = new_vol;
        current_pt += dir*dist;
      } else {
        vol = graveyard;
      }
      
    }

    for ( unsigned int slab_num=0; slab_num<slab_length.size(); slab_num++) {
      std::cout << slab_mat_name[slab_num] << "\t" << slab_density[slab_num] << "\t" 
                << slab_length[slab_num] << std::endl;
    }
    
  }

}
