#include "DagMC.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"

#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>

#define DEFAULT_NUMRAYS 1000
#define DAG moab::DagMC::instance()

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

moab::ErrorCode get_mat_rho(moab::EntityHandle vol,
                            std::string &mat_name,
                            double &density) {

  moab::ErrorCode rval;
  std::vector< std::string > tmp_properties;

  mat_name = "";
  if (DAG->has_prop(vol,"mat")) {
    rval = DAG->prop_values(vol,"mat",tmp_properties);
    mat_name = tmp_properties[0];
  } else {
    mat_name = "Void";
  }


  density = 0;
  if (DAG->has_prop(vol,"rho")) {
    rval = DAG->prop_values(vol,"rho",tmp_properties);
    density = atof(tmp_properties[0].c_str());
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode find_start_grave_vols(moab::CartVect &ref_point, 
                                      moab::EntityHandle &graveyard, 
                                      moab::EntityHandle &start_vol) {

  moab::ErrorCode rval;

  // for volume in list of volumes check point_in_volume for reference point
  //   & find graveyard
  for (int vol_num=DAG->num_entities(3); vol_num > 0; vol_num--) {
    int inside;
    moab::EntityHandle test_vol = DAG->entity_by_index(3,vol_num);

    std::vector< std::string > tmp_properties;
    
    // find graveyard
    if (DAG->has_prop(test_vol,"mat")){
      rval = DAG->prop_values(test_vol,"mat",tmp_properties);
      for ( int mat_num=tmp_properties.size()-1; mat_num >= 0; mat_num--) {
        if ( !tmp_properties[mat_num].compare("Graveyard") ) {
          graveyard = test_vol;
        }
      }
    }

    rval = DAG->point_in_volume(test_vol,&ref_point[0],inside);
    if (moab::MB_SUCCESS != rval) {
      return rval;
    }
    
    // find start volume
    if (inside) {
      start_vol = test_vol;
      std::cout << "Found start volume: idx=" << DAG->index_by_handle(start_vol) 
                << " id=" << DAG->get_entity_id(start_vol) << std::endl;
    }
    
  }
  
  if (0 == start_vol) {
    std::cerr << "Error: Can't find starting volume" << std::endl;
    return moab::MB_FAILURE;
  }
  
  if (0 == graveyard) {
    std::cerr << "Warning: No graveryard found.  Using lost ray for boundary." << std::endl;
  }

  return moab::MB_SUCCESS;
  
}




int main(int argc, char* argv[]) {

  moab::ErrorCode rval;

  std::string geom_file, dir_filename;
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
  po.addOpt<std::string>("direction_file","List of directions");
  po.addOpt<double>("ref_point_x,x","X location of point to which rays should be fired from boundary");
  po.addOpt<double>("ref_point_y,y","Y location of point to which rays should be fired from boundary");
  po.addOpt<double>("ref_point_z,z","Z location of point to which rays should be fired from boundary");
  po.addOpt<int>("num_rays,n","Number of rays to fire (default=1000)");

  po.parseCommandLine( argc, argv );
  
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

  // read dir_file
  std::vector< moab::CartVect > dir_list;

  if (!po.getOpt("direction_file",&dir_filename)) {
    // read the number of rays from command line
    if (!po.getOpt("num_rays",&num_rays)) {
      num_rays = DEFAULT_NUMRAYS;
    }
    for (int ray_num = 0; ray_num < num_rays; ray_num++) {
      moab::CartVect dir(3);
      get_rand_dir(dir);
      dir_list.push_back(dir);
    }
  } else {
    std::ifstream dir_file(dir_filename.c_str());
    while (!dir_file.eof()) {
      moab::CartVect dir(3);
      dir_file >> dir[0] >> dir[1] >> dir[2];
      dir_list.push_back(dir);
      dir_file.get();
    }
  }
      
  // load geometry
  rval = DAG->load_file(geom_file.c_str());
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "Could not load file " << geom_file << std::endl;
    exit(rval);
  }

  // initialize OBB tree
  rval = DAG->init_OBBTree();

  // populate keywords
  std::vector< std::string > prop_keywords;
  std::map< std::string, std::string > prop_keyword_synonyms;
  std::string delimiters = ":/";

  prop_keywords.push_back("mat");
  prop_keywords.push_back("rho");
  rval = DAG->parse_properties(prop_keywords,prop_keyword_synonyms, delimiters.c_str());

  // function to get graveyard & start vol
  rval = find_start_grave_vols(ref_point,graveyard,start_vol);

  // for each ray requested in the input
  for (unsigned int ray_num = 0; ray_num < dir_list.size(); ray_num++) {

    std::cout << "Ray " << ray_num << std::endl;

    // initialize the new list
    moab::EntityHandle vol = start_vol;
    moab::EntityHandle surf = 0;
    moab::CartVect current_pt = ref_point;
    moab::CartVect dir = dir_list[ray_num];
    double dist = 0;
    slab_length.clear();
    slab_density.clear();
    slab_mat_name.clear();
    
    // while not at the graveyard
    while (vol != graveyard) {
      // std::cout << current_pt << "\t" << dir << std::endl;
      rval = DAG->ray_fire(vol,current_pt.array(),dir.array(),surf,dist);
      if (dist < huge && surf != 0) {
        // std::cout << dist << "\t" << surf << std::endl;
        slab_length.push_front(dist);
        double density;
        std::string mat_name;
        rval = get_mat_rho(vol,mat_name,density);
        rval = DAG->prop_value(vol,"mat",mat_name);
        // std::cout << mat_name << std::endl;
        slab_mat_name.push_front(mat_name);
        slab_density.push_front(density);
        moab::EntityHandle new_vol;
        rval = DAG->next_vol(surf,vol,new_vol);
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
  
