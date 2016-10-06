// std includes
#include <iostream>
#include <string>
// moab includes
#include "moab/ProgOptions.hpp"
#include "moab/ScdInterface.hpp"
#include "DagMC.hpp"

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

ErrorCode test_sphere();
ErrorCode test_cylinder();

//dagmc includes
int main() {

  ErrorCode rval;
  
  rval = test_sphere();
  MB_CHK_SET_ERR(rval, "DAGMC Preconditioner Sphere Test Failed");

  rval = test_cylinder();
  MB_CHK_SET_ERR(rval, "DAGMC Preconditioner Sphere Test Failed");
  
				    
  return 0;
}

ErrorCode dag_init_file(char *filename, DagMC* &dagmc, ScdBox* &box) {
  //create a new dagmc instance and load the file
  ErrorCode rval;
  dagmc = new DagMC();
  rval = dagmc->load_file(filename);
  MB_CHK_SET_ERR(rval,"Could not load the preconditioner test file");

  //get the ScdBox of the first (and only) volume
  Range vols;
  const int three = 3;
  const void* const three_val[] = {&three};
  Tag geomTag = dagmc->geom_tag();
  rval = dagmc->moab_instance()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, three_val, 1, vols );
  MB_CHK_SET_ERR(rval,"Could not get the model volumes");

  //there shoul only be one volume before we init the OBB tree
  if(1 != vols.size()) {
    MB_CHK_SET_ERR(MB_FAILURE, "Incorrect number of volumes found after loading the model");
  }
  //initalize the OBBTree (and the preconditioning datastructs)
  rval = dagmc->init_OBBTree();
  MB_CHK_SET_ERR(rval,"Could not initialize the OBBTree");
  //get the preconditioning box of interest
  box = dagmc->get_preconditioning_box(vols[0]);

  return rval;
};

// returns the nearest point on a sphere centered on the origin
CartVect nearest_on_sphere(CartVect point, double radius) {
  CartVect nearest_location = point;
  //get unit direction vector from origin
  nearest_location.normalize();
  //extend vector to sphere radius
  nearest_location *= radius;
  return nearest_location;
};

// returns the nearest point on a z-axis aligned cylinder centered on the origin
CartVect nearest_on_cylinder(CartVect point, double height, double radius) {
  CartVect nearest_location;

  // find point on nearest lid plane
  CartVect nearest_to_lids;
  if(point[2] > height/2) {
    nearest_to_lids = CartVect(point[0],point[1],height/2);
  }
  else {
    nearest_to_lids = CartVect(point[0],point[1],-height/2);
  }

  //  std::cout << nearest_to_lids << std::endl;
  
  // find nearest point on an infinite barrel
  CartVect nearest_to_inf_barrel = CartVect(point[0],point[1],0);
  nearest_to_inf_barrel.normalize();
  nearest_to_inf_barrel*=5;
  nearest_to_inf_barrel[2] = point[2];
  //  std::cout << nearest_to_inf_barrel << std::endl;
  // check if point is inside or outside the barrel
  if(CartVect(point[0],point[1],0).length() < radius) {
    // point is inside the barrel, then check for the minimum between the barrel and lids
    if((nearest_to_inf_barrel-point).length() < (nearest_to_lids-point).length()) {
      nearest_location = nearest_to_inf_barrel;
    }
    else {
      nearest_location = nearest_to_lids;
    }
  }
  else {
    // point is outside the barrel, return the inf barrel location
    nearest_location = nearest_to_inf_barrel;
  }
  return nearest_location;
};

  
ErrorCode test_sphere(){

  ErrorCode rval;
  DagMC* dagmc;
  ScdBox* box = NULL;

  char *filename = STRINGIFY(MESHDIR) "/sphere_rad5.h5m";
  rval = dag_init_file(filename, dagmc, box);
  MB_CHK_SET_ERR(rval,"Could not initalize DAGMC instance and retrieve preconditioner box");

  //get the signed distance field tag
  Tag sdfTag = dagmc->sdf_tag();
  int xints, yints, zints;
  box->box_size(xints,yints,zints);
  for(unsigned int i = 0 ; i < xints; i++){
    for(unsigned int j = 0 ; j < yints; j++){
      for(unsigned int k = 0 ; k < zints; k++){
	//retrieve vertex handle from box
	EntityHandle vert = box->get_vertex(i,j,k);
	//get the vertex coordinates
	CartVect vert_coords;
	rval = dagmc->moab_instance()->get_coords(&vert, 1, vert_coords.array());
	MB_CHK_SET_ERR(rval,"Could not get SCD vertex coords");

	//get the distance to the nearest point on the sphere
	double sphere_radius = 5;
	CartVect expected_location = nearest_on_sphere(vert_coords, sphere_radius);
	double expected_distance = fabs((vert_coords-expected_location).length());	
	
	//retrieve the stored distance value of this vertex
	double distance;
	void *ptr = &distance;
	rval = dagmc->moab_instance()->tag_get_data( sdfTag, &vert, 1, ptr);
	MB_CHK_SET_ERR(rval,"Could not retrieve signed distance field tag value");
	distance = fabs(distance);
	
	//use facet tolerance as maximal error
	double facet_tol = dagmc->faceting_tolerance();
	//compare the values - they should be off by no more than the faceting tolerance
	if(fabs(expected_distance-distance) > facet_tol) {
	  MB_CHK_SET_ERR(MB_FAILURE, "Incorrect distance value found");
	  std::cout << vert_coords << std::endl;
	  std::cout << expected_distance << std::endl;
	  std::cout << distance << std::endl;
	}
	
      }
    }
  }
  return rval;
};

ErrorCode test_cylinder(){

  ErrorCode rval;
  DagMC* dagmc;
  ScdBox* box = NULL;

  char *filename = STRINGIFY(MESHDIR) "/cyl_h5_r5.h5m";
  rval = dag_init_file(filename, dagmc, box);
  MB_CHK_SET_ERR(rval,"Could not initalize DAGMC instance and retrieve preconditioner box");

  //get the signed distance field tag
  Tag sdfTag = dagmc->sdf_tag();
  int xints, yints, zints;
  box->box_size(xints,yints,zints);
  for(unsigned int i = 0 ; i < xints; i++){
    for(unsigned int j = 0 ; j < yints; j++){
      for(unsigned int k = 0 ; k < zints; k++){
	//retrieve vertex handle from box
	EntityHandle vert = box->get_vertex(i,j,k);
	//get the vertex coordinates
	CartVect vert_coords;
	rval = dagmc->moab_instance()->get_coords(&vert, 1, vert_coords.array());
	MB_CHK_SET_ERR(rval,"Could not get SCD vertex coords");

	//get the distance to the nearest point on the sphere
	double cylinder_radius = 5, cylinder_height = 5;
	CartVect expected_location = nearest_on_cylinder(vert_coords, cylinder_radius, cylinder_radius);
	double expected_distance = fabs((vert_coords-expected_location).length());
	
	//retrieve the stored distance value of this vertex
	double distance;
	void *ptr = &distance;
	rval = dagmc->moab_instance()->tag_get_data( sdfTag, &vert, 1, ptr);
	MB_CHK_SET_ERR(rval,"Could not retrieve signed distance field tag value");
	distance = fabs(distance);
	
	//use facet tolerance as maximal error
	double facet_tol = dagmc->faceting_tolerance();
	//compare the values - they should be off by no more than the faceting tolerance
	if(fabs(expected_distance-distance) > facet_tol) {
	  MB_CHK_SET_ERR(MB_FAILURE, "Incorrect distance value found");
	  std::cout << vert_coords << std::endl;
	  std::cout << expected_distance << std::endl;
	  std::cout << distance << std::endl;
	}
	
      }
    }
  }
  return rval;
};
