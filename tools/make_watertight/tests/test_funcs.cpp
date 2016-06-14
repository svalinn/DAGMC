
#include "test_funcs.hpp"

moab::Interface* MBI() {
 static moab::Core instance;
 return &instance;
}

moab::ErrorCode move_vert(moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose) {

 moab::ErrorCode result;
 
 //get coordinates from the mesh
 double coords[3];
 result= MBI()->get_coords(&vertex, 1, coords);
 if(gen::error(moab::MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

 if(verbose) {
   //original coordinates
   std::cout << std::endl << "Original Coordinates" << std::endl;
   std::cout << "x = " << coords[0] << std::endl;
   std::cout << "y = " << coords[1] << std::endl;
   std::cout << "z = " << coords[2] << std::endl;
 }
  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;
  
 if(verbose)
 {
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result = MBI()->set_coords(&vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode rand_vert_move(moab::EntityHandle vertex, double tol, bool verbose) {

 moab::ErrorCode result;
 //get coordinates from the mesh
 double coords[3];
 result = MBI()->get_coords(&vertex, 1, coords);
 if(gen::error(moab::MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

 // get random values for the changes in x,y,z
 double dx,dy,dz;
 dx = rand(); dy = rand(); dz = rand();

 double mag = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

 // set the change in the vertex to be a unit vector
 dx/=mag; dy/=mag; dz/=mag;

 // set the change in the vertex to be something slightly less than the facet tolerance
 dx*=tol*0.9; dy*=tol*0.9; dz*=tol*0.9;
 
 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;

 if(verbose)
 {
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result = MBI()->set_coords(&vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return moab::MB_SUCCESS;
}

moab::ErrorCode single_vert_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex = verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert(vertex, bump_dist_x, bump_dist_y, bump_dist_z);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, std::string root_name, bool verbose) {
  
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode locked_pair_bump_rand(moab::Range verts, double facet_tol, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode rand_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
 
  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1 = (verts.back()-index);
  moab::EntityHandle vertex2 = (verts.back()-index-1);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode rand_locked_pair_bump_rand(moab::Range verts, double facet_tol, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
 
  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1 = (verts.back()-index);
  moab::EntityHandle vertex2 = (verts.back()-index-1);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode rand_vert_bump(moab::Range verts, double facet_tol, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex = verts.back();
  //move the desired vertex by the allotted distance
  result = rand_vert_move(vertex, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;

  return moab::MB_SUCCESS;
}

/// moves the third to last and the last vertices in the model the same distance in x, y, and z
moab::ErrorCode adjplone_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode adjplone_locked_pair_bump_rand(moab::Range verts, double facet_tol, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode nonadj_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode nonadj_locked_pair_bump_rand(moab::Range verts, double facet_tol, std::string root_name, bool verbose) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode write_mod_file(std::string filename) {

  moab::ErrorCode result;
  std::string output_filename = filename + "_mod.h5m";
  //write file
  result = MBI()->write_mesh(output_filename.c_str());
  if(gen::error(moab::MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 
  
  return moab::MB_SUCCESS;
}

moab::ErrorCode reload_mesh(const char* filename,  moab::EntityHandle &meshset, bool debug) {
  
  moab::ErrorCode result;
  // delete meshset
  result = MBI()->delete_mesh();
  if(gen::error(moab::MB_SUCCESS!=result, "could not delete the mesh set")) return result;

  // re-initialize meshset
  result = MBI()->create_meshset(moab::MESHSET_SET, meshset);
  if(gen::error(moab::MB_SUCCESS!=result, "could not create the new meshset")) return result; 
 
  //reload the file
  result = MBI()->load_file(filename, &meshset);
  if(gen::error(moab::MB_SUCCESS!=result, "could not re-load the file into the mesh set")) return result; 

  //check that something was loaded into the meshset
  int num_meshsets;   
  result = MBI()->num_contained_meshsets (meshset, &num_meshsets);
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the number of meshsets")) return result; 
  
  if(debug) std::cout << "number of mesh sets after reload = " << num_meshsets << std::endl;

  if(num_meshsets == 0) return moab::MB_FAILURE;

  return result;
}

void single_test_output(bool test_result) {
  if (test_result) {
    std::cout << " PASS" << std::endl;
  }
  else {
    std::cout << " FAIL" << std::endl;
  }
}

void test_set_output(std::string& set_title, bool set_result) {
  std::cout << set_title;
  if (set_result) {
    std::cout << " PASSED" << std::endl;
  }
  else {
    std::cout << " FAILED" << std::endl;
    exit(0);
  }
  std::cout << std::endl; //extra line to separate test set output
}

