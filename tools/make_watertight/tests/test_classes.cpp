
#include "test_classes.hpp"

moab::Interface* MBI()
{
  static moab::Core instance;
  return &instance;
}

void MakeWatertightTest::SetUp()
{
  // set this classes test filename
  setFilename();

  // setup the MakeWatertight class
  mw = new MakeWatertight(MBI());
  cw = new CheckWatertight(MBI());

  // delete meshset
  result = MBI()->delete_mesh();
  EXPECT_EQ(result,moab::MB_SUCCESS);

  // re-initialize meshset
  result = MBI()->create_meshset(moab::MESHSET_SET, input_fileset);
  EXPECT_EQ(result,moab::MB_SUCCESS);

  //reload the file
  result = MBI()->load_file(filename.c_str(), &input_fileset);
  EXPECT_EQ(result,moab::MB_SUCCESS);

  //// get faceting tolerance ////
  moab::Tag faceting_tol_tag;
  //get faceting tolerance handle from file
  result = MBI()->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
                                 faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  EXPECT_EQ(result,moab::MB_SUCCESS);

  //get the faceting tolerance of any entity
  moab::Range file_set;
  result = MBI()->get_entities_by_type_and_tag(0, moab::MBENTITYSET,
           &faceting_tol_tag, NULL, 1, file_set);

  //get facetint tolerance value
  result = MBI()->tag_get_data(faceting_tol_tag, &file_set.front(), 1, &facet_tol);
  EXPECT_EQ(result,moab::MB_SUCCESS);

  //check that something was loaded into the meshset
  int num_meshsets;
  result = MBI()->num_contained_meshsets (input_fileset, &num_meshsets);
  EXPECT_EQ(result,moab::MB_SUCCESS);
  EXPECT_NE(0, num_meshsets);

  //retrieve the verticies again so the model can be broken
  int dim = 0;
  result = MBI()->get_entities_by_dimension(input_fileset, dim, verts, false);
  EXPECT_EQ(result,moab::MB_SUCCESS);
};

void MakeWatertightTest::TearDown()
{
  result = MBI()->delete_mesh();
  EXPECT_EQ(result,moab::MB_SUCCESS);
};

moab::ErrorCode MakeWatertightTest::check_num_ents(int ent_dimension, int expected_num)
{
  moab::ErrorCode result;
  moab::Range entities;
  moab::Tag geom_tag;
  result = MBI()->tag_get_handle( "GEOM_DIMENSION", geom_tag);
  MB_CHK_SET_ERR(result, "could not get the geometry dimension tag");
  EXPECT_EQ(result, moab::MB_SUCCESS);
  void *tag_ptr = &geom_tag;
  const void *val_ptr = &ent_dimension;
  result = MBI()->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &geom_tag, &val_ptr, 1, entities, moab::Interface::INTERSECT, true);
  MB_CHK_SET_ERR(result, "could not get the number of entities by dimension");
  EXPECT_EQ(result, moab::MB_SUCCESS);

  EXPECT_EQ(expected_num, entities.size());

  return result;
}

moab::ErrorCode MakeWatertightTest::move_vert(moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose)
{

//get coordinates from the mesh
  double coords[3];
  result= MBI()->get_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not get the vertex coordinates");

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

  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }

  //write new coordinates to the mesh
  result = MBI()->set_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not set the vertex coordinates");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::rand_vert_move(moab::EntityHandle vertex, double tol, bool verbose)
{

  //get coordinates from the mesh
  double coords[3];
  result = MBI()->get_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not get the vertex coordinates");

  // get random values for the changes in x,y,z
  double dx,dy,dz;
  dx = rand();
  dy = rand();
  dz = rand();

  double mag = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

  // set the change in the vertex to be a unit vector
  dx/=mag;
  dy/=mag;
  dz/=mag;

  // set the change in the vertex to be something slightly less than the facet tolerance
  dx*=tol*0.9;
  dy*=tol*0.9;
  dz*=tol*0.9;

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

  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }

  //write new coordinates to the mesh
  result = MBI()->set_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not set the vertex coordinates");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::single_vert_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose)
{

  moab::EntityHandle vertex = verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert(vertex, bump_dist_x, bump_dist_y, bump_dist_z);
  MB_CHK_SET_ERR(result, "could not move single vert");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);

  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);

  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = rand_vert_move(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::rand_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);

  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1 = (verts.back()-index);
  moab::EntityHandle vertex2 = (verts.back()-index-1);

  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::rand_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);

  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1 = (verts.back()-index);
  moab::EntityHandle vertex2 = (verts.back()-index-1);

  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = rand_vert_move(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::rand_vert_bump(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex = verts.back();
  //move the desired vertex by the allotted distance
  result = rand_vert_move(vertex, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move single vert");

  return moab::MB_SUCCESS;
}

/// moves the third to last and the last vertices in the model the same distance in x, y, and z
moab::ErrorCode MakeWatertightTest::adjplone_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);

  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::adjplone_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);

  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = rand_vert_move(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::nonadj_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);

  //move the desired verticies by the allotted distance(s)
  result = move_vert(vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert(vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::nonadj_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);

  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = rand_vert_move(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightTest::write_mod_file(std::string filename)
{

  std::string output_filename = filename + "_mod.h5m";
  //write file
  result = MBI()->write_mesh(output_filename.c_str());
  MB_CHK_SET_ERR(result,"could not write the mesh to output_filename");

  return moab::MB_SUCCESS;
}

bool MakeWatertightTest::seal_and_check(moab::EntityHandle input_set, double facet_tolerance, bool verbose)
{
  // seal the model using make_watertight
  result = mw->make_mesh_watertight(input_set, facet_tolerance, verbose);
  MB_CHK_SET_ERR(result, "could not make the mesh watertight");

  // Check to see if make_watertight fixed the model by topology
  bool sealed, test = true, topo_test = true;
  result = cw->check_mesh_for_watertightness(input_set, facet_tolerance, sealed, test, verbose, topo_test);
  MB_CHK_SET_ERR(result, "could not check model for watertightness");
  if(!sealed) return sealed;

  // Check to see if make_watertight fixed the model via proximity tolerance
  topo_test = false;
  result = cw->check_mesh_for_watertightness(input_set, facet_tolerance, sealed, test, verbose, topo_test);
  MB_CHK_SET_ERR(result, "could not check model for watertightness");

  return sealed;
}


/////////////////////////////////////////////
//        CYLINDER TEST FUNCTIONS          //
/////////////////////////////////////////////

/// moves a vertex along the rim of the cylinder in the theta direction a distance equal to the faceting_tolerance
moab::ErrorCode MakeWatertightCylinderTest::move_vert_theta(moab::EntityHandle vertex, double tolerance, bool verbose)
{


  //get vertex coordinates
  double coords[3];
  result = MBI()->get_coords(&vertex, 1, coords);

  // determine radius
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));

  // get the current theta value
  // need both because of the oddness/evenness of the inverse functions
  double theta_x = acos(coords[0]/radius);
  double theta_y = asin(coords[1]/radius);

  // set the vertex bump distance
  double dtheta = tolerance/(radius);

  if(verbose) {
    //original coordinates
    std::cout << std::endl << "Original Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  // create new x and y values
  coords[0] = radius*cos(theta_x+dtheta);
  coords[1] = radius*sin(theta_y+dtheta);

  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  //set new vertex coordinates
  result = MBI()->set_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not set the vertex coordinates");

  return moab::MB_SUCCESS;
}

/// moves the vertex in R some distance less than tol
moab::ErrorCode MakeWatertightCylinderTest::move_vert_R(moab::EntityHandle vertex, double tol, bool verbose)
{


  //get vertex coordinates
  double coords[3];
  result = MBI()->get_coords(&vertex, 1, coords);

  if(verbose) {
    //original coordinates
    std::cout << std::endl << "Original Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));
  //get unit vector in x-y plane
  coords[0]/=radius;
  coords[1]/=radius;

  //alter radius to new value of radius+tol
  radius-=tol;
  coords[0]*=radius;
  coords[1]*=radius;

  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }

  //set new vertex coordinates
  result = MBI()->set_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not set the vertex coordinates");

  return moab::MB_SUCCESS;
}

/// bumps the last vertex in the cylinder model in the R direction
moab::ErrorCode MakeWatertightCylinderTest::single_vert_bump_R(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert_R(vertex, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move single vert");

  return moab::MB_SUCCESS;
}


/// selects a random pair of verticies and moves them along theta a distance less than the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::rand_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose)
{


  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);

  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1=(verts.back()-index);
  moab::EntityHandle vertex2=(verts.back()-index-1);

  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_theta(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

// FOR CYLINDER TESTING ONLY
/// moves the last vertex in the model along the curve of the cylinder some distance bump distance theta
moab::ErrorCode MakeWatertightCylinderTest::theta_vert_bump(moab::Range verts, double bump_dist_theta, bool verbose)
{


  //get vertex coordinates
  double coords[3];
  moab::EntityHandle vertex = verts.back();
  result = MBI()->get_coords(&vertex, 1, coords);

  // determine radius
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));

  // get the current theta value
  double theta = asin(coords[1]/radius);

  // set the vertex bump distance
  double dtheta = 0.5*bump_dist_theta/(radius);

  if(verbose) {
    //original coordinates
    std::cout << std::endl << "Original Coordinates HEre I am" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  // create new x and y values
  coords[0] = radius*cos(theta+dtheta);
  coords[1] = radius*sin(theta+dtheta);
  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  //write new coordinates to the mesh
  // might not be necesarry any longer as we move to doing tests on a moab-instance basis
  result = MBI()->set_coords(&vertex, 1, coords);
  MB_CHK_SET_ERR(result, "could not set the vertex coordinates");
  // alter output filename

  return moab::MB_SUCCESS;
}

/// moves two adjacent vertices along theta a distance equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::locked_pair_bump_theta(moab::Range verts, double tolerance, bool verbose)
{


  //get vertex coordinates
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = verts.back()-1;

  result = move_vert_theta(vertex1, tolerance, verbose);
  MB_CHK_SET_ERR(result,"could not move vertex1 along theta");
  result = move_vert_theta(vertex2, tolerance, verbose);
  MB_CHK_SET_ERR(result,"could not move vertex1 along theta");

  return moab::MB_SUCCESS;
}


/// moves the third to last and the last verticies in the model in theta the same distance along theta equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::adjplone_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);

  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_theta(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertightCylinderTest::locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);

  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_R(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

/// selects random verticies from verts and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::rand_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);

  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1=(verts.back()-index);
  moab::EntityHandle vertex2=(verts.back()-index-1);
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_R(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

/// selects a the last vertex and third to last vertex in the model and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::adjplone_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose)
{

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);

  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_R(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

/// moves the last vertex in the model and a randomly selected, non-adjacent vertex and moves them both in R a distance equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::nonadj_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_R(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}

/// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
moab::ErrorCode MakeWatertightCylinderTest::nonadj_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose)
{

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);

  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex1");
  result = move_vert_theta(vertex2, facet_tol, verbose);
  MB_CHK_SET_ERR(result, "could not move vertex2");

  return moab::MB_SUCCESS;
}
