
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/Types.hpp"

#include "gtest/gtest.h"
#include "MakeWatertight.hpp"
#include "CheckWatertight.hpp"

moab::Interface *MBI();

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class MakeWatertightTest : public ::testing::Test
{
 protected:
  virtual void SetUp();
  void reload_mesh();
  virtual void TearDown();
  virtual void setFilename() {};

  // make sure the expected number of entities with dimension are present
  moab::ErrorCode check_num_ents(int ent_dimension, int expected_num);

  // moves the vertex by dx, dy, dz
  moab::ErrorCode move_vert(moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose = false);

  // moves the vertex by a random dx, dy, dz
  moab::ErrorCode rand_vert_move(moab::EntityHandle vertex, double tol, bool verbose = false);

  /// bumps the last vertex in the model by the x,y,z values given to the problem
  moab::ErrorCode single_vert_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false);

  /// bumps the last vertex in the cylinder model in the R direction
  moab::ErrorCode locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false);

  /// moves the last two verticies in the model in the same direction some distance less than the faceting tolerance
  moab::ErrorCode locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose = false);

  /// selects a random pair of adjacent verticies and bumps them in x, y, and z
  moab::ErrorCode rand_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false);

  /// selects a random pair of verticies from verts and moves them in random directions some distance less than the faceting tolerance
  moab::ErrorCode rand_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose = false);

  /// moves the last vertex in the model in a random direction by a distance less than the faceting tolerance
  moab::ErrorCode rand_vert_bump(moab::Range verts, double facet_tol, bool verbose = false);

  /// moves the third to last and the last vertices in the model the same distance in x, y, and z
  moab::ErrorCode adjplone_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false);

  /// moves the third to last and the last verticies in the model in rand directions some distance less than the facet_tolerance
  moab::ErrorCode adjplone_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose = false);

  /// selects a random pair of adjacent verticies and bumps them the same distance in x, y, and z
  moab::ErrorCode nonadj_locked_pair_bump(moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false);

  /// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
  moab::ErrorCode nonadj_locked_pair_bump_rand(moab::Range verts, double facet_tol, bool verbose = false);

  /// appends "_mod" to the original file name and writes to a new .h5m
  moab::ErrorCode write_mod_file(std::string filename);

  /// runs the make_watertight algorithm and checks the watertightness of the model
  bool seal_and_check(moab::EntityHandle input_set, double facet_tolerance, bool verbose = false);

 protected:
  std::string filename;
  MakeWatertight* mw;
  CheckWatertight* cw;
  moab::ErrorCode result;
  moab::EntityHandle input_fileset;
  moab::Range verts;
  double facet_tol;
};

// Rename of the general test class
class MakeWatertightConeTest : public MakeWatertightTest
{
 protected:
  // set test file name
  virtual void setFilename() {
    filename = "cones.h5m";
  };
};


// FOR CYLINDER TESTING ONLY
// Extension of the general test class to include cylinder-specific tests
class MakeWatertightCylinderTest : public MakeWatertightTest
{

 protected:
  /// set test file name
  virtual void setFilename() {
    filename = "cyl.h5m";
  };

  /// selects random verticies from verts and moves them in R a distance equal to the faceting tolerance
  moab::ErrorCode rand_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false);

  /// moves a vertex along the rim of the cylinder in the theta direction a distance equal to the faceting_tolerance
  moab::ErrorCode move_vert_theta(moab::EntityHandle vertex, double tolerance, bool verbose = false);

  /// moves the vertex in R some distance less than tol
  moab::ErrorCode move_vert_R(moab::EntityHandle vertex, double tol, bool verbose = false);

  /// bumps the last vertex in the cylinder model in the R direction
  moab::ErrorCode single_vert_bump_R(moab::Range verts, double facet_tol, bool verbose = false);

  /// selects a random pair of verticies and moves them along theta a distance less than the faceting tolerance
  moab::ErrorCode rand_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false);

  /// moves the last vertex in the model along the curve of the cylinder some distance bump distance theta
  moab::ErrorCode theta_vert_bump(moab::Range verts, double bump_dist_theta, bool verbose = false);

  /// moves two adjacent vertices along theta a distance equal to the faceting tolerance
  moab::ErrorCode locked_pair_bump_theta(moab::Range verts, double tolerance, bool verbose = false);

  /// moves the third to last and the last verticies in the model in theta the same distance along theta equal to the faceting tolerance
  moab::ErrorCode adjplone_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false);

  //moves a locked pair of vertices in the radial direction
  moab::ErrorCode locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false);

  /// selects a the last vertex and third to last vertex in the model and moves them in R a distance equal to the faceting tolerance
  moab::ErrorCode adjplone_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false);

  /// moves the last vertex in the model and a randomly selected, non-adjacent vertex and moves them both in R a distance equal to the faceting tolerance
  moab::ErrorCode nonadj_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false);

  /// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
  moab::ErrorCode nonadj_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false);

};

// Rename of the general test class
class MakeWatertightNoCurveSphereTest : public MakeWatertightTest
{
 protected:
  // set test file name
  virtual void setFilename() {
    filename = "no_curve_sphere.h5m";
  };

  moab::ErrorCode sphere_deletion_test(moab::EntityHandle input_set, double facet_tolerance, bool verbose = false);
};


