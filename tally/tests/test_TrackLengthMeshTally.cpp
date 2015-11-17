// test_TrackLengthMeshTally.cpp
#include "gtest/gtest.h"

#include "moab/CartVect.hpp"

#include "../Tally.hpp"
#include "../TallyEvent.hpp"
#include "../TrackLengthMeshTally.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class TrackLengthMeshTallyTest : public ::testing::Test
{
 protected:
  // initialize variables for each test
  virtual void SetUp() {
    std::multimap<std::string, std::string> options;
    // Minimum required member settings
    input.energy_bin_bounds.push_back(0.0);
    input.energy_bin_bounds.push_back(10.0);
    input.tally_id = 4;

  }

  // deallocate memory resources
  virtual void TearDown() {
    delete mesh_tally;
  }

  void make_event(TallyEvent &event) {
    event.type = TallyEvent::TRACK; // track length tally
    event.particle = 1;
    event.current_cell = 1;
    event.track_length = 1.0;
    event.position = {0.0,0.0,0.0};
    event.direction = {1.0,0.0,0.0};
    event.total_cross_section = 1.0;
    event.particle_energy = 5.0;
    event.particle_weight = 1.0;
    event.multipliers.push_back(1.0);
  }

  void mod_event(TallyEvent &event, double x_coord, double x_dir,
                 double track_length) {

    event.position[0] = x_coord;
    event.direction[0] = x_dir;
    event.track_length = track_length;
  }

  void mod_event_3d(TallyEvent &event, double coord[3], double dir[3],
                    double track_length) {

    event.position[0] = coord[0];
    event.position[1] = coord[1];
    event.position[2] = coord[2];

    event.direction[0] = dir[0];
    event.direction[1] = dir[1];
    event.direction[2] = dir[2];

    event.track_length = track_length;
  }


 protected:
  // data needed for each test
  Tally *mesh_tally;
  TallyInput input;
};

// Tests Tally constructor for default number of energy bins
TEST_F(TrackLengthMeshTallyTest, DefaultNumEnergyBins)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);

  const TallyData& data = mesh_tally->getTallyData();
  EXPECT_EQ(1, data.get_num_energy_bins());
}

// Tests Tally constructor for non default number of energy bins
TEST_F(TrackLengthMeshTallyTest, NonDefaultNumEnergyBins)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));

  // increase the number of bins
  input.energy_bin_bounds.push_back(12.0);
  input.energy_bin_bounds.push_back(15.0);
  input.energy_bin_bounds.push_back(21.0);

  mesh_tally = Tally::create_tally(input);

  const TallyData& data = mesh_tally->getTallyData();

  EXPECT_EQ(5, data.get_num_energy_bins());
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, CreateTrackLengthMeshTally)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());
}

//---------------------------------------------------------------------------//
// A common theme in the following tests is that we introduce a particle at some
// point and follow it as it traverses a mesh. This is done by setting an
// initial direction and position and using mod_event_3d to move it specific
// distances. This does not update the position, so it must be updated manually
// for each track.
//
// At the end of each test, the total track length through the mesh is
// calculated. This is done by using get_tally_data to get the track length
// through every tet in the entire mesh, and then summing these track lengths.
//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScore1Ray)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;

  make_event(event);
  mod_event(event,1.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_EQ(total, event.track_length);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScore2Ray)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);
  mod_event(event,0.0,1.0,2.5);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // compute score to halfway
  // now do the rest

  mod_event(event,2.5,1.0,2.5);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 5.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScore1RaySplit)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstr_mesh_split.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);
  mod_event(event,0.0,1.0,5.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScore4RaySplit)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstr_mesh_split.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);
  // 2nd track 1-2 cm
  mod_event(event,1.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 3rd track 2-3 cm
  mod_event(event,2.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 4th track 3-4 cm
  mod_event(event,3.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 5th track
  mod_event(event,4.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 2.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScore5RaySplit)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstr_mesh_split.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;

  make_event(event);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 2nd track 1-2 cm
  mod_event(event,1.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 3rd track 2-3 cm
  mod_event(event,2.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 4th track 3-4 cm
  mod_event(event,3.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 5th track
  mod_event(event,4.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScorePointOnBoundary)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstr_mesh_split.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);
  // 2nd track 1-2 cm
  mod_event(event,1.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 3rd track 2-3 cm
  mod_event(event,2.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 4th track 3-4 cm
  mod_event(event,3.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();
  // 5th track
  mod_event(event,4.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 2.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ComputeScorePointOnAdjacentBoundary)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstr_mesh_split.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;

  // 3rd track 2-3 cm
  make_event(event);
  mod_event(event,2.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 3-4 cm
  mod_event(event,3.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track
  mod_event(event,4.0,1.0,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 2.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, PointOutside5of5RayReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  // first track 0-1 cm
  make_event(event);
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {-1.0,4.0,0.5};
  mod_event_3d(event,position,direction,2.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}


//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,4.0,0.5};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 4of5RayReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,4.0,0.5};

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();


  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 2.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 3of5RayReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,4.0,0.5};

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 2.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 2of5RayReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,4.0,0.5};

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 1.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, ReEntrantMeshRayOffCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-11 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,4.0,0.5};
  mod_event_3d(event,position,direction,11.0);

  mesh_tally->compute_score(event);
  mesh_tally->end_history();


  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOnCenter)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,0.0,0.0};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 11.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenterY)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,0.5,0.0};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 11.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenterY2)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,0.51,0.0};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenterY2Z)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,0.51,2.0};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenterY2Z2)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,5.0,5.0};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 3.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, 5of5RayReEntrantMeshRayOffCenterY3Z3)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "rune_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // first track 0-1 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {0.0,5.01,5.01};
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track 1-5 cm
  position[0]=1.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3th track 5-6 cm
  position[0]=5.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 6-10
  position[0]=6.0;
  mod_event_3d(event,position,direction,4.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 10-11
  position[0]=10.0;
  mod_event_3d(event,position,direction,1.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ ) {
    total += track_data[i];
  }

  EXPECT_DOUBLE_EQ(total, 0.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, HashtagMeshReentrantSeparateTracks_1)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "hashtag_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // 1st track -50 to -30 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {-50.0,18.0,2.0};
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track -30 to -20 cm
  position[0]=-30.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3rd track -20 to 20 cm
  position[0]=-20.0;
  mod_event_3d(event,position,direction,40.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 20 to 30 cm
  position[0]=20.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 30 to 50 cm
  position[0]=30.0;
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length);

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    total += track_data[i];

  EXPECT_DOUBLE_EQ(total, 20.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, HashtagMeshReentrantSeparateTracks_2)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "hashtag_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // 1st track -50 to -30 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {-50.0,-7.0,0.0};
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track -30 to -20 cm
  position[0]=-30.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3rd track -20 to 20 cm
  position[0]=-20.0;
  mod_event_3d(event,position,direction,40.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 20 to 30 cm
  position[0]=20.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 30 to 50 cm
  position[0]=30.0;
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length);

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    total += track_data[i];

  EXPECT_DOUBLE_EQ(total, 20.0);
}

//---------------------------------------------------------------------------//
TEST_F(TrackLengthMeshTallyTest, HashtagMeshReentrantSeparateTracks_3)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "hashtag_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  make_event(event);

  // 1st track -50 to -30 cm
  double direction[3]= {1.0,0.0,0.0};
  double position[3]= {-50.0,11.0,3.0};
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 2nd track -30 to -20 cm
  position[0]=-30.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 3rd track -20 to 20 cm
  position[0]=-20.0;
  mod_event_3d(event,position,direction,40.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 4th track 20 to 30 cm
  position[0]=20.0;
  mod_event_3d(event,position,direction,10.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // 5th track 30 to 50 cm
  position[0]=30.0;
  mod_event_3d(event,position,direction,20.0);
  mesh_tally->compute_score(event);
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length);

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    total += track_data[i];

  EXPECT_DOUBLE_EQ(total, 20.0);
}
