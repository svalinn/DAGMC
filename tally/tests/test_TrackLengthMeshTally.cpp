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
    virtual void SetUp()
    {
      std::multimap<std::string, std::string> options;
      // Minimum required member settings
      input.energy_bin_bounds.push_back(0.0);
      input.energy_bin_bounds.push_back(10.0);
      input.tally_id = 4;

    }

    // deallocate memory resources
    virtual void TearDown()
    {
      delete mesh_tally;
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
TEST_F(TrackLengthMeshTallyTest, ComputeScore1Ray)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "unstructured_mesh.h5m"));
  mesh_tally = Tally::create_tally(input);
  EXPECT_TRUE(mesh_tally != NULL);
  EXPECT_EQ("unstr_track", mesh_tally->get_tally_type());

  TallyEvent event;
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 5.0;
  event.position = {0.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

 
  mesh_tally->compute_score(event);
  
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    {
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
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 2.5;
  event.position = {0.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();

  // compute score to halfway
  // now do the rest

  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 2.5;
  event.position = {2.5,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    {
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
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 5.0;
  event.position = {0.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();

  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    {
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
  // 2nd track 1-2 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {1.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 3rd track 2-3 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {2.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 4th track 3-4 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {3.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 5th track
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {4.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    {
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
  // 1st track 0-1 cm
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

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 2nd track 1-2 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {1.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 3rd track 2-3 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {2.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 4th track 3-4 cm
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {3.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();
  // 5th track
  event.type = TallyEvent::TRACK; // track length tally
  event.particle = 1;
  event.current_cell = 1;
  event.track_length = 1.0;
  event.position = {4.0,0.0,0.0};
  event.direction = {1.0,0.0,0.0};
  event.total_cross_section = 1.0;
  event.particle_energy = 5.0;
  event.particle_weight = 1.0;
  event.multipliers.push_back(1.0);

  mesh_tally->compute_score(event); 
  mesh_tally->end_history();

  // now figure out the track
  TallyData data = mesh_tally->getTallyData();
  int length;
  double *track_data = data.TallyData::get_tally_data(length); // << std::endl;

  double total = 0.0;
  for ( int i = 0 ; i < length ; i++ )
    {
      total += track_data[i];
    }

  EXPECT_DOUBLE_EQ(total, 3.0);
}
