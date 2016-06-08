
#include <gtest/gtest.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include <unistd.h>
#include <stdio.h>

#include "uwuw_preprocessor.hpp"

namespace
{

class UWUWTest : public ::testing::Test
{
 protected:

  UWUWTest() {}
  virtual ~UWUWTest() {}

  virtual void SetUp() {
  }

};


/*
 * Empty common setup function
 */
TEST_F(UWUWTest,SetUp)
{
}

/*
 * Test to make sure the name is only uppercased
 */
TEST_F(UWUWTest,name8)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "SpecialSteel";
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIALS");
  return;
}

// test to make sure the space is removed
TEST_F(UWUWTest,name8RemoveSpace)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Special Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIALS");
  return;
}

// test to make sure the , is removed
TEST_F(UWUWTest,name8RemoveCommaAndSpace)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Special, Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIALS");
  return;
}

// test to make sure the number is incremented
TEST_F(UWUWTest,name8IncrementID)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Special, Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIALS");
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIAL1");
  return;
}

// test to make sure the number is incremented
TEST_F(UWUWTest,name8lessthan8)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name),"  AIRDRY");
  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace
TEST_F(UWUWTest,name8lessthan8preserve)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name),"  AIRDRY");
  EXPECT_EQ(ncr->make_name_8bytes(name)," AIRDRY1");

  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace and the number is good
TEST_F(UWUWTest,name8lessthan8preserve100)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name),"  AIRDRY");

  for ( int i = 1 ; i < 10 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name)," AIRDRY"+str);
  }
  for ( int i = 10 ; i < 100 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name),"AIRDRY"+str);
  }
  for ( int i = 100 ; i < 999 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name),"AIRDR"+str);
  }

  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace and the number is good
TEST_F(UWUWTest,name8lessthan8preserve100Steel)
{
  // make new class
  name_concatenator *ncr = new name_concatenator();

  std::string name = "Special Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIALS");

  for ( int i = 1 ; i < 10 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIAL"+str);
  }
  for ( int i = 10 ; i < 100 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name),"SPECIA"+str);
  }
  for ( int i = 100 ; i < 999 ; i++ ) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name),"SPECI"+str);
  }

  return;
}


};
