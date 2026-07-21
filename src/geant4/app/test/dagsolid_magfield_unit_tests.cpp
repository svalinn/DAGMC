//  DagSolid_test.cpp
#include <gtest/gtest.h>

#include <cassert>
#include <cmath>
#include <iostream>

#include "DagSolidMagneticField.hh"

/*
 * Notes: DAGMC geometries implicitly have units of cm's due to the history of
 * development with centimetre based physics engines. Please note that Geant4
 * uses units of mm internally (and can course be changed) for the purposes of
 * this test note that mm's are being passed into the various function calls,
 * e.g. a call to the PointInVolume(10,0,0) refers to (10mm,0,0).
 */

class DagSolidMagneticFieldTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    magneticField = new MagneticField();
    magneticField->LoadFile("test_field.txt");
  }

 protected:
  MagneticField* magneticField;
};

class DagSolidWireMagneticFieldTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    // note wire radius is 10mm, mu is 1, current is 1e6 A
    magneticField = new WireMagneticField(10,1,1e6);
  }

 protected:
  MagneticField* magneticField;
};

/*
 * Empty common setup function
 */
TEST_F(DagSolidMagneticWireFieldTest, SetUp) {}

//
TEST_F(DagSolidMagneticFieldTest, sample_test) {
  G4double point[4] = {5.0, 0.0, 0.0, 0.0};  // mm
  G4double field[3];
  magneticField->GetFieldValue(point,field);

  G4double b = 1*1e6/2*CLHEP::pi*5.0;
  // expect field x to be 1.0
  G4double field_mag = std::sqrt(field[0]*field[0] + field[1]*field[1] + field[2]*field[2]);
  EXPECT_NEAR(b, field_mag, 1e-6);
  EXPECT_EQ(0.0, field[0]);
  EXPECT_NEAR(b, field[1],1e-6);
  EXPECT_EQ(0.0, field[2]);

  return;
}
