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

/*
 * Empty common setup function
 */
TEST_F(DagSolidMagneticFieldTest, SetUp) {}

// Interpolation test
TEST_F(DagSolidMagneticFieldTest, interpolation_test) {
  const std::array<double,2> x_stencil = {0.0, 1000};
  const std::array<double,2> y_stencil = {-1000, 1000};
  const std::array<double,4> f_values = {1.0, 1.0, 1.0, 1.0};

  EXPECT_EQ(1.0, BilinearInterpolation(x_stencil, y_stencil, f_values, 0., 0.));
  EXPECT_EQ(1.0, BilinearInterpolation(x_stencil, y_stencil, f_values, 0., 1000.));
  EXPECT_EQ(1.0, BilinearInterpolation(x_stencil, y_stencil, f_values, 1000., 1000.));
  EXPECT_EQ(1.0, BilinearInterpolation(x_stencil, y_stencil, f_values, 1000., -1000.));

  const std::array<double,2> x_stencil2 = {0.0, 1000};
  const std::array<double,2> y_stencil2 = {-1000, 1000};
  const std::array<double,4> f_values2 = {1.0, 2.0, 3.0, 4.0};

  // check the 4 corners
  EXPECT_EQ(1.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 0., -1000.));
  EXPECT_EQ(2.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 0., 1000.));
  EXPECT_EQ(3.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 1000., -1000.));
  EXPECT_EQ(4.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 1000., 1000.));
  // check the middle
  EXPECT_EQ(2.5, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 500., 0.));
  // check middle left  
  EXPECT_EQ(1.5, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 0., 0.));
  // check middle right
  EXPECT_EQ(3.5, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 1000., 0.));
  // check middle top
  EXPECT_EQ(3.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 500., 1000.));
  // check middle bottom
  EXPECT_EQ(2.0, BilinearInterpolation(x_stencil2, y_stencil2, f_values2, 500., -1000.));

  return;
}

// 

//
TEST_F(DagSolidMagneticFieldTest, sample_test) {
  G4double point[4] = {500.0, 0.0, 0.0, 0.0};  // mm
  G4double field[3];
  magneticField->GetFieldValue(point,field);

  // expect field x to be 1.0
  EXPECT_EQ(1.0, field[0]);
  EXPECT_EQ(0.0, field[1]);
  EXPECT_EQ(0.0, field[2]);

  point[0] = 0.0;
  magneticField->GetFieldValue(point,field);
  // expect field x to be 1.0
  EXPECT_EQ(1.0, field[0]);
  EXPECT_EQ(0.0, field[1]);
  EXPECT_EQ(0.0, field[2]);

  point[0] = 1000.0;
  point[1] = 0.0;
  point[2] = 0.0;

  magneticField->GetFieldValue(point,field);
  // expect field x to be 1.0
  EXPECT_EQ(1.0, field[0]);
  EXPECT_EQ(0.0, field[1]);
  EXPECT_EQ(0.0, field[2]);

  return;
}
