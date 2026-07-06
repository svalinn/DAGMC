#ifndef DAG_SOLID_MAGNETIC_FIELD_HH
#define DAG_SOLID_MAGNETIC_FIELD_HH

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <array>
#include <vector>
#include <string>
#include <set>

double LinearInterpolation(const double x1, const double x2,
  const double f1, const double f2, const double x) {
  
  if (x2 == x1) {
    std::cout << "x1: " << x1 << " x2: " << x2 << std::endl;
    std::cout << "f1: " << f1 << " f2: " << f2 << std::endl;
    throw std::invalid_argument("x0 and x1 cannot be the same (division by zero).");
  } 
  return f1 + ( (f2 - f1) / (x2 - x1) ) * (x - x1);
}

double BilinearInterpolation(const std::array<double,2> x_stencil,
  const std::array<double,2> y_stencil,
  const std::array<double,4> f_values,
  const double x,
  const double y) {

  double x1 = x_stencil[0];
  double x2 = x_stencil[1];
  double y1 = y_stencil[0];
  double y2 = y_stencil[1];
  
  double denom = (x2 - x1) * (y1 - y2);
  if (denom == 0.0) {
    if (y1 - y2 == 0.0)
      return LinearInterpolation(x1, x2, f_values[0], f_values[1], x);
    else
      return LinearInterpolation(y1, y2, f_values[2], f_values[3], y);
  }

  // order of f_values is: p[0,0], p[0,1], p[0,1], p[1,1];
  double q11 = f_values[0];
  double q12 = f_values[1];
  double q21 = f_values[2];
  double q22 = f_values[3];

  double fxy1 = ((x2 - x) * q11 + (x - x1) * q21) / (x2 - x1);
  double fxy2 = ((x2 - x) * q12 + (x - x1) * q22) / (x2 - x1);
  double fxy = ((y2 - y) * fxy1 + (y - y1) * fxy2) / (y2 - y1);

  return fxy;
}

enum Sampling { NEAREST, BILINEAR_INTERPOLATION };

// Class to provide a simple magnetic field in and around
// a wire following amperes law, it is assumed that the field
// is centred on 0,0 in the x,y plane and that the wire points along
// the z axis. 
class WireMagneticField : public G4MagneticField {
  public:
  // Constructor
  WireMagneticField(const G4double radius = 5*cm,
		    const G4double mu = 1.256*e-6*m,
		    const G4double current = 1e6);

  // Destructor
  ~WireMagneticField();

  // main lookup function
  void GetFieldValue(const G4double Point[4],
		     double *field) const override;
  private:
  G4double wireRadius;
  G4double current;
  
}

class MagneticField : public G4MagneticField {
  public:
  // Constructor 
  MagneticField(Sampling samplingMode = BILINEAR_INTERPOLATION);

  // Destructor
  ~MagneticField() override;

  // main lookup function that Geant4 will call
  void GetFieldValue(const G4double Point[4],
                    double *field) const override;                     

  private:
  // load the file containing white space delimited
  // r, z, Br, Bz, Btoroial
  void LoadFile(std::string filename);

  // main function called by GetFieldValue
  // to interpolate the 3d B field as a function 
  // of Cartesian coordinates
  void InterpolateField(const double x,
                        const double y,
                        const double z,
                        double &field_x,
                        double &field_y,
                        double &field_z) const;

  // process the field data which has been initialised from
  // file. Should not be called before LoadFile                         
  void ProcessField();

  // get the magnetic field by lookup to nearest point
  void FieldByNearest(const double r,
    const double z,
    double &final_r,
    double &field_z,
    double &field_t) const;

  // perform a bilinear interpolation of the magnetic field
  // at r,z using the 4 nearest points
  void FieldByInterpolation(const double r, 
    const double z,
    double &field_x,
    double &field_y,
    double &field_z) const;

  // helper function for FieldByInterpolation to return the Br,Bz,Bt at given
  // r,z indices
  std::array<std::array<double,4>, 3> GetFourFieldValuesByIndex(const int r_index,
    const int z_index) const;

  // helper function for FieldByInterpolation to return the Br,Bz,Bt at given
  // r,z values
  std::array<std::array<double,4>, 3> GetFourFieldValuesByValue(const double r,
      const double z) const;

  // indexing functions
  int FindVerticalIndex(const double z) const;
  int FindRadialIndex(const double r) const;
  void FindRadialandVerticalIndices(const double r, const double z,
    int &r_ind, int &z_ind) const;
  int GlobalIndexByRadialandVerticalIndex(const int r_index, 
    const int z_index) const;
  int GlobalIndexByRadialandVerticalPosition(const double, 
    const double z) const;
  std::array<double,2> LookupRPtsByIndex(const int r_index) const;
  std::array<double,2> LookupRPtsByValue(const double r) const;
  std::array<double,2> LookupZPtsByIndex(const int z_index) const;
  std::array<double,2> LookupZPtsByValue(const double z) const;
  
  private:
  Sampling mode;
  int n_r; // number of r points
  int n_z; // number of z points
  std::vector<double> r_point; // the unique grid points in r
  std::vector<double> z_point; // the unique grid points in z

  public:
  std::vector<double> radial_pos; /// temp
  std::vector<double> vertical_pos; /// temp
  std::vector<std::array<double,3> > b_field;

};
#endif 
				
