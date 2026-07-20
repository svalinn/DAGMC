#include "DagSolidMagneticField.hh"

#include <iostream>
#include <fstream>

WireMagneticField::WireMagneticField(const G4double radius,
				     const G4double mu,
				     const G4double current) {
  wireRadius = radius;
  wireCurrent = current;
  wireMu = mu;
}

WireMagneticField::~WireMagneticField() {
}

// set the magnetic field value given the position
void WireMagneticField::GetFieldValue(const G4double Point[4],
				 double *field) {
  // radius squared
  G4double r2 = Point[0]*Point[0] + Point[1]*Point[1];
  G4double B = 0;
  
  // calculate the field strength
  if (r2 < wireRadius) {
    B = wireMu*std::sqrt(r2)*current/2*pi*wireRadius*wireRadius;
  } else if ( r2 == wireRadius*wireRadius ) {
    B  = wireMu*current/(twopi*wireRadius);
  } else {
    B = (CLHEP::mu0*current)/(twopi*std::sqrt(r2));
  }
  
  // calculate field direction, always tangent to circle
  // if current +ve anticlockwise, if -ve clockwise
  G4double i = 0;
  G4double j = 0;
  if ( wireCurrent > 0 ) {
    i = -Point[1];
    j = Point[0];
  } else {
    i = Point[1];
    j = -Point[0];
  }

  // dont forget to normalise vector
  Field[0] = B*i/std::sqrt(r2);
  Field[1] = B*j/std::sqrt(r2);
  Field[2] = 0;
  return;
}

MagneticField::MagneticField(const Sampling samplingMode) {
  // set the sampling mode
  mode = samplingMode;
}

MagneticField::~MagneticField() {
}

void MagneticField::LoadFile(const std::string filename) {
  // open the file
  std::ifstream input;
  input.open(filename, std::ifstream::in);

  // if the input file can't be found
  if(!input.good()) {
    std::cout << "Failed to open file" << std::endl;
    exit(10);
  }

  double r,z,br,bz,btheta;
  // while we can read
  while ( input >> r >> z >> br >> bz >> btheta ) {
    // push_back to the storage array 
    // note using G4units so that data in natively in mm;
    radial_pos.push_back(r*m);
    vertical_pos.push_back(z*m);
    // note the order br bz btheta
    b_field.push_back({br,bz,btheta});
  }

  // process the field
  return ProcessField();
}

// return the magnetic field components
void MagneticField::GetFieldValue(const G4double Point[4],
  double *field) const {
  // interpolate the field given the point
  InterpolateField(Point[0],Point[1],Point[2],
    field[0],field[1],field[2]);
  return; 
}

// process the field into a form which 
// is usable for interpolation
void MagneticField::ProcessField() {
  // get the unique r and z values
  std::set<double> unique_r(radial_pos.begin(), radial_pos.end());
  std::set<double> unique_z(vertical_pos.begin(), vertical_pos.end());

  // set the size of the data
  n_r = unique_r.size();
  n_z = unique_z.size();
  
  // set the size of the vectors
  // copy in the set data
  std::copy(unique_r.begin(), unique_r.end(),
	    std::back_inserter(r_point));
  std::copy(unique_z.begin(), unique_z.end(),
	    std::back_inserter(z_point));

  // sort the vectors
  std::sort(r_point.begin(), r_point.end());
  std::sort(z_point.begin(), z_point.end());

  return;
}

// given the radial position return its index
int MagneticField::FindRadialIndex(const double r) const {
  int radialIndex = 0;
  // given r find the owning radial index
  const std::vector<double>::const_iterator it_r = std::upper_bound(r_point.begin(),
    r_point.end(), r);

  // if we are pointing outside return the last element
  if (it_r - r_point.begin() > r_point.size() - 1)
    radialIndex = r_point.size() - 1;
  else if (it_r - r_point.begin() <= 0 )
    radialIndex = 0;
  else
    radialIndex = it_r - r_point.begin();

  return radialIndex;
}

// given the vertical position return its index
int MagneticField::FindVerticalIndex(const double z)  const {
  int verticalIndex = 0;
  
  // given z find the owning vertical index
  std::vector<double>::const_iterator it_z = std::lower_bound(z_point.begin(),
    z_point.end(), z);

  // if we are pointing outside return the last element
  if (it_z - z_point.begin() > z_point.size() - 1)
    verticalIndex = z_point.size() - 1;
  else if (it_z - z_point.begin() <= 0 )
     verticalIndex = 0;
  else
     verticalIndex = it_z - z_point.begin();

  return verticalIndex;
}

// find out which radial and vertical bins are nearest
void MagneticField::FindRadialandVerticalIndices(const double r, const double z,
				 int &r_ind, int &z_ind) const {

  r_ind = FindRadialIndex(r);
  z_ind = FindVerticalIndex(z);
  return; 
}

// array index given radial and vertical index
int MagneticField::GlobalIndexByRadialandVerticalIndex(const int r_index, 
  const int z_index) const {
  int index = (r_index)*n_z + (z_index) ;
  return index;
}

// array index given radial and vertical index
int MagneticField::GlobalIndexByRadialandVerticalPosition(const double r, 
  const double z) const {
  int r_index = FindRadialIndex(r);
  int z_index = FindVerticalIndex(z);
  return GlobalIndexByRadialandVerticalIndex(1,1);
}

std::array<double,2> MagneticField::LookupRPtsByIndex(const int r_index) const {
  std::array<double,2> r_pts;
  // set the r coordinate values
  // if radial value returned is the 1st bin
  // resort to 2d interpolation
  if ( r_index == 0 ) {
    r_pts[0] = r_point[r_index];
    r_pts[1] = r_point[r_index];
  } else if ( r_index > n_r ) {
    // if the r_index is greater than the
    // number of indices
    r_pts[0] = r_point[r_index];
    r_pts[1] = r_point[r_index];
  } else {
    // standard case
    r_pts[0] = r_point[r_index-1];
    r_pts[1] = r_point[r_index];
  }
  return r_pts;
}

std::array<double,2> MagneticField::LookupRPtsByValue(const double r) const {
  int r_index = FindRadialIndex(r);
  return LookupRPtsByIndex(r_index);
}

// lookup the values of the z coordinates bounded by the
// z index
std::array<double,2> MagneticField::LookupZPtsByIndex(const int z_index) const {
  std::array<double,2> z_pts;
  // set the z coordinate values
  // if radial value returned is the 1st bin
  // resort to 2d interpolation
  if ( z_index == 0 ) {
    z_pts[0] = z_point[z_index];
    z_pts[1] = z_point[z_index];
  } else if ( z_index > n_z ) {
    // if the r_index is greater than the
    // number of indices
    z_pts[0] = z_point[z_index];
    z_pts[1] = z_point[z_index];
  } else {
    // standard case
    z_pts[0] = z_point[z_index-1];
    z_pts[1] = z_point[z_index];
  }
  return z_pts;
}

std::array<double,2> MagneticField::LookupZPtsByValue(const double z) const {
  int z_index = FindVerticalIndex(z);
  return LookupZPtsByIndex(z_index);
}

std::array<std::array<double,4>, 3> MagneticField::GetFourFieldValuesByIndex(const int r_index,
									    const int z_index) const {
  std::array<std::array<double,4>, 3> values;
  // interpolation array
  
  // order of values needs to be: p[0,0], p[0,1], p[1,0], p[1,1];
  for ( int dir = 0 ; dir < 3 ; dir++ ) {
    int index = r_index*n_z + z_index; // q11
    values[dir][0] = b_field[index][dir];
    index = (r_index-1)*n_z + z_index; // q12
    values[dir][1] = b_field[index][dir];
    index = (r_index*n_z) + z_index-1; // q21
    values[dir][2] = b_field[index][dir];
    index = (r_index*n_z) + z_index; // q22
    values[dir][3] = b_field[index][dir];
  }

  return values;
}

std::array<std::array<double,4>, 3> MagneticField::GetFourFieldValuesByValue(const double r,
  const double z) const {

  int r_index = FindRadialIndex(r);
  int z_index = FindVerticalIndex(z);

  return GetFourFieldValuesByIndex(r_index,z_index);
}

// get the magnetic field by lookup to nearest point
void MagneticField::FieldByNearest(const double r,
  const double z,
  double &field_r,
  double &field_z,
  double &field_t) const {
  
  // find the global index
  int index = GlobalIndexByRadialandVerticalPosition(r,z);

  // given the index return the field  
  field_r = b_field[index][0];
  field_z = b_field[index][1];
  field_t = b_field[index][2];
  return;
}

// get the magnetic field at r,z by interpolating using the 4 nearest
// values
void MagneticField::FieldByInterpolation(const double r,
						  const double z,
						  double &field_x,
						  double &field_y,
						  double &field_z) const {

  // get the radial and vertical bounding points
  std::array<double,2> r_pts = LookupRPtsByValue(r);
  std::array<double,2> z_pts = LookupZPtsByValue(z);

  // get the 4 values nearest us
  std::array<std::array<double,4>, 3> values = GetFourFieldValuesByValue(r,z);

  // b field 
  field_x = BilinearInterpolation(r_pts,z_pts,
					  values[0],r,z);
  field_y = BilinearInterpolation(r_pts,z_pts,
					  values[1],r,z);
  field_z = BilinearInterpolation(r_pts,z_pts,
					  values[2],r,z);

  return;
}

// given some position, iterpolate the field from the nearest points
void MagneticField::InterpolateField(const double x,
                                     const double y,
                                     const double z,
                                     double &field_x,
                                     double &field_y,
                                     double &field_z) const {
  // determine the radial value
  double r = std::sqrt(x*x + y*y);

  // early return ensure r and z are within bounds else return 0
  if ( r < r_point.front() || r > r_point.back() ) {
    field_x = 0.;
    field_y = 0.;
    field_z = 0.;
    return;
  }

  // early return ensure r and z are within bounds else return 0
  if ( z < z_point.front() || z > z_point.back() ) {
    field_x = 0.;
    field_y = 0.;
    field_z = 0.;
    return;
  }


  double field_r, field_t; // radial field and toroidal dal
 
  // note passing field_r and field_t below so we can rotate the field
  if( mode == NEAREST ) {
    FieldByNearest(r, z, field_r, field_z, field_t);
  }
  // 
  if( mode == BILINEAR_INTERPOLATION ) {
    FieldByInterpolation(r, z, field_r, field_z, field_t);
  }
  
  // need to rotate the b-field appropraitely
  // determine theta angle
  double theta = std::atan2(y,x);
  //
  field_x = field_r*std::cos(theta) - field_t*std::sin(theta);
  field_y = field_r*std::sin(theta) + field_t*std::cos(theta);

  return;
}



