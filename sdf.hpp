
#include <iostream>
#include <vector>
#include <assert.h>

namespace moab {
  
class SignedDistanceField {

public:
  //  SignedDistanceField();
  
  // data storage of signed distance values
  std::vector<double> signed_distance_values;
  // reference location 
  double lower_left_corner[3];
  // size of data in each dimension
  int field_dims[3];
  // distance between data points
  double step_size;

  inline void set_dims(int x_steps, int y_steps, int z_steps) {
    field_dims[0] = x_steps; field_dims[1] = y_steps; field_dims[2] = z_steps;
  }

  inline void set_step(double step) { step_size = step; }

  inline void set_data(std::vector<double> data) {
    assert(data.size() == field_dims[0]*field_dims[1]*field_dims[2]);
    signed_distance_values = data;
  }

  double find_sdv(const double pnt[3]);

  // gets the mesh element for this point
  void get_ijk(const double pnt[3], int &i, int &j, int &k);
  void get_ijk(const double pnt[3], int &i, int &j, int &k, bool &outside);

};

}
