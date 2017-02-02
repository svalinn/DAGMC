
#include <iostream>
#include <vector>
#include <assert.h>
#include "moab/CartVect.hpp"

namespace moab {
  
class SignedDistanceField {

public:
  //  SignedDistanceField();
  
  // data storage of signed distance values
  std::vector<double> signed_distance_values;
  // reference location 
  double lower_left_corner[3];
  // size of data in each dimension
  int dims[3];
  // distance between data points
  double step_size;

  inline void get_dims(int &x_steps, int &y_steps, int &z_steps) {
    x_steps = dims[0]; y_steps = dims[1]; z_steps = dims[2];    
  }
  
  inline void set_dims(int x_steps, int y_steps, int z_steps) {
    dims[0] = x_steps; dims[1] = y_steps; dims[2] = z_steps;
  }

  inline void set_corner(double x_min, double y_min, double z_min) {
    lower_left_corner[0] = x_min;
    lower_left_corner[1] = y_min;
    lower_left_corner[2] = z_min;
  };
  
  inline void set_step(double step) { step_size = step; }

  inline void set_data(std::vector<double> data) {
    assert(data.size() == dims[0]*dims[1]*dims[2]);
    signed_distance_values = data;
  }

  inline CartVect get_vert_coords(int i, int j, int k){
    return CartVect(lower_left_corner)+((i)*step_size,(j)*step_size,(k)*step_size);
  }

  inline double get_data_ijk(int i, int j, int k) {
       return signed_distance_values[(i)+(j)*dims[0]+(k)*dims[0]*dims[1]];
  }

  double find_sdv(const double pnt[3]);

  // gets the mesh element for this point
  void get_element_ijk(const double pnt[3], int &i, int &j, int &k);
  void get_element_ijk(const double pnt[3], int &i, int &j, int &k, bool &outside);


};

}
