
#include <iostream>
#include <vector>
#include <assert.h>
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/ScdInterface.hpp"

namespace moab {
  
class SignedDistanceField {

public:
  // Constructor
  SignedDistanceField(double x_min, double y_min, double z_min, double step_size, int num_x_points, int num_y_points, int num_z_points);
  
private:
  // data storage of signed distance values
  std::vector<double> signed_distance_values;
  // reference location 
  double lower_left_corner[3];
  // size of data in each dimension
  int dims[3];
  // distance between data points
  double step_size;
  // signed distance value tag name
  const std::string sdf_tag_name = "SIGNED_DISTANCE_FIELD";
  
public:
  // method for getting the lower left corner of the signed distance field
  inline void get_llc(double &x, double &y, double &z) {
    x = lower_left_corner[0]; y = lower_left_corner[1]; z = lower_left_corner[2];
  }

  // method for getting the dimensions of the field
  inline void get_dims(int &x_steps, int &y_steps, int &z_steps) {
    x_steps = dims[0]; y_steps = dims[1]; z_steps = dims[2];    
  }

  // get the step size between points for this structured field
  inline void get_step(double &step) { step = step_size; }

  // sets the data for the signed distance field
  inline void set_data(std::vector<double> data) {
    assert(data.size() == dims[0]*dims[1]*dims[2]);
    signed_distance_values = data;
  }

  // get the coordinates of a point at the i,j,k index
  inline CartVect get_coords(int i, int j, int k){
    return CartVect(lower_left_corner)+((i)*step_size,(j)*step_size,(k)*step_size);
  }

  // get the data at the i,j,k index
  inline double get_data(int i, int j, int k) {
       return signed_distance_values[(i)+(j)*dims[0]+(k)*dims[0]*dims[1]];
  }

  //returns the trilinear interpolated signed distance value for the point
  double find_sdv(const double pnt[3]);

  // gets the mesh element for this point
  void get_element_indices_for_pnt(const double pnt[3], int &i, int &j, int &k);
  
  // gets the mesh element for this point and indicates whether or not the point
  // is inside of the data structure
  void get_element_indices_for_pnt(const double pnt[3], int &i, int &j, int &k, bool &outside);

  // generates an scdBox for this field and tags the vertices with its data
  ErrorCode create_scdBox(ScdBox* sdfBox, Interface* mbi = NULL);
  
};

}
