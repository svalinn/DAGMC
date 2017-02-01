

#include "sdf.hpp"
#include "moab/CartVect.hpp"

namespace moab {

// Constructor
//  SignedDistanceField::SignedDistanceField() {};

  double SignedDistanceField::find_sdv(const double pnt[3]) {

    int i,j,k;
    get_ijk(pnt,i,j,k);
    
    
  };

  
  void SignedDistanceField::get_ijk(const double pnt[3], int &i, int &j, int &k) {
    CartVect this_pnt(pnt);
    CartVect llc(lower_left_corner);
    CartVect vec = this_pnt-llc;
    i = vec[0]/step_size;
    j = vec[1]/step_size;
    k = vec[2]/step_size;
  };


  void SignedDistanceField::get_ijk(const double pnt[3], int &i, int &j, int &k, bool &outside) {
    get_ijk(pnt, i, j, k);
    if( i < 0 || j < 0 || k < 0 || i > field_dims[0]-2 || j > field_dims[1]-2 || k > field_dims[2]-2) {
      outside = true;
     }
  };

}
