

#include "sdf.hpp"
#include "moab/CartVect.hpp"

namespace moab {

// Constructor
//  SignedDistanceField::SignedDistanceField() {};

  double SignedDistanceField::find_sdv(const double pnt[3]) {

    //first figure out what element we're in
    int i,j,k;
    get_element_ijk(pnt,i,j,k);

    //now get the data for that element and construct associated points
    double sdvs[8];
    sdvs[0] = signed_distance_values[(i)+(j)*dims[0]+(k)*dims[0]*dims[1]];
    sdvs[1] = signed_distance_values[(i+1)+(j)*dims[0]+(k)*dims[0]*dims[1]];
    sdvs[2] = signed_distance_values[(i+1)+(j+1)*dims[0]+(k)*dims[0]*dims[1]];
    sdvs[3] = signed_distance_values[(i)+(j+1)*dims[0]+(k)*dims[0]*dims[1]];
    sdvs[4] = signed_distance_values[(i)+(j)*dims[0]+(k+1)*dims[0]*dims[1]];
    sdvs[5] = signed_distance_values[(i+1)+(j)*dims[0]+(k+1)*dims[0]*dims[1]];
    sdvs[6] = signed_distance_values[(i+1)+(j+1)*dims[0]+(k+1)*dims[0]*dims[1]];
    sdvs[7] = signed_distance_values[(i)+(j+1)*dims[0]+(k+1)*dims[0]*dims[1]];

    CartVect coords[8];
    coords[0] = CartVect((i)*step_size,(j)*step_size,(k)*step_size);
    coords[1] = CartVect((i+1)*step_size,(j)*step_size,(k)*step_size);
    coords[2] = CartVect((i+1)*step_size,(j+1)*step_size,(k)*step_size);
    coords[3] = CartVect((i)*step_size,(j+1)*step_size,(k)*step_size);
    coords[4] = CartVect((i)*step_size,(j)*step_size,(k+1)*step_size);
    coords[5] = CartVect((i+1)*step_size,(j)*step_size,(k+1)*step_size);
    coords[6] = CartVect((i+1)*step_size,(j+1)*step_size,(k+1)*step_size);
    coords[7] = CartVect((i)*step_size,(j+1)*step_size,(k+1)*step_size);

    //interpolate in x
    double a,b,c,d;
    a = sdvs[0] + ((sdvs[1]-sdvs[0])/(coords[1][0]-coords[0][0]))*(pnt[0]-coords[0][0]);
    b = sdvs[2] + ((sdvs[3]-sdvs[2])/(coords[3][0]-coords[2][0]))*(pnt[0]-coords[2][0]);
    c = sdvs[4] + ((sdvs[5]-sdvs[4])/(coords[5][0]-coords[4][0]))*(pnt[0]-coords[4][0]);
    d = sdvs[6] + ((sdvs[7]-sdvs[6])/(coords[7][0]-coords[6][0]))*(pnt[0]-coords[6][0]);
  
    // interpolate in y
    double aa,bb;
    aa = a+((b-a)/(coords[3][1]-coords[1][1]))*(pnt[1]-coords[1][1]);
    // aa = a*(fabs(pnt[1]-coords[3][1])/precondStep)+b*(fabs(pnt[1]-coords[1][1])/precondStep);
    bb = c+((d-c)/(coords[7][1]-coords[5][1]))*(pnt[1]-coords[5][1]);
    // bb = c*(fabs(pnt[1]-coords[7][1])/precondStep)+d*(fabs(pnt[1]-coords[5][1])/precondStep);

    //finally interpolate in z
    double result = aa+((bb-aa)/(coords[5][2]-coords[1][2]))*(pnt[2]-coords[1][2]);

    return result;
  };

  
  void SignedDistanceField::get_element_ijk(const double pnt[3], int &i, int &j, int &k) {
    CartVect this_pnt(pnt);
    CartVect llc(lower_left_corner);
    CartVect vec = this_pnt-llc;
    i = vec[0]/step_size;
    j = vec[1]/step_size;
    k = vec[2]/step_size;
  };


  void SignedDistanceField::get_element_ijk(const double pnt[3], int &i, int &j, int &k, bool &outside) {
    get_element_ijk(pnt, i, j, k);
    if( i < 0 || j < 0 || k < 0 || i > dims[0]-2 || j > dims[1]-2 || k > dims[2]-2) {
    outside = true;
     }
  };

}
