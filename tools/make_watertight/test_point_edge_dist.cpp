#include <iostream>
#include <math.h>

double dot_product( double vect0[], double vect1[] ) {                                        
  return vect0[0]*vect1[0] + vect0[1]*vect1[1] + vect0[2]*vect1[2];                           
}                                                                                             
                                                                                              
void cross_product( double vect0[], double vect1[], double vect2[] ) {                        
  vect2[0] = vect0[1]*vect1[2] - vect0[2]*vect1[1];                                           
  vect2[1] = vect0[2]*vect1[0] - vect0[0]*vect1[2];                                           
  vect2[2] = vect0[0]*vect1[1] - vect0[1]*vect1[0];                                           
}

double dist_between_verts( double coords0[], double coords1[] ) {                             
  return sqrt( (coords0[0]-coords1[0])*(coords0[0]-coords1[0]) +                              
               (coords0[1]-coords1[1])*(coords0[1]-coords1[1]) +                              
               (coords0[2]-coords1[2])*(coords0[2]-coords1[2]) );                             
}

// from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1                    
double point_edge_dist( double a[], double b[], double c[] ) {                                
  double ab[3], bc[3], ba[3], ac[3];                                                          
  for(int i=0; i<3; i++) {                                                                    
    ab[i] = b[i] - a[i];                                                                      
    bc[i] = c[i] - b[i];                                                                      
    ba[i] = a[i] - b[i];                                                                      
    ac[i] = c[i] - a[i];                                                                      
  }                                                                                           
                                                                                              
  // find the magnitude of the cross product and test the line                                
  double cross_prod[3];                                                                       
  cross_product(ab, ac, cross_prod);                                                          
  double dist = sqrt(dot_product(cross_prod,cross_prod)) / dist_between_verts(a,b);           
                                                                                              
  // test endpoint1                                                                           
  double dot1 = dot_product( ab, bc );                                                        
  if (dot1 > 0) {                                                                             
    std::cout << "point_edge_dist=" << dist_between_verts(b,c) << " at endpt1" << std::endl;  
    return dist_between_verts(b,c);                                                           
  }                                                                                           
                                                                                              
  // test endpoint0                                                                           
  double dot2 = dot_product( ba, ac );                                                        
  if (dot2 > 0) {                                                                             
    std::cout << "point_edge_dist=" << dist_between_verts(a,c) << " at endpt0" << std::endl;  
    return dist_between_verts(a,c);                                                           
  }                                                                                           
                                                                                              
  std::cout << "point_edge_dist=" << fabs(dist) << " at middle" << std::endl;                 
  return fabs(dist);                                                                          
}     

int main() {
  double a[] = { 859.5510374 ,90.1085829 ,1.929114349 };
  double b[] = { 859.5571858, 90.1057543, 2.56411115 };
  double c[] = { 859.556885,    90.105728,     2.564118 };
  std::cout << point_edge_dist( a,b,c) << std::endl;
}
