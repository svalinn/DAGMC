
#include "sdf.hpp"

namespace moab {

  // Constructors
  SignedDistanceField::SignedDistanceField(double x_min, double y_min, double z_min, double mesh_step_size, int num_x_points, int num_y_points, int num_z_points) {
    lower_left_corner[0] = x_min; lower_left_corner[1] = y_min; lower_left_corner[2] = z_min;
    step_size = mesh_step_size;
    dims[0] = num_x_points; dims[1] = num_y_points; dims[2] = num_z_points;
  }

  double SignedDistanceField::find_sdv(const double pnt[3]) {

    //first figure out what element we're in
    int i,j,k;
    bool outside;
    get_element_indices_for_pnt(pnt,i,j,k,outside);

    //if point is outside the data structure, return big negative value
    if(outside) {
    return -1.e37;
    }

    //now get the data for that element and construct associated points
    double sdvs[8];
    sdvs[0] = get_data(i,j,k);
    sdvs[1] = get_data(i+1,j,k);
    sdvs[2] = get_data(i+1,j+1,k);
    sdvs[3] = get_data(i,j+1,k);
    sdvs[4] = get_data(i,j,k+1);
    sdvs[5] = get_data(i+1,j,k+1);
    sdvs[6] = get_data(i+1,j+1,k+1);
    sdvs[7] = get_data(i,j+1,k+1);

    CartVect coords[8];
    coords[0] = get_coords(i,j,k);
    coords[1] = get_coords(i+1,j,k);
    coords[2] = get_coords(i+1,j+1,k);
    coords[3] = get_coords(i,j+1,k);
    coords[4] = get_coords(i,j,k+1);
    coords[5] = get_coords(i+1,j,k+1);
    coords[6] = get_coords(i+1,j+1,k+1);
    coords[7] = get_coords(i,j+1,k+1);

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
  }

  
  void SignedDistanceField::get_element_indices_for_pnt(const double pnt[3], int &i, int &j, int &k) {
    CartVect this_pnt(pnt);
    CartVect llc(lower_left_corner);
    CartVect vec = this_pnt-llc;
    i = vec[0]/step_size;
    j = vec[1]/step_size;
    k = vec[2]/step_size;
  }


  void SignedDistanceField::get_element_indices_for_pnt(const double pnt[3], int &i, int &j, int &k, bool &outside) {
    get_element_indices_for_pnt(pnt, i, j, k);
    if( i < 0 || j < 0 || k < 0 || i > dims[0]-2 || j > dims[1]-2 || k > dims[2]-2) {
      outside = true;
     }
    else {
      outside = false;
    }
  }

  ErrorCode SignedDistanceField::create_scdBox(ScdBox *sdfBox, Interface* mbi) {

    //create new moab instance if one doesn't already exist
    bool created_new_moab = false;
    if( !mbi ) {
      Interface* mbi = new Core();
      created_new_moab = true;
    }
    
    ScdInterface* scdi  = new ScdInterface(mbi);

    int num_x_pnts,num_y_pnts,num_z_pnts;
    get_dims(num_x_pnts,num_y_pnts,num_z_pnts);

    double x_min,y_min,z_min;
    get_llc(x_min,y_min,z_min);

    double step_size;
    get_step(step_size);
 
    //create the relevant locations and 
    std::vector<double> pnts;
    for(unsigned int k = 0; k < num_z_pnts; k++){
      for(unsigned int j = 0; j < num_y_pnts; j++){
    	for(unsigned int i = 0; i < num_x_pnts; i++){
	  CartVect coords = get_coords(i,j,k);
    	  pnts.push_back(coords[0]);
    	  pnts.push_back(coords[1]);
    	  pnts.push_back(coords[2]);
    	}
      }
    }

    // now create an SCD box to track all of this info
    HomCoord l = HomCoord(0,0,0);
    HomCoord h = HomCoord(num_x_pnts-1, num_y_pnts-1, num_z_pnts-1);
    ErrorCode rval = scdi->construct_box(l, h, &(pnts[0]), pnts.size(), sdfBox);

    Tag sdfTag;
    rval = mbi->tag_get_handle( sdf_tag_name.c_str(), 1, MB_TYPE_DOUBLE, sdfTag, MB_TAG_DENSE|MB_TAG_CREAT );
    MB_CHK_SET_ERR(rval, "Could not create the signed distance field tag.");
    for(unsigned int k = 0 ; k < num_z_pnts; k++){
       for(unsigned int j = 0 ; j < num_y_pnts; j++){  
	 for(unsigned int i = 0 ; i < num_x_pnts; i++){
	   EntityHandle vert = sdfBox->get_vertex(i,j,k);
	   double sdv = get_data(i,j,k);
	   void *ptr = &sdv;
	   rval = mbi->tag_set_data(sdfTag, &vert, 1, ptr);
	   MB_CHK_SET_ERR(rval, "Could not tag vert with signed distance value");
	 }
       }
    }
  
    //if we had to create a new moab instance, clean it out and the ptr
    if(created_new_moab) {
      rval = mbi->delete_mesh();
      MB_CHK_SET_ERR(rval,"Could not delete the mesh");
      delete mbi;
    }    

    delete scdi;

    return MB_SUCCESS;
  }
  
}
