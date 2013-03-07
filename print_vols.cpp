// Load in a h5m file and then find a given surface
#include <iostream>
#include <fstream>

#include "MBCore.hpp"
#include "MBSkinner.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"

#include "functions.h"

int main(int argc, char *argv[])
{  
  std::string filename;
  filename = "iter.h5m";
  // instantiate & load a mesh from a file
  MBCore *mb = new MBCore();
  MBErrorCode rval = mb->load_mesh(filename.c_str());
  if (moab::MB_SUCCESS != rval) 
  {
    std::cout << "Failed to read file, " << filename << std::endl;
    return 1;
  }

  std::cout << "File loaded" << std::endl;
  get_num_ents(mb); // get numer of entities

  return 0;
}

// Get entity information
void get_num_ents( MBInterface *mb)
{
  MBRange ents;
  MBErrorCode rval;

  // iterate over dimensions
  for (int d = 0 ; d <= 3; d++) 
    {
    ents.clear();
    // clear entity information
    rval = mb->get_entities_by_dimension(0, d, ents); 
    // get entity information based on dimension

    if ( MB_SUCCESS != rval) //if we fail to 
      {
	return; // return to caller
      }

    std::cout << "Found " << ents.size() << " " << d 
	                  << "-dimensional entities:" << std::endl;

    if ( d == 0 )  // if vertices
      { 
	std::ofstream vtk_output;
	double *xyz;
	xyz = new double[3*ents.size()];
	MBErrorCode error_code = mb->get_coords(ents,xyz);
	if ( error_code != MB_SUCCESS )
	  {
	    std::cout << "Failed to get vertex positions" << std::endl;
	  }

	unsigned int i,j;

	vtk_output.open("points.vtk");

	vtk_output << "# vtk DataFile Version 2.0" << std::endl;
	vtk_output << "Bunch of points from h5m file" << std::endl;
	vtk_output << "ASCII" << std::endl;
	vtk_output << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtk_output << "POINTS " << ents.size() << " FLOAT" << std::endl;

	for ( i = 0 ; i <= (ents.size())-1 ; i++ )
	  {
	    j = 3*i;
	    vtk_output << xyz[j] << " " 
		      << xyz[j+1] << " " << xyz[j+2] << std::endl;
	  }
	
          vtk_output << " " << std::endl;
	  vtk_output << "CELLS " << ents.size() << " " << 2*ents.size() << std::endl;
	  for ( i = 0 ; i <= (ents.size())-1 ; i++ )
	    {
	      vtk_output << "1 " << i << std::endl;
	    }

	  vtk_output << " " << std::endl;
	  vtk_output << "CELL_TYPES " << ents.size() << std::endl;
	  for ( i = 0 ; i <= (ents.size())-1 ; i++ )
	    {
	      vtk_output << "1" << std::endl;
	    }	  
	  delete [] xyz;
	  vtk_output.close();
      }
   }  
  return;
}

