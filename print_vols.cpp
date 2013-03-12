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
  int i;

  if( argc <= 1 ) // only the executable name provided
    {
      std::cout << "Usage is print_vols <mesh_file> " << std::endl;
      return 1;
    }
  else
    {
      for ( i = 1 ; i <= argc-1 ; i++ ) // loop over input args
	{
	  if ( std::string::npos != std::string(argv[i]).find(std::string(".h5m")) ) //if input file contains .h5m
	    {
	      filename = argv[i]; //set the filename
	    }
	  else
	    {
	      std::cout << "Cannot open h5m file, " << argv[i] << std::endl;
	      return 1;
	    }
	}
    }



  //filename = "iter.h5m";
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

  get_line_data(mb); // get line data

  get_tri_data(mb); // get facet data

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

// get entities of dimension 1, 
void get_line_data (MBInterface *mb)
{
  int dimension = 1;
  MBRange entities;
  MBErrorCode rval;
  std::ofstream vtk_output;

  vtk_output.open("lines.vtk");

  entities.clear();  // clear entity information

  rval = mb->get_entities_by_dimension(0, dimension, entities); // get 1 d entities

  std::cout << "number of 1 d entities = " << entities.size() << std::endl;
  
  double *xyz;
  xyz = new double[3*entities.size()];
  MBErrorCode error_code = mb->get_coords(entities,xyz);


  vtk_output << "# vtk DataFile Version 2.0" << std::endl;
  vtk_output << "Bunch of points from h5m file" << std::endl;
  vtk_output << "ASCII" << std::endl;
  vtk_output << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtk_output << "POINTS " << entities.size() << " FLOAT" << std::endl;

  unsigned int i,j;

  for ( i = 0 ; i <= entities.size()-1 ; i++ )
    {
      j = 3*i;
      vtk_output << xyz[j] << " " << xyz[j+1] << " " << xyz[j+2] << std::endl; 
    }

  vtk_output << " " << std::endl;
  vtk_output << "CELLS " << entities.size() << " " << 3*entities.size() << std::endl;
  for ( i = 0 ; i <= (entities.size())-1 ; i++ )
    {
      j = i*2;
      vtk_output << "2 " << j  << " " << j+1 << std::endl;
    }
	  
  vtk_output << " " << std::endl;
  vtk_output << "CELL_TYPES " << entities.size() << std::endl;
  for ( i = 0 ; i <= (entities.size())-1 ; i++ )
    {
      vtk_output << "3" << std::endl;
    }	  
  delete [] xyz;
  vtk_output.close();

  return;

}

// get entities of dimension 1, 
void get_tri_data (MBInterface *mb)
{
  int dimension = 2;
  MBRange entities;
  MBErrorCode rval;
  std::ofstream vtk_output;

  vtk_output.open("tris.vtk");

  entities.clear();  // clear entity information

  rval = mb->get_entities_by_dimension(0, dimension, entities); // get 1 d entities

  std::cout << "number of 2 d entities = " << entities.size() << std::endl;
  
  double *xyz;
  xyz = new double[3*entities.size()];
  MBErrorCode error_code = mb->get_coords(entities,xyz);


  vtk_output << "# vtk DataFile Version 2.0" << std::endl;
  vtk_output << "Bunch of points from h5m file" << std::endl;
  vtk_output << "ASCII" << std::endl;
  vtk_output << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtk_output << "POINTS " << entities.size() << " FLOAT" << std::endl;

  unsigned int i,j;

  for ( i = 0 ; i <= entities.size()-1 ; i++ )
    {
      j = 3*i;
      vtk_output << xyz[j] << " " << xyz[j+1] << " " << xyz[j+2] << std::endl; 
    }

  vtk_output << " " << std::endl;
  vtk_output << "CELLS " << entities.size() << " " << 4*entities.size() << std::endl;
  for ( i = 0 ; i <= (entities.size())-1 ; i++ )
    {
      j = i*3;
      vtk_output << "3 " << j  << " " << j+1 << " " << j+2 << std::endl;
    }
	  
  vtk_output << " " << std::endl;
  vtk_output << "CELL_TYPES " << entities.size() << std::endl;
  for ( i = 0 ; i <= (entities.size())-1 ; i++ )
    {
      vtk_output << "5" << std::endl;
    }	  
  delete [] xyz;
  vtk_output.close();

  return;

}
