#include <iostream>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "MBTypes.h"

#include "MBSkinner.hpp"
#include "MBCartVect.hpp"

#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

#include "gen.hpp"
#include "geometry.hpp"

MBInterface *MOAB()
{
    static MBCore instance;
    return &instance;
}

MBErrorCode calculate_geometry_size_before_sealing(const moab::Range geometry_sets[], 
						   const MBTag geom_tag, const MBTag size_tag)
{
  MBErrorCode rval;
  for ( int d = 1 ; d <= 3 ; d++ )
    // dont bother with vertices (d=0) since they 
    // have no area/volume
    {
      for ( moab::Range::iterator i = geometry_sets[d].begin() ; i != geometry_sets[d].end() ; i++ )
	{
	  double size;
	  rval = geometry::measure( *i, geom_tag, size );
	  if(gen::error(MB_SUCCESS!=rval,"could not measure")) 
	    {
	      return rval;
	    } 
	  if(gen::error(MB_SUCCESS!=rval,"could not set size tag")) 
	    {
	      return rval;
	    }
	}
    }

  return MB_SUCCESS;
}

int main(int argc, char *argv[])
{
  // moab::Core *moab = new moab::Core(); // new moab instance

  std::string filename = "../test/cones.h5m"; // set the fname

  MBErrorCode rval = MOAB()->load_file(filename.c_str()); // load the h5m file
  if(rval != MB_SUCCESS)
    {
      std::cout << "Could not open mesh file";
      exit(0);
    }
  
  MBTag geom_tag, id_tag, true_int, orig_curve_tag, size_tag,
    geometry_resabs_tag,faceting_tol_tag,merge_tag,normal_tag;

  
  rval = MOAB()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,MB_TYPE_INTEGER, geom_tag, 
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(rval!=MB_SUCCESS,"could not set geom_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1,MB_TYPE_INTEGER, id_tag, 
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
  if(gen::error(rval!=MB_SUCCESS,"could not set id_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "NORMAL", sizeof(MBCartVect), MB_TYPE_OPAQUE,normal_tag, 
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
  if(gen::error(rval!=MB_SUCCESS,"could not set normal_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "MERGE", 1, MB_TYPE_HANDLE, merge_tag, 
			       moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  if(gen::error(rval!=MB_SUCCESS,"could not set merge_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "FACETING_TOL", 1, MB_TYPE_DOUBLE,faceting_tol_tag , 
			       moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  if(gen::error(rval!=MB_SUCCESS,"could not set faceting_tol_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "GEOMETRY_RESABS", 1,     MB_TYPE_DOUBLE,geometry_resabs_tag, 
			       moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT  );
  if(gen::error(rval!=MB_SUCCESS,"could not set geometry_resabs_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "GEOM_SIZE", 1, MB_TYPE_DOUBLE,size_tag, 
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT  );
  if(gen::error(rval!=MB_SUCCESS,"could not set size_tag tag"))
    {
      return rval;
    }

  rval = MOAB()->tag_get_handle( "ORIG_CURVE", 1, MB_TYPE_INTEGER, orig_curve_tag, 
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT, &true_int );
  if(gen::error(rval!=MB_SUCCESS,"could not set orig_curve_tag tag"))
    {
      return rval;
    }

  //moab::Range file_set;
  MBRange file_set;
  rval = MOAB()->get_entities_by_type_and_tag( 0, MBENTITYSET, 
		    &faceting_tol_tag, NULL, 1, file_set, false);
  if(gen::error(MB_SUCCESS!=rval,"could not get faceting_tol_tag")) 
    {
      return rval;
    }

  if(gen::error(file_set.empty(),"file set not found"))
    {
      return rval;
    }

  if(gen::error(1!=file_set.size(),"Refacet with newer version of ReadCGM.")) 
    {
      return MB_FAILURE;
    }

  // get the faceting tolerance from the mesh
  double facet_tol, sme_resabs_tol=1e-6;
  rval = MOAB()->tag_get_data( faceting_tol_tag, &file_set.front(), 1,  
                                  &facet_tol );
  if (rval != MB_SUCCESS )
    {
      std::cout << "Could  not get faceting tolerance from file" << std::endl;
    }

  rval = MOAB()->tag_get_data( geometry_resabs_tag, &file_set.front(), 1,  
                                  &sme_resabs_tol );
  if (rval != MB_SUCCESS )
    {
      std::cout << "Could not get absolute tolerance from file" << std::endl;
    }

  std::cout << "Faceting tolerance = " << facet_tol << " Absolute tolerance = " << sme_resabs_tol << std::endl;
  
  // get geometry set information
  //

  moab::Range geometry_sets[4]; // moab range of entities
  for ( int d = 0 ; d <= 3 ; d++) // loop over dimensionality
    {
      void *val[] = {&d};
      rval = MOAB()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, val , 1, geometry_sets[d] );
      if(gen::error(MB_SUCCESS!=rval,"could not get geometry information")) 
	{
	  return rval;
	}

      // looop over the members of the set of dimensionality
      // d
      for ( moab::Range::iterator i=geometry_sets[d].begin() ; i != geometry_sets[d].end() ; i++ )
	{
	  unsigned int options;
	  rval = MOAB()->get_meshset_options(*i, options);
	  if( rval != MB_SUCCESS)
	    {
	      std::cout << "could not get mesh options" << std::endl;
	      return rval;
	    }

	    // if options are wrong change them
	    if(d==1) 
	      {
		if( !(MESHSET_TRACK_OWNER &options) || !(MESHSET_ORDERED &options) ) 
		  {
		    rval = MOAB()->set_meshset_options(*i, MESHSET_TRACK_OWNER | MESHSET_ORDERED);
		  }
	      } 
	    else 
	      {
		if( !(MESHSET_TRACK_OWNER &options) ) 
		  {        
		    rval = MOAB()->set_meshset_options(*i, MESHSET_TRACK_OWNER);
		  }
	      }

	}
    }

  // check the geometry size before sealing the geometry
  // find each entities size (vol,area etc)

  rval = calculate_geometry_size_before_sealing(geometry_sets, geom_tag, size_tag); 
  if(gen::error(MB_SUCCESS!=rval,"Calculation of geometry size failed"))
    {
      return rval;
    }

  std::cout << "Get entity count before sealing" << std::endl;
  // Get entity count before sealing.
  int orig_n_tris;
  rval = MOAB()->get_number_entities_by_type( 0, MBTRI, orig_n_tris );
  if(gen::error(rval!=MB_SUCCESS,"could not calculate the number of entities"))
    {
      exit(rval);
    }
  std::cout << "  input faceted geometry contains " << geometry_sets[3].size() << " volumes, " 
	    << geometry_sets[2].size() << " surfaces, " << geometry_sets[1].size() 
	    << " curves, and " << orig_n_tris << " triangles" << std::endl;  
  
  std::cout << "Finding degenerate triangles " << std::endl;



  return(0);
}


