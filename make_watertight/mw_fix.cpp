// ********************************************************************
// Patrick Shriwise
// August 2013
/* _curve_to_be_tested_for_watertightness_
      vert1 X X vert1
            | |
      vert2 X |
  surf1     | |    surf2
            | |
      vert3 X X vert2
            | |
      vert4 X X vert3                   */

// input:  h5m filename, tolerance
// output: watertight h5m

#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"



const char GEOM_SENSE_2_TAG_NAME[] = "GEOM_SENSE_2";
const char GEOM_SENSE_N_ENTS_TAG_NAME[] = "GEOM_SENSE_N_ENTS";
const char GEOM_SENSE_N_SENSES_TAG_NAME[] = "GEOM_SENSE_N_SENSES"; 

moab::Interface *MBI();


char* sense_printer( int sense){

  if ( sense == 1 ) return "FORWARD (1)";
  if ( sense == -1 ) return  "REVERSE (-1)";
  if ( sense == 0 ) return "UNKNOWN (0)";
}

void moab_printer(moab::ErrorCode error_code)
{
  if ( error_code == moab::MB_INDEX_OUT_OF_RANGE )
    {
      std::cerr << "ERROR: moab::MB_INDEX_OUT_OF_RANGE" << std::endl;
    }
  if ( error_code == moab::MB_MEMORY_ALLOCATION_FAILED )
    {
      std::cerr << "ERROR: moab::MB_MEMORY_ALLOCATION_FAILED" << std::endl;
    }
  if ( error_code == moab::MB_ENTITY_NOT_FOUND )
    {
      std::cerr << "ERROR: moab::MB_ENTITY_NOT_FOUND" << std::endl;
    }
  if ( error_code == moab::MB_MULTIPLE_ENTITIES_FOUND )
    {
      std::cerr << "ERROR: moab::MB_MULTIPLE_ENTITIES_FOUND" << std::endl;
    }
  if ( error_code == moab::MB_TAG_NOT_FOUND )
    {
      std::cerr << "ERROR: moab::MB_TAG_NOT_FOUND" << std::endl;
    }
  if ( error_code == moab::MB_FILE_DOES_NOT_EXIST )
    {
      std::cerr << "ERROR: moab::MB_FILE_DOES_NOT_EXIST" << std::endl;
    }    
  if ( error_code == moab::MB_FILE_WRITE_ERROR )
    {
      std::cerr << "ERROR: moab::MB_FILE_WRITE_ERROR" << std::endl;
    }    
  if ( error_code == moab::MB_ALREADY_ALLOCATED )
    {
      std::cerr << "ERROR: moab::MB_ALREADY_ALLOCATED" << std::endl;
    }    
  if ( error_code == moab::MB_VARIABLE_DATA_LENGTH )
    {
      std::cerr << "ERROR: moab::MB_VARIABLE_DATA_LENGTH" << std::endl;
    }  
  if ( error_code == moab::MB_INVALID_SIZE )
    {
      std::cerr << "ERROR: moab::MB_INVALID_SIZE" << std::endl;
    }  
  if ( error_code == moab::MB_UNSUPPORTED_OPERATION )
    {
      std::cerr << "ERROR: moab::MB_UNSUPPORTED_OPERATION" << std::endl;
    }  
  if ( error_code == moab::MB_UNHANDLED_OPTION )
    {
      std::cerr << "ERROR: moab::MB_UNHANDLED_OPTION" << std::endl;
    }  
  if ( error_code == moab::MB_FAILURE )
    {
      std::cerr << "ERROR: moab::MB_FAILURE" << std::endl;
    }  
  return;
}


moab::ErrorCode get_geom_size_before_sealing( const moab::Range geom_sets[], 
                                          const moab::Tag geom_tag,
                                          const moab::Tag size_tag,
                                          bool verbose ) {
  moab::ErrorCode rval;
  for(int dim=1; dim <= 3 ; dim++) {
    std::cout << "dim = " << dim << std::endl;
    for(moab::Range::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
      double size = 0;

      rval = gen::measure( *i, geom_tag, size, false, verbose );
      if(gen::error(moab::MB_SUCCESS!=rval,"could not measure")) return rval;
      rval = MBI()->tag_set_data( size_tag, &(*i), 1, &size );
      if(gen::error(moab::MB_SUCCESS!=rval,"could not set size tag")) return rval;

    }
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode get_senses(moab::EntityHandle entity,
    std::vector<moab::EntityHandle> &wrt_entities, std::vector<int> &senses)
{
  int edim = 1;

  if (-1 == edim)
    return moab::MB_FAILURE;// not geometry entity

  moab::ErrorCode rval;
  wrt_entities.clear();
  senses.clear();

  if (1 == edim)// edge
  {
    
    const void *dum_ptr;
    int num_ents;
    unsigned flags = moab::MB_TAG_SPARSE;
  
    moab::Tag senseNEntsTag;
    rval = MBI() -> tag_get_handle(GEOM_SENSE_N_ENTS_TAG_NAME, 0 , moab::MB_TYPE_HANDLE, senseNEntsTag, flags);
    if (gen::error(moab::MB_SUCCESS!=rval, "could not get ent sense handles")) return rval;
    rval = MBI()->tag_get_by_ptr(senseNEntsTag, &entity, 1, &dum_ptr, &num_ents);
    if (gen::error(moab::MB_SUCCESS!=rval, "could not get ent sense data")) return rval;
    const moab::EntityHandle *ents_data = static_cast<const moab::EntityHandle*> (dum_ptr);
    std::copy(ents_data, ents_data + num_ents, std::back_inserter(wrt_entities));

    moab::Tag senseNSensesTag;
    rval = MBI()->tag_get_handle(GEOM_SENSE_N_SENSES_TAG_NAME, 0 , moab::MB_TYPE_INTEGER, senseNSensesTag, flags);
    if (gen::error(moab::MB_SUCCESS!=rval, "could not get senses handle")) return rval;
    rval = MBI()->tag_get_by_ptr(senseNSensesTag, &entity, 1, &dum_ptr,
        &num_ents);
    if (gen::error(moab::MB_SUCCESS!=rval, "could not get senses data")) return rval;

    const int *senses_data = static_cast<const int*> (dum_ptr);
    std::copy(senses_data, senses_data + num_ents, std::back_inserter(senses));

  }

  return moab::MB_SUCCESS;
}
  int main(int argc, char **argv) {

   // ******************************************************************
    // Load the h5m file and create tags.
    // ******************************************************************

    clock_t start_time = clock();
    const bool debug = false;
    const bool check_geom_size = true;
    bool verbose = true;

    // check input args
    if( 2 > argc || 3 < argc ) 
      {
	std::cout << "To zip a faceted h5m file:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.h5m>" << std::endl;
	std::cout << "To facet and zip an ACIS file using the default facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat>" << std::endl;
	std::cout << "To facet and zip an ACIS file using a specified facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat> <facet_tolerance>" << std::endl;
	return moab::MB_FAILURE;
      }

    // The root name does not have an extension
    std::string input_name = argv[1];
    std::string root_name = argv[1];
    int len = root_name.length();
    root_name.erase(len - 4);
    bool is_acis;

    // load the input file
    moab::ErrorCode result, rval;
    moab::EntityHandle input_set;

    rval = MBI()->create_meshset( moab::MESHSET_SET, input_set );

    if(gen::error(moab::MB_SUCCESS!=rval,"failed to create_meshset"))
      {
	return rval;
      }

    std::cout << "Loading input file..." << std::endl;

    // If reading an h5m file, the facet tolerance has already been determined.
    // Read the facet_tol from the file_set. There should only be one input
    // argument.

    if(std::string::npos!=input_name.find("h5m") && (2==argc)) 
      {
	rval = MBI()->load_file( input_name.c_str(), &input_set );
	if(gen::error(moab::MB_SUCCESS!=rval,"failed to load_file 0")) 
	  {
	    return rval;      
	  }
	
	is_acis = false;

      } 

    // create tags
    clock_t load_time = clock();    
    moab::Tag geom_tag, id_tag, normal_tag, merge_tag, faceting_tol_tag, 
      geometry_resabs_tag, size_tag, orig_curve_tag;
  
    result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
				moab::MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
    assert( moab::MB_SUCCESS == result );
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1,
				moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert( moab::MB_SUCCESS == result );
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "NORMAL", sizeof(moab::CartVect), moab::MB_TYPE_OPAQUE,
        normal_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert( moab::MB_SUCCESS == result );
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "MERGE", 1, moab::MB_TYPE_HANDLE,
        merge_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
    assert( moab::MB_SUCCESS == result ); 
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      } 
    result = MBI()->tag_get_handle( "FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
        faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
    assert( moab::MB_SUCCESS == result );  
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "GEOMETRY_RESABS", 1,     moab::MB_TYPE_DOUBLE,
                             geometry_resabs_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT  );
    assert( moab::MB_SUCCESS == result );  
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "GEOM_SIZE", 1, moab::MB_TYPE_DOUBLE,
				    size_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT  );
    assert( (moab::MB_SUCCESS == result) );
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    int true_int = 1;    
    result = MBI()->tag_get_handle( "ORIG_CURVE", 1,
				moab::MB_TYPE_INTEGER, orig_curve_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT, &true_int );
    assert( moab::MB_SUCCESS == result );
    if ( result != moab::MB_SUCCESS )
      {
	moab_printer(result);
      }
    // PROBLEM: MOAB is not consistent with file_set behavior. The tag may not be
    // on the file_set.
    moab::Range file_set;
    result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &faceting_tol_tag,
                                                  NULL, 1, file_set );

    if(gen::error(moab::MB_SUCCESS!=result,"could not get faceting_tol_tag")) 
      {
	return result;
      }

    gen::error(file_set.empty(),"file set not found");

    if(gen::error(1!=file_set.size(),"Refacet with newer version of ReadCGM.")) 
      {
	return moab::MB_FAILURE;
      }

    double facet_tol, sme_resabs_tol=1e-6;
    result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1,  
                                  &facet_tol );
    assert(moab::MB_SUCCESS == result);
    result = MBI()->tag_get_data( geometry_resabs_tag, &file_set.front(), 1,  
                                  &sme_resabs_tol );
    if(moab::MB_SUCCESS != result) 
      {
	std::cout <<  "absolute tolerance could not be read from file" << std::endl;
      }

    // In practice, use 2*facet_tol because we are always comparing 2 faceted
    // entities. If instead we were comparing a faceted entity and a geometric
    // entitiy, then 1*facet_tol is correct.

    const double SME_RESABS_TOL = sme_resabs_tol; // from ACIS through CGM
    const double FACET_TOL = facet_tol; // specified by user when faceting cad
    std::cout << "  faceting tolerance=" << facet_tol << " cm" << std::endl;
    std::cout << "  absolute tolerance=" << sme_resabs_tol << " cm" << std::endl;
    
    
    // get all geometry sets
    moab::Range geom_sets[4];
    for(unsigned dim=0; dim<4; dim++) 
      {
	void *val[] = {&dim};
	result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
	std::cout << "Get entities by type and tag" << std::endl;

	assert(moab::MB_SUCCESS == result);

	// make sure that sets TRACK membership and curves are ordered
	for(moab::Range::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) 
	  {
	    unsigned int options;
	    result = MBI()->get_meshset_options(*i, options );
	    assert(moab::MB_SUCCESS == result);
    
	    // if options are wrong change them
	    if(dim==1) 
	      {
		if( !(moab::MESHSET_TRACK_OWNER&options) || !(moab::MESHSET_ORDERED&options) ) 
		  {
		    result = MBI()->set_meshset_options(*i, moab::MESHSET_TRACK_OWNER|moab::MESHSET_ORDERED);
		    assert(moab::MB_SUCCESS == result);
		  }
	      } 
	    else 
	      {
		if( !(moab::MESHSET_TRACK_OWNER&options) ) 
		  {        
		    result = MBI()->set_meshset_options(*i, moab::MESHSET_TRACK_OWNER);
		    assert(moab::MB_SUCCESS == result);
		  }
	      }
	  }
      }

    std::cout << "I am here" << std::endl;

    // this could be related to when there are sat files rather than mesh?
    // If desired, find each entity's size before sealing.
    if(check_geom_size) 
      {
	std::cout << "I am checking the geometry size" << std::endl;
	result = get_geom_size_before_sealing( geom_sets, geom_tag, size_tag, verbose );
	if(gen::error(moab::MB_SUCCESS!=result,"measuring geom size failed"))
	  {
	    return result;
	  }
      }
    std::cout << "Get entity count before sealing" << std::endl;
    // Get entity count before sealing.
    int orig_n_tris;
    result = MBI()->get_number_entities_by_type( 0, moab::MBTRI, orig_n_tris );
    std::cout << result << std::endl;

    assert(moab::MB_SUCCESS == result);

    std::cout << "==================================" << std::endl;
    std::cout << "  Input faceted geometry contains: " << std::endl;
    std::cout << "==================================" << std::endl;

    std::cout << geom_sets[3].size() << " volumes, " 
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size() 
              << " curves, " << orig_n_tris << " triangles, and " 
              << geom_sets[0].size() << " vertices" << std::endl;  

    std::cout << "==================================" << std::endl;
    

// Get all curve senses
    std::cout << "=====================================" << std::endl;
    std::cout << " CURVE SENSES " << std::endl;
    std::cout << "=====================================" << std::endl;

    moab::GeomTopoTool gt(MBI(), false);
    
    for( unsigned int i=0; i<geom_sets[1].size(); i++)
    {
    moab::EntityHandle curve = geom_sets[1][i];
    std::vector<moab::EntityHandle> surfs;
    std::vector<int> senses;
    rval = gt.get_senses( curve, surfs, senses);
    if(gen::error(moab::MB_SUCCESS!=result,"could not get curve senses")) return result;
    std::cout << "Number of senses for curve " << gen::geom_id_by_handle(curve) << " = " << senses.size() << std::endl;
    for (unsigned int index=0; index<senses.size() ; index++)
    { 
     std::cout << "surf = " << gen::geom_id_by_handle(surfs[index]) << std::endl;
     std::cout << "sense = " << sense_printer( senses[index] ) << std::endl;
    }
    std::cout << std::endl;
    }

    

//Get all surface senses

    std::cout << "=====================================" << std::endl;
    std::cout << " SURFACE SENSES " << std::endl;
    std::cout << "=====================================" << std::endl;


    for( unsigned int i=0; i<geom_sets[2].size(); i++)
    {
    moab::EntityHandle surf = geom_sets[2][i];
    std::vector<moab::EntityHandle> vols;
    std::vector<int> surf_senses;
    rval = gt.get_senses( surf, vols, surf_senses);
    if(gen::error(moab::MB_SUCCESS!=result,"could not get surface senses")) return result;
    std::cout << "Number of senses for surface " << gen::geom_id_by_handle(surf) << " = " << surf_senses.size() << std::endl;
    for (unsigned int index=0; index<surf_senses.size() ; index++)
    { 
     std::cout << "vol = " << gen::geom_id_by_handle(vols[index]) << std::endl;
     std::cout << "sense = " << sense_printer( surf_senses[index] ) << std::endl;
    }
    std::cout << std::endl;
    }

// Print all Vertex Coordinates

    std::cout << "=====================================" << std::endl;
    std::cout << " Vertex Coordinates " << std::endl;
    std::cout << "=====================================" << std::endl;

    moab::Range verts;
    rval = MBI()->get_entities_by_type(0, moab::MBVERTEX, verts);
    if(gen::error(moab::MB_SUCCESS!=result,"could not get vertex handles")) return result;
    double x[geom_sets[0].size()];
    double y[geom_sets[0].size()];
    double z[geom_sets[0].size()];


    rval = MBI()-> get_coords( verts, &x[0], &y[0], &z[0]);
    if(gen::error(moab::MB_SUCCESS!=result,"could not get coordinates of the vertices")) return result;


    int j=0;
    for (moab::Range::const_iterator i = verts.begin(); i!=verts.end(); i++)
    {
        
     std::cout << "Vertex ID = " << *i << std::endl;
     std::cout << "X = " << x[j] << std::endl;
     std::cout << "Y = " << y[j] << std::endl;
     std::cout << "Z = " << z[j] << std::endl;

     j++;

   }
}
//==========EOL=============//

moab::Interface *MBI() {
  static moab::Core instance;
    return &instance;
  }
