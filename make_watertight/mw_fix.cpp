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

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

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

MBInterface *MBI();


void moab_printer(MBErrorCode error_code)
{
  if ( error_code == MB_INDEX_OUT_OF_RANGE )
    {
      std::cerr << "ERROR: MB_INDEX_OUT_OF_RANGE" << std::endl;
    }
  if ( error_code == MB_MEMORY_ALLOCATION_FAILED )
    {
      std::cerr << "ERROR: MB_MEMORY_ALLOCATION_FAILED" << std::endl;
    }
  if ( error_code == MB_ENTITY_NOT_FOUND )
    {
      std::cerr << "ERROR: MB_ENTITY_NOT_FOUND" << std::endl;
    }
  if ( error_code == MB_MULTIPLE_ENTITIES_FOUND )
    {
      std::cerr << "ERROR: MB_MULTIPLE_ENTITIES_FOUND" << std::endl;
    }
  if ( error_code == MB_TAG_NOT_FOUND )
    {
      std::cerr << "ERROR: MB_TAG_NOT_FOUND" << std::endl;
    }
  if ( error_code == MB_FILE_DOES_NOT_EXIST )
    {
      std::cerr << "ERROR: MB_FILE_DOES_NOT_EXIST" << std::endl;
    }    
  if ( error_code == MB_FILE_WRITE_ERROR )
    {
      std::cerr << "ERROR: MB_FILE_WRITE_ERROR" << std::endl;
    }    
  if ( error_code == MB_ALREADY_ALLOCATED )
    {
      std::cerr << "ERROR: MB_ALREADY_ALLOCATED" << std::endl;
    }    
  if ( error_code == MB_VARIABLE_DATA_LENGTH )
    {
      std::cerr << "ERROR: MB_VARIABLE_DATA_LENGTH" << std::endl;
    }  
  if ( error_code == MB_INVALID_SIZE )
    {
      std::cerr << "ERROR: MB_INVALID_SIZE" << std::endl;
    }  
  if ( error_code == MB_UNSUPPORTED_OPERATION )
    {
      std::cerr << "ERROR: MB_UNSUPPORTED_OPERATION" << std::endl;
    }  
  if ( error_code == MB_UNHANDLED_OPTION )
    {
      std::cerr << "ERROR: MB_UNHANDLED_OPTION" << std::endl;
    }  
  if ( error_code == MB_FAILURE )
    {
      std::cerr << "ERROR: MB_FAILURE" << std::endl;
    }  
  return;
}


MBErrorCode get_geom_size_before_sealing( const MBRange geom_sets[], 
                                          const MBTag geom_tag,
                                          const MBTag size_tag ) {
  MBErrorCode rval;
  for(int dim=1; dim <= 3 ; dim++) {
    std::cout << "dim = " << dim << std::endl;
    for(MBRange::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
      double size = 0;
	std::cout << "*i =" << *i << std::endl;
	std::cout << "geom_tag =" << geom_tag << std::endl;
	std::cout << "size =" << size << std::endl;


      rval = gen::measure( *i, geom_tag, size );
      if(gen::error(MB_SUCCESS!=rval,"could not measure")) return rval;
      rval = MBI()->tag_set_data( size_tag, &(*i), 1, &size );
      if(gen::error(MB_SUCCESS!=rval,"could not set size tag")) return rval;

	std::cout << "*i =" << *i << std::endl;
	std::cout << "geom_tag =" << geom_tag << std::endl;
	std::cout << "size =" << size << std::endl;

    }
  }
  return MB_SUCCESS;
}

MBErrorCode get_senses(MBEntityHandle entity,
    std::vector<MBEntityHandle> &wrt_entities, std::vector<int> &senses)
{
  //
  // the question here is: the wrt_entities is supplied or not?
  // I assume not, we will obtain it !!
  int edim = 1;

  if (-1 == edim)
    return MB_FAILURE;// not geometry entity

  MBErrorCode rval;
  wrt_entities.clear();
  senses.clear();

  if (1 == edim)// edge
  {
    
    const void *dum_ptr;
    int num_ents;
    unsigned flags = MB_TAG_SPARSE;
  
    MBTag senseNEntsTag;
    rval = MBI() -> tag_get_handle(GEOM_SENSE_N_ENTS_TAG_NAME, 0 , MB_TYPE_HANDLE, senseNEntsTag, flags);
    if (gen::error(MB_SUCCESS!=rval, "could not get ent sense handles")) return rval;
    rval = MBI()->tag_get_by_ptr(senseNEntsTag, &entity, 1, &dum_ptr, &num_ents);
    if (gen::error(MB_SUCCESS!=rval, "could not get ent sense data")) return rval;
    const MBEntityHandle *ents_data = static_cast<const MBEntityHandle*> (dum_ptr);
    std::copy(ents_data, ents_data + num_ents, std::back_inserter(wrt_entities));

    MBTag senseNSensesTag;
    rval = MBI()->tag_get_handle(GEOM_SENSE_N_SENSES_TAG_NAME, 0 , MB_TYPE_INTEGER, senseNSensesTag, flags);
    if (gen::error(MB_SUCCESS!=rval, "could not get senses handle")) return rval;
    rval = MBI()->tag_get_by_ptr(senseNSensesTag, &entity, 1, &dum_ptr,
        &num_ents);
    if (gen::error(MB_SUCCESS!=rval, "could not get senses data")) return rval;

    const int *senses_data = static_cast<const int*> (dum_ptr);
    std::copy(senses_data, senses_data + num_ents, std::back_inserter(senses));

  }/* else // face in volume, edim == 2
  {
    
    MBEntityHandle sense_data[2] = { 0, 0 };
    rval = MBI()->tag_get_data(GEOM_SENSE_2_TAG_NAME, &entity, 1, sense_data);
    if (MB_SUCCESS != rval)
      return rval;
    if (sense_data[0] != 0 && sense_data[1] == sense_data[0]) {
      wrt_entities.push_back(sense_data[0]);
      senses.push_back(0);// both
    } else {
      if (sense_data[0] != 0) {
        wrt_entities.push_back(sense_data[0]);
        senses.push_back(1);
      }
      if (sense_data[1] != 0) {
        wrt_entities.push_back(sense_data[1]);
        senses.push_back(-1);
      }

    }

  }
  */
  // filter the results with the sets that are in the model at this time
  // this was introduced because extracting some sets (e.g. neumann set, with mbconvert)
  //   from a model would leave some sense tags not defined correctly
  // also, the geom ent set really needs to be part of the current model set
 /*
  unsigned int currentSize =0;

  for (unsigned int index=0; index<wrt_entities.size(); index++)
  {
    MBEntityHandle wrt_ent=wrt_entities[index];
    if (wrt_ent )
    {
      if (MBI()->contains_entities(modelSet, &wrt_ent, 1))
      {
        wrt_entities[currentSize] = wrt_entities[index];
        senses[currentSize] = senses[index];
        currentSize++;
      }
    }
  }
  wrt_entities.resize(currentSize);
  senses.resize(currentSize);
  //
  */
  return MB_SUCCESS;
}
  int main(int argc, char **argv) {

   // ******************************************************************
    // Load the h5m file and create tags.
    // ******************************************************************

    clock_t start_time = clock();
    const bool debug = false;
    const bool check_geom_size = true;

    // check input args
    if( 2 > argc || 3 < argc ) 
      {
	std::cout << "To zip a faceted h5m file:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.h5m>" << std::endl;
	std::cout << "To facet and zip an ACIS file using the default facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat>" << std::endl;
	std::cout << "To facet and zip an ACIS file using a specified facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat> <facet_tolerance>" << std::endl;
	return MB_FAILURE;
      }

    // The root name does not have an extension
    std::string input_name = argv[1];
    std::string root_name = argv[1];
    int len = root_name.length();
    root_name.erase(len - 4);
    bool is_acis;

    // load the input file
    MBErrorCode result, rval;
    MBEntityHandle input_set;

    rval = MBI()->create_meshset( MESHSET_SET, input_set );

    if(gen::error(MB_SUCCESS!=rval,"failed to create_meshset"))
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
	if(gen::error(MB_SUCCESS!=rval,"failed to load_file 0")) 
	  {
	    return rval;      
	  }
	
	is_acis = false;

    // If reading a sat file, the facet toleance will default to 1e-3 if it is
    // not specified. If the user does not specify a facet_tol, default to 1e-3.
    // This is the same as what ReadCGM uses.
      } 

    /*
     // recreate to only perform these operations on h5m meshes  
    else if(std::string::npos!=input_name.find("sat") && 
	      ((2==argc) || (3==argc)) ) 
      {
	double facet_tol;
	if(3 == argc) 
	  {
	    facet_tol = atof(argv[2]);
	  }
	else 
	  {
	    facet_tol = 1e-3;
	  }

	std::string options;
	options += "FACET_DISTANCE_TOLERANCE=";
	std::stringstream facet_tol_ss;
	facet_tol_ss << facet_tol; 
	options += facet_tol_ss.str();
	if(debug) std::cout << "  options=" << options << std::endl;
	rval = MBI()->load_file( input_name.c_str(), &input_set, options.c_str() );
	if(gen::error(MB_SUCCESS!=rval,"failed to load_file 1")) return rval;      

      // write an HDF5 file of facets with known tolerance   
	std::string facet_tol_filename = root_name + "_" + facet_tol_ss.str() + ".h5m";
	rval = MBI()->write_mesh( facet_tol_filename.c_str() );
	if(gen::error(MB_SUCCESS!=rval,"failed to write_mesh 0")) return rval;      
	is_acis = true;
      } 
    else 
      {
	std::cout << "incorrect input arguments" << std::endl;
	return MB_FAILURE;
      }
     //not required if  only doing this with h5m files
     */

    // create tags
    clock_t load_time = clock();    
    MBTag geom_tag, id_tag, normal_tag, merge_tag, faceting_tol_tag, 
      geometry_resabs_tag, size_tag, orig_curve_tag;
  
    result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
				MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
    assert( MB_SUCCESS == result );
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1,
				MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert( MB_SUCCESS == result );
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "NORMAL", sizeof(MBCartVect), MB_TYPE_OPAQUE,
        normal_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert( MB_SUCCESS == result );
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "MERGE", 1, MB_TYPE_HANDLE,
        merge_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
    assert( MB_SUCCESS == result ); 
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      } 
    result = MBI()->tag_get_handle( "FACETING_TOL", 1, MB_TYPE_DOUBLE,
        faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
    assert( MB_SUCCESS == result );  
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "GEOMETRY_RESABS", 1,     MB_TYPE_DOUBLE,
                             geometry_resabs_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT  );
    assert( MB_SUCCESS == result );  
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    result = MBI()->tag_get_handle( "GEOM_SIZE", 1, MB_TYPE_DOUBLE,
				    size_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT  );
    assert( (MB_SUCCESS == result) );
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    int true_int = 1;    
    result = MBI()->tag_get_handle( "ORIG_CURVE", 1,
				MB_TYPE_INTEGER, orig_curve_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT, &true_int );
    assert( MB_SUCCESS == result );
    if ( result != MB_SUCCESS )
      {
	moab_printer(result);
      }
    // PROBLEM: MOAB is not consistent with file_set behavior. The tag may not be
    // on the file_set.
    MBRange file_set;
    result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &faceting_tol_tag,
                                                  NULL, 1, file_set );

    if(gen::error(MB_SUCCESS!=result,"could not get faceting_tol_tag")) 
      {
	return result;
      }

    gen::error(file_set.empty(),"file set not found");

    if(gen::error(1!=file_set.size(),"Refacet with newer version of ReadCGM.")) 
      {
	return MB_FAILURE;
      }

    double facet_tol, sme_resabs_tol=1e-6;
    result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1,  
                                  &facet_tol );
    assert(MB_SUCCESS == result);
    result = MBI()->tag_get_data( geometry_resabs_tag, &file_set.front(), 1,  
                                  &sme_resabs_tol );
    if(MB_SUCCESS != result) 
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
 
    
    // get all geometry sets
    MBRange geom_sets[4];
    for(unsigned dim=0; dim<4; dim++) 
      {
	void *val[] = {&dim};
	result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
	std::cout << "Get entities by type and tag" << std::endl;

	assert(MB_SUCCESS == result);

	// make sure that sets TRACK membership and curves are ordered
	// MESHSET_TRACK_OWNER=0x1, MESHSET_SET=0x2, MESHSET_ORDERED=0x4
	for(MBRange::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) 
	  {
	    unsigned int options;
	    result = MBI()->get_meshset_options(*i, options );
	    assert(MB_SUCCESS == result);
    
	    // if options are wrong change them
	    if(dim==1) 
	      {
		if( !(MESHSET_TRACK_OWNER&options) || !(MESHSET_ORDERED&options) ) 
		  {
		    result = MBI()->set_meshset_options(*i, MESHSET_TRACK_OWNER|MESHSET_ORDERED);
		    assert(MB_SUCCESS == result);
		  }
	      } 
	    else 
	      {
		if( !(MESHSET_TRACK_OWNER&options) ) 
		  {        
		    result = MBI()->set_meshset_options(*i, MESHSET_TRACK_OWNER);
		    assert(MB_SUCCESS == result);
		  }
	      }
	  }
      }

    std::cout << "I am here" << std::endl;

    // this could be related to when there are sat files rather than mesh?
    // If desired, find each entity's size before sealing.
    ///*
    if(check_geom_size) 
      {
	std::cout << "I am checking the geometry size" << std::endl;
	result = get_geom_size_before_sealing( geom_sets, geom_tag, size_tag );
	if(gen::error(MB_SUCCESS!=result,"measuring geom size failed"))
	  {
	    return result;
	  }
      }
    
    //*/

    std::cout << "Get entity count before sealing" << std::endl;
    // Get entity count before sealing.
    int orig_n_tris;
    result = MBI()->get_number_entities_by_type( 0, MBTRI, orig_n_tris );
    std::cout << result << std::endl;

    assert(MB_SUCCESS == result);

    std::cout << "  input faceted geometry contains " << geom_sets[3].size() << " volumes, " 
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size() 
              << " curves, and " << orig_n_tris << " triangles" << std::endl;  


    //Print all geometry entities
    std::cout << "Surfaces" << std::endl;
    for (unsigned int index=0; index < geom_sets[2].size(); index++)
    {
      std::cout << "surface handle = " << geom_sets[2][index] << std::endl;
      std::cout << "surface id = " << gen::geom_id_by_handle(geom_sets[2][index]) << std::endl;
    }
    
    std::cout << "Curves" << std::endl;
    for (unsigned int index=0; index < geom_sets[1].size(); index++)
    {
      std::cout << "curve handle = " << geom_sets[1][index] << std::endl;
      std::cout << "curve id = " << gen::geom_id_by_handle(geom_sets[1][index]) << std::endl;
    }
    
    std::cout << "Volumes" << std::endl;
    for (unsigned int index=0; index < geom_sets[3].size(); index++)
    {
      std::cout << "volume handle = " << geom_sets[3][index] << std::endl;
      std::cout << "volume id = " << gen::geom_id_by_handle(geom_sets[3][index]) << std::endl;
    }

// Trying to get sense data

    //Sense Handle
    std::cout << std::endl << "Getting sense_handle directly..." << std::endl;
    MBTag sense_handle; 
    unsigned flags = MB_TAG_SPARSE;
    result = MBI()->tag_get_handle(GEOM_SENSE_2_TAG_NAME, 2 , MB_TYPE_HANDLE, sense_handle, flags );
    if (gen::error(MB_SUCCESS!=result, "could not get sense tag handle")) return result;
    std::cout << "sense_handle = " << sense_handle << std::endl;

   //Sense Name

    std::string sense_name;
    result = MBI()-> tag_get_name(sense_handle, sense_name);
    if(gen::error(MB_SUCCESS!=result, "could not get sense_name")) return result;
    std::cout << "sense name = " << sense_name << std::endl;


   //Get Sense Data
    MBEntityHandle test_surf= geom_sets[2][0];
    const void *dum_ptr;
    int size;
    MBEntityHandle sense_data[2]= {0 , 0};
    
    //result = MBI() -> tag_get_data( sense_handle, &test_surf, 1, 


    //USING POINTER METHOD
    result = MBI() -> tag_get_by_ptr( sense_handle, &test_surf , 1, &dum_ptr, &size);
    if (gen::error(MB_SUCCESS!=result, "could not retrieve sense data")) return result;

    const int *senses_data = static_cast<const int*> (dum_ptr);
   
    std::cout << "senses_data = " << *senses_data << std::endl;
    
    std::vector<int> senses;

    std::copy(senses_data, senses_data + size, std::back_inserter(senses));

    for ( unsigned int index=0; index < senses.size() ; index ++)
    {
     std::cout << "senses[" << index << "] = " << senses[index] << std::endl;
    }
    
    MBEntityHandle test_edge = geom_sets[1][0];
    std::vector<MBEntityHandle> wrt_entities;
    std::vector<int> sens;   
    result = get_senses( test_edge, wrt_entities, sens);
    if (gen::error( MB_SUCCESS != result, "could not get_senses")) return result;
    
    std::cout << std::endl << "size of data returned = " << sens.size() << std::endl;
    
    for(unsigned int index=0; index< sens.size() ; index ++)
    {
     std::cout << "wrt_entities[" << index << "] = " << wrt_entities[index] << std::endl;
     std::cout << "sense[" << index << "] = " << sens[index] << std::endl;
    }
      

    

  
}
//==========EOL=============//

MBInterface *MBI() {
    static MBCore instance;
    return &instance;
  }
