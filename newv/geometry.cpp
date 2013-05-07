#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/GeomTopoTool.hpp"

#include "MBSkinner.hpp"
#include "DagMC.hpp"

#include "geometry.hpp"

namespace geometry 
{
  // measure the property of the geometric entity pointed to
  MBErrorCode measure(const MBEntityHandle set, const MBTag geom_tag, double &size ) 
  {
    MBErrorCode rval;
    int dimension;
    rval = MOAB()->tag_get_data( geom_tag, &set, 1, &dimension );  // get the tag data geom_tag for the set pointed to               
    if( dimension == 0 )  // we dont do anything with vertices
      {
	std::cout << "measure: cannot measure vertex" << std::endl;
	return MB_FAILURE;
      } 
    else if( 1 == dimension ) // curves
      {
	std::vector<MBEntityHandle> vctr;
	rval = MOAB()->get_entities_by_handle( set, vctr );
	if (rval != MB_SUCCESS )
	  {
	    return rval;
	  }
	size = length(vctr);
      } 
    else if(2 == dimension) // facets
      {
	MBRange tris;
	rval = MOAB()->get_entities_by_type( set, MBTRI, tris );
	assert(MB_SUCCESS == rval);
	size = triangle_area( tris );
      } 
    else if(3 == dimension) // volumes
      {
	//	moab::DagMC &dagmc = *moab::DagMC::instance( MOAB() );
	//rval = dagmc.measure_volume( set, size );
	//	std::cout << "in measure volume" << std::endl;
	rval = measure_volume( set, size );
	//	std::cout << "in measure volume" << std::endl;
	if(MB_SUCCESS != rval) 
	  {
	    std::cout << "rval=" << rval <<  std::endl;
	  }
      } 
    else 
      {
	std::cout << "measure: incorrect dimension" << std::endl;
	return MB_FAILURE;
      }
    return MB_SUCCESS;

  }

  double length( std::vector<MBEntityHandle> edges ) 
  {
    MBErrorCode result;
    std::vector<MBEntityHandle>::iterator i; 
    double dist = 0;
    MBEntityType type = MOAB()->type_from_handle( edges[0] ); 
    
    // if there is no work return
    if(edges.empty()) 
      {
	return 0;
      }
    
    // if vector has both edges and verts, only use edges
    // NOTE: The curve sets from ReadCGM do not contain duplicate endpoints for loops!
    MBEntityType end_type = MOAB()->type_from_handle( edges.back() );
    if(type != end_type) 
      {
	for(std::vector<MBEntityHandle>::iterator i=edges.begin(); i!=edges.end(); i++) 
	  {
	    if(MBVERTEX == MOAB()->type_from_handle( *i )) 
	      {
		i = edges.erase(i) - 1;
	      }
	  }
      }
    // determine if vector defines an arc by edges of verts
    type = MOAB()->type_from_handle( edges[0] ); 
    if (MBEDGE == type) 
      {
	if(edges.empty()) 
	  {
	    return 0.0; 
	  }
	for( i=edges.begin(); i!=edges.end(); i++ ) 
	  {
	    int n_verts;
	    const MBEntityHandle *conn;
	    result = MOAB()->get_connectivity( *i, conn, n_verts );
	    if( MB_SUCCESS!=result ) 
	      {
		std::cout << "result=" << result << std::endl; 
	      }
	    assert(MB_SUCCESS == result);
	    assert( 2 == n_verts );
	    if(conn[0] == conn[1])
	      {
		continue;
	      }
	    dist += dist_between_verts( conn[0], conn[1] );
	    //std::cout << "length: " << dist << std::endl;
	  }
      } 
    else if (MBVERTEX == type) 
      {
	if(2 > edges.size()) 
	  {
	    return 0.0;
	  }
	MBEntityHandle front_vert = edges.front();
	for( i=edges.begin()+1; i!=edges.end(); i++) 
	  {
	    dist += dist_between_verts( front_vert, *i );
	    front_vert = *i;
	  }
      } 
    else 
      {
	return MB_FAILURE;
      }

    return dist;
  }

  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 ) 
  {
    MBCartVect v2 = v0 - v1;
    return v2.length();
  }

  MBErrorCode dist_between_verts( const MBEntityHandle v0, const MBEntityHandle v1, double &d) 
  {
    MBErrorCode result;
    MBCartVect coords0, coords1;
    result = MOAB()->get_coords( &v0, 1, coords0.array() );
    if(MB_SUCCESS != result) 
      {
	std::cout << "dist_between_verts: get_coords on v0=" << v0 << " result=" 
		  << result << std::endl; 
	return result;
      }
    result = MOAB()->get_coords( &v1, 1, coords1.array() );
    if(MB_SUCCESS != result) 
      {
	std::cout << "dist_between_verts: get_coords on v1=" << v1 << " result=" 
		  << result << std::endl; 
	return result;
      }
    d = dist_between_verts( coords0, coords1 );
    return MB_SUCCESS;
  }

  double dist_between_verts( double coords0[], double coords1[] ) 
  {
    return sqrt( (coords0[0]-coords1[0])*(coords0[0]-coords1[0]) +
		 (coords0[1]-coords1[1])*(coords0[1]-coords1[1]) +
		 (coords0[2]-coords1[2])*(coords0[2]-coords1[2]) );
  }

  double dist_between_verts( MBEntityHandle vert0, MBEntityHandle vert1 ) 
  {
    double coords0[3], coords1[3];
    MBErrorCode result;
    result = MOAB()->get_coords( &vert0, 1, coords0 );
    if(MB_SUCCESS!=result)
      {
	std::cout << "result=" << result << " vert=" 
		  << vert0 << std::endl;
      }
    assert(MB_SUCCESS == result);
    result = MOAB()->get_coords( &vert1, 1, coords1 );
    if(MB_SUCCESS!=result) 
      {
	std::cout << "result=" << result << " vert=" 
		  << vert1 << std::endl;
      }
    assert(MB_SUCCESS == result);
    return dist_between_verts( coords0, coords1 );
  }

double triangle_area( const MBCartVect a, const MBCartVect b, 
                        const MBCartVect c) {
    MBCartVect d = c - a;
    MBCartVect e = c - b;
    MBCartVect f = d*e;
    return 0.5*f.length();
  }
  MBErrorCode triangle_area( const MBEntityHandle conn[], double &area ) {
    MBCartVect coords[3];
    MBErrorCode result = MOAB()->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == result);
    area = triangle_area( coords[0], coords[1], coords[2] );
    return MB_SUCCESS;
  }
  MBErrorCode triangle_area( const MBEntityHandle tri, double &area ) {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;
    result = MOAB()->get_connectivity( tri, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    
    result = triangle_area( conn, area );
    assert(MB_SUCCESS == result);
    return MB_SUCCESS;
  }
  double triangle_area( const MBRange tris ) {
    double a, area = 0;
    MBErrorCode result;
    for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
      result = triangle_area( *i, a);
      assert(MB_SUCCESS == result);
      area += a;
    }
    return area;
  }

  // calculate volume of polyhedron
  MBErrorCode measure_volume( MBEntityHandle volume, double& result )
  {
    MBErrorCode rval;
    std::vector<MBEntityHandle> surfaces, surf_volumes;
    result = 0.0;
  
    // get surfaces from volume
    rval = MOAB()->get_child_meshsets( volume, surfaces );
    if (MB_SUCCESS != rval)
      {
	return rval;
      }
  
    // get surface senses
    std::vector<int> senses( surfaces.size() );
    rval = surface_sense ( volume, surfaces.size(), &surfaces[0], &senses[0] );
    if (MB_SUCCESS != rval) 
      {
	std::cerr << "ERROR: Surface-Volume relative sense not available. "
		  << "Cannot calculate volume." << std::endl;
	return rval;
      }
  
    for (unsigned i = 0; i < surfaces.size(); ++i) 
      {
	// skip non-manifold surfaces
	if (!senses[i])
	  {
	    continue;
	  }
    
	// get triangles in surface
	MBRange triangles;
	rval = MOAB()->get_entities_by_dimension( surfaces[i], 2, triangles );
	if (MB_SUCCESS != rval) 
	  {
	    return rval;
	  }
	if (!triangles.all_of_type(MBTRI)) 
	  {
	    std::cout << "WARNING: Surface " << surfaces[i] 
		      << " contains non-triangle elements. Volume calculation may be incorrect." 
		      << std::endl;
	    triangles.clear();
	    rval = MOAB()->get_entities_by_type( surfaces[i], MBTRI, triangles );
	    if (MB_SUCCESS != rval) return rval;
	  }
    
	// calculate signed volume beneath surface (x 6.0)
	double surf_sum = 0.0;
	const MBEntityHandle *conn;
	int len;
	MBCartVect coords[3];
	for (MBRange::iterator j = triangles.begin(); j != triangles.end(); ++j) 
	  {
	    rval = MOAB()->get_connectivity( *j, conn, len, true );
	    if (MB_SUCCESS != rval)
	      {
		return rval;
	      }
	    assert(3 == len);
	    rval = MOAB()->get_coords( conn, 3, coords[0].array() );
	    if (MB_SUCCESS != rval) 
	      {
		return rval;
	      }
	    coords[1] -= coords[0];
	    coords[2] -= coords[0];
	    surf_sum += (coords[0] % (coords[1] * coords[2]));
	  }
	result += senses[i] * surf_sum;
      }
  
    result /= 6.0;
    return MB_SUCCESS;
  }

  MBErrorCode surface_sense( MBEntityHandle volume, 
                           int num_surfaces,
                           const MBEntityHandle* surfaces,
                           int* senses_out )
  {
    std::vector<MBEntityHandle> surf_volumes( 2*num_surfaces );
    MBTag senseTag = get_tag( "GEOM_SENSE_2", 2, MB_TAG_SPARSE, MB_TYPE_HANDLE, NULL, false );
    MBErrorCode rval = MOAB()->tag_get_data( senseTag , surfaces, num_surfaces, &surf_volumes[0] );
    if (MB_SUCCESS != rval)
      {
	return rval;
      }
  
    const MBEntityHandle* end = surfaces + num_surfaces;
    std::vector<MBEntityHandle>::const_iterator surf_vols = surf_volumes.begin();
    while (surfaces != end) 
      {
	MBEntityHandle forward = *surf_vols; ++surf_vols;
	MBEntityHandle reverse = *surf_vols; ++surf_vols;
	if (volume == forward) 
	  {
	    *senses_out = (volume != reverse); // zero if both, otherwise 1
	  }
	else if (volume == reverse)
	  {
	    *senses_out = -1;
	  }
	else 
	  {
	    return MB_ENTITY_NOT_FOUND;
	  }
    
	++surfaces;
	++senses_out;
      }
  
    return MB_SUCCESS;
  }

  // get sense of surface(s) wrt volume
  MBErrorCode surface_sense( MBEntityHandle volume, 
                             MBEntityHandle surface,
                             int& sense_out )
  {
    // get sense of surfaces wrt volumes
    MBEntityHandle surf_volumes[2];
    MBTag senseTag = get_tag( "GEOM_SENSE_2", 2, MB_TAG_SPARSE, MB_TYPE_HANDLE, NULL, false );
    MBErrorCode rval = MOAB()->tag_get_data( senseTag , &surface, 1, surf_volumes );
    if (MB_SUCCESS != rval) 
      {
	return rval;
      }
  
    if (surf_volumes[0] == volume)
      {
	sense_out = (surf_volumes[1] != volume); // zero if both, otherwise 1
      }
    else if (surf_volumes[1] == volume)
      {
	sense_out = -1;
      }
    else
      {
	return MB_ENTITY_NOT_FOUND;
      }
  
  return MB_SUCCESS;
 }

 MBTag get_tag( const char* name, int size, MBTagType store,
		     MBDataType type, const void* def_value,
		     bool create_if_missing)
 {
   MBTag retval = 0;
   unsigned flags = store|moab::MB_TAG_CREAT;
   if (!create_if_missing)
     {
       flags |= moab::MB_TAG_EXCL;
     }
   MBErrorCode result = MOAB()->tag_get_handle(name, size, type, retval, flags, def_value);
   if (create_if_missing && MB_SUCCESS != result)
     {
       std::cerr << "Couldn't find nor create tag named " << name << std::endl;
     }

   return retval;
 }

  bool triangle_degenerate( const MBEntityHandle tri ) 
  {
    MBErrorCode result;
    const MBEntityHandle *conn;
    int n_verts;
    result = MOAB()->get_connectivity( tri, conn, n_verts );
    assert(MB_SUCCESS == result);
    assert(3 == n_verts);
    return triangle_degenerate( conn[0], conn[1], conn[2] );
  }

  bool triangle_degenerate( const MBEntityHandle v0, const MBEntityHandle v1,
			    const MBEntityHandle v2 ) 
  { 
    if(v0==v1 || v1==v2 || v2==v0) return true;
    return false;
  }
}

