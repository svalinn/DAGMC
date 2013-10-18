// ********************************************************************
// Patrick Shriwise
// October, 2013

// functions needed to seal unwatertight models in make_watertight


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

#include "mw_func.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"


const char GEOM_SENSE_2_TAG_NAME[] = "GEOM_SENSE_2";
const char GEOM_SENSE_N_ENTS_TAG_NAME[] = "GEOM_SENSE_N_ENTS";
const char GEOM_SENSE_N_SENSES_TAG_NAME[] = "GEOM_SENSE_N_SENSES"; 


MBInterface *MOAB();

namespace mw_func {

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

MBErrorCode delete_all_edges() {
  // delete all of the edges. Never keep edges. They are too hard to track and use
  // due to orientation and multiple entities errors when merging.
  MBErrorCode result;
  MBRange edges;
  result = MBI()->get_entities_by_type( 0, MBEDGE, edges );
  if(gen::error(MB_SUCCESS!=result,"could not get edges")) return result;
  assert(MB_SUCCESS == result);
  result = MBI()->delete_entities( edges );
  if(gen::error(MB_SUCCESS!=result,"could not delete edges")) return result;
  assert(MB_SUCCESS == result);
  return MB_SUCCESS;
}

MBErrorCode find_degenerate_tris() {
  MBErrorCode result;
  MBRange tris;
  result = MBI()->get_entities_by_type( 0, MBTRI, tris );
  if(gen::error(MB_SUCCESS!=result,"could not get tris")) return result;
  assert(MB_SUCCESS == result);
  int counter = 0;
  for(MBRange::const_iterator i=tris.begin(); i!=tris.end(); ++i) {
    if( gen::triangle_degenerate(*i) ) {
      result = MBI()->list_entity(*i);
      if(gen::error(MB_SUCCESS!=result,"found degenerate tri")) return result;
      assert(MB_SUCCESS == result);
      ++counter;
    }
  }
  if(counter != 0)
  {
  std::cout << "Found " << counter << " degenerate triangles. " << std::endl;
  }
  return MB_SUCCESS;
}

// Input: Possibly unordered sets of curves that do track ownership. Curves contain
//        edges and vertices. Parents sets are surfaces. Child sets are endpoint
//        vertices.
// Output: Ordered sets of verts that do track ownership. All edges are deleted.
MBErrorCode prepare_curves(MBRange &curve_sets, 
                           MBTag geom_tag, MBTag id_tag, MBTag merge_tag, 
                           const double FACET_TOL, const bool debug, bool verbose ) {
  MBErrorCode result;
  if (verbose) std::cout << "Modifying faceted curve representation and removing small curves..." 
            << std::endl;

  // process each curve
  for(MBRange::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++ ) {
    // get the curve id of the curve meshset
    int id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &id );
    if(gen::error(MB_SUCCESS!=result,"could not get id tag")) return result;
    if(debug) std::cout << "curve " << id << std::endl;

    // get the range of edges of the curve meshset
    std::vector<MBEntityHandle> curve_edges;
    result = MBI()->get_entities_by_type( *i, MBEDGE, curve_edges );
    if(gen::error(MB_SUCCESS!=result,"could not get curve_edges")) return result;

    /* Merge the endpoints of the curve and remove its edges if it is too small.
    Use the MERGE_TOL because these edges will be merged with the MERGE_TOL 
    during zipping anyhow. Doing this now removes small curves from zipping and
    reduces ambiguity. */
    if(FACET_TOL > gen::length(curve_edges)) {
      std::cout << "  deleted curve " << id << ", length=" << gen::length(curve_edges) 
                << " cm, n_verts=" << curve_edges.size()+1 << std::endl;

      // get the endpoints of the curve
      MBRange endpt_sets;
      result = MBI()->get_child_meshsets( *i, endpt_sets );
      if(gen::error(MB_SUCCESS!=result,"could not get curve child sets")) return result;

      if(endpt_sets.empty()) {
        if(gen::error(true,"curve has no child sets")) return result;
      } else if(1 == endpt_sets.size()) {
        // The edges are no longer needed. Remove them before altering the range
        // by deleting degenerate edges below.
        result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
        if(gen::error(MB_SUCCESS!=result,"could not delete edges")) return result;

      } else if(2 == endpt_sets.size()) {
        // The edges are no longer needed. Remove them before altering the range
        // by deleting degenerate edges below.
        result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
        if(gen::error(MB_SUCCESS!=result,"could not delete edges")) return result;

        MBRange front_endpt, back_endpt;
        result = MBI()->get_entities_by_type( endpt_sets.front(), MBVERTEX, front_endpt);
        if(gen::error(MB_SUCCESS!=result,"could not get vert from front endpt set")) return result;
        if(gen::error(1!=front_endpt.size(),"front endpt set does not have 1 vert")) return result;

        result = MBI()->get_entities_by_type( endpt_sets.back(), MBVERTEX, back_endpt);
        if(gen::error(MB_SUCCESS!=result,"could not get vert from back endpt set")) return result;
        if(gen::error(1!=back_endpt.size(),"back endpt set does not have 1 vert")) return result;

        // merge the endpoints-ALWAYS CHECK TO AVOID MERGING THE SAME ENTITY!!!
        if(front_endpt[0] != back_endpt[0]) {
	  std::vector<MBEntityHandle> temp;
          result = zip::merge_verts( front_endpt.front(), back_endpt.front(), temp, temp );
          if(gen::error(MB_SUCCESS!=result,"could not merge verts")) return result;

          // check for and remove degenerate edges caused by the merge
          MBRange edges;
          MBEntityHandle temp_pt = front_endpt[0];
          result = MBI()->get_adjacencies( &temp_pt, 1, 1, false, edges);
          if(gen::error(MB_SUCCESS!=result,"could not get adj edges")) return result;

          for(MBRange::iterator j=edges.begin(); j!=edges.end(); j++) {
            const MBEntityHandle *conn;
            int n_verts;
            result = MBI()->get_connectivity( *j, conn, n_verts); 
            if(gen::error(MB_SUCCESS!=result,"could not get edge conn")) return result;

            if(conn[0] == conn[1]) {
              result = MBI()->delete_entities( &(*j), 1 );
              if(gen::error(MB_SUCCESS!=result,"could not delete degenerate edge")) return result;
            }
          }
	}
      } else {
        assert(false);
      }
      // It is possible that the endpoints themselves are orphaned. Should these
      // be deleted?

      // Remove the curve set. This also removes parent-child relationships.
      result = MBI()->delete_entities( &(*i), 1);
      if(gen::error(MB_SUCCESS!=result,"could not delete curve set")) return result;
      i = curve_sets.erase(i) - 1;
    } else {

      // convert the curve of edges into a curve of verts 
      std::vector<MBEntityHandle> ordered_verts;
      result = gen::ordered_verts_from_ordered_edges( curve_edges, ordered_verts);
      if(gen::error(MB_SUCCESS!=result,"could not order_verts_by_edge")) return result;

      // replace the unordered edges with the ordered verts
      result = arc::set_meshset( *i, ordered_verts );
      if(gen::error(MB_SUCCESS!=result,"could not set_meshset")) return result;
      
      // The edges are no longer needed.
      result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
      if(gen::error(MB_SUCCESS!=result,"could not delete edges")) return result;
    }
  }

  // merge curves that are the same within facet_tol 
  if (verbose) std::cout << "Identifying coincident curves to be merged..." << std::endl;
  result = arc::merge_curves(curve_sets, FACET_TOL, id_tag, merge_tag, debug );
  if(gen::error(MB_SUCCESS!=result,"could not merge_curves")) return result;

  return MB_SUCCESS;
}

/* The chosen curve may need to be reversed to be in the same direction as the skin.
   The skin itself is not changed. */
    /* PROBLEM: Some sliver surfaces are smaller in every dimension than the 5x 
     facet_tol that the skin was found to vary within.   
     This method finds the average distance between a curve and skin arc. The
     smallest distance is the next curve. */
// Cut an arc out of the skin. Return a corresponding curve, the curve's set, and
// is the curve has ben reversed. The returned skin has the arc cut away. The
// returned vector of curve sets has the curve removed.
MBErrorCode create_arc_pair(  const double FACET_TOL,
                              const MBEntityHandle surf_set,
			      std::vector<MBEntityHandle> &skin_loop,
			      std::vector<MBEntityHandle> &curve_sets,
			      const MBEntityHandle front_endpt,
                              const bool debug,
			      MBEntityHandle &curve_set,
			      bool &curve_is_reversed,
			      std::vector<MBEntityHandle> &curve,
			      std::vector<MBEntityHandle> &skin_arc ) {

  /* Now we have a topological connection between the curve and skin. Find other
     curves that share this vert. Be aware if the curve is reversed so the
     original direction can be preserved when done zipping. */
  MBErrorCode rval;
  if(debug) {
    std::cout << curve_sets.size() << " curves remain to be zipped" 
              << " skin_loop.size()=" << skin_loop.size() << " front_endpt="
              << front_endpt << " skin:" << std::endl; 
    //gen::print_loop(skin_loop);
  }

  // initialize stuff
  double min_dist = std::numeric_limits<double>::max();
  curve_set = 0;
  double curve_set_idx, skin_pos;;
  curve.clear();
  skin_arc.clear();
  skin_arc.reserve( skin_loop.size() );
  curve.reserve(    skin_loop.size() );
  
  // Compare all curves, keeping the best pair
  for(unsigned i=0; i<curve_sets.size(); ++i) {
    // get geometric vertex sets
    MBRange endpt_sets;
    rval = MBI()->get_child_meshsets(curve_sets[i], endpt_sets );
    if(gen::error(MB_SUCCESS!=rval,"could not get endpt_sets")) return rval;
    if(gen::error(endpt_sets.empty() || 2<endpt_sets.size(),
                  "too many endpt_sets")) return MB_FAILURE;
    // get the vertex handles
    std::vector<MBEntityHandle> endpts;
    for(unsigned j=0; j<endpt_sets.size(); ++j) {
      MBRange endpt;
      rval = MBI()->get_entities_by_type( endpt_sets[j], MBVERTEX, endpt );
      if(gen::error(MB_SUCCESS!=rval,"could not get endpt")) return rval;
      if(gen::error(1!=endpt.size(),"not one endpt")) return MB_FAILURE;
      endpts.push_back( endpt.front() );
      if(debug) std::cout << "curve " << gen::geom_id_by_handle(curve_sets[i]) 
                          << " endpt=" << endpt.front() << std::endl;
    }

    // if an endpt is not the front_endpt, then this curve isn't adjacent
    if(front_endpt!=endpts.front() && front_endpt!=endpts.back()) continue;

    // get the point representation
    std::vector<MBEntityHandle> temp_curve;
    rval = arc::get_meshset( curve_sets[i], temp_curve );
    if(gen::error(MB_SUCCESS!=rval,"could not get curve set")) return rval;
    //if(gen::error(2>temp_curve.size(),"curve is degenerate")) return MB_FAILURE;
    if(2>temp_curve.size()) std::cout << "warning11: curve is degenerate" << std::endl;
    if(debug) {
      std::cout << "  adj curve " << gen::geom_id_by_handle(curve_sets[i]) << ":" << std::endl; 
      //gen::print_loop( temp_curve );
    }

    // Get the curve-surface relative sense. From CGMA/builds/dbg/include/CubitDefines,
    
    int sense;
    if(debug)
    {
    std::cout << "surf_set = " << gen::geom_id_by_handle(surf_set) << std::endl;
    std::cout << "curve_set = " << gen::geom_id_by_handle(curve_sets[i]) << std::endl;
    }
    rval = gen::get_curve_surf_sense( surf_set, curve_sets[i], sense );
    if(gen::error(MB_SUCCESS!=rval,"could not get_curve_surf_sense")) return rval;

    // get the curve length, for efficient find_closest_vert, plus a tolerance.
    // This also helps to find the correct skin arc in special (~1D surfs) cases
    // where the closest skin pt is not the correct skin.
    const double temp_curve_len = gen::length(temp_curve);
    const double extra = 1.0;

    // is it forward oriented? To allow for ambiguous cases, instead of accepting
    // only forward-oriented curves, only reject reverse-oriented curves.
    if(front_endpt==temp_curve.front() && SENSE_REVERSE!=sense) {
      // find the closest skin pt to the curve endpt
      unsigned pos;
      // if only one curve is left, take the entire skin. Otherwise, due to 
      // coordinate drift via merging the closest skin pt may not be the correct one.
      if(1==curve_sets.size()) {
        pos = skin_loop.size()-1;
      } else {
        rval = gen::find_closest_vert( temp_curve.back(), skin_loop, pos, temp_curve_len+extra );
        if(gen::error(MB_SUCCESS!=rval,"could not find_closest_vert")) return rval;
      }
      if(debug) std::cout << "  end of skin arc=" << skin_loop[pos] << std::endl; 
      // SPECIAL CASE: If the skin is a circle, create an arc out of the circle
      if(skin_loop[pos] == skin_loop.front()) pos = skin_loop.size()-1;

      // create a skin arc to test against the curve
      std::vector<MBEntityHandle> temp_skin(skin_loop.begin(), skin_loop.begin()+pos+1); 
      double d;
      rval = gen::dist_between_arcs( debug, temp_skin, temp_curve, d );
      if(gen::error(MB_SUCCESS!=rval,"could not get dist_between_arcs")) return rval;
      if(debug) std::cout << " curve-skin dist=" << d << std::endl;
 
      // if less than the min_dist, this curve is the best (thus far)
      if(d < min_dist) {
	min_dist          = d;
	curve_set         = curve_sets[i];
	curve_set_idx     = i;
	curve_is_reversed = false;
	curve             = temp_curve;
	skin_arc          = temp_skin;
	skin_pos          = pos;
      }
    }
    // is it reverse oriented?
    if(front_endpt==temp_curve.back() && SENSE_FORWARD!=sense) {
      reverse( temp_curve.begin(), temp_curve.end() );
      // find the closest skin pt to the curve endpt
      unsigned pos;
      // if only one curve is left, take the entire skin. Otherwise, due to 
      // coordinate drift via merging the closest skin pt may not be the correct one.
      if(1==curve_sets.size()) {
        pos = skin_loop.size()-1;
      } else {
        rval = gen::find_closest_vert( temp_curve.back(), skin_loop, pos, temp_curve_len+extra );
        if(gen::error(MB_SUCCESS!=rval,"could not find_closest_vert")) return rval;
      }
      if(debug) std::cout << "  end of skin arc=" << skin_loop[pos] << std::endl; 
      // SPECIAL CASE: If the skin is a circle, create an arc out of the circle
      if(skin_loop[pos] == skin_loop.front()) pos = skin_loop.size()-1;

      // create a skin arc to test against the curve
      std::vector<MBEntityHandle> temp_skin(skin_loop.begin(), skin_loop.begin()+pos+1); 
      double d;
      rval = gen::dist_between_arcs( debug, temp_skin, temp_curve, d );
      if(gen::error(MB_SUCCESS!=rval,"could not get dist_between_arcs")) return rval;
      if(debug) std::cout << " curve-skin dist=" << d << std::endl;

      // if less than the min_dist, this curve is the best (thus far)
      if(d < min_dist) {
	min_dist          = d;
	curve_set         = curve_sets[i];
	curve_set_idx     = i;
	curve_is_reversed = true;
	curve             = temp_curve;
	skin_arc          = temp_skin;
	skin_pos          = pos;
      }
    }
  } // loop over all remaining curves of the surface

  // If no adjacent curves were found, something is wrong
  if(0==curve_set) {
    std::cout << "      no adjacent curve found" << std::endl;
    for(unsigned i=0; i<curve_sets.size(); ++i) {
      std::cout << "  curve " << gen::geom_id_by_handle(curve_sets[i]) 
                << " is unsealed" << std::endl;
    }
    return MB_FAILURE;
  }

  // If the average distance between the closest skin and curve is too far...
  if(100*FACET_TOL<=min_dist) {
    std::cout << "  warning8: curve too far from skin, average_dist=" 
              << min_dist << std::endl;
    if(1.0<=min_dist) return MB_FAILURE;
  }
     
  // remove the chosen curve the set of unsealed curves
  curve_sets.erase( curve_sets.begin()+curve_set_idx );

  // remove the chosen skin arc from the rest of the skin
  skin_loop.erase( skin_loop.begin(), skin_loop.begin()+skin_pos );  

  // If the entire skin loop has been sectioned into arcs, a single point remains.
  // This is because the skin_arc needs to have the same back point as the front
  // of the remaining skin_loop. If a single point remains, remove it.
  if(1==skin_loop.size()) skin_loop.clear();

  if(debug) std::cout << "  curve " << gen::geom_id_by_handle(curve_set) 
                      << " paired with skin, min_dist =" << min_dist << std::endl;
  return MB_SUCCESS;
}

//  -Instead of altering the skin and curve vectors, make them const. Put the
// sealed curve in a new vector, using only push_back. This
// would avoid all of the really slow inserts and erases in the curve and skin vectors.
MBErrorCode seal_arc_pair( const bool debug,
                           const double FACET_TOL,
                           const MBTag normal_tag,
                           std::vector<MBEntityHandle> &edge, /* in */
                           std::vector<MBEntityHandle> &skin /* in/out */,
                           const int surf_id ) {
  if(debug) {
    std::cout << "edge before sealing:" << std::endl;
    gen::print_loop(edge);
    std::cout << "skin before sealing:" << std::endl;
    gen::print_loop(skin);
   }

  MBErrorCode rval;
  const double TOL_SQR = FACET_TOL*FACET_TOL;
  if(gen::error(edge.empty() || skin.empty(),"edge or skin has no verts")) 
    return MB_FAILURE;
  //  sealed_edge.clear();
  //sealed_edge.reserve( 2*edge.size() );
  //std::vector<MBEntityHandle>::iterator e, s;
  //e = edge.begin(); 
  //s = skin.begin();
  //MBEntityHandle curr_pt = *e;
  //sealed_edge.push_back( curr_pt );
  //++e;

  //**************************************************************************
  // Merge the front of the skin to the front of the curve
  //**************************************************************************
  {
    MBEntityHandle keep_vert   = edge.front();
    MBEntityHandle delete_vert = skin.front();
    if(keep_vert != delete_vert) {
      double merge_dist;
      rval = gen::dist_between_verts( keep_vert, delete_vert, merge_dist );
      if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist g")) return rval;
      if(debug) {
	std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" << keep_vert
		  << " merge_dist=" << merge_dist << std::endl;
      }  
      if(FACET_TOL < merge_dist) {
	std::cout << "  warning0: front pt merge_dist=" << merge_dist << std::endl;
      }  
      rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
      if(gen::error(MB_SUCCESS!=rval,"could not merge verts g")) return rval;
    }
  }

  //**************************************************************************
  // Merge the back of the skin to the back of the curve
  //**************************************************************************
  {
    MBEntityHandle keep_vert   = edge.back();
    MBEntityHandle delete_vert = skin.back();
    if(keep_vert != delete_vert) {
      double merge_dist;
      rval = gen::dist_between_verts( keep_vert, delete_vert, merge_dist );
      if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist h")) return rval;
      if(debug) {
	std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" << keep_vert
		  << " merge_dist=" << merge_dist << std::endl;
      }  
      if(FACET_TOL < merge_dist) {
	std::cout << "  warning0: surf " << surf_id << " back pt merge_dist=" 
                  << merge_dist << std::endl;
        if(1000*FACET_TOL < merge_dist) return MB_FAILURE;
      }  
      rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
      if(gen::error(MB_SUCCESS!=rval,"could not merge verts g")) return rval;
    }
  }

  // ************************************************************************
  // zip the skin to the curve until the end of the curve is reached.
  // ************************************************************************
  // Zip the edge until we reach the end of the curve.
  unsigned int e_pos = 1, s_pos = 1; 
  bool edge_is_next;
  double dist, e_dist, s_dist, es_dist;
  while(true) {
 
    // Special Case: Through merging, a curve may have only a single vertex
    if(e_pos==edge.size() || s_pos==skin.size()) break;

    rval = gen::squared_dist_between_verts( edge[e_pos-1], edge[e_pos], e_dist );
    if(gen::error(MB_SUCCESS!=rval,"could not get e_dist")) return rval;
    rval = gen::squared_dist_between_verts( edge[e_pos-1], skin[s_pos], s_dist );
    if(gen::error(MB_SUCCESS!=rval,"could not get s_dist")) return rval;
    if(debug) {
      std::cout << " e_pos=" << e_pos << " e_dist=" 
		<< e_dist << " vert=" << edge[e_pos] << " size="
		<< edge.size() << std::endl;
      std::cout << " s_pos=" << s_pos << " s_dist=" 
		<< s_dist << " vert=" << skin[s_pos] << " size="
		<< skin.size() << std::endl;
    }

    // If skin is next, move the skin pt to the line extending from the next 
    // curve edge.
    if(e_dist > s_dist) {
      edge_is_next = false;
      double move_dist;
      rval = gen::line_point_dist( edge[e_pos-1], edge[e_pos], skin[s_pos], move_dist );
      if(gen::error(MB_SUCCESS!=rval,"could not get line_point_dist")) return rval;
      if(10*FACET_TOL < move_dist) {
	std::cout << "  warning5: surf " << surf_id << " vertex move_dist=" 
                  << move_dist << std::endl;
      }
      rval = gen::point_line_projection( edge[e_pos-1], edge[e_pos], skin[s_pos]);
      if(gen::error(MB_SUCCESS!=rval,"could not get point_line_projection")) return rval;
      rval = gen::squared_dist_between_verts( edge[e_pos-1], skin[s_pos], s_dist );
      if(gen::error(MB_SUCCESS!=rval,"could not get s_dist b")) return rval;
      dist = s_dist;
      if(debug) std::cout << "skin is next, projected dist=" << dist << std::endl;
    } else {
      edge_is_next = true;
      dist = e_dist;
      if(debug) std::cout << "edge is next" << std::endl;
    }

    // find the cs_dist after moving the skin to the curve (if skin_is_next)
    rval = gen::squared_dist_between_verts( edge[e_pos], skin[s_pos], es_dist );
    if(gen::error(MB_SUCCESS!=rval,"could not get es_dist")) return rval;
  
    // **************************************************************************
    // Merge with previous vert if it is too close
    // **************************************************************************
    if(dist < TOL_SQR) {
      MBEntityHandle keep_vert = edge[e_pos-1];
      if(edge_is_next) {
	MBEntityHandle delete_vert = edge[e_pos];
	if(keep_vert != delete_vert) { // cannot merge the same vert
	  if(debug) {
	    double merge_dist;
            rval = gen::dist_between_verts( keep_vert, delete_vert, merge_dist );
            if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist")) return rval;
	    std::cout << "  merged edge_vert=" << delete_vert << " to edge_vert=" 
		      << keep_vert << " merge_dist=" << merge_dist <<std::endl;
	  }  
	  rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
          if(gen::error(MB_SUCCESS!=rval,"could not merge_verts a")) return rval;
	}
	if(edge.size() < e_pos+1) {
	  std::cout << "edge.size()=" << edge.size() << " e_pos=" << e_pos << std::endl;
	}
	edge.erase( edge.begin() + e_pos );
      } else {
	MBEntityHandle delete_vert = skin[s_pos];
	if(keep_vert != delete_vert) {
	  if(debug) {
	    double merge_dist;
            rval  = gen::dist_between_verts( keep_vert, delete_vert, merge_dist );
            if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist b")) return rval;
	    std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" 
		      << keep_vert << " merge_dist=" << merge_dist << std::endl;  
	  } 
	  rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
          if(gen::error(MB_SUCCESS!=rval,"could not merge_verts b")) return rval;
	}
	if(skin.size() < s_pos+1) {
	  std::cout << "skin.size()=" << skin.size() << " s_pos=" << s_pos << std::endl;
	}
	skin.erase( skin.begin() + s_pos );      
      }
      // **************************************************************************
      // merge with next vert if it is too close
      // **************************************************************************
      // If the next hits are at the same distance or within merge tol, 
      //   advance the curve ahead. We need to merge all the verts between the
      //   current skin vert and the curve vert that we are merging to. This
      //   could be more than one! 
      //} else if(FACET_TOL > fabs(c_dist-s_dist)) {
    } else if(TOL_SQR > es_dist) {
      // merge the verts if they are not the same
      MBEntityHandle keep_vert  = edge[e_pos];
      if(skin[s_pos] != keep_vert) {
	MBEntityHandle delete_vert= skin[s_pos]; 
	if(debug) {
	  double merge_dist;
          rval  = gen::dist_between_verts( keep_vert, delete_vert, merge_dist );
          if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist c")) return rval;
	  std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" 
		    << keep_vert << " merge_dist=" << merge_dist << std::endl;
	}
	rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
        if(gen::error(MB_SUCCESS!=rval,"could not merge_verts b")) return rval;
      }
      s_pos++;
      e_pos++;
      // Are there more skin verts in proximity to the merged pt?
      while(true) {
	// check to see if skin still exists
	if(s_pos == skin.size()) break;
	// check if the distance is less than merge tol
	MBEntityHandle delete_vert= skin[s_pos]; 
	double merge_dist;
        rval = gen::dist_between_verts( keep_vert, delete_vert, merge_dist); 
        if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist d")) return rval;
	if(FACET_TOL < merge_dist) break;
	// merge the verts if they are not the same
	if(keep_vert != delete_vert) {
	  if(debug) {
	    std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" 
		      << keep_vert << " merge_dist=" << merge_dist << std::endl;
	  }
	  rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
          if(gen::error(MB_SUCCESS!=rval,"could not merge_verts d")) return rval;
	}
	skin.erase( skin.begin() + s_pos );      
      }
      // Are there more curve verst in proximity to the merged pt?
      while(true) {
	// check to see if curve still exists
	if(e_pos == edge.size()) break;
	// check if the distance is less than merge tol
	MBEntityHandle delete_vert= edge[e_pos]; 
	double merge_dist;
        rval = gen::dist_between_verts( keep_vert, delete_vert, merge_dist ); 
        if(gen::error(MB_SUCCESS!=rval,"could not get merge_dist e")) return rval;
	if(FACET_TOL < merge_dist) break;
	// merge the verts if they are not the same
	if(keep_vert != delete_vert) {
	  if(debug) {
	    std::cout << "  merged edge_vert=" << delete_vert << " to edge_vert=" 
		      << keep_vert << " merge_dist=" << merge_dist << std::endl;
	  }
	  rval = zip::merge_verts( keep_vert, delete_vert, skin, edge );
          if(gen::error(MB_SUCCESS!=rval,"could not merge_verts e")) return rval;
	}
	edge.erase( edge.begin() + e_pos );
      }
	  
      // **************************************************************************
      // Otherwise do a t_joint
      // **************************************************************************
    } else {
      if(edge_is_next) {
	if(debug) {
	  std::cout << "  zip phase: t_joint is inserting edge vert " 
		    << edge[e_pos] << " between edge vert " << edge[e_pos-1] 
		    << " and skin vert " << skin[s_pos] << std::endl;
	}
	double move_dist;
	rval = gen::line_point_dist( edge[e_pos-1], skin[s_pos], edge[e_pos], move_dist );
	if(gen::error(MB_SUCCESS!=rval,"could not get line_point_dist")) return rval;
	if(10*FACET_TOL < move_dist) {
	  std::cout << "  warning6: surf " << surf_id << " vertex move_dist=" 
                    << move_dist << std::endl;
	}
	rval = zip::t_joint( normal_tag, edge[e_pos-1], edge[e_pos], skin[s_pos] );
        if(gen::error(MB_SUCCESS!=rval,"tjoint failed a")) return rval;
	skin.insert( skin.begin()+s_pos, edge[e_pos] );
	e_pos++;
	s_pos++;
      } else { // skin is next: move the t to the curve
	if(debug) {
	  std::cout << "  zip phase: inserting projected skin vert " 
		    << skin[s_pos] << " between edge verts " 
		    << edge[e_pos-1] << " and " << edge[e_pos] << std::endl;
	}
        double move_dist;
	rval = gen::line_point_dist( edge[e_pos-1], edge[e_pos], skin[s_pos], move_dist );
	if(gen::error(MB_SUCCESS!=rval,"could not get line_point_dist")) return rval;
	if(10*FACET_TOL < move_dist) {
	  std::cout << "  warning6: surf " << surf_id << " vertex move_dist=" 
                    << move_dist << std::endl;
	}
	rval = zip::t_joint( normal_tag, edge[e_pos-1], skin[s_pos], edge[e_pos] );
        if(gen::error(MB_SUCCESS!=rval,"tjoint failed b")) return rval;
	edge.insert( edge.begin() + e_pos, skin[s_pos] );
	e_pos++;
	s_pos++;
      }
    }


    // exit if at the end of the curve
    if(e_pos==edge.size() || s_pos==skin.size()) break;

    // The smallest edge length should be no less than MERGE_TOL.
    if(2 <= e_pos) {
      double d;
      rval = gen::squared_dist_between_verts( edge[e_pos-1], edge[e_pos-2], d );
      if(gen::error(MB_SUCCESS!=rval,"could not get dist")) return rval;
      if(TOL_SQR > d) { 
	std::cout << "zip_loop: d=" << d << std::endl;
	gen::print_vertex_coords(edge[e_pos-1]);
	gen::print_vertex_coords(edge[e_pos-2]);
	std::cout << "warning7: surf " << surf_id 
                  << " adjacent edge points are closer than FACET_TOL" << std::endl;
      }
    }
    // The position should be the same. Do not exceed array bounds when checking.
    if(gen::error(e_pos!=s_pos,"skin and edge positions do not match")) return rval;
    if(edge[e_pos-1] != skin[s_pos-1]) {
      std::cout << "edge[" << e_pos-1 << "]=" << edge[e_pos-1]
		<< " skin[" << s_pos-1 << "]=" << skin[s_pos-1] << std::endl;
    }
    if(gen::error(edge[e_pos-1]!=skin[s_pos-1],"skin and edge vert does not match")) 
      return rval;
  }

  // The skin and curve should be the same size
  if(edge.size()!=skin.size()) {
    std::cout << "    surf " << surf_id 
              << " sealed skin and curve are not the same size" << std::endl;
    if(debug) {
      std::cout << "edge:" << std::endl;
      gen::print_loop(edge);
      std::cout << "skin:" << std::endl;
      gen::print_loop(skin);
    }
    //return MB_FAILURE;
  }

  if(debug) {
    std::vector< std::vector<MBEntityHandle> > temp;
    temp.push_back(edge);
    temp.push_back(skin);
    rval = zip::test_zipping(FACET_TOL, temp);
    if(gen::error(MB_SUCCESS!=rval,"sealing test failed")) return rval;
  }

  return MB_SUCCESS;

}

// 
MBErrorCode seal_loop( bool debug,
                       const double FACET_TOL,
                       const MBTag normal_tag,
                       const MBTag orig_curve_tag,
                       const MBEntityHandle surf_set,
                       std::vector<MBEntityHandle> &curve_sets,
                       std::vector<MBEntityHandle> &skin_loop ) {
  MBErrorCode rval;
  debug = false;
  // Find a curve that corresponds to the skin. Note that because there is more
  // then one skin loop, not all curves of the surface correspond to each skin
  // Loop.

  // To establish an inital connection, find the curve with endpoint closest
  // to a skin point.
  if(debug) std::cout << "seal_loop: new loop contains " << skin_loop.size() 
                      << " skin pts" << std::endl;
   
  // If skin remains but all the curves are zipped an error has occured.
  if( curve_sets.empty() ) {
    std::cout << "seal_loop: no curves are left, but skin remains" << std::endl;
    gen::print_loop(skin_loop);
    skin_loop.clear();
    return MB_FAILURE;
  }

  //**************************************************************************
  // choose the closest curve endpoint and align the loop with it
  //**************************************************************************
  /* Find the closest skin pt to the first curve's front endpt.
     MAJOR DESIGN PUSHER: The curves cannot be assembled into loops without
     knowing orientation (which we don't know). Instead this information is
     pulled from the skin loops. This is why we do not first create loops from
     curves. 

     If zipping fails, the failed loop and curve are deleted from the candidate
     set of things to be zipped. It is likely that curves of the failed zip
     still exist even though their corresponding loop has been removed. Not all
     curves in this list will ever be zipped. Select a curve that can be zipped. */
    unsigned pos, curve_idx;
    MBEntityHandle closest_skin_pt, closest_front_curve_endpt;
    double min_dist = std::numeric_limits<double>::max();
    for(unsigned i=0; i<curve_sets.size(); ++i) {
      // get geometric vertices
      MBRange endpt_sets;
      rval = MBI()->get_child_meshsets(curve_sets[i], endpt_sets );
      if(gen::error(MB_SUCCESS!=rval,"could not get endpt_sets")) return rval;
      if(gen::error(endpt_sets.empty() || 2<endpt_sets.size(),
                    "too many endpt_sets")) return MB_FAILURE;

      std::vector<MBEntityHandle> endpts;
      for(unsigned j=0; j<endpt_sets.size(); ++j) {
        MBRange endpt;
        rval = MBI()->get_entities_by_type( endpt_sets[j], MBVERTEX, endpt );
        if(gen::error(MB_SUCCESS!=rval,"could not get endpt")) return rval;
        if(gen::error(1!=endpt.size(),"not one endpt")) return MB_FAILURE;
        endpts.push_back( endpt.front() );
      }

      // check to ensure that the endpt sets aren't degenerate
      if(2==endpt_sets.size() && endpts.front()==endpts.back()) {
        //if(gen::error(endpts.front()==endpts.back(),
        //  "geometric endpts degenerate")) return MB_FAILURE;
	std::cout << "  warning9: curve " << gen::geom_id_by_handle(curve_sets[i]) 
                  << " geometric endpoints degenerate" << std::endl;
      }
      
      // check to ensure that geometric verts are the curve endpts
      std::vector<MBEntityHandle> curve;
      rval = arc::get_meshset( curve_sets[i], curve );
      if(gen::error(MB_SUCCESS!=rval,"could not get_meshset")) return rval;
      if(1==endpt_sets.size()) {
        if(gen::error(curve.front()!=curve.back(),"endpt discrepancy")) return MB_FAILURE;
        if(gen::error(curve.front()!=endpts.front(),
          "geometric verts inconsistent with curve")) return MB_FAILURE;
      } else {
	//        if(gen::error(curve.front()==curve.back(),"degenerate curve endpts")) return MB_FAILURE;
        if(curve.front()==curve.back()) 
          std::cout << "  warning10: degenerate curve endpts" << std::endl;
        if(gen::error(curve.front()!=endpts.front() && curve.front()!=endpts.back(),
		      "endpts not consistent")) return MB_FAILURE;
        if(gen::error(curve.back()!=endpts.front() && curve.back()!=endpts.back(),
		      "endpts not consistent")) return MB_FAILURE;
      }

      // determine the orientation of the curve wrt the surf.
      int sense;
      if(debug)
      {
      std::cout << "surf_set = " << gen::geom_id_by_handle(surf_set) << std::endl;
      std::cout << "curve_set = " << gen::geom_id_by_handle(curve_sets[i]) << std::endl;
      }
      rval = gen::get_curve_surf_sense( surf_set, curve_sets[i], sense );  
      if(gen::error(MB_SUCCESS!=rval,"could not get_curve_surf_sense")) return rval;
      // select the front wrt the skin.
      MBEntityHandle curve_endpt = (SENSE_FORWARD==sense) ? curve.front() : curve.back();

      // find closest skin vert to front of curve
      std::vector<double> d;
      std::vector<unsigned> p;
      rval = gen::find_closest_vert( 0, curve_endpt, skin_loop, p, d);
      if(gen::error(MB_SUCCESS!=rval,"could not find_closest_vert")) return rval;
      if(debug) std::cout << "zip_loop: loop-curve endpt dist=" << d.front() << " skin_vert="
                          << skin_loop[p.front()] << " curve="
                          << gen::geom_id_by_handle(curve_sets[i]) << " front_endpt="
                          << curve_endpt << std::endl;
      if(d.front() < min_dist) {
        min_dist = d.front();
        curve_idx = i;
        pos = p.front();
        // select_next_curve uses the front_endpt to choose the next curve
        closest_front_curve_endpt = curve_endpt;
        closest_skin_pt = skin_loop[p.front()];
      }
    }

    MBEntityHandle front_endpt = closest_front_curve_endpt;

    // The closest points should be within facet tolerance
    if(100*FACET_TOL<min_dist) {
      std::cout << "closest skin pt farther than 100*FACET_TOL from curve vert (" 
                << min_dist << ")" << std::endl; 
      if(true) {
	std::cout << "  skin pt:" << std::endl;
        rval = MBI()->list_entity(closest_skin_pt);
        if(gen::error(MB_SUCCESS!=rval,"error listing skin pt")) return rval;
	std::cout << "  curve vert:" << std::endl;
        rval = MBI()->list_entity(front_endpt);
        if(gen::error(MB_SUCCESS!=rval,"error listing curve_vert")) return rval;
      }
      return rval;
    }

    if(debug) {
      std::cout << "closest skin vert=" << skin_loop[pos] << " pos=" << pos 
                << " min_dist=" << min_dist << " curve=" 
                << gen::geom_id_by_handle(curve_sets[curve_idx]) 
                << " front_endpt=" << closest_front_curve_endpt << std::endl;
    }

    /* Take the skin loops and advance it such that our common point is at 
       location 0. Ensure that that both endpoints are the same. */
    std::vector<MBEntityHandle> temp_loop;
    temp_loop.reserve(skin_loop.size());
    std::vector<MBEntityHandle>::iterator j = temp_loop.begin();
    std::vector<MBEntityHandle>::iterator k = skin_loop.begin();
    // Do not insert the back endpoint (avoid duplicate).   
    temp_loop.insert( j, k+pos, k+skin_loop.size()-1 );
    j = temp_loop.begin(); // j became invalid because temp_loop resized
    j += skin_loop.size() - pos - 1;
    temp_loop.insert( j, k, k+pos+1 );
    if(gen::error(temp_loop.size()!=skin_loop.size(),"loop size not conserved")) return MB_FAILURE;
    assert(temp_loop.size() == skin_loop.size()); // same size
    if(gen::error(temp_loop[0]!=temp_loop[temp_loop.size()-1],
      "loop endpts not continuous")) return MB_FAILURE;
    assert(temp_loop[0] == temp_loop[temp_loop.size()-1]); // same endpoint
    skin_loop = temp_loop;
    if(debug) {
      std::cout << "  skin:" << std::endl;
      //gen::print_loop(skin_loop);
    }

    // Create a set to store skin_loop so that handles are updated during merging.
    // This prevents stale handles, and avoids O(n) search through every handle in
    // the skin loop if manually updating the skin_loop vector.
    MBEntityHandle skin_loop_set;
    rval = MBI()->create_meshset( MESHSET_TRACK_OWNER|MESHSET_ORDERED, skin_loop_set );
    if(gen::error(MB_SUCCESS!=rval,"creating skin_loop_set failed")) return rval;      

    while(!skin_loop.empty()) {
      //**************************************************************************
      // Select the curve adjacent to the common point that best matches the skin 
      // using average distance. Chop the corresponding arc of skin, to form a
      // curve-skin_arc pair.
      //**************************************************************************

      bool curve_is_reversed;
      MBEntityHandle curve_set;
      std::vector<MBEntityHandle> curve, skin_arc;
      rval = create_arc_pair( FACET_TOL, surf_set, skin_loop, curve_sets, front_endpt,
			      debug, curve_set, curve_is_reversed, curve, skin_arc );
      if(gen::error(MB_SUCCESS!=rval,"  pair creation failed")) return rval;

      // Let moab store skin loop to avoid stale vert handles from merging.
      rval = arc::set_meshset( skin_loop_set, skin_loop );    
      if(gen::error(MB_SUCCESS!=rval,"setting skin_loop_set failed")) return rval;      

      // The original faceted curves are never used. Instead they are replaced by
      // skin. This reduces the number of new triangles created.
      int orig_curve;
      rval = MBI()->tag_get_data( orig_curve_tag, &curve_set, 1, &orig_curve );
      if(gen::error(MB_SUCCESS!=rval,"can't get tag")) return rval;

      // If the tag is non-zero, the facet edge has already been replaced.
      if(orig_curve) {
        // this tag is used to mark the edge as updated
        int false_int = 0;
        rval = MBI()->tag_set_data( orig_curve_tag, &curve_set, 1, &false_int );
        if(gen::error(MB_SUCCESS!=rval,"can't set tag")) return rval;

        // merge new endpoints to old endpoints  
        if(curve.front()!=skin_arc.front()) {
          rval = zip::merge_verts( curve.front(), skin_arc.front(), curve, skin_arc );
	  if(gen::error(MB_SUCCESS!=rval,"merge verts failed")) return rval;
        }
        if(curve.back()!=skin_arc.back()) {
          rval = zip::merge_verts( curve.back(), skin_arc.back(), curve, skin_arc );
	  if(gen::error(MB_SUCCESS!=rval,"merge verts failed")) return rval;
        }

	// replace the faceted edge with a skin arc
        curve = skin_arc;
        if(debug) std::cout << "  curve " << gen::geom_id_by_handle(curve_set) 
                            << " has been replaced by skin:" << std::endl;
        //if(debug) gen::print_loop(curve);

      } else {
	// seal the pair together
	std::vector<MBEntityHandle> sealed_curve;
	rval = seal_arc_pair( debug, FACET_TOL, normal_tag, curve, skin_arc,
                              gen::geom_id_by_handle(surf_set) );
        if(gen::error(MB_SUCCESS!=rval,"    can't seal pair")) return rval;
      }

      // get new front_endpt to guide selection of next curve
      front_endpt = curve.back();
     
      // To preserve the original rotation, reverse the curve if it has been reversed
      if(curve_is_reversed) reverse( curve.begin(), curve.end() );        

      // set the sealed edge
      rval = arc::set_meshset( curve_set, curve );    
      if(gen::error(MB_SUCCESS!=rval,"setting curve set failed")) return rval;

      // Get skin_loop to cut an arc from it
      rval = arc::get_meshset( skin_loop_set, skin_loop );    
      if(gen::error(MB_SUCCESS!=rval,"getting skin_loop set failed")) return rval;      

      
    }

    // The skin_loop_set is no longer needed.
    rval = MBI()->delete_entities( &skin_loop_set, 1 );
    if(gen::error(MB_SUCCESS!=rval,"deleting skin_loop_set failed")) return rval;      

    return MB_SUCCESS;
  }

  // input: surface sets, ordered curve sets,
  // output: skin arcs corresponding to curves are added to parent surface sets
MBErrorCode prepare_surfaces(MBRange &surface_sets,
                             MBTag geom_tag, MBTag id_tag, MBTag normal_tag, MBTag merge_tag,
                             MBTag orig_curve_tag,
                             const double SME_RESABS_TOL, const double FACET_TOL, 
                               const bool debug, bool verbose) 
{
   
    MBErrorCode result;
   
    // loop over each surface meshset
    for(MBRange::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) 
      {

	// get the surf id of the surface meshset
	int surf_id;
	result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );

	if(gen::error(MB_SUCCESS!=result,"could not get id tag"))
	  {
	    return result;
	  }

	assert(MB_SUCCESS == result);

	if(debug) std::cout << "  surf id= " << surf_id << std::endl;

	// get the 2D entities in the surface set
	MBRange dim2_ents;
	result = MBI()->get_entities_by_dimension( *i, 2, dim2_ents );
	if(gen::error(MB_SUCCESS!=result,"could not get 3D entities")) 
	  {
	    return result;
	  }
	assert(MB_SUCCESS == result);

	// get facets of the surface meshset
	MBRange tris;
	result = MBI()->get_entities_by_type( *i, MBTRI, tris );

	if(gen::error(MB_SUCCESS!=result,"could not get tris")) return result;

	assert(MB_SUCCESS == result);

      // Remove and 2D entities that are not triangles. This is needed because
      // ReadCGM will add quads and polygons to the surface set. This code only
      // works with triangles.
      MBRange not_tris = subtract( dim2_ents, tris );
      if(!not_tris.empty()) {
        result = MBI()->delete_entities( not_tris );
        if(gen::error(MB_SUCCESS!=result,"could not delete not_tris")) return result;
        assert(MB_SUCCESS == result);
        std::cout << "  removed " << not_tris.size() 
                  << " 2D elements that were not triangles from surface " 
                  << surf_id << std::endl;
      }

      // Get the curves sets
      std::vector<MBEntityHandle> curve_sets, unmerged_curve_sets;
      result = MBI()->get_child_meshsets( *i, curve_sets );
      if(gen::error(MB_SUCCESS!=result,"could not get child sets")) return result;
      assert(MB_SUCCESS==result);

      // Update the curve_sets with that contain entity_to_delete curves with their
      // entity_to_keep curves. Unmerged_curve_sets will end up holding the curves
      // of this surface that are not merged with another curve in this surface.
      for(std::vector<MBEntityHandle>::iterator j=curve_sets.begin();
	  j!=curve_sets.end(); j++) {
        MBEntityHandle merged_curve, curve;
        result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );
        assert(MB_TAG_NOT_FOUND==result || MB_SUCCESS==result);     
        if(MB_TAG_NOT_FOUND==result) {
          curve = *j;
        } else if(MB_SUCCESS == result) {
	  if(1) {
            std::cout << "  curve_id=" << gen::geom_id_by_handle(*j) 
                      << " is entity_to_delete" << std::endl;
          }
          curve = merged_curve;
          // should parent-childs be updated for the entity_to_keep?
        } else {
	  std::cout << "prepare_surfaces: result=" << result << std::endl;
          return result;        
        }
      
        // If the curve is in unmerged_curve_sets, then remove it. Otherwise add it.
	std::vector<MBEntityHandle>::iterator k=find(unmerged_curve_sets.begin(), 
	  unmerged_curve_sets.end(), curve);
        if(unmerged_curve_sets.end() == k) {
  	  //std::cout << "  curve " << gen::geom_id_by_handle(*k) 
          //          << " is entity_to_keep" << std::endl;
          unmerged_curve_sets.push_back(curve);
        } else {
          unmerged_curve_sets.erase(k);
        }
      }

      // If all of the curves are merged, remove the surfaces facets.
      // THIS ALSO REMOVES SURFACES THAT HAVE ALL CURVES DELETED?
      if(unmerged_curve_sets.empty()) {
        //measure area of the surface  
        double area;
        result = gen::measure( *i, geom_tag, area, verbose );
        if(gen::error(MB_SUCCESS!=result,"could not measure area")) return result;
        assert(MB_SUCCESS == result);

        //remove triagngles from the surface
        result = MBI()->remove_entities( *i, tris);                          
        if(gen::error(MB_SUCCESS!=result,"could not remove tris")) return result;
	assert(MB_SUCCESS == result);                                                        

	std::cout << "  deleted surface " << surf_id
                  << ", area=" << area << " cm^2, n_facets=" << tris.size() << std::endl;

        // delete triangles from mesh data
	result = MBI()->delete_entities( tris );                               
        if(gen::error(MB_SUCCESS!=result,"could not delete tris")) return result;
	assert(MB_SUCCESS == result);

        //remove the sense data for the surface from its child curves
        result = remove_surf_sense_data(*i);
        if(gen::error(MB_SUCCESS!=result, "could not remove surface's sense data")) return result;

        // remove the surface set itself
        result = MBI()->delete_entities( &(*i), 1);
        if(gen::error(MB_SUCCESS!=result,"could not delete surface set")) return result;
        assert(MB_SUCCESS == result);
        i = surface_sets.erase(i) - 1;
        continue;
      }


      // Save the normals of the facets. These will later be used to determine if
      // the tri became inverted.
      result = gen::save_normals( tris, normal_tag );
      if(gen::error(MB_SUCCESS!=result,"could not save_normals")) return result;
      assert(MB_SUCCESS == result);
  
      // do edges already exist?
      int n_edges;
      result = MBI()->get_number_entities_by_type(0, MBEDGE, n_edges );
      if(gen::error(MB_SUCCESS!=result,"could not get number of edges")) return result;
      assert(MB_SUCCESS == result);
      if(gen::error(0!=n_edges,"edges exist")) return result;
      assert(0 == n_edges); //*** Why can't we have edges? (Also, this assertion is never used)

      // get the range of skin edges from the range of facets
      MBSkinner tool(MBI());
      MBRange skin_edges, skin_edges2;
      if(tris.empty()) continue; // nothing to zip
      result = tool.find_skin( tris, 1, skin_edges, false );
      if(gen::error(MB_SUCCESS!=result,"could not find_skin")) return result;
      assert(MB_SUCCESS == result);
      //std::cout << "my skinner size=" << skin_edges2.size() << std::endl;
      //gen::print_range_of_edges(skin_edges2);
      //result = MBI()->delete_entities( skin_edges2 );
      //assert(MB_SUCCESS == result);
      //result = tool.find_skin( tris, 1, skin_edges, false );
      //MBRange temp_verts;
      //result = tool.find_skin_vertices( tris, temp_verts, &skin_edges, true );
      //assert(MB_SUCCESS == result);
      //std::cout << "moab skinner size=" << skin_edges.size() << std::endl;
      //gen::print_range_of_edges(skin_edges);
      //assert(skin_edges.size() == skin_edges2.size());

      //find_degenerate_tris();
      //if(881 == surf_id) gen::print_triangles( tris);

      // merge the vertices of the skin
      // BRANDON: For some reason cgm2moab does not do this? This was the 
      // problem with mod13 surf 881. Two skin verts were coincident. A tol=1e-10
      // found the verts, but tol=0 did not.
      MBRange skin_verts;
      result = MBI()->get_adjacencies( skin_edges, 0, false, skin_verts, 
                                       MBInterface::UNION );
      if(gen::error(MB_SUCCESS!=result,"could not get adj verts")) return result;
      assert(MB_SUCCESS == result);
      result = gen::merge_vertices( skin_verts, SME_RESABS_TOL );         
      if (MB_SUCCESS != result) {                                             
	if(debug) std::cout << "result= " << result << std::endl;                  
	std::cout << "  surface " << surf_id << " failed to zip: could not merge vertices"   
		  << surf_id << std::endl;  
        result = MBI()->delete_entities(skin_edges);
        assert(MB_SUCCESS == result);
	continue;
      }  

      // Merging vertices create degenerate edges.
      result = arc::remove_degenerate_edges( skin_edges, debug );
      if(MB_SUCCESS!=result) {
	std::cout << "  surface " << surf_id 
                  << " failed to zip: could not remove degenerate edges" << std::endl;
        result = MBI()->delete_entities(skin_edges);
        if(gen::error(MB_SUCCESS!=result,"could not delete skin edges")) return result;
        assert(MB_SUCCESS == result);
	continue;
      }  

 
      /* Remove pairs of edges that are geometrically the same but in opposite 
        order due to a faceting error ACIS. In other words, merge (a,b) and (b,a). 
        Examples are mod13surf280 and tbm_surf1605. */
      result = arc::remove_opposite_pairs_of_edges_fast( skin_edges, debug );
      if (MB_SUCCESS != result) {                                             
	std::cout << "  surface " << surf_id << " failed to zip: could not remove opposite edges"   
		  << surf_id << std::endl;  
        result = MBI()->delete_entities(skin_edges);
        if(gen::error(MB_SUCCESS!=result,"could not delete skin edges")) return result;
        assert(MB_SUCCESS == result);
	continue;
      }  

      /* Order the edges so that the triangle is on their left side. Skin
	 edges are adjacent to only one triangle. */
      bool catch_error = false;
      for(MBRange::iterator j=skin_edges.begin(); j!=skin_edges.end(); j++) {
	MBRange adj_tris;
	result = MBI()->get_adjacencies( &(*j), 1, 2, false, adj_tris );
        if(gen::error(MB_SUCCESS!=result,"could not get adj tris")) return result;
	assert(MB_SUCCESS == result);
	MBRange skin_tri = intersect( adj_tris, tris );
        if(1 != skin_tri.size()) {
          std::cout << "skin_tri.size()=" << skin_tri.size() << std::endl;
          catch_error = true;
          break;
        }
	result = arc::orient_edge_with_tri( *j, skin_tri.front() );
        if(gen::error(MB_SUCCESS!=result,"could not orient_edge_with_tri")) return result;
	assert(MB_SUCCESS == result);
      }
      // I NEED TO ADD BETTER CLEANUP AFTER THESE FAILURE CONDITIONS
      if(catch_error) {
	std::cout << "  surface " << surf_id << " failed to zip: could not orient edge"   
                  << std::endl;  
        result = MBI()->delete_entities(skin_edges);
        if(gen::error(MB_SUCCESS!=result,"could not delete skin edges")) return result;
        assert(MB_SUCCESS == result);
        continue;
      }
						      
      // Create loops with the skin edges.  
      std::vector< std::vector<MBEntityHandle> > skin_loops_of_edges;
      result = arc::create_loops_from_oriented_edges( skin_edges, skin_loops_of_edges, debug );
      if(MB_SUCCESS != result) {
	std::cout << "  surface " << surf_id << " failed to zip: could not create loops" 
		  << std::endl;
        result = MBI()->delete_entities(skin_edges);
        assert(MB_SUCCESS == result);
	continue;
      }
      if(debug) std::cout << skin_loops_of_edges.size() << " skin loop(s)" << std::endl;

      // Convert the loops of skin edges to loops of skin verts.
      std::vector< std::vector<MBEntityHandle> > skin(skin_loops_of_edges.size());
      for(unsigned int j=0; j<skin_loops_of_edges.size(); j++) {
        result = gen::ordered_verts_from_ordered_edges( skin_loops_of_edges[j], skin[j] );
        assert(MB_SUCCESS == result);
	// check to make sure that the loop is closed
	assert(skin[j].front() == skin[j].back());
      }

      // edges are no longer needed       
      result = MBI()->delete_entities(skin_edges);
      if(gen::error(MB_SUCCESS!=result,"could not delete skin_edges")) return result;
      assert(MB_SUCCESS == result);

      /* Get the curves that are part of the surface. Use vectors to store all curve
	 stuff so that we can remove curves from the set as they are zipped. */
      //std::vector<int> curve_ids;
      //int curve_id;
      //std::vector<std::vector<MBEntityHandle> > curves;
      for(std::vector<MBEntityHandle>::iterator j=curve_sets.begin(); 
        j!=curve_sets.end(); j++) {

        // If a delete_curve, replace it with the keep_curve. This approach allows
        // for duplicates because we are using vectors instead of ranges. Note that
        // parent-child links also cannot store duplicate handles.
        MBEntityHandle merged_curve;
        result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );     
        assert(MB_TAG_NOT_FOUND==result || MB_SUCCESS==result);
        if(MB_SUCCESS == result) {
          // add senses from curve_to_delete to curve_to_keep
	  std::vector<MBEntityHandle> curve_to_keep_surfs, curve_to_delete_surfs, combined_surfs;
	  std::vector<int> curve_to_keep_senses, curve_to_delete_senses, combined_senses;
	  moab::GeomTopoTool gt(MBI(), false);
          result = gt.get_senses( *j, curve_to_delete_surfs, curve_to_delete_senses );
          if(gen::error(MB_SUCCESS!=result,"failed to get senses 0")) return result;
          result = gt.get_senses( merged_curve, curve_to_keep_surfs, curve_to_keep_senses );
          if(gen::error(MB_SUCCESS!=result,"failed to get senses 1")) return result;
          std::cout << "curve to keep id = " << gen::geom_id_by_handle(merged_curve) << std::endl;
          std::cout << "curve to delete id = " << gen::geom_id_by_handle(*j) << std::endl;
          for(unsigned int index=0; index < curve_to_keep_surfs.size(); index++)
          {
           std::cout << "curve_to_keep_surf " << index << " id = " << gen::geom_id_by_handle(curve_to_keep_surfs[index]) << std::endl;
           std::cout << "curve_to_keep_sense " << index << " = " << curve_to_keep_senses[index] << std::endl;
           
          }
          for(unsigned int index=0; index < curve_to_keep_surfs.size(); index++)
          {
          
          std::cout << "curve_to_delete_surf " << index << " id = " << gen::geom_id_by_handle(curve_to_delete_surfs[index]) << std::endl;
          std::cout << "curve_to_delete_sense " << index << " = " << curve_to_delete_senses[index] << std::endl;
          }

          combined_surfs.insert( combined_surfs.end(), curve_to_keep_surfs.begin(), 
                                                       curve_to_keep_surfs.end() );
          combined_surfs.insert( combined_surfs.end(), curve_to_delete_surfs.begin(), 
                                                       curve_to_delete_surfs.end() );
          combined_senses.insert(combined_senses.end(),curve_to_keep_senses.begin(),
		                		       curve_to_keep_senses.end() );
          combined_senses.insert(combined_senses.end(),curve_to_delete_senses.begin(),
		                		       curve_to_delete_senses.end() );
          std::cout << combined_surfs.size() << std::endl;
          std::cout << combined_senses.size() << std::endl;
          int edim;
           
          //don't send unknown senses to set_senses
          /*
          for(unsigned int index=0; index < combined_senses.size() ; index++)
          {
          if (combined_senses[index]==0)
            {
             combined_surfs.erase(combined_surfs.begin() + index);
  
             combined_senses.erase(combined_senses.begin() + index);
             index = -1;
            }
          }
          */
          for(unsigned int index=0; index < combined_senses.size(); index++)
          {

           std::cout << "combined_surfs{"<< index << "] = " << gen::geom_id_by_handle(combined_surfs[index]) << std::endl;
           std::cout << "combined_sense["<< index << "] = " << combined_senses[index] << std::endl;
          }
          result = gt.set_senses( merged_curve, combined_surfs, combined_senses );
          
          if(gen::error(MB_SUCCESS!=result && MB_MULTIPLE_ENTITIES_FOUND!=result,"failed to set senses: "))
          {
           moab_printer(result);
           return result;
          }
          
          // add the duplicate curve_to_keep to the vector of curves
          *j = merged_curve;
          
          
        }
	//	result = MBI()->tag_get_data( id_tag, &(*j), 1, &curve_id );
	//assert(MB_SUCCESS == result);
	//std::cout << "  curve_id=" << curve_id << " handle=" << *j << std::endl;
	//curve_ids.push_back(curve_id);
	//std::vector<MBEntityHandle> curve;
	//result = arc::get_meshset( *j, curve );
	//assert(MB_SUCCESS == result);
	//curves.push_back( curve );
      }

      // Place skin loops in meshsets to allow moab to update merged entities.
      // This prevents stale handles. I suspect moab can use upward adjacencies
      // to do this more efficiently than a manual O(n) search through an 
      // unsorted vector.
      MBErrorCode rval;
      MBEntityHandle skin_loop_sets[skin.size()];
      for(unsigned j=0; j<skin.size(); ++j) {
        rval = MBI()->create_meshset( MESHSET_TRACK_OWNER|MESHSET_ORDERED, skin_loop_sets[j] );
        if(gen::error(MB_SUCCESS!=rval,"failed to zip: creating skin_loop_set failed"))
          return rval;   

        rval = arc::set_meshset( skin_loop_sets[j], skin[j] );    
        if(gen::error(MB_SUCCESS!=rval,"failed ot zip: setting skin_loop_set failed")) 
          return rval;      
      }

      // Keep zipping loops until each is either zipped or failed. This function
      // returns only after all loops are zipped or a failure occurs.
      for(unsigned j=0; j<skin.size(); ++j) {
	std::vector<MBEntityHandle> skin_loop;
        rval = arc::get_meshset( skin_loop_sets[j], skin_loop );    
        if(gen::error(MB_SUCCESS!=rval,"failed to zip: setting skin_loop_set failed")) 
          return rval;      

        rval = seal_loop( debug, FACET_TOL, normal_tag, orig_curve_tag, *i, 
                          curve_sets, skin_loop );
        if(MB_SUCCESS != rval) {
	  std::cout << "failed to zip: surface " << surf_id << ": failed to seal a loop" 
                    << std::endl;
          //assert(false);
        }
      }

      // Remove the sets of skin loops
      rval = MBI()->delete_entities( &skin_loop_sets[0], skin.size() );
      if(gen::error(MB_SUCCESS!=rval,"failed to zip: deleting skin_loop_sets failed")) 
        return rval;

      // mod13surf2996, 3028 and 2997 are adjacent to the same bad geometry (figure 8 loop)
      //assert(MB_SUCCESS==result || 2996==surf_id || 2997==surf_id || 3028==surf_id);
    } // loop over each surface
    return MB_SUCCESS;
  }

// removes sense data from all curves associated with the surface given to the function

MBErrorCode remove_surf_sense_data(MBEntityHandle del_surf) {
 
  MBErrorCode result;
  moab::GeomTopoTool gt(MBI(), false);
    int edim = gt.dimension(del_surf);

    if(gen::error(edim!=2,"could not remove sense data: entity is of the wrong dimension")) return MB_FAILURE;

 // get the curves of the surface
        MBRange del_surf_curves;
        result = MBI() -> get_child_meshsets( del_surf, del_surf_curves);
        if(gen::error(MB_SUCCESS!=result,"could not get the curves of the surface to delete")) return result;
        std::cout << "got the curves" << std::endl;
       
  

        std::cout << "number of curves to the deleted surface = " << del_surf_curves.size() << std::endl;
        for(unsigned int index =0 ; index < del_surf_curves.size() ; index++)
        {
          std::cout << "deleted surface's child curve id " << index << " = " << gen::geom_id_by_handle(del_surf_curves[index]) << std::endl;
        }        
  
        //get the sense data for each curve

        //get sense_tag handles from MOAB
        MBTag senseEnts, senseSenses;
        unsigned flags = MB_TAG_SPARSE;
     
        //get tag for the entities with sense data associated with a given moab entity
        result = MBI()-> tag_get_handle(GEOM_SENSE_N_ENTS_TAG_NAME, 0, MB_TYPE_HANDLE, senseEnts, flags);
        if(gen::error(MB_SUCCESS!=result, "could not get senseEnts tag")) return result;
        
        //get tag for the sense data associated with the senseEnts entities for a given moab entity
        result = MBI()-> tag_get_handle(GEOM_SENSE_N_SENSES_TAG_NAME, 0, MB_TYPE_INTEGER, senseSenses, flags);
        if(gen::error(MB_SUCCESS!=result,"could not get senseSenses tag")) return result;
        
        //initialize vectors for entities and sense data
        std::vector<MBEntityHandle> surfaces;
        std::vector<int> senses;
        for(MBRange::iterator i=del_surf_curves.begin(); i!=del_surf_curves.end(); i++ ) 
        {
       
        result = gt.get_senses(*i, surfaces, senses);
        if(gen::error(MB_SUCCESS!=result, "could not get the senses for the del_surf_curve")) return result;
          // if the surface to be deleted (del_surf) exists in the sense data (which it should), then remove it
          for(unsigned int index = 0; index < senses.size() ; index++)
          {
            if(surfaces[index]==del_surf)
            {
             surfaces.erase(surfaces.begin() + index);
             senses.erase(senses.begin() +index);
             index=-1;
            }
          }
          //remove existing sense entity data for the curve
          result= MBI()-> tag_delete_data( senseEnts, &*i, 1);
          if(gen::error(MB_SUCCESS!=result, "could not delete sense entity data")) return result;
       
          //remove existing sense data for the curve
          result = MBI()-> tag_delete_data(senseSenses, &*i, 1);
          if(gen::error(MB_SUCCESS!=result, "could not delete sense data")) return result;

          //reset the sense data for each curve 
          result = gt.set_senses( *i, surfaces, senses);
          if(gen::error(MB_SUCCESS!=result, "could not update sense data for surface deletion")) return result;

        }

return MB_SUCCESS;
} 

MBErrorCode fix_normals(MBRange surface_sets, MBTag id_tag, MBTag normal_tag, const bool debug, const bool verbose) {
    MBErrorCode result;
    if(debug) std::cout<< "number of surfaces=" << surface_sets.size() << std::endl;
    int inverted_tri_counter = 0;

    // loop over each surface meshset
    for(MBRange::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) {

      // get the surf id of the surface meshset
      int surf_id;
      result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );
      assert(MB_SUCCESS == result);
      if(debug) std::cout << "fix_normals surf id=" << surf_id << std::endl;

      // get facets from the surface meshset
      MBRange tris;
      result = MBI()->get_entities_by_type( *i, MBTRI, tris );
      assert(MB_SUCCESS == result);

      // get the normals, pre zipping
      std::vector<MBCartVect> old_normals(tris.size()), new_normals(tris.size());
      result = MBI()->tag_get_data( normal_tag, tris, &old_normals[0]);
      assert(MB_SUCCESS == result);

      // get the normals, post zipping
      result = gen::triangle_normals( tris, new_normals );
      assert(MB_SUCCESS == result);

      // test the normals, finding the inverted tris
      std::vector<int> inverted_tri_indices;
      result = zip::test_normals( old_normals, new_normals, inverted_tri_indices);
      assert(MB_SUCCESS == result);

      // insert the inverted tris into a range
      MBRange inverted_tris;
      for(unsigned int j=0; j<inverted_tri_indices.size(); j++) {
        inverted_tris.insert( tris[inverted_tri_indices[j]] );
        if(debug) gen::print_triangle( tris[inverted_tri_indices[j]], false );
      }

      // do edges exist?
      int n_edges;
      result = MBI()->get_number_entities_by_type( 0, MBEDGE, n_edges );
      assert(MB_SUCCESS == result);
      assert(0 == n_edges); // *** Why can't we have edges?
  
      // fix the inverted tris
      inverted_tri_counter += inverted_tris.size();
      result = zip::remove_inverted_tris(normal_tag, inverted_tris, debug );
      if(MB_SUCCESS != result) 
        std::cout << "  failed to fix inverted triangles in surface " << surf_id << std::endl;

      // if fix_normals exits on an error, we still need to remove its edges
      result = delete_all_edges();
      assert(MB_SUCCESS == result);
    }
    if(verbose)
    {
    std::cout << "  Before fixing, " << inverted_tri_counter 
              << " inverted triangles were found." << std::endl;
    }
    return MB_SUCCESS;
  }

  MBErrorCode restore_moab_curve_representation( const MBRange curve_sets ) {
    MBErrorCode result;
    for(MBRange::const_iterator i=curve_sets.begin(); i!=curve_sets.end(); ++i) {
      // get the ordered verts
      std::vector<MBEntityHandle> ordered_verts;
      result = arc::get_meshset( *i, ordered_verts );
      assert(MB_SUCCESS==result);
      if(MB_SUCCESS != result) return result;

      // Check for duplicate verts. This should not happen, but could if line
      // surfaces exist. This happens when feature size is violated and skin
      // from curves on the other side of the 1D line surf gets sealed. 
      if( 1<ordered_verts.size() ) {
        for(std::vector<MBEntityHandle>::iterator j=ordered_verts.begin()+1;
	    j!=ordered_verts.end(); ++j) {
          if( *j == *(j-1) ) {
	    std::cout << "duplicate vertex found in curve " 
                      << gen::geom_id_by_handle(*i) << std::endl;
            j = ordered_verts.erase(j) - 1;
          }
        }
      }

      // Check for a point curve (should never happen),
      // a single degenerate edge (should never happen), or
      // a degenerate loop of two edges (should never happen)
      if( 4>ordered_verts.size() && 
          ordered_verts.front()==ordered_verts.back() ) {
	std::cout << "warning: curve " << gen::geom_id_by_handle(*i) 
                  << " is one degenerate edge" << std::endl;
        //return MB_FAILURE;
      }

      // Determine if the curve is a loop or point curve. At least 4 verts are
      // needed to form a loop.
      bool is_loop = ( ordered_verts.front()==ordered_verts.back() &&
                       3<ordered_verts.size() );

      // get geometric endpoint sets of the curve (may be stale)
      MBRange endpt_sets;
      result = MBI()->get_child_meshsets( *i, endpt_sets );
      assert(MB_SUCCESS==result);
      if(MB_SUCCESS != result) return result;

      // do the correct number of endpt sets exist?
      const unsigned int n_endpts = (is_loop) ? 1 : 2;
      if(n_endpts != endpt_sets.size()) {
	std::cout << "curve " << gen::geom_id_by_handle(*i) << " has " << n_endpts 
                  << " endpoints, but " << endpt_sets.size()
                  << " endpoint sets exist" << std::endl;
      }

      // do they match the current endpoints?
      for(MBRange::iterator j=endpt_sets.begin(); j!=endpt_sets.end(); ++j) {
        MBRange endpt_vert;
        result = MBI()->get_entities_by_handle( *j, endpt_vert );
        if(MB_SUCCESS != result) return result;
        assert(MB_SUCCESS==result);
        if(1 != endpt_vert.size()) {
	  std::cout << "curve " << gen::geom_id_by_handle(*i) 
                    << " has" << endpt_vert.size() 
                    << " endpoint vertices in the geometric vertex set" << std::endl;     
          return MB_INVALID_SIZE;
        }
        if(endpt_vert.front()!=ordered_verts.front() &&
          endpt_vert.front()!=ordered_verts.back() ) {
	  std::cout << "curve " << gen::geom_id_by_handle(*i) 
                    << " endpt sets do not match" << std::endl;
        }
      }

      // create the edges of the curve
      std::vector<MBEntityHandle> ordered_edges(ordered_verts.size()-1);
      for(unsigned int j=0; j<ordered_verts.size()-1; ++j) {
        MBEntityHandle edge;
        result = MBI()->create_element( MBEDGE, &ordered_verts[j], 2, edge );
        assert(MB_SUCCESS==result);
        if(MB_SUCCESS != result) return result;
        ordered_edges[j] = edge;
      }

      // If the curve is a loop, remove the duplicate endpoint.
      if(is_loop) ordered_verts.pop_back();

      // clear the set
      result = MBI()->clear_meshset( &(*i), 1 );
      assert(MB_SUCCESS==result);
      if(MB_SUCCESS != result) return result;

      // add the verts then edges to the curve set
      result = MBI()->add_entities( *i, &ordered_verts[0], ordered_verts.size() );
      assert(MB_SUCCESS==result);
      if(MB_SUCCESS != result) return result;
      result = MBI()->add_entities( *i, &ordered_edges[0], ordered_edges.size() );
      assert(MB_SUCCESS==result);
      if(MB_SUCCESS != result) return result;
    }
    return MB_SUCCESS;
  }

MBErrorCode get_geom_size_before_sealing( const MBRange geom_sets[], 
                                          const MBTag geom_tag,
                                          const MBTag size_tag,
                                          bool verbose ) 
{
  MBErrorCode rval;
  for( int dim = 1 ; dim < 4 ; ++dim ) 
    {
    for(MBRange::iterator i=geom_sets[dim].begin() ; i != geom_sets[dim].end() ; i++) 
      {
	double size;
	//std::cout << "dim = " << dim << " *i =" << *i << std::endl;
	rval = gen::measure( *i, geom_tag, size, verbose );
	//std::cout << " here in gen mesaure" << std::endl;
	if(gen::error(MB_SUCCESS!=rval,"could not measure")) 
	  {
	    return rval;
	  }

	//std::cout <<  " *i = " << *i << " size = " << size << std::endl;

	rval = MBI()->tag_set_data( size_tag, &(*i), 1, &size );
	//std::cout << " here in set tag data" << std::endl;
	if(gen::error(MB_SUCCESS!=rval,"could not set size tag")) 
	  {
	    return rval;
	  }
      }
    }
  if (verbose)
  { 
  std::cout << "finished in get_geom_size_before_sealing" << std::endl;
  }
  return MB_SUCCESS;
}

MBErrorCode get_geom_size_after_sealing( const MBRange geom_sets[], 
                                         const MBTag geom_tag,
                                         const MBTag size_tag,
                                         const double FACET_TOL,
                                         bool verbose ) {
  const bool debug = false;
      // save the largest difference for each dimension
      struct size_data {
        double orig_size, new_size, diff, percent;
        int id;
      };
      size_data largest_diff[3];
      size_data largest_percent[3];
      size_data smallest_size[3];

      MBErrorCode rval;
      for(unsigned dim=1; dim<4; dim++) {
        largest_diff[dim-1].diff       = 0; 
        largest_percent[dim-1].percent = 0; 
        smallest_size[dim-1].orig_size = std::numeric_limits<int>::max(); 
	for(MBRange::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
	  double orig_size = 0, new_size = 0;
	  rval = MBI()->tag_get_data( size_tag, &(*i), 1, &orig_size );
	  if(MB_SUCCESS != rval) {
	    std::cout << "rval=" << rval << " id=" << gen::geom_id_by_handle(*i) << std::endl;
	  }
	  assert(MB_SUCCESS == rval);
	  rval = gen::measure( *i, geom_tag, new_size, verbose );
	  assert(MB_SUCCESS == rval);

          // Remember the largest difference and associated percent difference
          double diff = fabs(new_size - orig_size);
          double percent_diff = 100.0*diff/orig_size;
          if(diff > largest_diff[dim-1].diff) {
            largest_diff[dim-1].orig_size = orig_size;
            largest_diff[dim-1].new_size  = new_size;
            largest_diff[dim-1].diff    = diff;
            largest_diff[dim-1].percent = percent_diff;
            largest_diff[dim-1].id      = gen::geom_id_by_handle(*i);
          }
          if(orig_size < smallest_size[dim-1].orig_size) {
            smallest_size[dim-1].orig_size = orig_size;
            smallest_size[dim-1].new_size  = new_size;
            smallest_size[dim-1].diff      = diff;
            smallest_size[dim-1].percent   = percent_diff;
            smallest_size[dim-1].id        = gen::geom_id_by_handle(*i);
          }
          if(percent_diff > largest_percent[dim-1].percent) {
            largest_percent[dim-1].orig_size = orig_size;
            largest_percent[dim-1].new_size  = new_size;
            largest_percent[dim-1].diff      = diff;
            largest_percent[dim-1].percent   = percent_diff;
            largest_percent[dim-1].id        = gen::geom_id_by_handle(*i);
          }

	  bool print_warning = false;
          // PROBLEM: There is no analytical "maximum" change for this. There are special
          // cases in each that could be infinitely large. For example, a curve can oscillate
          // up and down infinitely and still be within FACET_TOL of the geometric curve.
	  // The curve is changed only by merging vertices within FACET_TOL.
	  if(1 == dim) {
            if(FACET_TOL < diff) print_warning = true;

	  // The surface cannot change more than its zipped curves.
	  } else if(2 == dim) {
	    MBRange curve_sets;
	    rval = MBI()->get_child_meshsets( *i, curve_sets );
	    double total_length = 0;
	    for(MBRange::iterator j=curve_sets.begin(); j!=curve_sets.end(); ++j) {
	      double length;
	      rval = MBI()->tag_get_data( size_tag, &(*j), 1, &length );
	      total_length += length;
	    }
	    if(total_length*FACET_TOL < fabs(new_size - orig_size)) print_warning = true;
	    // The volume cannot change more than the "cylinder" of error around each curve.
	  } else if(3 == dim) {
	    MBRange surf_sets;
	    rval = MBI()->get_child_meshsets( *i, surf_sets );
	    double total_length = 0;
	    for(MBRange::iterator j=surf_sets.begin(); j!=surf_sets.end(); ++j) {          
	      MBRange curve_sets;
	      rval = MBI()->get_child_meshsets( *j, curve_sets );
	      for(MBRange::iterator k=curve_sets.begin(); k!=curve_sets.end(); ++k) {
		double length;
		rval = MBI()->tag_get_data( size_tag, &(*j), 1, &length );
		total_length += length;
	      }
	    }
	    if(total_length*FACET_TOL*FACET_TOL*3.14 < fabs(new_size - orig_size)) {
              print_warning = true;
            }
	  } else {
	    return MB_FAILURE;
	  }
	  if(print_warning && debug) {
	    std::cout << "  dim=" << dim << " id=" << gen::geom_id_by_handle(*i)
		      << " orig_size=" << orig_size << " new_size=" << new_size << std::endl;
	  }
	}
      }
      // print largest size change among each dimension
      std::cout << "Summary of deformation due to sealing: largest absolute change" << std::endl;
      std::cout.width(6);
      for(unsigned dim=1; dim<4; dim++) {
	std::cout << "  dim="            << dim
                  << ", id="             << largest_diff[dim-1].id
                  << ", orig_size="      << largest_diff[dim-1].orig_size
                  << ", new_size="       << largest_diff[dim-1].new_size
                  << ", abs_change="     << largest_diff[dim-1].diff 
                  << ", percent_change=" << largest_diff[dim-1].percent << std::endl;
      }
      std::cout << "Summary of deformation due to sealing: largest percent change" << std::endl;
      for(unsigned dim=1; dim<4; dim++) {
	std::cout << "  dim="            << dim
                  << ", id="             << largest_percent[dim-1].id
                  << ", orig_size="      << largest_percent[dim-1].orig_size
                  << ", new_size="       << largest_percent[dim-1].new_size
                  << ", abs_change="     << largest_percent[dim-1].diff
  	          << ", percent_change=" << largest_percent[dim-1].percent << std::endl;
      }
      std::cout << "Summary of deformation due to sealing: smallest size" << std::endl;
      for(unsigned dim=1; dim<4; dim++) {
	std::cout << "  dim="            << dim 
                  << ", id="             << smallest_size[dim-1].id
                  << ", orig_size="      << smallest_size[dim-1].orig_size
                  << ", new_size="       << smallest_size[dim-1].new_size
                  << ", abs_change="     << smallest_size[dim-1].diff
  	          << ", percent_change=" << smallest_size[dim-1].percent << std::endl;
      }
      std::cout.unsetf(std::ios::scientific|std::ios::showpos);
      return MB_SUCCESS;
}

MBErrorCode make_mesh_watertight(MBEntityHandle input_set, double &facet_tol, bool verbose) 
  {



    MBErrorCode result;

    //added to this function because they are called in make_watertight but aren't used until here
    const bool debug = false;
    const bool check_geom_size = true;
    // duplicated here from make_watertight, seemed to always be set to false
    bool is_acis=false;

    // create tags
   
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

    double sme_resabs_tol=1e-6;
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
    if(verbose)
    {
    std::cout << "  faceting tolerance=" << facet_tol << " cm" << std::endl;
    std::cout << "  absolute tolerance=" << sme_resabs_tol << " cm" << std::endl;
    }

    // get all geometry sets
    MBRange geom_sets[4];
    for(unsigned dim=0; dim<4; dim++) 
      {
	void *val[] = {&dim};
	result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
        if(verbose)
        {
	std::cout << "Get entities by type and tag" << std::endl;
        }
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

    

    // this could be related to when there are sat files rather than mesh?
    // If desired, find each entity's size before sealing.
    ///*
    if(check_geom_size) 
      {
	//std::cout << "I am checking the geometry size" << std::endl;
	result = get_geom_size_before_sealing( geom_sets, geom_tag, size_tag, verbose );
	if(gen::error(MB_SUCCESS!=result,"measuring geom size failed"))
	  {
	    return result;
	  }
      }
    
    //*/
    
    if (verbose) std::cout << "Get entity count before sealing" << std::endl;
    // Get entity count before sealing.
    int orig_n_tris;
    result = MBI()->get_number_entities_by_type( 0, MBTRI, orig_n_tris );
    std::cout << result << std::endl;

    assert(MB_SUCCESS == result);
    if(verbose)
    {
    std::cout << "  input faceted geometry contains " << geom_sets[3].size() << " volumes, " 
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size() 
              << " curves, and " << orig_n_tris << " triangles" << std::endl;  
    
    std::cout << "Finding degenerate triangles " << std::endl;
    }
    result = find_degenerate_tris();
    if(gen::error(result!=MB_SUCCESS,"could not determine if triangles were degenerate or not"))
      {
	return(result);
      }

    result = prepare_curves(geom_sets[1], geom_tag, id_tag, merge_tag, FACET_TOL, debug, verbose);
    if(gen::error(result!=MB_SUCCESS,"could not prepare the curves"))
      {
	return(result);
      }

    if (verbose) 
    {
    std::cout << "Zipping loops and removing small surfaces whose curves were all merged as pairs..." << std::endl;
    std::cout << "SME_RESABS_TOL = " << SME_RESABS_TOL << " FACET_TOL = " << FACET_TOL << std::endl;
    }
    result = prepare_surfaces(geom_sets[2], geom_tag, id_tag, normal_tag, merge_tag,
                              orig_curve_tag,SME_RESABS_TOL, FACET_TOL, debug);
    if ( result != MB_SUCCESS ) 
      {
	std::cout << "I have failed to zip" << std::endl;
	return result;
      }

    // After zipping surfaces, merged curve entity_to_deletes are no longer needed.
    // Swap surface parent-childs for curves to delete with curve to keep. 
    // ARE THEIR ORPHANED CHILD VERTEX SETS STILL AROUND? 
    if(verbose) std::cout << "Adjusting parent-child links then removing merged curves..." << std::endl;
    MBRange curves_to_delete;
    result = MBI()->get_entities_by_type_and_tag(0, MBENTITYSET, &merge_tag, NULL,
                                                 1, curves_to_delete);
    assert(MB_SUCCESS == result);
    // loop over the curves to delete 
    for(MBRange::const_iterator i=curves_to_delete.begin(); i!=curves_to_delete.end(); ++i) {
      // get the curve_to_keep
      MBEntityHandle curve_to_keep;
      result = MBI()->tag_get_data( merge_tag, &(*i), 1, &curve_to_keep );
      if(MB_SUCCESS != result) return result;
      // get parent surface of the curve_to_delete
      MBRange parent_surfs;
      result = MBI()->get_parent_meshsets( *i, parent_surfs );
      if(MB_SUCCESS != result) return result;
      // remove the curve_to_delete and replace with curve_to_keep
      for(MBRange::iterator j=parent_surfs.begin(); j!=parent_surfs.end(); ++j) {
        result = MBI()->remove_parent_child( *j, *i );
        if(MB_SUCCESS != result) return result;
        result = MBI()->add_parent_child( *j, curve_to_keep );
        if(MB_SUCCESS != result) return result;
      }
    }
    result = MBI()->delete_entities( curves_to_delete );
    assert(MB_SUCCESS == result);
    if ( result != MB_SUCCESS ) 
      {
	std::cout << "Houston, we have a problem" << std::endl;
      }
    geom_sets[1] = subtract(geom_sets[1], curves_to_delete );
    if(debug) std::cout << "deleted " << curves_to_delete.size() << " curves." << std::endl;

    // SHOULD COINCIDENT SURFACES ALSO BE MERGED?
    // 20100205 Jason: If child curves are the same, check to make sure each
    // facet pt is within 2x tol of the opposing surface.

    // This function is still screwed up, but 99% right.
    if (verbose) std::cout << "Fixing inverted triangles..." << std::endl;
    result = fix_normals(geom_sets[2], id_tag, normal_tag, debug, verbose);
    assert(MB_SUCCESS == result);

    /*
    // As sanity check, did zipping drastically change the entity's size?
    if(check_geom_size) {
      std::cout << "Checking size change of zipped entities..." << std::endl;
      result = get_geom_size_after_sealing( geom_sets, geom_tag, size_tag, FACET_TOL );
      if(gen::error(MB_SUCCESS!=result,"measuring geom size failed")) return result;
    }
    */

    // This tool stores curves as a set of ordered edges. MOAB stores curves as
    // ordered vertices and edges. MOAB represents curve endpoints as geometry
    // sets containing a singe vertex. This function restores MOAB's curve
    // representation.
    if (verbose) std::cout << "Restoring faceted curve representation..." << std::endl;
    result = restore_moab_curve_representation( geom_sets[1] );
    if(gen::error(MB_SUCCESS!=result,"restore_moab_curve_representation failed")) return result;    
    // If all of a volume's surfaces have been deleted, delete the volume.
    if (verbose) std::cout << "Removing small volumes if all surfaces have been removed..." << std::endl;
    for(MBRange::iterator i=geom_sets[3].begin(); i!=geom_sets[3].end(); ++i) {
      int n_surfs;
      result = MBI()->num_child_meshsets( *i, &n_surfs );
      assert(MB_SUCCESS == result);
      if(0 == n_surfs) {
        // Remove the volume set. This also removes parent-child relationships.
	std::cout << "  deleted volume " << gen::geom_id_by_handle(*i)  << std::endl;
        result = MBI()->delete_entities( &(*i), 1);
        assert(MB_SUCCESS == result);
        i = geom_sets[3].erase(i) - 1;
      } 
    }

    // The obbtrees are no longer valid because the triangles have been altered.
    // Surface and volume sets are tagged with tags holding the obb tree
    // root handles. Somehow, delete the old tree without deleting the
    // surface and volume sets, then build a new tree.    
    
    // removing this for now, has something to do with an interaction with DAGMC 
    // which doesn't actually occur any more
    //if (verbose) std::cout << "Removing stale OBB trees..." << std::endl;
    //result = cleanup::remove_obb_tree();
    //assert(MB_SUCCESS == result);

    //std::cout << "INSERT FUNCTION HERE TO REMOVE STALE VERTS, EDGES, TRIS, VERT SETS, ETC"
    //        << std::endl;

    // Resetting meshsets so that they no longer track.
    if (verbose) std::cout << "Restoring original meshset options and tags..." << std::endl;
    for(unsigned dim=0; dim<4; dim++) {
      for(MBRange::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
        result = MBI()->set_meshset_options(*i, 1==dim ? MESHSET_ORDERED : MESHSET_SET);
        if(MB_SUCCESS != result) return result;
      }
    }    

    // Tags for merging curves and checking the change in geometry size were added.
    // Delete these because they are no longer needed.
    result = MBI()->tag_delete( normal_tag );
    if(MB_SUCCESS != result) return result;
    result = MBI()->tag_delete( merge_tag );
    if(MB_SUCCESS != result) return result;
    result = MBI()->tag_delete( size_tag );
    if(MB_SUCCESS != result) return result;
    result = MBI()->tag_delete( orig_curve_tag );
    if(MB_SUCCESS != result) return result;

    // Write output file
    int sealed_n_tris;
    if (verbose) std::cout << "Writing zipped file..." << std::endl;
    result = MBI()->get_number_entities_by_type( 0, MBTRI, sealed_n_tris );
    assert(MB_SUCCESS == result);
    if (verbose)
    {
    std::cout << "  output file contains " << geom_sets[3].size() << " volumes, " 
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size() 
              << " curves, and " << sealed_n_tris << " triangles" << std::endl;  
    std::cout << "  triangle count changed " << (double)sealed_n_tris/orig_n_tris
              << "x (sealed/unsealed)" << std::endl;
    }
    return MB_SUCCESS;  
  }

}

