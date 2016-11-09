// ********************************************************************
// Patrick Shriwise
// October, 2013

// functions needed to seal unwatertight models in make_watertight


// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "MakeWatertight.hpp"
#include "Arc.hpp"
#include "Zip.hpp"
#include "Cleanup.hpp"


moab::ErrorCode MakeWatertight::delete_all_edges()
{
  // delete all of the edges. Never keep edges. They are too hard to track and use
  // due to orientation and multiple entities errors when merging.
  moab::ErrorCode result;
  moab::Range edges;
  result = MBI()->get_entities_by_type( 0, moab::MBEDGE, edges );
  MB_CHK_SET_ERR(result,"could not get edges");
  assert(moab::MB_SUCCESS == result);
  result = MBI()->delete_entities( edges );
  MB_CHK_SET_ERR(result,"could not delete edges");
  assert(moab::MB_SUCCESS == result);
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::find_degenerate_tris()
{
  moab::ErrorCode result;
  moab::Range tris;
  result = MBI()->get_entities_by_type( 0, moab::MBTRI, tris );
  MB_CHK_SET_ERR(result,"could not get tris");
  assert(moab::MB_SUCCESS == result);
  int counter = 0;
  for(moab::Range::const_iterator i=tris.begin(); i!=tris.end(); ++i) {
    if( gen->triangle_degenerate(*i) ) {
      result = MBI()->list_entity(*i);
      MB_CHK_SET_ERR(result,"found degenerate tri");
      assert(moab::MB_SUCCESS == result);
      ++counter;
    }
  }
  if(counter != 0) {
    std::cout << "Found " << counter << " degenerate triangles. " << std::endl;
  }
  return moab::MB_SUCCESS;
}

// Input: Possibly unordered sets of curves that do track ownership. Curves contain
//        edges and vertices. Parents sets are surfaces. Child sets are endpoint
//        vertices.
// Output: Ordered sets of verts that do track ownership. All edges are deleted.
moab::ErrorCode MakeWatertight::prepare_curves(moab::Range &curve_sets,
    moab::Tag geom_tag, moab::Tag id_tag, moab::Tag merge_tag,
    const double FACET_TOL, const bool debug, bool verbose )
{
  moab::ErrorCode result;
  if (verbose) std::cout << "Modifying faceted curve representation and removing small curves..."
                           << std::endl;

  // process each curve
  for(moab::Range::iterator i=curve_sets.begin(); i!=curve_sets.end(); i++ ) {
    // get the curve id of the curve meshset
    int id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &id );
    MB_CHK_SET_ERR(result,"could not get id tag");
    if(debug) std::cout << "curve " << id << std::endl;

    // get the range of edges of the curve meshset
    std::vector<moab::EntityHandle> curve_edges;
    result = MBI()->get_entities_by_type( *i, moab::MBEDGE, curve_edges );
    MB_CHK_SET_ERR(result,"could not get curve_edges");

    /* Merge the endpoints of the curve and remove its edges if it is too small.
    Use the MERGE_TOL because these edges will be merged with the MERGE_TOL
    during zipping anyhow. Doing this now removes small curves from zipping and
    reduces ambiguity. */
    if(FACET_TOL > gen->length(curve_edges)) {
      if (debug) {
        std::cout << "  deleted curve " << id << ", length=" << gen->length(curve_edges)
                  << " cm, n_verts=" << curve_edges.size()+1 << std::endl;
      }
      // get the endpoints of the curve
      moab::Range endpt_sets;
      result = MBI()->get_child_meshsets( *i, endpt_sets );
      MB_CHK_SET_ERR(result,"could not get curve child sets");

      if(endpt_sets.empty()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"curve has no child sets");
      } else if(1 == endpt_sets.size()) {
        // The edges are no longer needed. Remove them before altering the range
        // by deleting degenerate edges below.
        result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
        MB_CHK_SET_ERR(result,"could not delete edges");

      } else if(2 == endpt_sets.size()) {
        // The edges are no longer needed. Remove them before altering the range
        // by deleting degenerate edges below.
        result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
        MB_CHK_SET_ERR(result,"could not delete edges");

        moab::Range front_endpt, back_endpt;
        result = MBI()->get_entities_by_type( endpt_sets.front(), moab::MBVERTEX, front_endpt);
        MB_CHK_SET_ERR(result,"could not get vert from front endpt set");
        if(1!=front_endpt.size()) {
          MB_CHK_SET_ERR(moab::MB_FAILURE,"front endpt set does not have 1 vert");
        }

        result = MBI()->get_entities_by_type( endpt_sets.back(), moab::MBVERTEX, back_endpt);
        MB_CHK_SET_ERR(result,"could not get vert from back endpt set");
        if(1!=back_endpt.size()) {
          MB_CHK_SET_ERR(moab::MB_FAILURE,"back endpt set does not have 1 vert");
        }

        // merge the endpoints-ALWAYS CHECK TO AVOID MERGING THE SAME ENTITY!!!
        if(front_endpt[0] != back_endpt[0]) {
          std::vector<moab::EntityHandle> temp;
          result = zip->merge_verts( front_endpt.front(), back_endpt.front(), temp, temp );
          MB_CHK_SET_ERR(result,"could not merge verts");

          // check for and remove degenerate edges caused by the merge
          moab::Range edges;
          moab::EntityHandle temp_pt = front_endpt[0];
          result = MBI()->get_adjacencies( &temp_pt, 1, 1, false, edges);
          MB_CHK_SET_ERR(result,"could not get adj edges");

          for(moab::Range::iterator j=edges.begin(); j!=edges.end(); j++) {
            const moab::EntityHandle *conn;
            int n_verts;
            result = MBI()->get_connectivity( *j, conn, n_verts);
            MB_CHK_SET_ERR(result,"could not get edge conn");

            if(conn[0] == conn[1]) {
              result = MBI()->delete_entities( &(*j), 1 );
              MB_CHK_SET_ERR(result,"could not delete degenerate edge");
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
      MB_CHK_SET_ERR(result,"could not delete curve set");
      i = curve_sets.erase(i) - 1;
    } else {

      // convert the curve of edges into a curve of verts
      std::vector<moab::EntityHandle> ordered_verts;
      result = gen->ordered_verts_from_ordered_edges( curve_edges, ordered_verts);
      MB_CHK_SET_ERR(result,"could not order_verts_by_edge");

      // replace the unordered edges with the ordered verts
      result = arc->set_meshset( *i, ordered_verts );
      MB_CHK_SET_ERR(result,"could not set_meshset");

      // The edges are no longer needed.
      result = MBI()->delete_entities( &curve_edges[0], curve_edges.size() );
      MB_CHK_SET_ERR(result,"could not delete edges");
    }
  }

  // merge curves that are the same within facet_tol
  if (verbose) std::cout << "Identifying coincident curves to be merged..." << std::endl;
  result = arc->merge_curves(curve_sets, FACET_TOL, id_tag, merge_tag, debug );
  MB_CHK_SET_ERR(result,"could not merge_curves");

  return moab::MB_SUCCESS;
}

/* The chosen curve may need to be reversed to be in the same direction as the skin.
   The skin itself is not changed. */
/* PROBLEM: Some sliver surfaces are smaller in every dimension than the 5x
 facet_tol that the skin was found to vary within.
 This method finds the average distance between a curve and skin arc. The
 smallest distance is the next curve. */
/// Cut an arc out of the skin. Return a corresponding curve, the curve's set, and
/// is the curve has ben reversed. The returned skin has the arc cut away. The
/// returned vector of curve sets has the curve removed.
moab::ErrorCode MakeWatertight::create_arc_pair(  const double FACET_TOL,
    const moab::EntityHandle surf_set,
    std::vector<moab::EntityHandle> &skin_loop,
    std::vector<moab::EntityHandle> &curve_sets,
    const moab::EntityHandle front_endpt,
    const bool debug,
    moab::EntityHandle &curve_set,
    bool &curve_is_reversed,
    std::vector<moab::EntityHandle> &curve,
    std::vector<moab::EntityHandle> &skin_arc )
{

  /* Now we have a topological connection between the curve and skin. Find other
     curves that share this vert. Be aware if the curve is reversed so the
     original direction can be preserved when done zipping. */
  moab::ErrorCode rval;
  if(debug) {
    std::cout << curve_sets.size() << " curves remain to be zipped"
              << " skin_loop.size()=" << skin_loop.size() << " front_endpt="
              << front_endpt << " skin:" << std::endl;
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
    moab::Range endpt_sets;
    rval = MBI()->get_child_meshsets(curve_sets[i], endpt_sets );
    MB_CHK_SET_ERR(rval,"could not get endpt_sets");
    if(endpt_sets.empty() || 2<endpt_sets.size()) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"too many endpt_sets");
    }
    // get the vertex handles
    std::vector<moab::EntityHandle> endpts;
    for(unsigned j=0; j<endpt_sets.size(); ++j) {
      moab::Range endpt;
      rval = MBI()->get_entities_by_type( endpt_sets[j], moab::MBVERTEX, endpt );
      MB_CHK_SET_ERR(rval,"could not get endpt");
      if(1!=endpt.size()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"not one endpt");
      }
      endpts.push_back( endpt.front() );
      if(debug) std::cout << "curve " << gen->geom_id_by_handle(curve_sets[i])
                            << " endpt=" << endpt.front() << std::endl;
    }

    // if an endpt is not the front_endpt, then this curve isn't adjacent
    if(front_endpt!=endpts.front() && front_endpt!=endpts.back()) continue;

    // get the point representation
    std::vector<moab::EntityHandle> temp_curve;
    rval = arc->get_meshset( curve_sets[i], temp_curve );
    MB_CHK_SET_ERR(rval,"could not get curve set");
    if(2>temp_curve.size()) std::cout << "warning11: curve is degenerate" << std::endl;
    if(debug) {
      std::cout << "  adj curve " << gen->geom_id_by_handle(curve_sets[i]) << ":" << std::endl;
      //gen->print_loop( temp_curve );
    }

    // Get the curve-surface relative sense. From CGMA/builds/dbg/include/CubitDefines,

    int sense;
    if(debug) {
      std::cout << "surf_set = " << gen->geom_id_by_handle(surf_set) << std::endl;
      std::cout << "curve_set = " << gen->geom_id_by_handle(curve_sets[i]) << std::endl;
    }
    rval = gen->get_curve_surf_sense( surf_set, curve_sets[i], sense );
    MB_CHK_SET_ERR(rval,"could not get_curve_surf_sense");

    // get the curve length, for efficient find_closest_vert, plus a tolerance.
    // This also helps to find the correct skin arc in special (~1D surfs) cases
    // where the closest skin pt is not the correct skin.
    const double temp_curve_len = gen->length(temp_curve);
    const double extra = 1.0;

    // is it forward oriented? To allow for ambiguous cases, instead of accepting
    // only forward-oriented curves, only reject reverse-oriented curves.
    if(front_endpt==temp_curve.front() && moab::SENSE_REVERSE!=sense) {
      // find the closest skin pt to the curve endpt
      unsigned pos;
      // if only one curve is left, take the entire skin. Otherwise, due to
      // coordinate drift via merging the closest skin pt may not be the correct one.
      if(1==curve_sets.size()) {
        pos = skin_loop.size()-1;
      } else {
        rval = gen->find_closest_vert( temp_curve.back(), skin_loop, pos, temp_curve_len+extra );
        MB_CHK_SET_ERR(rval,"could not find_closest_vert");
      }
      if(debug) std::cout << "  end of skin arc=" << skin_loop[pos] << std::endl;
      // SPECIAL CASE: If the skin is a circle, create an arc out of the circle
      if(skin_loop[pos] == skin_loop.front()) pos = skin_loop.size()-1;

      // create a skin arc to test against the curve
      std::vector<moab::EntityHandle> temp_skin(skin_loop.begin(), skin_loop.begin()+pos+1);
      double d;
      rval = gen->dist_between_arcs( debug, temp_skin, temp_curve, d );
      MB_CHK_SET_ERR(rval,"could not get dist_between_arcs");
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
    if(front_endpt==temp_curve.back() && moab::SENSE_FORWARD!=sense) {
      reverse( temp_curve.begin(), temp_curve.end() );
      // find the closest skin pt to the curve endpt
      unsigned pos;
      // if only one curve is left, take the entire skin. Otherwise, due to
      // coordinate drift via merging the closest skin pt may not be the correct one.
      if(1==curve_sets.size()) {
        pos = skin_loop.size()-1;
      } else {
        rval = gen->find_closest_vert( temp_curve.back(), skin_loop, pos, temp_curve_len+extra );
        MB_CHK_SET_ERR(rval,"could not find_closest_vert");
      }
      if(debug) std::cout << "  end of skin arc=" << skin_loop[pos] << std::endl;
      // SPECIAL CASE: If the skin is a circle, create an arc out of the circle
      if(skin_loop[pos] == skin_loop.front()) pos = skin_loop.size()-1;

      // create a skin arc to test against the curve
      std::vector<moab::EntityHandle> temp_skin(skin_loop.begin(), skin_loop.begin()+pos+1);
      double d;
      rval = gen->dist_between_arcs( debug, temp_skin, temp_curve, d );
      MB_CHK_SET_ERR(rval,"could not get dist_between_arcs");
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
      std::cout << "  curve " << gen->geom_id_by_handle(curve_sets[i])
                << " is unsealed" << std::endl;
    }
    return moab::MB_FAILURE;
  }

  // If the average distance between the closest skin and curve is too far...
  if(100*FACET_TOL<=min_dist) {
    std::cout << "  warning8: curve too far from skin, average_dist="
              << min_dist << std::endl;
    if(1.0<=min_dist) return moab::MB_FAILURE;
  }

  // remove the chosen curve the set of unsealed curves
  curve_sets.erase( curve_sets.begin()+curve_set_idx );

  // remove the chosen skin arc from the rest of the skin
  skin_loop.erase( skin_loop.begin(), skin_loop.begin()+skin_pos );

  // If the entire skin loop has been sectioned into arcs, a single point remains.
  // This is because the skin_arc needs to have the same back point as the front
  // of the remaining skin_loop. If a single point remains, remove it.
  if(1==skin_loop.size()) skin_loop.clear();

  if(debug) std::cout << "  curve " << gen->geom_id_by_handle(curve_set)
                        << " paired with skin, min_dist =" << min_dist << std::endl;
  return moab::MB_SUCCESS;
}


///Runs the make_watertight algorithm on the edge (curve) and skin (arc).
/// See PHD Thesis of Brandon Smith University of Wisconsin - Madison 2011
/// for details


//  -Instead of altering the skin and curve vectors, make them const. Put the
// sealed curve in a new vector, using only push_back. This
// would avoid all of the really slow inserts and erases in the curve and skin vectors.
moab::ErrorCode MakeWatertight::seal_arc_pair( const bool debug,
    const double FACET_TOL,
    const moab::Tag normal_tag,
    std::vector<moab::EntityHandle> &edge, /* in */
    std::vector<moab::EntityHandle> &skin /* in/out */,
    const int surf_id )
{
  if(debug) {
    std::cout << "edge before sealing:" << std::endl;
    gen->print_loop(edge);
    std::cout << "skin before sealing:" << std::endl;
    gen->print_loop(skin);
  }

  moab::ErrorCode rval;
  const double TOL_SQR = FACET_TOL*FACET_TOL;
  if(edge.empty() || skin.empty()) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"edge or skin has no verts");
  }

  //**************************************************************************
  // Merge the front of the skin to the front of the curve
  //**************************************************************************
  {
    moab::EntityHandle keep_vert   = edge.front();
    moab::EntityHandle delete_vert = skin.front();
    if(keep_vert != delete_vert) {
      double merge_dist;
      rval = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
      MB_CHK_SET_ERR(rval,"could not get merge_dist g");
      if(debug) {
        std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" << keep_vert
                  << " merge_dist=" << merge_dist << std::endl;
      }
      if(FACET_TOL < merge_dist) {
        std::cout << "  warning0: front pt merge_dist=" << merge_dist << std::endl;
      }
      rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
      MB_CHK_SET_ERR(rval,"could not merge verts g");
    }
  }

  //**************************************************************************
  // Merge the back of the skin to the back of the curve
  //**************************************************************************
  {
    moab::EntityHandle keep_vert   = edge.back();
    moab::EntityHandle delete_vert = skin.back();
    if(keep_vert != delete_vert) {
      double merge_dist;
      rval = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
      MB_CHK_SET_ERR(rval,"could not get merge_dist h");
      if(debug) {
        std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert=" << keep_vert
                  << " merge_dist=" << merge_dist << std::endl;
      }
      if(FACET_TOL < merge_dist) {
        std::cout << "  warning0: surf " << surf_id << " back pt merge_dist="
                  << merge_dist << std::endl;
        if(1000*FACET_TOL < merge_dist) return moab::MB_FAILURE;
      }
      rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
      MB_CHK_SET_ERR(rval,"could not merge verts g");
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

    rval = gen->squared_dist_between_verts( edge[e_pos-1], edge[e_pos], e_dist );
    MB_CHK_SET_ERR(rval,"could not get e_dist");
    rval = gen->squared_dist_between_verts( edge[e_pos-1], skin[s_pos], s_dist );
    MB_CHK_SET_ERR(rval,"could not get s_dist");
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
      rval = gen->line_point_dist( edge[e_pos-1], edge[e_pos], skin[s_pos], move_dist );
      MB_CHK_SET_ERR(rval,"could not get line_point_dist");
      if(10*FACET_TOL < move_dist) {
        std::cout << "  warning5: surf " << surf_id << " vertex move_dist="
                  << move_dist << std::endl;
      }
      rval = gen->point_line_projection( edge[e_pos-1], edge[e_pos], skin[s_pos]);
      MB_CHK_SET_ERR(rval,"could not get point_line_projection");
      rval = gen->squared_dist_between_verts( edge[e_pos-1], skin[s_pos], s_dist );
      MB_CHK_SET_ERR(rval,"could not get s_dist b");
      dist = s_dist;
      if(debug) std::cout << "skin is next, projected dist=" << dist << std::endl;
    } else {
      edge_is_next = true;
      dist = e_dist;
      if(debug) std::cout << "edge is next" << std::endl;
    }

    // find the cs_dist after moving the skin to the curve (if skin_is_next)
    rval = gen->squared_dist_between_verts( edge[e_pos], skin[s_pos], es_dist );
    MB_CHK_SET_ERR(rval,"could not get es_dist");

    // **************************************************************************
    // Merge with previous vert if it is too close
    // **************************************************************************
    if(dist < TOL_SQR) {
      moab::EntityHandle keep_vert = edge[e_pos-1];
      if(edge_is_next) {
        moab::EntityHandle delete_vert = edge[e_pos];
        if(keep_vert != delete_vert) { // cannot merge the same vert
          if(debug) {
            double merge_dist;
            rval = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
            MB_CHK_SET_ERR(rval,"could not get merge_dist");
            std::cout << "  merged edge_vert=" << delete_vert << " to edge_vert="
                      << keep_vert << " merge_dist=" << merge_dist <<std::endl;
          }
          rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
          MB_CHK_SET_ERR(rval,"could not merge_verts a");
        }
        if(edge.size() < e_pos+1) {
          std::cout << "edge.size()=" << edge.size() << " e_pos=" << e_pos << std::endl;
        }
        edge.erase( edge.begin() + e_pos );
      } else {
        moab::EntityHandle delete_vert = skin[s_pos];
        if(keep_vert != delete_vert) {
          if(debug) {
            double merge_dist;
            rval  = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
            MB_CHK_SET_ERR(rval,"could not get merge_dist b");
            std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert="
                      << keep_vert << " merge_dist=" << merge_dist << std::endl;
          }
          rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
          MB_CHK_SET_ERR(rval,"could not merge_verts b");
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
      moab::EntityHandle keep_vert  = edge[e_pos];
      if(skin[s_pos] != keep_vert) {
        moab::EntityHandle delete_vert= skin[s_pos];
        if(debug) {
          double merge_dist;
          rval  = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
          MB_CHK_SET_ERR(rval,"could not get merge_dist c");
          std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert="
                    << keep_vert << " merge_dist=" << merge_dist << std::endl;
        }
        rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
        MB_CHK_SET_ERR(rval,"could not merge_verts b");
      }
      s_pos++;
      e_pos++;
      // Are there more skin verts in proximity to the merged pt?
      while(true) {
        // check to see if skin still exists
        if(s_pos == skin.size()) break;
        // check if the distance is less than merge tol
        moab::EntityHandle delete_vert= skin[s_pos];
        double merge_dist;
        rval = gen->dist_between_verts( keep_vert, delete_vert, merge_dist);
        MB_CHK_SET_ERR(rval,"could not get merge_dist d");
        if(FACET_TOL < merge_dist) break;
        // merge the verts if they are not the same
        if(keep_vert != delete_vert) {
          if(debug) {
            std::cout << "  merged skin_vert=" << delete_vert << " to edge_vert="
                      << keep_vert << " merge_dist=" << merge_dist << std::endl;
          }
          rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
          MB_CHK_SET_ERR(rval,"could not merge_verts d");
        }
        skin.erase( skin.begin() + s_pos );
      }
      // Are there more curve verst in proximity to the merged pt?
      while(true) {
        // check to see if curve still exists
        if(e_pos == edge.size()) break;
        // check if the distance is less than merge tol
        moab::EntityHandle delete_vert= edge[e_pos];
        double merge_dist;
        rval = gen->dist_between_verts( keep_vert, delete_vert, merge_dist );
        MB_CHK_SET_ERR(rval,"could not get merge_dist e");
        if(FACET_TOL < merge_dist) break;
        // merge the verts if they are not the same
        if(keep_vert != delete_vert) {
          if(debug) {
            std::cout << "  merged edge_vert=" << delete_vert << " to edge_vert="
                      << keep_vert << " merge_dist=" << merge_dist << std::endl;
          }
          rval = zip->merge_verts( keep_vert, delete_vert, skin, edge );
          MB_CHK_SET_ERR(rval,"could not merge_verts e");
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
        rval = gen->line_point_dist( edge[e_pos-1], skin[s_pos], edge[e_pos], move_dist );
        MB_CHK_SET_ERR(rval,"could not get line_point_dist");
        if(10*FACET_TOL < move_dist) {
          std::cout << "  warning6: surf " << surf_id << " vertex move_dist="
                    << move_dist << std::endl;
        }
        rval = zip->t_joint( normal_tag, edge[e_pos-1], edge[e_pos], skin[s_pos], debug );
        MB_CHK_SET_ERR(rval,"tjoint failed a");
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
        rval = gen->line_point_dist( edge[e_pos-1], edge[e_pos], skin[s_pos], move_dist );
        MB_CHK_SET_ERR(rval,"could not get line_point_dist");
        if(10*FACET_TOL < move_dist) {
          std::cout << "  warning6: surf " << surf_id << " vertex move_dist="
                    << move_dist << std::endl;
        }
        rval = zip->t_joint( normal_tag, edge[e_pos-1], skin[s_pos], edge[e_pos], debug );
        MB_CHK_SET_ERR(rval,"tjoint failed b");
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
      rval = gen->squared_dist_between_verts( edge[e_pos-1], edge[e_pos-2], d );
      MB_CHK_SET_ERR(rval,"could not get dist");
      if(TOL_SQR > d) {
        std::cout << "zip_loop: d=" << d << std::endl;
        gen->print_vertex_coords(edge[e_pos-1]);
        gen->print_vertex_coords(edge[e_pos-2]);
        std::cout << "warning7: surf " << surf_id
                  << " adjacent edge points are closer than FACET_TOL" << std::endl;
      }
    }
    // The position should be the same. Do not exceed array bounds when checking.
    if(e_pos!=s_pos) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"skin and edge positions do not match");
    }
    if(edge[e_pos-1] != skin[s_pos-1]) {
      std::cout << "edge[" << e_pos-1 << "]=" << edge[e_pos-1]
                << " skin[" << s_pos-1 << "]=" << skin[s_pos-1] << std::endl;
    }
    if(edge[e_pos-1]!=skin[s_pos-1]) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"skin and edge vert does not match");
    }
  }

  // The skin and curve should be the same size
  if(edge.size()!=skin.size()) {
    std::cout << "    surf " << surf_id
              << " sealed skin and curve are not the same size" << std::endl;
    if(debug) {
      std::cout << "edge:" << std::endl;
      gen->print_loop(edge);
      std::cout << "skin:" << std::endl;
      gen->print_loop(skin);
    }
    //return moab::MB_FAILURE;
  }

  if(debug) {
    std::vector< std::vector<moab::EntityHandle> > temp;
    temp.push_back(edge);
    temp.push_back(skin);
    rval = zip->test_zipping(FACET_TOL, temp);
    MB_CHK_SET_ERR(rval,"sealing test failed");
  }

  return moab::MB_SUCCESS;

}

/// seals the skin_loop to the closest curves in curve sets in a watertight fashion
moab::ErrorCode MakeWatertight::seal_loop( bool debug,
    const double FACET_TOL,
    const moab::Tag normal_tag,
    const moab::Tag orig_curve_tag,
    const moab::EntityHandle surf_set,
    std::vector<moab::EntityHandle> &curve_sets,
    std::vector<moab::EntityHandle> &skin_loop,
    bool verbose)
{
  moab::ErrorCode rval;
  debug = false;
  // Find a curve that corresponds to the skin. Note that because there is more
  // than one skin loop, not all curves of the surface correspond to each skin
  // Loop.

  // To establish an inital connection, find the curve with endpoint closest
  // to a skin point.
  if(debug) std::cout << "seal_loop: new loop contains " << skin_loop.size()
                        << " skin pts" << std::endl;

  // If skin remains but all the curves are zipped an error has occured.
  if( curve_sets.empty() ) {
    std::cout << "seal_loop: no curves are left, but skin remains" << std::endl;
    gen->print_loop(skin_loop);
    skin_loop.clear();
    return moab::MB_FAILURE;
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
  moab::EntityHandle closest_skin_pt, closest_front_curve_endpt;
  double min_dist = std::numeric_limits<double>::max();
  for(unsigned i=0; i<curve_sets.size(); ++i) {
    // get geometric vertices
    moab::Range endpt_sets;
    rval = MBI()->get_child_meshsets(curve_sets[i], endpt_sets );
    MB_CHK_SET_ERR(rval,"could not get endpt_sets");
    if(endpt_sets.empty() || 2<endpt_sets.size()) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"too many endpt_sets");
    }

    std::vector<moab::EntityHandle> endpts;
    for(unsigned j=0; j<endpt_sets.size(); ++j) {
      moab::Range endpt;
      rval = MBI()->get_entities_by_type( endpt_sets[j], moab::MBVERTEX, endpt );
      MB_CHK_SET_ERR(rval,"could not get endpt");
      if(1!=endpt.size()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"not one endpt");
      }
      endpts.push_back( endpt.front() );
    }

    // check to ensure that the endpt sets aren't degenerate

    if(2==endpt_sets.size() && endpts.front()==endpts.back() && debug ) {

      std::cout << "  warning9: curve " << gen->geom_id_by_handle(curve_sets[i])
                << " geometric endpoints degenerate" << std::endl;
    }

    // check to ensure that geometric verts are the curve endpts
    std::vector<moab::EntityHandle> curve;
    rval = arc->get_meshset( curve_sets[i], curve );
    MB_CHK_SET_ERR(rval,"could not get_meshset");
    if(1==endpt_sets.size()) {
      if(curve.front()!=curve.back()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"endpt discrepancy");
      }
      if(curve.front()!=endpts.front()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"geometric verts inconsistent with curve");
      }
    } else {
      if(curve.front()==curve.back())
        if(debug) std::cout << "  warning10: degenerate curve endpts" << std::endl;
      if(curve.front()!=endpts.front() && curve.front()!=endpts.back()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"endpts not consistent");
      }
      if(curve.back()!=endpts.front() && curve.back()!=endpts.back()) {
        MB_CHK_SET_ERR(moab::MB_FAILURE,"endpts not consistent");
      }
    }

    // determine the orientation of the curve wrt the surf.
    int sense;
    if(debug) {
      std::cout << "surf_set = " << gen->geom_id_by_handle(surf_set) << std::endl;
      std::cout << "curve_set = " << gen->geom_id_by_handle(curve_sets[i]) << std::endl;
    }
    rval = gen->get_curve_surf_sense( surf_set, curve_sets[i], sense );
    MB_CHK_SET_ERR(rval,"could not get_curve_surf_sense");
    // select the front wrt the skin.
    moab::EntityHandle curve_endpt = (moab::SENSE_FORWARD==sense) ? curve.front() : curve.back();

    // find closest skin vert to front of curve
    std::vector<double> d;
    std::vector<unsigned> p;
    rval = gen->find_closest_vert( 0, curve_endpt, skin_loop, p, d);
    MB_CHK_SET_ERR(rval,"could not find_closest_vert");
    if(debug) std::cout << "zip_loop: loop-curve endpt dist=" << d.front() << " skin_vert="
                          << skin_loop[p.front()] << " curve="
                          << gen->geom_id_by_handle(curve_sets[i]) << " front_endpt="
                          << curve_endpt << std::endl;
    if(d.front() < min_dist) {
      min_dist = d.front();
      curve_idx = i;
      pos = p.front();
      closest_front_curve_endpt = curve_endpt;
      closest_skin_pt = skin_loop[p.front()];
    }
  }

  moab::EntityHandle front_endpt = closest_front_curve_endpt;

  // The closest points should be within facet tolerance
  if(100*FACET_TOL<min_dist && verbose) {
    std::cout << "closest skin pt farther than 100*FACET_TOL from curve vert ("
              << min_dist << ")" << std::endl;
    if(true) {
      std::cout << "  skin pt:" << std::endl;
      rval = MBI()->list_entity(closest_skin_pt);
      MB_CHK_SET_ERR(rval,"error listing skin pt");
      std::cout << "  curve vert:" << std::endl;
      rval = MBI()->list_entity(front_endpt);
      MB_CHK_SET_ERR(rval,"error listing curve_vert");
    }
    return rval;
  }

  if(debug) {
    std::cout << "closest skin vert=" << skin_loop[pos] << " pos=" << pos
              << " min_dist=" << min_dist << " curve="
              << gen->geom_id_by_handle(curve_sets[curve_idx])
              << " front_endpt=" << closest_front_curve_endpt << std::endl;
  }

  /* Take the skin loops and advance it such that our common point is at
     location 0. Ensure that that both endpoints are the same. */
  std::vector<moab::EntityHandle> temp_loop;
  temp_loop.reserve(skin_loop.size());
  std::vector<moab::EntityHandle>::iterator j = temp_loop.begin();
  std::vector<moab::EntityHandle>::iterator k = skin_loop.begin();
  // Do not insert the back endpoint (avoid duplicate).
  temp_loop.insert( j, k+pos, k+skin_loop.size()-1 );
  j = temp_loop.begin(); // j became invalid because temp_loop resized
  j += skin_loop.size() - pos - 1;
  temp_loop.insert( j, k, k+pos+1 );
  if(temp_loop.size()!=skin_loop.size()) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"loop size not conserved");
  }
  assert(temp_loop.size() == skin_loop.size()); // same size
  if(temp_loop[0]!=temp_loop[temp_loop.size()-1]) {
    MB_CHK_SET_ERR(moab::MB_FAILURE,"loop endpts not continuous");
  }
  assert(temp_loop[0] == temp_loop[temp_loop.size()-1]); // same endpoint
  skin_loop = temp_loop;
  if(debug) {
    std::cout << "  skin:" << std::endl;
  }

  // Create a set to store skin_loop so that handles are updated during merging.
  // This prevents stale handles, and avoids O(n) search through every handle in
  // the skin loop if manually updating the skin_loop vector.
  moab::EntityHandle skin_loop_set;
  rval = MBI()->create_meshset( moab::MESHSET_TRACK_OWNER|moab::MESHSET_ORDERED, skin_loop_set );
  MB_CHK_SET_ERR(rval,"creating skin_loop_set failed");

  while(!skin_loop.empty()) {
    //**************************************************************************
    // Select the curve adjacent to the common point that best matches the skin
    // using average distance. Chop the corresponding arc of skin, to form a
    // curve-skin_arc pair.
    //**************************************************************************

    bool curve_is_reversed;
    moab::EntityHandle curve_set;
    std::vector<moab::EntityHandle> curve, skin_arc;
    rval = create_arc_pair( FACET_TOL, surf_set, skin_loop, curve_sets, front_endpt,
                            debug, curve_set, curve_is_reversed, curve, skin_arc );
    MB_CHK_SET_ERR(rval,"  pair creation failed");

    // Let moab store skin loop to avoid stale vert handles from merging.
    rval = arc->set_meshset( skin_loop_set, skin_loop );
    MB_CHK_SET_ERR(rval,"setting skin_loop_set failed");

    // The original faceted curves are never used. Instead they are replaced by
    // skin. This reduces the number of new triangles created.
    int orig_curve;
    rval = MBI()->tag_get_data( orig_curve_tag, &curve_set, 1, &orig_curve );
    MB_CHK_SET_ERR(rval,"can't get tag");

    // If the tag is non-zero, the facet edge has already been replaced.
    if(orig_curve) {
      // this tag is used to mark the edge as updated
      int false_int = 0;
      rval = MBI()->tag_set_data( orig_curve_tag, &curve_set, 1, &false_int );
      MB_CHK_SET_ERR(rval,"can't set tag");

      // merge new endpoints to old endpoints
      if(curve.front()!=skin_arc.front()) {
        rval = zip->merge_verts( curve.front(), skin_arc.front(), curve, skin_arc );
        MB_CHK_SET_ERR(rval,"merge verts failed");
      }
      if(curve.back()!=skin_arc.back()) {
        rval = zip->merge_verts( curve.back(), skin_arc.back(), curve, skin_arc );
        MB_CHK_SET_ERR(rval,"merge verts failed");
      }

      // replace the faceted edge with a skin arc
      curve = skin_arc;
      if(debug) std::cout << "  curve " << gen->geom_id_by_handle(curve_set)
                            << " has been replaced by skin:" << std::endl;

    } else {
      // seal the pair together
      std::vector<moab::EntityHandle> sealed_curve;
      rval = seal_arc_pair( debug, FACET_TOL, normal_tag, curve, skin_arc,
                            gen->geom_id_by_handle(surf_set) );
      MB_CHK_SET_ERR(rval,"can't seal pair");
    }

    // get new front_endpt to guide selection of next curve
    front_endpt = curve.back();

    // To preserve the original rotation, reverse the curve if it has been reversed
    if(curve_is_reversed) reverse( curve.begin(), curve.end() );

    // set the sealed edge
    rval = arc->set_meshset( curve_set, curve );
    MB_CHK_SET_ERR(rval,"setting curve set failed");

    // Get skin_loop to cut an arc from it
    rval = arc->get_meshset( skin_loop_set, skin_loop );
    MB_CHK_SET_ERR(rval,"getting skin_loop set failed");


  }

  // The skin_loop_set is no longer needed.
  rval = MBI()->delete_entities( &skin_loop_set, 1 );
  MB_CHK_SET_ERR(rval,"deleting skin_loop_set failed");

  return moab::MB_SUCCESS;
}

// input: surface sets, ordered curve sets,
// output: skin arcs corresponding to curves are added to parent surface sets
moab::ErrorCode MakeWatertight::prepare_surfaces(moab::Range &surface_sets,
    moab::Tag geom_tag, moab::Tag id_tag, moab::Tag normal_tag, moab::Tag merge_tag,
    moab::Tag orig_curve_tag,
    const double SME_RESABS_TOL, const double FACET_TOL,
    const bool debug, bool verbose)
{

  moab::ErrorCode result;
  moab::Skinner tool(MBI());
  // loop over each surface meshset
  for(moab::Range::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) {

    // get the surf id of the surface meshset
    int surf_id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );
    MB_CHK_SET_ERR(result,"could not get id tag");
    assert(moab::MB_SUCCESS == result);
    if(debug) std::cout << "  surf id= " << surf_id << std::endl;

    // get the 2D entities in the surface set
    moab::Range dim2_ents;
    result = MBI()->get_entities_by_dimension( *i, 2, dim2_ents );
    MB_CHK_SET_ERR(result,"could not get 3D entities");
    assert(moab::MB_SUCCESS == result);

    // get facets of the surface meshset
    moab::Range tris;
    result = MBI()->get_entities_by_type( *i, moab::MBTRI, tris );
    MB_CHK_SET_ERR(result,"could not get tris");
    assert(moab::MB_SUCCESS == result);

    // Remove any 2D entities that are not triangles. This is needed because
    // ReadCGM will add quads and polygons to the surface set. This code only
    // works with triangles.
    moab::Range not_tris = subtract( dim2_ents, tris );
    if(!not_tris.empty()) {
      result = MBI()->delete_entities( not_tris );
      MB_CHK_SET_ERR(result,"could not delete not_tris");
      assert(moab::MB_SUCCESS == result);
      std::cout << "  removed " << not_tris.size()
                << " 2D elements that were not triangles from surface "
                << surf_id << std::endl;
    }

    // Get the curves and determine the number of unmerged curves
    std::vector<moab::EntityHandle> curve_sets, unmerged_curve_sets;
    result = get_unmerged_curves( *i , curve_sets, unmerged_curve_sets, merge_tag, verbose, debug);
    MB_CHK_SET_ERR(result, "could not get the curves and unmerged curves");

    // Save the normals of the facets. These will later be used to determine if
    // the tri became inverted.
    result = gen->save_normals( tris, normal_tag );
    MB_CHK_SET_ERR(result,"could not save_normals");
    assert(moab::MB_SUCCESS == result);

    // If all of the curves are merged, remove the surfaces facets.
    if(unmerged_curve_sets.empty()) {

      // if this surface is closed and contains triangles
      // then we will leave it alone but continue as we normally would
      //verify that the surface is closed
      moab::Range temp_skin_edges;
      if( tris.size() >= 4 && curve_sets.empty() ) {
        result = tool.find_skin(0, tris, 1, temp_skin_edges, false);
        MB_CHK_SET_ERR(result, "could not skin the triangles");
        // if the surface is closed, continue to next surface
        if(temp_skin_edges.empty()) {
          continue;
        }
      }

      result = gen->delete_surface( *i , geom_tag, tris, surf_id, debug, verbose);
      MB_CHK_SET_ERR(result, "could not delete surface");
      // adjust iterator so *i is still the same surface
      i = surface_sets.erase(i) - 1;
      continue;
    }

    // combine merged curve's surface senses
    result = gen->combine_merged_curve_senses( curve_sets, merge_tag, debug );
    MB_CHK_SET_ERR(result,"could not combine the merged curve sets");



    // Check if edges exist
    int n_edges;
    result = MBI()->get_number_entities_by_type(0, moab::MBEDGE, n_edges );
    MB_CHK_SET_ERR(result,"could not get number of edges");
    assert(moab::MB_SUCCESS == result);
    if(0!=n_edges) {
      MB_CHK_SET_ERR(moab::MB_FAILURE,"edges exist");
    }
    assert(0 == n_edges); //*** Why can't we have edges? (Also, this assertion is never used)

    // get the range of skin edges from the range of facets
    moab::Range skin_edges, skin_edges2;
    if(tris.empty()) continue; // nothing to zip
    // The MOAB skinner is not used here currently as it doesn't allow
    // make_watertight to close loops. The local version of find_skin is used instead.
    // This should be ok as the local find_skin fundtion should only be avoided when checking meshes for watertightness
    // to keep from altering the data set when checking.
    result = tool.find_skin( 0, tris, 1, skin_edges, false);
    MB_CHK_SET_ERR(result,"could not find_skin");
    assert(moab::MB_SUCCESS == result);


    // merge the vertices of the skin
    // BRANDON: For some reason cgm2moab does not do this? This was the
    // problem with mod13 surf 881. Two skin verts were coincident. A tol=1e-10
    // found the verts, but tol=0 did not.
    moab::Range skin_verts;
    bool cont = false;
    result = merge_skin_verts( skin_verts, skin_edges, SME_RESABS_TOL, surf_id, cont, debug);
    MB_CHK_SET_ERR(result,"could not merge the skin verts");
    if(cont) continue;


    // take skin edges and create loops of vertices
    std::vector < std::vector <moab::EntityHandle> > skin;
    cont = false;
    result = create_skin_vert_loops ( skin_edges, tris, skin , surf_id, cont, debug);
    MB_CHK_SET_ERR(result, " could not create skin loops of vertices");
    if(cont) continue;


// separate the remainder into a new function seal surface??

// separate this part from prepare surfaces into make_mesh_watertight??


    moab::EntityHandle skin_loop_sets[skin.size()];
    result = seal_surface_loops ( *i , skin_loop_sets , skin,  curve_sets, normal_tag, orig_curve_tag, FACET_TOL, surf_id, debug);
    MB_CHK_SET_ERR(result,"could not seal the surface loops");


    // Remove the sets of skin loops
    result = MBI()->delete_entities( &skin_loop_sets[0], skin.size() );
    MB_CHK_SET_ERR(result,"failed to zip: deleting skin_loop_sets failed");


  } // loop over each surface
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::fix_normals(moab::Range surface_sets,
    moab::Tag id_tag,
    moab::Tag normal_tag,
    const bool debug,
    const bool verbose)
{
  moab::ErrorCode result;
  if(debug) std::cout<< "number of surfaces=" << surface_sets.size() << std::endl;
  int inverted_tri_counter = 0;

  // loop over each surface meshset
  for(moab::Range::iterator i=surface_sets.begin(); i!=surface_sets.end(); i++ ) {

    // get the surf id of the surface meshset
    int surf_id;
    result = MBI()->tag_get_data( id_tag, &(*i), 1, &surf_id );
    assert(moab::MB_SUCCESS == result);
    if(debug) std::cout << "fix_normals surf id=" << surf_id << std::endl;

    // get facets from the surface meshset
    moab::Range tris;
    result = MBI()->get_entities_by_type( *i, moab::MBTRI, tris );
    assert(moab::MB_SUCCESS == result);

    // get the normals, pre zipping
    std::vector<moab::CartVect> old_normals(tris.size()), new_normals(tris.size());
    result = MBI()->tag_get_data( normal_tag, tris, &old_normals[0]);
    assert(moab::MB_SUCCESS == result);

    // get the normals, post zipping
    result = gen->triangle_normals( tris, new_normals );
    assert(moab::MB_SUCCESS == result);

    // test the normals, finding the inverted tris
    std::vector<int> inverted_tri_indices;
    result = zip->test_normals( old_normals, new_normals, inverted_tri_indices);
    assert(moab::MB_SUCCESS == result);

    // insert the inverted tris into a range
    moab::Range inverted_tris;
    for(unsigned int j=0; j<inverted_tri_indices.size(); j++) {
      inverted_tris.insert( tris[inverted_tri_indices[j]] );
      if(debug) gen->print_triangle( tris[inverted_tri_indices[j]], false );
    }

    // do edges exist?
    int n_edges;
    result = MBI()->get_number_entities_by_type( 0, moab::MBEDGE, n_edges );
    assert(moab::MB_SUCCESS == result);
    assert(0 == n_edges); // *** Why can't we have edges?

    // fix the inverted tris
    inverted_tri_counter += inverted_tris.size();
    result = zip->remove_inverted_tris(normal_tag, inverted_tris, debug );
    if(moab::MB_SUCCESS != result)
      std::cout << "  failed to fix inverted triangles in surface " << surf_id << std::endl;

    // if fix_normals exits on an error, we still need to remove its edges
    result = delete_all_edges();
    assert(moab::MB_SUCCESS == result);
  }
  if(verbose) {
    std::cout << "  Before fixing, " << inverted_tri_counter
              << " inverted triangles were found." << std::endl;
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::restore_moab_curve_representation( const moab::Range curve_sets )
{
  moab::ErrorCode result;
  for(moab::Range::const_iterator i=curve_sets.begin(); i!=curve_sets.end(); ++i) {
    // get the ordered verts
    std::vector<moab::EntityHandle> ordered_verts;
    result = arc->get_meshset( *i, ordered_verts );
    assert(moab::MB_SUCCESS==result);
    if(moab::MB_SUCCESS != result) return result;

    // Check for duplicate verts. This should not happen, but could if line
    // surfaces exist. This happens when feature size is violated and skin
    // from curves on the other side of the 1D line surf gets sealed.
    if( 1<ordered_verts.size() ) {
      for(std::vector<moab::EntityHandle>::iterator j=ordered_verts.begin()+1;
          j!=ordered_verts.end(); ++j) {
        if( *j == *(j-1) ) {
          std::cout << "duplicate vertex found in curve "
                    << gen->geom_id_by_handle(*i) << std::endl;
          j = ordered_verts.erase(j) - 1;
        }
      }
    }

    // Check for a point curve (should never happen),
    // a single degenerate edge (should never happen), or
    // a degenerate loop of two edges (should never happen)
    if( 4>ordered_verts.size() &&
        ordered_verts.front()==ordered_verts.back() ) {
      std::cout << "warning: curve " << gen->geom_id_by_handle(*i)
                << " is one degenerate edge" << std::endl;
      //return moab::MB_FAILURE; **** when this is uncommented, problems occur
    }

    // Determine if the curve is a loop or point curve. At least 4 verts are
    // needed to form a loop.
    bool is_loop = ( ordered_verts.front()==ordered_verts.back() &&
                     3<ordered_verts.size() );

    // get geometric endpoint sets of the curve (may be stale)
    moab::Range endpt_sets;
    result = MBI()->get_child_meshsets( *i, endpt_sets );
    assert(moab::MB_SUCCESS==result);
    if(moab::MB_SUCCESS != result) return result;

    // do the correct number of endpt sets exist?
    const unsigned int n_endpts = (is_loop) ? 1 : 2;
    if(n_endpts != endpt_sets.size()) {
      std::cout << "curve " << gen->geom_id_by_handle(*i) << " has " << n_endpts
                << " endpoints, but " << endpt_sets.size()
                << " endpoint sets exist" << std::endl;
    }

    // do they match the current endpoints?
    for(moab::Range::iterator j=endpt_sets.begin(); j!=endpt_sets.end(); ++j) {
      moab::Range endpt_vert;
      result = MBI()->get_entities_by_handle( *j, endpt_vert );
      if(moab::MB_SUCCESS != result) return result;
      assert(moab::MB_SUCCESS==result);
      if(1 != endpt_vert.size()) {
        std::cout << "curve " << gen->geom_id_by_handle(*i)
                  << " has" << endpt_vert.size()
                  << " endpoint vertices in the geometric vertex set" << std::endl;
        return moab::MB_INVALID_SIZE;
      }
      if(endpt_vert.front()!=ordered_verts.front() &&
         endpt_vert.front()!=ordered_verts.back() ) {
        std::cout << "curve " << gen->geom_id_by_handle(*i)
                  << " endpt sets do not match" << std::endl;
      }
    }

    // create the edges of the curve
    std::vector<moab::EntityHandle> ordered_edges(ordered_verts.size()-1);
    for(unsigned int j=0; j<ordered_verts.size()-1; ++j) {
      moab::EntityHandle edge;
      result = MBI()->create_element( moab::MBEDGE, &ordered_verts[j], 2, edge );
      assert(moab::MB_SUCCESS==result);
      if(moab::MB_SUCCESS != result) return result;
      ordered_edges[j] = edge;
    }

    // If the curve is a loop, remove the duplicate endpoint.
    if(is_loop) ordered_verts.pop_back();

    // clear the set
    result = MBI()->clear_meshset( &(*i), 1 );
    assert(moab::MB_SUCCESS==result);
    if(moab::MB_SUCCESS != result) return result;

    // add the verts then edges to the curve set
    result = MBI()->add_entities( *i, &ordered_verts[0], ordered_verts.size() );
    assert(moab::MB_SUCCESS==result);
    if(moab::MB_SUCCESS != result) return result;
    result = MBI()->add_entities( *i, &ordered_edges[0], ordered_edges.size() );
    assert(moab::MB_SUCCESS==result);
    if(moab::MB_SUCCESS != result) return result;
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::get_geom_size_before_sealing( const moab::Range geom_sets[],
    const moab::Tag geom_tag,
    const moab::Tag size_tag,
    bool debug,
    bool verbose )
{
  moab::ErrorCode rval;
  for( int dim = 1 ; dim < 4 ; ++dim ) {
    for(moab::Range::iterator i=geom_sets[dim].begin() ; i != geom_sets[dim].end() ; i++) {
      double size;
      //std::cout << "dim = " << dim << " *i =" << *i << std::endl;
      rval = gen->measure( *i, geom_tag, size, debug, verbose );
      //std::cout << " here in gen mesaure" << std::endl;
      MB_CHK_SET_ERR(rval,"could not measure");

      rval = MBI()->tag_set_data( size_tag, &(*i), 1, &size );
      //std::cout << " here in set tag data" << std::endl;
      MB_CHK_SET_ERR(rval,"could not set size tag");
    }
  }
  if (verbose) {
    std::cout << "finished in get_geom_size_before_sealing" << std::endl;
  }
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::get_geom_size_after_sealing( const moab::Range geom_sets[],
    const moab::Tag geom_tag,
    const moab::Tag size_tag,
    const double FACET_TOL,
    bool debug,
    bool verbose )
{

  // save the largest difference for each dimension
  struct size_data {
    double orig_size, new_size, diff, percent;
    int id;
  };
  size_data largest_diff[3];
  size_data largest_percent[3];
  size_data smallest_size[3];

  moab::ErrorCode rval;
  for(unsigned dim=1; dim<4; dim++) {
    largest_diff[dim-1].diff       = 0;
    largest_percent[dim-1].percent = 0;
    smallest_size[dim-1].orig_size = std::numeric_limits<int>::max();
    for(moab::Range::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
      double orig_size = 0, new_size = 0;
      rval = MBI()->tag_get_data( size_tag, &(*i), 1, &orig_size );
      if(moab::MB_SUCCESS != rval) {
        std::cout << "rval=" << rval << " id=" << gen->geom_id_by_handle(*i) << std::endl;
      }
      assert(moab::MB_SUCCESS == rval);
      rval = gen->measure( *i, geom_tag, new_size, debug, verbose );
      assert(moab::MB_SUCCESS == rval);

      // Remember the largest difference and associated percent difference
      double diff = fabs(new_size - orig_size);
      double percent_diff = 100.0*diff/orig_size;
      if(diff > largest_diff[dim-1].diff) {
        largest_diff[dim-1].orig_size = orig_size;
        largest_diff[dim-1].new_size  = new_size;
        largest_diff[dim-1].diff    = diff;
        largest_diff[dim-1].percent = percent_diff;
        largest_diff[dim-1].id      = gen->geom_id_by_handle(*i);
      }
      if(orig_size < smallest_size[dim-1].orig_size) {
        smallest_size[dim-1].orig_size = orig_size;
        smallest_size[dim-1].new_size  = new_size;
        smallest_size[dim-1].diff      = diff;
        smallest_size[dim-1].percent   = percent_diff;
        smallest_size[dim-1].id        = gen->geom_id_by_handle(*i);
      }
      if(percent_diff > largest_percent[dim-1].percent) {
        largest_percent[dim-1].orig_size = orig_size;
        largest_percent[dim-1].new_size  = new_size;
        largest_percent[dim-1].diff      = diff;
        largest_percent[dim-1].percent   = percent_diff;
        largest_percent[dim-1].id        = gen->geom_id_by_handle(*i);
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
        moab::Range curve_sets;
        rval = MBI()->get_child_meshsets( *i, curve_sets );
        double total_length = 0;
        for(moab::Range::iterator j=curve_sets.begin(); j!=curve_sets.end(); ++j) {
          double length;
          rval = MBI()->tag_get_data( size_tag, &(*j), 1, &length );
          total_length += length;
        }
        if(total_length*FACET_TOL < fabs(new_size - orig_size)) print_warning = true;
        // The volume cannot change more than the "cylinder" of error around each curve.
      } else if(3 == dim) {
        moab::Range surf_sets;
        rval = MBI()->get_child_meshsets( *i, surf_sets );
        double total_length = 0;
        for(moab::Range::iterator j=surf_sets.begin(); j!=surf_sets.end(); ++j) {
          moab::Range curve_sets;
          rval = MBI()->get_child_meshsets( *j, curve_sets );
          for(moab::Range::iterator k=curve_sets.begin(); k!=curve_sets.end(); ++k) {
            double length;
            rval = MBI()->tag_get_data( size_tag, &(*j), 1, &length );
            total_length += length;
          }
        }
        if(total_length*FACET_TOL*FACET_TOL*3.14 < fabs(new_size - orig_size)) {
          print_warning = true;
        }
      } else {
        return moab::MB_FAILURE;
      }
      if(print_warning && debug) {
        std::cout << "  dim=" << dim << " id=" << gen->geom_id_by_handle(*i)
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
  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::delete_merged_curves( moab::Range &existing_curve_sets,
    moab::Tag merge_tag,
    bool debug)
{

  moab::ErrorCode result;

  moab::Range curves_to_delete;
  result = MBI()->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &merge_tag, NULL,
           1, curves_to_delete);
  assert(moab::MB_SUCCESS == result);
  // loop over the curves to delete
  for(moab::Range::const_iterator i=curves_to_delete.begin(); i!=curves_to_delete.end(); ++i) {
    // get the curve_to_keep
    moab::EntityHandle curve_to_keep;
    result = MBI()->tag_get_data( merge_tag, &(*i), 1, &curve_to_keep );
    if(moab::MB_SUCCESS != result) return result;
    // get parent surface of the curve_to_delete
    moab::Range parent_surfs;
    result = MBI()->get_parent_meshsets( *i, parent_surfs );
    if(moab::MB_SUCCESS != result) return result;
    // remove the curve_to_delete and replace with curve_to_keep
    for(moab::Range::iterator j=parent_surfs.begin(); j!=parent_surfs.end(); ++j) {
      result = MBI()->remove_parent_child( *j, *i );
      if(moab::MB_SUCCESS != result) return result;
      result = MBI()->add_parent_child( *j, curve_to_keep );
      if(moab::MB_SUCCESS != result) return result;
    }
  }
  result = MBI()->delete_entities( curves_to_delete );
  assert(moab::MB_SUCCESS == result);
  if ( result != moab::MB_SUCCESS ) {
    std::cout << "Houston, we have a problem" << std::endl;
  }
  existing_curve_sets = subtract(existing_curve_sets, curves_to_delete );
  if(debug) std::cout << "deleted " << curves_to_delete.size() << " curves." << std::endl;


  return result;

}

moab::ErrorCode MakeWatertight::delete_sealing_tags( moab::Tag normal_tag,
    moab::Tag merge_tag,
    moab::Tag size_tag,
    moab::Tag orig_curve_tag)
{

  moab::ErrorCode result;

  result = MBI()->tag_delete( normal_tag );
  if(moab::MB_SUCCESS != result) return result;
  result = MBI()->tag_delete( merge_tag );
  if(moab::MB_SUCCESS != result) return result;
  result = MBI()->tag_delete( size_tag );
  if(moab::MB_SUCCESS != result) return result;
  result = MBI()->tag_delete( orig_curve_tag );
  if(moab::MB_SUCCESS != result) return result;

  return result;
}

moab::ErrorCode MakeWatertight::get_unmerged_curves( moab::EntityHandle surface, std::vector<moab::EntityHandle> &curves,
    std::vector<moab::EntityHandle> &unmerged_curves,
    moab::Tag merge_tag,
    bool verbose,
    bool debug)
{

  moab::ErrorCode result;

  result = MBI()->get_child_meshsets( surface, curves );
  MB_CHK_SET_ERR(result,"could not get child sets");
  assert(moab::MB_SUCCESS==result);

  // Update the curve_sets with that contain entity_to_delete curves with their
  // entity_to_keep curves. Unmerged_curve_sets will end up holding the curves
  // of this surface that are not merged with another curve in this surface.
  for(std::vector<moab::EntityHandle>::iterator j=curves.begin();
      j!=curves.end(); j++) {
    moab::EntityHandle merged_curve, curve;
    result = MBI()->tag_get_data( merge_tag, &(*j), 1, &merged_curve );
    assert(moab::MB_TAG_NOT_FOUND==result || moab::MB_SUCCESS==result);
    if(moab::MB_TAG_NOT_FOUND == result) {
      curve = *j;
    } else if(moab::MB_SUCCESS == result) {
      if(debug) {
        std::cout << "  curve_id=" << gen->geom_id_by_handle(*j)
                  << " is entity_to_delete" << std::endl;
      }
      curve = merged_curve;
      // should parent-childs be updated for the entity_to_keep?
    } else {
      std::cout << "prepare_surfaces: result=" << result << std::endl;
      return result;
    }

    // Add a curve (whether it is merged or not) if it is not in unmerged_merged_curve_sets.
    // If it is merged, then the curve handle will appear later and we will remove it.
    std::vector<moab::EntityHandle>::iterator k=find(unmerged_curves.begin(),
        unmerged_curves.end(), curve);
    if(unmerged_curves.end() == k) {

      unmerged_curves.push_back(curve);
    } else {

      unmerged_curves.erase(k);

    }
    // If all curves have a partner to be merged to, then unmerged_curve_sets will be empty.
  }


  return moab::MB_SUCCESS;

}


moab::ErrorCode MakeWatertight::create_skin_vert_loops( moab::Range &skin_edges,
    moab::Range tris,
    std::vector < std::vector <moab::EntityHandle> > &skin ,
    int surf_id,
    bool &cont,
    bool debug)
{

  moab::ErrorCode result;

  /* Remove pairs of edges that are geometrically the same but in opposite
     order due to a faceting error ACIS. In other words, merge (a,b) and (b,a).
     Examples are mod13surf280 and tbm_surf1605. */
  result = arc->remove_opposite_pairs_of_edges_fast( skin_edges, debug );
  if (moab::MB_SUCCESS != result) {
    std::cout << "  surface " << surf_id << " failed to zip: could not remove opposite edges"
              << surf_id << std::endl;
    result = MBI()->delete_entities(skin_edges);
    MB_CHK_SET_ERR(result,"could not delete skin edges");
    assert(moab::MB_SUCCESS == result);
    cont = true;
    return result;
  }




  /* Order the edges so that the triangle is on their left side. Skin
  edges are adjacent to only one triangle. */
  bool catch_error = false;
  for(moab::Range::iterator j=skin_edges.begin(); j!=skin_edges.end(); j++) {
    moab::Range adj_tris;
    result = MBI()->get_adjacencies( &(*j), 1, 2, false, adj_tris );
    MB_CHK_SET_ERR(result,"could not get adj tris");
    assert(moab::MB_SUCCESS == result);
    moab::Range skin_tri = intersect( adj_tris, tris );
    if(1 != skin_tri.size()) {
      std::cout << "skin_tri.size()=" << skin_tri.size() << std::endl;
      catch_error = true;
      break;
    }
    result = arc->orient_edge_with_tri( *j, skin_tri.front() );
    MB_CHK_SET_ERR(result,"could not orient_edge_with_tri");
    assert(moab::MB_SUCCESS == result);
  }
  // I NEED TO ADD BETTER CLEANUP AFTER THESE FAILURE CONDITIONS
  if(catch_error) {
    std::cout << "  surface " << surf_id << " failed to zip: could not orient edge"
              << std::endl;
    result = MBI()->delete_entities(skin_edges);
    MB_CHK_SET_ERR(result,"could not delete skin edges");
    assert(moab::MB_SUCCESS == result);
    cont = true;
    return result;
  }



  // Create loops with the skin edges.
  std::vector< std::vector<moab::EntityHandle> > skin_loops_of_edges;
  result = arc->create_loops_from_oriented_edges( skin_edges, skin_loops_of_edges, debug );
  if(moab::MB_SUCCESS != result) {
    std::cout << "  surface " << surf_id << " failed to zip: could not create loops"
              << std::endl;
    result = MBI()->delete_entities(skin_edges);
    assert(moab::MB_SUCCESS == result);
    cont = true;
    return result;
  }
  if(debug) std::cout << skin_loops_of_edges.size() << " skin loop(s)" << std::endl;

  // Convert the loops of skin edges to loops of skin verts.
  std::vector< std::vector<moab::EntityHandle> > skin_temp(skin_loops_of_edges.size());
  for(unsigned int j=0; j<skin_loops_of_edges.size(); j++) {
    result = gen->ordered_verts_from_ordered_edges( skin_loops_of_edges[j], skin_temp[j] );
    assert(moab::MB_SUCCESS == result);
    // check to make sure that the loop is closed
    assert(skin_temp[j].front() == skin_temp[j].back());
  }

  skin=skin_temp;

  // edges are no longer needed
  result = MBI()->delete_entities(skin_edges);
  MB_CHK_SET_ERR(result,"could not delete skin_edges");
  assert(moab::MB_SUCCESS == result);



  return moab::MB_SUCCESS;

}

moab::ErrorCode MakeWatertight::merge_skin_verts ( moab::Range &skin_verts,
    moab::Range &skin_edges,
    double SME_RESABS_TOL,
    int surf_id,
    bool cont,
    bool debug)
{

  moab::ErrorCode result;

  result = MBI()->get_adjacencies( skin_edges, 0, false, skin_verts,
                                   moab::Interface::UNION );
  MB_CHK_SET_ERR(result,"could not get adj verts");
  assert(moab::MB_SUCCESS == result);
  result = gen->merge_vertices( skin_verts, SME_RESABS_TOL );
  if (moab::MB_SUCCESS != result) {
    if(debug) std::cout << "result= " << result << std::endl;
    std::cout << "  surface " << surf_id << " failed to zip: could not merge vertices"
              << surf_id << std::endl;
    result = MBI()->delete_entities(skin_edges);
    assert(moab::MB_SUCCESS == result);
    cont = true;
    return result;
  }

  // Merging vertices create degenerate edges.
  result = arc->remove_degenerate_edges( skin_edges, debug );
  if(moab::MB_SUCCESS!=result) {
    std::cout << "  surface " << surf_id
              << " failed to zip: could not remove degenerate edges" << std::endl;
    result = MBI()->delete_entities(skin_edges);
    MB_CHK_SET_ERR(result,"could not delete skin edges");
    assert(moab::MB_SUCCESS == result);
    cont = true;
    return result;
  }


  return moab::MB_SUCCESS;
}

moab::ErrorCode MakeWatertight::seal_surface_loops ( moab::EntityHandle surf,
    moab::EntityHandle skin_loops[] ,
    std::vector < std::vector<moab::EntityHandle> > skin,
    std::vector<moab::EntityHandle> curves,
    moab::Tag normal_tag,
    moab::Tag orig_curve_tag,
    double FACET_TOL,
    int surf_id,
    bool debug)
{


  /* Get the curves that are part of the surface. Use vectors to store all curve
  stuff so that we can remove curves from the set as they are zipped. */

  // Place skin loops in meshsets to allow moab to update merged entities.
  // This prevents stale handles. I suspect moab can use upward adjacencies
  // to do this more efficiently than a manual O(n) search through an
  // unsorted vector.

  moab::ErrorCode rval;
  for(unsigned j=0; j<skin.size(); ++j) {
    rval = MBI()->create_meshset( moab::MESHSET_TRACK_OWNER|moab::MESHSET_ORDERED, skin_loops[j] );
    MB_CHK_SET_ERR(rval,"failed to zip: creating skin_loop_set failed");


    rval = arc->set_meshset( skin_loops[j], skin[j] );
    MB_CHK_SET_ERR(rval,"failed ot zip: setting skin_loop_set failed");

  }

  // Keep zipping loops until each is either zipped or failed. This function
  // returns only after all loops are zipped or a failure occurs.
  for(unsigned j=0; j<skin.size(); ++j) {
    std::vector<moab::EntityHandle> skin_loop;
    rval = arc->get_meshset( skin_loops[j], skin_loop );
    MB_CHK_SET_ERR(rval,"failed to zip: setting skin_loop_set failed");

    rval = seal_loop( debug, FACET_TOL, normal_tag, orig_curve_tag, surf,
                      curves, skin_loop );
    if(moab::MB_SUCCESS != rval) {
      std::cout << "failed to zip: surface " << surf_id << ": failed to seal a loop"
                << std::endl;
    }
  }

  return moab::MB_SUCCESS;

}

moab::ErrorCode MakeWatertight::make_mesh_watertight(moab::EntityHandle input_set, double &facet_tol, bool verbose)
{



  moab::ErrorCode result;

  //added to this function because they are called in make_watertight but aren't used until here
  const bool debug = false;
  const bool check_geom_size = true;
  // duplicated here from make_watertight, seemed to always be set to false
  bool is_acis=false;

  // initialize mesh tags
  moab::Tag geom_tag, id_tag, normal_tag, merge_tag, faceting_tol_tag,
       geometry_resabs_tag, size_tag, orig_curve_tag;
  double sme_resabs_tol=1e-6;

  // retrieve mesh tags necessary for sealing the mesh
  result = gen->get_sealing_mesh_tags( facet_tol, sme_resabs_tol, geom_tag, id_tag, normal_tag, merge_tag, faceting_tol_tag,
                                       geometry_resabs_tag, size_tag, orig_curve_tag);
  MB_CHK_SET_ERR(result, "could not get the mesh tags");


  // In practice, use 2*facet_tol because we are always comparing 2 faceted
  // entities. If instead we were comparing a faceted entity and a geometric
  // entitiy, then 1*facet_tol is correct.

  const double SME_RESABS_TOL = sme_resabs_tol; // from ACIS through CGM
  const double FACET_TOL = facet_tol; // specified by user when faceting cad

  if(verbose) {
    std::cout << "  faceting tolerance=" << facet_tol << " cm" << std::endl;
    std::cout << "  absolute tolerance=" << sme_resabs_tol << " cm" << std::endl;
  }

  // get all geometry sets and set tracking ordering options appropriately
  moab::Range geom_sets[4];
  result=gen->get_geometry_meshsets( geom_sets, geom_tag, verbose);
  MB_CHK_SET_ERR(result, "could not get the geometry meshsets");

  // If desired, find each entity's size before sealing.
  if(check_geom_size) {
    result = get_geom_size_before_sealing( geom_sets, geom_tag, size_tag, debug, verbose );
    MB_CHK_SET_ERR(result,"measuring geom size failed");
  }

  if (verbose) std::cout << "Getting entity count before sealing..." << std::endl;
  // Get entity count before sealing.
  int orig_n_tris;
  result = MBI()->get_number_entities_by_type( 0, moab::MBTRI, orig_n_tris );
  assert(moab::MB_SUCCESS == result);

  if(verbose) {
    std::cout << "  input faceted geometry contains " << geom_sets[3].size() << " volumes, "
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size()
              << " curves, and " << orig_n_tris << " triangles" << std::endl;
  }

  if(verbose) std::cout << "Finding degenerate triangles... " << std::endl;
  result = find_degenerate_tris();
  MB_CHK_SET_ERR(moab::MB_SUCCESS,"could not determine if triangles were degenerate or not");


  result = prepare_curves(geom_sets[1], geom_tag, id_tag, merge_tag, FACET_TOL, debug, verbose);
  MB_CHK_SET_ERR(moab::MB_SUCCESS,"could not prepare the curves");

  result = gen->check_for_geometry_sets(geom_tag, verbose);
  MB_CHK_SET_ERR(result,"no geometry sets exist in the model. Please check curve faceting.");


  if (verbose) {
    std::cout << "Zipping loops and removing small surfaces whose curves were all merged as pairs..." << std::endl;
    std::cout << "SME_RESABS_TOL = " << SME_RESABS_TOL << " FACET_TOL = " << FACET_TOL << std::endl;
  }
  result = prepare_surfaces(geom_sets[2], geom_tag, id_tag, normal_tag, merge_tag,
                            orig_curve_tag,SME_RESABS_TOL, FACET_TOL, debug);
  MB_CHK_SET_ERR( moab::MB_SUCCESS, "I have failed to zip");


  // After zipping surfaces, merged curve entity_to_deletes are no longer needed.
  // Swap surface parent-childs for curves to delete with curve to keep.
  // ARE THEIR ORPHANED CHILD VERTEX SETS STILL AROUND?
  if(verbose) std::cout << "Adjusting parent-child links then removing merged curves..." << std::endl;
  result = delete_merged_curves( geom_sets[1], merge_tag, debug);
  MB_CHK_SET_ERR(result, "could not delete the merged curves");


  // SHOULD COINCIDENT SURFACES ALSO BE MERGED?
  // 20100205 Jason: If child curves are the same, check to make sure each
  // facet pt is within 2x tol of the opposing surface.

  // This function is still screwed up, but 99% right.
  if (verbose) std::cout << "Fixing inverted triangles..." << std::endl;
  result = fix_normals(geom_sets[2], id_tag, normal_tag, debug, verbose);
  assert(moab::MB_SUCCESS == result);

  //clear out old geometry sets and repopulate (some may have been deleted)
  for(unsigned int i = 0; i < 4; i++) {
    geom_sets[i].clear();
  }
  result=gen->get_geometry_meshsets( geom_sets, geom_tag, verbose);
  MB_CHK_SET_ERR(result, "could not get the geometry meshsets");

  // As sanity check, did zipping drastically change the entity's size?
  if(check_geom_size && verbose) {
    std::cout << "Checking size change of zipped entities..." << std::endl;
    result = get_geom_size_after_sealing( geom_sets, geom_tag, size_tag, FACET_TOL, debug, verbose );
    MB_CHK_SET_ERR(result,"measuring geom size failed");
  }


  // This tool stores curves as a set of ordered edges. MOAB stores curves as
  // ordered vertices and edges. MOAB represents curve endpoints as geometry
  // sets containing a singe vertex. This function restores MOAB's curve
  // representation.
  if (verbose) std::cout << "Restoring faceted curve representation..." << std::endl;
  result = restore_moab_curve_representation( geom_sets[1] );
  MB_CHK_SET_ERR(result,"restore_moab_curve_representation failed");
  // If all of a volume's surfaces have been deleted, delete the volume.
  if (verbose) std::cout << "Removing small volumes if all surfaces have been removed..." << std::endl;
  for(moab::Range::iterator i=geom_sets[3].begin(); i!=geom_sets[3].end(); ++i) {
    int n_surfs;
    moab::EntityHandle vol = *i;
    result = MBI()->num_child_meshsets( vol, &n_surfs );
    assert(moab::MB_SUCCESS == result);
    if(0 == n_surfs) {
      // Remove the volume set. This also removes parent-child relationships.
      result = gen->delete_vol(vol);
      assert(moab::MB_SUCCESS == result);
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
  //assert(moab::MB_SUCCESS == result);

  //std::cout << "INSERT FUNCTION HERE TO REMOVE STALE VERTS, EDGES, TRIS, VERT SETS, ETC"
  //        << std::endl;

  // Resetting meshsets so that they no longer track.
  if (verbose) std::cout << "Restoring original meshset options and tags..." << std::endl;
  for(unsigned dim=0; dim<4; dim++) {
    for(moab::Range::iterator i=geom_sets[dim].begin(); i!=geom_sets[dim].end(); i++) {
      result = MBI()->set_meshset_options(*i, 1==dim ? moab::MESHSET_ORDERED : moab::MESHSET_SET);
      if(moab::MB_SUCCESS != result) return result;
    }
  }

  // Tags for merging curves and checking the change in geometry size were added.
  // Delete these because they are no longer needed.
  result = delete_sealing_tags( normal_tag, merge_tag, size_tag, orig_curve_tag);
  MB_CHK_SET_ERR(result, "could not delete sealing tags");

  // Print new size of the file to the user
  int sealed_n_tris;
  result = MBI()->get_number_entities_by_type( 0, moab::MBTRI, sealed_n_tris );
  assert(moab::MB_SUCCESS == result);
  if (verbose) {
    std::cout << "  output file contains " << geom_sets[3].size() << " volumes, "
              << geom_sets[2].size() << " surfaces, " << geom_sets[1].size()
              << " curves, and " << sealed_n_tris << " triangles" << std::endl;
    std::cout << "  triangle count changed " << (double)sealed_n_tris/orig_n_tris
              << "x (sealed/unsealed)" << std::endl;
  }
  return moab::MB_SUCCESS;
}



