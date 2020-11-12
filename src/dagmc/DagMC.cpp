#include "DagMC.hpp"
#include "util.hpp"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#ifndef M_PI /* windows */
#define M_PI 3.14159265358979323846
#endif

#ifdef DOUBLE_DOWN
#include "double_down/MOABRay.h"
#include "double_down/RTI.hpp"
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define FACETING_TOL_TAG_NAME "FACETING_TOL"
static const int null_delimiter_length = 1;

namespace moab {

/* Tolerance Summary

   Facet Tolerance:
   Maximum distance between continuous solid model surface and faceted surface.
     Performance:  increasing tolerance increased performance (fewer triangles)
     Robustness:   should not be affected
     Knowledge:    user must understand how coarser faceting influences accuracy
                   of results
*/

const bool counting = false; /* controls counts of ray casts and pt_in_vols */

// Empty synonym map for DagMC::parse_metadata()
const std::map<std::string, std::string> DagMC::no_synonyms;

// DagMC Constructor
DagMC::DagMC(std::shared_ptr<moab::Interface> mb_impl, double overlap_tolerance,
             double p_numerical_precision) {
#ifdef DOUBLE_DOWN
  std::cout << "Using the DOUBLE-DOWN interface to Embree." << std::endl;
#endif

  moab_instance_created = false;
  // if we arent handed a moab instance create one
  if (nullptr == mb_impl) {
    mb_impl = std::make_shared<Core>();
    moab_instance_created = true;
  }

  MBI_shared_ptr = mb_impl;
  // set the internal moab pointer
  MBI = MBI_shared_ptr.get();

  // make new GeomTopoTool and GeomQueryTool
  GTT = std::make_shared<GeomTopoTool>(MBI, false);
#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);
}

DagMC::DagMC(Interface* mb_impl, double overlap_tolerance,
             double p_numerical_precision) {
  moab_instance_created = false;
  // set the internal moab pointer
  MBI = mb_impl;
  MBI_shared_ptr = nullptr;

  // make new GeomTopoTool and GeomQueryTool
  GTT = std::make_shared<GeomTopoTool>(MBI, false);
#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);
}

// Destructor
DagMC::~DagMC() {
  // if we created the moab instance
  // clear it
  if (moab_instance_created) {
    MBI->delete_mesh();
  }
}

// get the float verision of dagmc version string
float DagMC::version(std::string* version_string) {
  if (NULL != version_string)
    *version_string =
        std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

/* SECTION I: Geometry Initialization and problem setup */

// the standard DAGMC load file method
ErrorCode DagMC::load_file(const char* cfile) {
  ErrorCode rval;
  std::string filename(cfile);
  std::cout << "Loading file " << cfile << std::endl;
  // load options
  char options[120] = {0};
  std::string file_ext = "";  // file extension

  // get the last 4 chars of file .i.e .h5m .sat etc
  int file_extension_size = 4;
  if (filename.size() > file_extension_size) {
    file_ext = filename.substr(filename.size() - file_extension_size);
  }
  EntityHandle file_set;
  rval = MBI->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != rval) return rval;

  rval = MBI->load_file(cfile, &file_set, options, NULL, 0, 0);

  if (MB_UNHANDLED_OPTION == rval) {
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end
    // in '.h5m'
    std::string filename(cfile);
    if (file_ext != ".h5m") {
      std::cerr << "DagMC warning: unhandled file loading options."
                << std::endl;
    }
  } else if (MB_SUCCESS != rval) {
    std::cerr << "DagMC Couldn't read file " << cfile << std::endl;
    std::string message;
    if (MB_SUCCESS == MBI->get_last_error(message) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;

    return rval;
  }

  return finish_loading();
}

// helper function to load the existing contents of a MOAB instance into DAGMC
ErrorCode DagMC::load_existing_contents() { return finish_loading(); }

// setup the implicit compliment
ErrorCode DagMC::setup_impl_compl() {
  // If it doesn't already exist, create implicit complement
  // Create data structures for implicit complement
  ErrorCode rval = GTT->setup_implicit_complement();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle."
              << std::endl;
    return rval;
  }
  return MB_SUCCESS;
}

// gets the entity sets tagged with geomtag 2 and 3
// surfaces and volumes respectively
ErrorCode DagMC::setup_geometry(Range& surfs, Range& vols) {
  ErrorCode rval;

  // get all surfaces
  rval = GTT->get_gsets_by_dimension(2, surfs);
  MB_CHK_SET_ERR(rval, "Could not get surfaces from GTT");

  // get all volumes
  rval = GTT->get_gsets_by_dimension(3, vols);
  MB_CHK_SET_ERR(rval, "Could not get volumes from GTT");

  return MB_SUCCESS;
}

// sets up the obb tree for the problem
ErrorCode DagMC::setup_obbs() {
  ErrorCode rval;

  // If we havent got an OBB Tree, build one.
  if (!GTT->have_obb_tree()) {
    std::cout << "Building acceleration data structures..." << std::endl;
#ifdef DOUBLE_DOWN
    rval = ray_tracer->init();
#else
    rval = GTT->construct_obb_trees();
#endif
    MB_CHK_SET_ERR(rval, "Failed to build obb trees");
  }
  return MB_SUCCESS;
}

// setups of the indices for the problem, builds a list of surface and volumes
// indices
ErrorCode DagMC::setup_indices() {
  Range surfs, vols;
  ErrorCode rval = setup_geometry(surfs, vols);

  // build the various index vectors used for efficiency
  rval = build_indices(surfs, vols);
  MB_CHK_SET_ERR(rval, "Failed to build surface/volume indices");
  return MB_SUCCESS;
}

bool DagMC::has_graveyard() {
  return get_graveyard_group() != 0;
}

EntityHandle DagMC::get_graveyard_group() {
  Range groups;
  ErrorCode rval = get_groups(groups);
  MB_CHK_SET_ERR_CONT(rval, "Failed to retrieve groups");

  EntityHandle graveyard_group = 0;

  // get the name of each group and check for the 'mat:graveyard' string
  for (auto group : groups) {
    std::string group_name;
    rval = get_group_name(group, group_name);
    MB_CHK_SET_ERR_CONT(rval, "Failed to get a group name");

    // convert name to lower case
    lowercase_str(group_name);

    // resize to match lengths
    group_name.resize(graveyard_name.size());

    // check for the graveyard string
    if (group_name == graveyard_name) {
      graveyard_group = group;
      break;
    }
  }

  return graveyard_group;
}

ErrorCode DagMC::remove_graveyard() {
  ErrorCode rval;

  EntityHandle graveyard_group = get_graveyard_group();

  if (graveyard_group == 0) { return MB_SUCCESS; }

  Range to_delete;

  to_delete.insert(graveyard_group);

  bool trees_exist = geom_tool()->have_obb_tree();

  // get the graveyard volume(s)
  Range graveyard_vols;
  rval = moab_instance()->get_entities_by_handle(graveyard_group, graveyard_vols);
  MB_CHK_SET_ERR(rval, "Failed to get the graveyard volume(s)");

  // get the implicit complement
  EntityHandle ic = 0;
  rval = geom_tool()->get_implicit_complement(ic);
  if (rval != MB_ENTITY_NOT_FOUND || rval != MB_SUCCESS) {
    MB_CHK_SET_ERR(rval, "Could not get the implicit complement");
  }

  // update the implicit complement tree if needed
  if (trees_exist && ic) {
    rval = geom_tool()->delete_obb_tree(ic, true);
    MB_CHK_SET_ERR(rval, "Failed to delete the implicit complement OBBTree/BVH");
  }


  for (auto vol : graveyard_vols) {
    if (trees_exist) {
      // will recursively delete the volume's surface trees as well
      rval = geom_tool()->delete_obb_tree(vol);
      MB_CHK_SET_ERR(rval, "Failed to delete the graveyard volume's tree");
    }

    Range surfs;
    rval = moab_instance()->get_child_meshsets(vol, surfs);
    MB_CHK_SET_ERR(rval, "Failed to get the surfaces of the graveyard volume");

    to_delete.merge(surfs);

    for (auto surf : surfs) {
      // add triangles, vertices to the list of entities to remove
      Range primitives;
      rval = moab_instance()->get_entities_by_handle(surf, primitives);
      MB_CHK_SET_ERR(rval, "Failed to get mesh entities of graveyard surface");
      //to_delete.merge(primitives);

      if (ic) {
        rval = moab_instance()->remove_child_meshset(ic, surf);
        MB_CHK_SET_ERR(rval, "Failed to remove graveyard surface relationship to implicit complement");
      }
    }
  }

  // delete accumulated entities
  rval = moab_instance()->delete_entities(to_delete);
  MB_CHK_SET_ERR(rval, "Failed to delete graveyard entities");

  if(trees_exist && ic) {
    rval = geom_tool()->construct_obb_tree(ic);
    MB_CHK_SET_ERR(rval, "Failed to re-create the implicit complement OBBTree/BVH");
  }

  rval = geom_tool()->find_geomsets();
  MB_CHK_SET_ERR(rval, "Failed to find geometry sets after removing the graveyard");

  return rval;
}

ErrorCode DagMC::create_graveyard(bool overwrite) {
  /* Methodology
    1) Determine the axis-aligned bounding box of the current model
    2) Create a new volume set for the graveyard
    3) Create a new group labeled as the graveyard and add volume to that group
    4) Create an inner surface from the bounding box
    5) Create an outer surface by extending the inner surface
    6) Add surfaces as children of the new graveyard volume and the implicit
       complement
    7) Set the surface senses w.r.t. the new graveyard volume and implicit
       complement
    8) Delete the implicit complement volume's out-of-date OBBTree/BVH
    9) Construct the new volume's OBBTree/BVH
    10) Construct the updated implicit complement OBBTree/BVH
  */
  ErrorCode rval;

  // if a graveyard already exists and we aren't overwriting it,
  // report an error
  if (has_graveyard() && !overwrite) {
    MB_CHK_SET_ERR(rval, "Graveyard already exists");
  }

  // currently relying on the BVH as looping over all vertices may
  // be prohibitive
  if (!geom_tool()->have_obb_tree()) {
    MB_CHK_SET_ERR(MB_FAILURE, "Graveyard creation attempted without BVH");
  }

  BBOX box;
  for(int i = 0; i < num_entities(3); i++) {
    moab::EntityHandle vol = this->entity_by_index(3, i+1);
    double vmin[3], vmax[3];
    rval = this->getobb(vol, vmin, vmax);
    MB_CHK_SET_ERR(rval, "Failed to get volume OBB");

    box.update(vmin);
    box.update(vmax);
  }

  if (!box.valid()) {
    MB_CHK_SET_ERR(rval, "Invalid model bounding box generated for graveyard volume");
  }

  // create a new volume meshset
  EntityHandle volume_set;
  rval = MBI->create_meshset(0, volume_set);
  MB_CHK_SET_ERR(rval, "Failed to create a graveyard volume set");

  // add volume set to the model
  rval = geom_tool()->add_geo_set(volume_set, 3);
  MB_CHK_SET_ERR(rval, "Failed to add the volume to the GeomTopoTool");

  // set the category tag
  std::string volume_str;
  volume_str.resize(CATEGORY_TAG_SIZE);
  volume_str = "Volume";
  rval = MBI->tag_set_data(category_tag(), &volume_set, 1, volume_str.c_str());
  MB_CHK_SET_ERR(rval, "Failed to set graveyard volume category");

  // create group set for the graveyard volume
  EntityHandle group_set;
  rval = MBI->create_meshset(0, group_set);
  MB_CHK_SET_ERR(rval, "Failed to create a new graveyard group set");

  rval = geom_tool()->add_geo_set(group_set, 4);
  MB_CHK_SET_ERR(rval, "Failed to add the graveyard group to the GeomTopoTool");

  // set the group category
  std::string group_str = "Group";
  group_str.resize(CATEGORY_TAG_SIZE);
  rval = MBI->tag_set_data(category_tag(), &group_set, 1, group_str.c_str());
  MB_CHK_SET_ERR(rval, "Failed to set the group category");

  // set the volume name tag data (material metadata)
  rval = MBI->tag_set_data(name_tag(), &group_set, 1, graveyard_name.c_str());
  MB_CHK_SET_ERR(rval, "Failed to set the graveyard name");

  // add the graveyard volume to this group
  rval = MBI->add_entities(group_set, &volume_set, 1);
  MB_CHK_SET_ERR(rval, "Failed to add the graveyard volume to the graveyard group");

  /// SURFACE CREATION ///

  // expand the box a bit
  for (int i = 0; i < 3; i++) {
    box.upper[i] += 10.0 * numerical_precision();
    box.lower[i] -= 10.0 * numerical_precision();
  }

  // tear down the implicit complement tree
  EntityHandle ic;
  rval = geom_tool()->get_implicit_complement(ic);
  MB_CHK_SET_ERR(rval, "Failed to get the implicit complement");

  EntityHandle inner_surface;
  rval = box_to_surf(box.lower, box.upper, inner_surface);

  // establish the volume-surface parent-child relationship with the inner surface
  rval = MBI->add_parent_child(volume_set, inner_surface);
  MB_CHK_SET_ERR(rval, "Failed to create the graveyard parent-child relationship");

  // establish the volume-surface parent-child relationship with the inner surface
  rval = MBI->add_parent_child(ic, inner_surface);
  MB_CHK_SET_ERR(rval, "Failed to create the graveyard parent-child relationship");

  // set the surface senses (all triangles have outward normals so this should
  // be REVERSE wrt the graveyard volume)
  EntityHandle inner_senses[2] = {ic, volume_set};
  rval = MBI->tag_set_data(sense_tag(), &inner_surface, 1, inner_senses);
  MB_CHK_SET_ERR(rval, "Failed to set graveyard surface senses");

  // expand the box a bit again for the outer surface
  for (int i = 0; i < 3; i++) {
    box.upper[i] += 10.0 * numerical_precision();
    box.lower[i] -= 10.0 * numerical_precision();
  }

  EntityHandle outer_surface;
  rval = box_to_surf(box.lower, box.upper, outer_surface);

  // establish the volume-surface parent-child relationship with tie outer surface
  rval = MBI->add_parent_child(volume_set, outer_surface);
  MB_CHK_SET_ERR(rval, "Failed to create the graveyard parent-child relationship");

  // establish the volume-surface parent-child relationship with tie outer surface
  rval = MBI->add_parent_child(ic, outer_surface);
  MB_CHK_SET_ERR(rval, "Failed to create the graveyard parent-child relationship");

  // set the surface senses (all triangles have outward normals so this should
  // be FORWARD wrt the graveyard volume)
  EntityHandle outer_senses[2] = {volume_set, ic};
  rval = MBI->tag_set_data(sense_tag(), &outer_surface, 1, outer_senses);
  MB_CHK_SET_ERR(rval, "Failed to set graveyard surface senses");

  // OBBTree/BVH update

  // delete the implicit complement tree (but not the surface trees)
  rval = geom_tool()->delete_obb_tree(ic, true);
  MB_CHK_SET_ERR(rval, "Failed to delete the implicit complement tree");

  // create BVH for both the new implicit complement and the new graveyard
  // volume
  rval = geom_tool()->construct_obb_tree(volume_set);
  MB_CHK_SET_ERR(rval, "Failed to build accel. data structure for the new graveyard volume");

  rval = geom_tool()->construct_obb_tree(ic);
  MB_CHK_SET_ERR(rval, "Failed to build accel. data structure for the new implicit complement");

  return rval;
}

ErrorCode DagMC::box_to_surf(double llc[3], double urc[3], EntityHandle& surface_set) {
  ErrorCode rval;

  //start with vertices
  std::vector<std::array<double,3>> vertex_coords;
  // vertex coordinates for the lower z face
  vertex_coords.push_back({llc[0], llc[1], llc[2]});
  vertex_coords.push_back({urc[0], llc[1], llc[2]});
  vertex_coords.push_back({urc[0], urc[1], llc[2]});
  vertex_coords.push_back({llc[0], urc[1], llc[2]});
  // vertex coordinate for the upper z face
  vertex_coords.push_back({llc[0], llc[1], urc[2]});
  vertex_coords.push_back({urc[0], llc[1], urc[2]});
  vertex_coords.push_back({urc[0], urc[1], urc[2]});
  vertex_coords.push_back({llc[0], urc[1], urc[2]});

  std::vector<moab::EntityHandle> box_verts;
  for(const auto& coords : vertex_coords) {
    EntityHandle new_vertex;
    rval = MBI->create_vertex(coords.data(), new_vertex);
    MB_CHK_SET_ERR(rval, "Failed to create graveyard vertex");
    box_verts.push_back(new_vertex);
  }

  // now we have 8 vertices to create triangles with
  std::vector<std::array<int, 3>> connectivity_indices;
  // lower z
  connectivity_indices.push_back({1, 0, 2});
  connectivity_indices.push_back({3, 2, 0});
  // upper z
  connectivity_indices.push_back({5, 6, 4});
  connectivity_indices.push_back({7, 4, 6});
  // lower x
  connectivity_indices.push_back({0, 1, 4});
  connectivity_indices.push_back({3, 4, 1});
  // upper x
  connectivity_indices.push_back({2, 3, 6});
  connectivity_indices.push_back({7, 6, 3});
  // lower y
  connectivity_indices.push_back({0, 4, 3});
  connectivity_indices.push_back({7, 3, 4});
  // upper y
  connectivity_indices.push_back({1, 2, 3});
  connectivity_indices.push_back({6, 3, 2});

  moab::Range new_tris;
  for(const auto& ind : connectivity_indices) {
    EntityHandle new_triangle;
    std::array<EntityHandle, 3> tri_conn = {box_verts[ind[0]], box_verts[ind[1]], box_verts[ind[2]]};
    rval = MBI->create_element(moab::MBTRI, tri_conn.data(), 3, new_triangle);
    MB_CHK_SET_ERR(rval, "Failed to create new graveyard triangle");
    new_tris.insert(new_triangle);
  }

  // create a surface set
  rval = MBI->create_meshset(0, surface_set);
  MB_CHK_SET_ERR(rval, "Failed to create a graveyard surface set");

  // add the triangles and vertices to the surface
  rval = MBI->add_entities(surface_set, new_tris);
  MB_CHK_SET_ERR(rval, "Failed to add triangles to the graveyard surface set");

  rval = MBI->add_entities(surface_set, box_verts.data(), box_verts.size());
  MB_CHK_SET_ERR(rval, "Failed to add vertices to the graveyard surface set");

  // tag the surface set with the appropriate info
  rval = geom_tool()->add_geo_set(surface_set, 2);
  MB_CHK_SET_ERR(rval, "Failed to add the surface to the GeomTopoTool");

  // set the category tag
  std::string surface_str;
  surface_str.resize(CATEGORY_TAG_SIZE);
  surface_str = "Surface";
  rval = MBI->tag_set_data(category_tag(), &surface_set, 1, surface_str.c_str());
  MB_CHK_SET_ERR(rval, "Failed to set graveyard volume category");

  return rval;
}

// initialise the obb tree
ErrorCode DagMC::init_OBBTree() {
  ErrorCode rval;

  // find all geometry sets
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "GeomTopoTool could not find the geometry sets");

  // implicit compliment
  rval = setup_impl_compl();
  MB_CHK_SET_ERR(rval, "Failed to setup the implicit compliment");

  // build obbs
  rval = setup_obbs();
  MB_CHK_SET_ERR(rval, "Failed to setup the OBBs");

  // setup indices
  rval = setup_indices();
  MB_CHK_SET_ERR(rval, "Failed to setup problem indices");

  return MB_SUCCESS;
}

// helper function to finish setting up required tags.
ErrorCode DagMC::finish_loading() {
  ErrorCode rval;

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE,
                    NULL, false);

  facetingTolTag =
      get_tag(FACETING_TOL_TAG_NAME, 1, MB_TAG_SPARSE, MB_TYPE_DOUBLE);

  // search for a tag that has the faceting tolerance
  Range tagged_sets;
  double facet_tol_tagvalue = 0;
  bool other_set_tagged = false, root_tagged = false;

  // get list of entity sets that are tagged with faceting tolerance
  // (possibly empty set)
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &facetingTolTag,
                                           NULL, 1, tagged_sets);
  // if NOT empty set
  if (MB_SUCCESS == rval && !tagged_sets.empty()) {
    rval = MBI->tag_get_data(facetingTolTag, &(*tagged_sets.begin()), 1,
                             &facet_tol_tagvalue);
    if (MB_SUCCESS != rval) return rval;
    other_set_tagged = true;
  } else if (MB_SUCCESS == rval) {
    // check to see if interface is tagged
    EntityHandle root = 0;
    rval = MBI->tag_get_data(facetingTolTag, &root, 1, &facet_tol_tagvalue);
    if (MB_SUCCESS == rval)
      root_tagged = true;
    else
      rval = MB_SUCCESS;
  }

  if ((root_tagged || other_set_tagged) && facet_tol_tagvalue > 0) {
    facetingTolerance = facet_tol_tagvalue;
  }

  // initialize ray_tracer
  std::cout << "Initializing the GeomQueryTool..." << std::endl;
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "Failed to find the geometry sets");

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;

  return MB_SUCCESS;
}

/* SECTION II: Fundamental Geometry Operations/Queries */

ErrorCode DagMC::ray_fire(const EntityHandle volume, const double point[3],
                          const double dir[3], EntityHandle& next_surf,
                          double& next_surf_dist, RayHistory* history,
                          double user_dist_limit, int ray_orientation,
                          OrientedBoxTreeTool::TrvStats* stats) {
  ErrorCode rval =
      ray_tracer->ray_fire(volume, point, dir, next_surf, next_surf_dist,
                           history, user_dist_limit, ray_orientation, stats);
  return rval;
}

ErrorCode DagMC::point_in_volume(const EntityHandle volume, const double xyz[3],
                                 int& result, const double* uvw,
                                 const RayHistory* history) {
  ErrorCode rval =
      ray_tracer->point_in_volume(volume, xyz, result, uvw, history);
  return rval;
}

ErrorCode DagMC::test_volume_boundary(const EntityHandle volume,
                                      const EntityHandle surface,
                                      const double xyz[3], const double uvw[3],
                                      int& result, const RayHistory* history) {
  ErrorCode rval = ray_tracer->test_volume_boundary(volume, surface, xyz, uvw,
                                                    result, history);
  return rval;
}

// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMC::point_in_volume_slow(EntityHandle volume, const double xyz[3],
                                      int& result) {
  ErrorCode rval = ray_tracer->point_in_volume_slow(volume, xyz, result);
  return rval;
}

// detemine distance to nearest surface
ErrorCode DagMC::closest_to_location(EntityHandle volume,
                                     const double coords[3], double& result,
                                     EntityHandle* surface) {
  ErrorCode rval =
      ray_tracer->closest_to_location(volume, coords, result, surface);
  return rval;
}

// calculate volume of polyhedron
ErrorCode DagMC::measure_volume(EntityHandle volume, double& result) {
  ErrorCode rval = ray_tracer->measure_volume(volume, result);
  return rval;
}

// sum area of elements in surface
ErrorCode DagMC::measure_area(EntityHandle surface, double& result) {
  ErrorCode rval = ray_tracer->measure_area(surface, result);
  return rval;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense(EntityHandle volume, int num_surfaces,
                               const EntityHandle* surfaces, int* senses_out) {
  ErrorCode rval =
      GTT->get_surface_senses(volume, num_surfaces, surfaces, senses_out);
  return rval;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense(EntityHandle volume, EntityHandle surface,
                               int& sense_out) {
  ErrorCode rval = GTT->get_sense(surface, volume, sense_out);
  return rval;
}

ErrorCode DagMC::get_angle(EntityHandle surf, const double in_pt[3],
                           double angle[3], const RayHistory* history) {
  ErrorCode rval = ray_tracer->get_normal(surf, in_pt, angle, history);
  return rval;
}

ErrorCode DagMC::next_vol(EntityHandle surface, EntityHandle old_volume,
                          EntityHandle& new_volume) {
  ErrorCode rval = GTT->next_vol(surface, old_volume, new_volume);
  return rval;
}

/* SECTION III: Indexing & Cross-referencing */

EntityHandle DagMC::entity_by_id(int dimension, int id) {
  return GTT->entity_by_id(dimension, id);
}

int DagMC::id_by_index(int dimension, int index) {
  EntityHandle h = entity_by_index(dimension, index);
  if (!h) return 0;

  int result = 0;
  MBI->tag_get_data(GTT->get_gid_tag(), &h, 1, &result);
  return result;
}

int DagMC::get_entity_id(EntityHandle this_ent) {
  return GTT->global_id(this_ent);
}

int DagMC::get_max_id(int dimension) {
    // determine what ID to give the volume
  int max_id = 0;
  for (int i = 0; i < num_entities(dimension); i++) {
    int id = id_by_index(dimension, i);
    max_id = std::max(id, max_id);
  }
  return max_id;
}

ErrorCode DagMC::build_indices(Range& surfs, Range& vols) {
  ErrorCode rval = MB_SUCCESS;

  if (surfs.size() == 0 || vols.size() == 0) {
    std::cout << "Volumes or Surfaces not found" << std::endl;
    return MB_ENTITY_NOT_FOUND;
  }
  setOffset = std::min(*surfs.begin(), *vols.begin());
  // surf/vol offsets are just first handles
  EntityHandle tmp_offset = std::max(surfs.back(), vols.back());

  // set size
  entIndices.resize(tmp_offset - setOffset + 1);

  // store surf/vol handles lists (surf/vol by index) and
  // index by handle lists
  surf_handles().resize(surfs.size() + 1);
  std::vector<EntityHandle>::iterator iter = surf_handles().begin();
  // MCNP wants a 1-based index but C++ has a 0-based index. So we need to set
  // the first value to 0 and then start at the next position in the vector
  // (iter++) thereafter.
  *(iter++) = 0;
  std::copy(surfs.begin(), surfs.end(), iter);
  int idx = 1;
  for (Range::iterator rit = surfs.begin(); rit != surfs.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  vol_handles().resize(vols.size() + 1);
  iter = vol_handles().begin();

  // MCNP wants a 1-based index but C++ has a 0-based index. So we need to set
  // the first value to 0 and then start at the next position in the vector
  // (iter++) thereafter.
  *(iter++) = 0;
  std::copy(vols.begin(), vols.end(), iter);
  idx = 1;
  for (Range::iterator rit = vols.begin(); rit != vols.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  // get group handles
  Range groups;
  rval = get_groups(groups);
  if (MB_SUCCESS != rval)
    return rval;

  group_handles().resize(groups.size() + 1);
  group_handles()[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles()[1]);

  return MB_SUCCESS;
}

ErrorCode DagMC::get_groups(Range& groups) {
  // get group handles
  Tag cat_tag = category_tag();
  std::string group_category;
  group_category.resize(CATEGORY_TAG_SIZE);
  group_category = "Group";
  const void* const group_val[] = {group_category.c_str()};
  ErrorCode rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &cat_tag,
                                                     group_val, 1, groups);
  MB_CHK_SET_ERR(rval, "Failed to retrieve groups from the MOAB instance");
  return rval;
}

Tag DagMC::category_tag() {
  return get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                 MB_TAG_SPARSE, MB_TYPE_OPAQUE);
}

/* SECTION IV: Handling DagMC settings */

double DagMC::overlap_thickness() {
  return ray_tracer->get_overlap_thickness();
}

double DagMC::numerical_precision() {
  return ray_tracer->get_numerical_precision();
}

void DagMC::set_overlap_thickness(double new_thickness) {
  ray_tracer->set_overlap_thickness(new_thickness);
}

void DagMC::set_numerical_precision(double new_precision) {
  ray_tracer->set_numerical_precision(new_precision);
}

ErrorCode DagMC::write_mesh(const char* ffile, const int flen) {
  ErrorCode rval;

  // write out a mesh file if requested
  if (ffile && 0 < flen) {
    rval = MBI->write_mesh(ffile);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return rval;
    }
  }

  return MB_SUCCESS;
}

/* SECTION V: Metadata handling */

ErrorCode DagMC::get_group_name(EntityHandle group_set, std::string& name) {
  ErrorCode rval;
  const void* v = NULL;
  int ignored;
  rval = MBI->tag_get_by_ptr(name_tag(), &group_set, 1, &v, &ignored);
  if (MB_SUCCESS != rval) return rval;
  name = static_cast<const char*>(v);
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_group_name(EntityHandle group_set, prop_map& result,
                                  const char* delimiters) {
  ErrorCode rval;
  std::string group_name;
  rval = get_group_name(group_set, group_name);
  if (rval != MB_SUCCESS) return rval;

  std::vector<std::string> group_tokens;
  tokenize(group_name, group_tokens, delimiters);

  // iterate over all the keyword positions
  // keywords are even indices, their values (optional) are odd indices
  for (unsigned int i = 0; i < group_tokens.size(); i += 2) {
    std::string groupkey = group_tokens[i];
    std::string groupval;
    if (i < group_tokens.size() - 1) groupval = group_tokens[i + 1];
    result[groupkey] = groupval;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::detect_available_props(std::vector<std::string>& keywords_list,
                                        const char* delimiters) {
  ErrorCode rval;
  std::set<std::string> keywords;
  for (std::vector<EntityHandle>::const_iterator grp = group_handles().begin();
       grp != group_handles().end(); ++grp) {
    std::map<std::string, std::string> properties;
    rval = parse_group_name(*grp, properties, delimiters);
    if (rval == MB_TAG_NOT_FOUND)
      continue;
    else if (rval != MB_SUCCESS)
      return rval;

    for (prop_map::iterator i = properties.begin(); i != properties.end();
         ++i) {
      keywords.insert((*i).first);
    }
  }
  keywords_list.assign(keywords.begin(), keywords.end());
  return MB_SUCCESS;
}

ErrorCode DagMC::append_packed_string(Tag tag, EntityHandle eh,
                                      std::string& new_string) {
  // When properties have multiple values, the values are tagged in a single
  // character array with the different values separated by null characters
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = MBI->tag_get_by_ptr(tag, &eh, 1, &p, &len);
  if (rval == MB_TAG_NOT_FOUND) {
    // This is the first entry, and can be set directly
    p = new_string.c_str();
    return MBI->tag_clear_data(tag, &eh, 1, p, new_string.length() + 1);
  } else if (rval != MB_SUCCESS)
    return rval;
  else {
    str = static_cast<const char*>(p);
  }

  // append a new value for the property to the existing property string
  unsigned int tail_len = new_string.length() + null_delimiter_length;
  int new_len = tail_len + len;

  char* new_packed_string = new char[new_len];
  memcpy(new_packed_string, str, len);
  memcpy(new_packed_string + len, new_string.c_str(), tail_len);

  p = new_packed_string;
  rval = MBI->tag_set_by_ptr(tag, &eh, 1, &p, &new_len);
  delete[] new_packed_string;
  return rval;
}

ErrorCode DagMC::unpack_packed_string(Tag tag, EntityHandle eh,
                                      std::vector<std::string>& values) {
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = MBI->tag_get_by_ptr(tag, &eh, 1, &p, &len);
  if (rval != MB_SUCCESS) return rval;
  str = static_cast<const char*>(p);
  int idx = 0;
  while (idx < len) {
    std::string item(str + idx);
    values.push_back(item);
    idx += item.length() + null_delimiter_length;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_properties(
    const std::vector<std::string>& keywords,
    const std::map<std::string, std::string>& keyword_synonyms,
    const char* delimiters) {
  ErrorCode rval;

  // master keyword map, mapping user-set words in cubit to canonical property
  // names
  std::map<std::string, std::string> keyword_map(keyword_synonyms);

  for (std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i) {
    keyword_map[*i] = *i;
  }

  // the set of all canonical property names
  std::set<std::string> prop_names;
  for (prop_map::iterator i = keyword_map.begin(); i != keyword_map.end();
       ++i) {
    prop_names.insert((*i).second);
  }

  // set up DagMC's property tags based on what's been requested
  for (std::set<std::string>::iterator i = prop_names.begin();
       i != prop_names.end(); ++i) {
    std::string tagname("DAGMCPROP_");
    tagname += (*i);

    Tag new_tag;
    rval = MBI->tag_get_handle(tagname.c_str(), 0, MB_TYPE_OPAQUE, new_tag,
                               MB_TAG_SPARSE | MB_TAG_VARLEN | MB_TAG_CREAT);
    if (MB_SUCCESS != rval) return rval;
    property_tagmap[(*i)] = new_tag;
  }

  // now that the keywords and tags are ready, iterate over all the actual
  // geometry groups
  for (std::vector<EntityHandle>::iterator grp = group_handles().begin();
       grp != group_handles().end(); ++grp) {
    prop_map properties;
    rval = parse_group_name(*grp, properties, delimiters);
    if (rval == MB_TAG_NOT_FOUND)
      continue;
    else if (rval != MB_SUCCESS)
      return rval;

    Range grp_sets;
    rval = MBI->get_entities_by_type(*grp, MBENTITYSET, grp_sets);
    if (MB_SUCCESS != rval) return rval;
    if (grp_sets.size() == 0) continue;

    for (prop_map::iterator i = properties.begin(); i != properties.end();
         ++i) {
      std::string groupkey = (*i).first;
      std::string groupval = (*i).second;

      if (property_tagmap.find(groupkey) != property_tagmap.end()) {
        Tag proptag = property_tagmap[groupkey];
        const unsigned int groupsize = grp_sets.size();
        for (unsigned int j = 0; j < groupsize; ++j) {
          rval = append_packed_string(proptag, grp_sets[j], groupval);
          if (MB_SUCCESS != rval) return rval;
        }
      }
    }
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_value(EntityHandle eh, const std::string& prop,
                            std::string& value) {
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr(proptag, &eh, 1, &data, &ignored);
  if (rval != MB_SUCCESS) return rval;
  value = static_cast<const char*>(data);
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_values(EntityHandle eh, const std::string& prop,
                             std::vector<std::string>& values) {
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  return unpack_packed_string(proptag, eh, values);
}

bool DagMC::has_prop(EntityHandle eh, const std::string& prop) {
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return false;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr(proptag, &eh, 1, &data, &ignored);
  return (rval == MB_SUCCESS);
}

ErrorCode DagMC::get_all_prop_values(const std::string& prop,
                                     std::vector<std::string>& return_list) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &proptag, NULL, 1,
                                           all_ents);
  if (MB_SUCCESS != rval) return rval;

  std::set<std::string> unique_values;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (MB_SUCCESS != rval) return rval;
    unique_values.insert(values.begin(), values.end());
  }

  return_list.assign(unique_values.begin(), unique_values.end());
  return MB_SUCCESS;
}

ErrorCode DagMC::entities_by_property(const std::string& prop,
                                      std::vector<EntityHandle>& return_list,
                                      int dimension, const std::string* value) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  // Note that we cannot specify values for proptag here-- the passed value,
  // if it exists, may be only a subset of the packed string representation
  // of this tag.
  Tag tags[2] = {proptag, GTT->get_geom_tag()};
  void* vals[2] = {NULL, (dimension != 0) ? &dimension : NULL};
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, tags, vals, 2,
                                           all_ents);
  if (MB_SUCCESS != rval) return rval;

  std::set<EntityHandle> handles;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (MB_SUCCESS != rval) return rval;
    if (value) {
      if (std::find(values.begin(), values.end(), *value) != values.end()) {
        handles.insert(*i);
      }
    } else {
      handles.insert(*i);
    }
  }

  return_list.assign(handles.begin(), handles.end());
  return MB_SUCCESS;
}

bool DagMC::is_implicit_complement(EntityHandle volume) {
  return GTT->is_implicit_complement(volume);
}

void DagMC::tokenize(const std::string& str, std::vector<std::string>& tokens,
                     const char* delimiters) const {
  std::string::size_type last = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos = str.find_first_of(delimiters, last);
  if (std::string::npos == pos)
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back(str.substr(last, pos - last));
      last = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, last);
      if (std::string::npos == pos) pos = str.size();
    }
}

Tag DagMC::get_tag(const char* name, int size, TagType store, DataType type,
                   const void* def_value, bool create_if_missing) {
  Tag retval = 0;
  unsigned flags = store | MB_TAG_CREAT;
  // NOTE: this function seems to be broken in that create_if_missing has
  // the opposite meaning from what its name implies.  However, changing the
  // behavior causes tests to fail, so I'm leaving the existing behavior
  // in place.  -- j.kraftcheck.
  if (!create_if_missing) flags |= MB_TAG_EXCL;
  ErrorCode result =
      MBI->tag_get_handle(name, size, type, retval, flags, def_value);
  if (create_if_missing && MB_SUCCESS != result)
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;

  return retval;
}

/* SECTION VI: Other */

ErrorCode DagMC::getobb(EntityHandle volume, double minPt[3], double maxPt[3]) {
#ifdef DOUBLE_DOWN
  ErrorCode rval = ray_tracer->get_bbox(volume, minPt, maxPt);
#else
  ErrorCode rval = GTT->get_bounding_coords(volume, minPt, maxPt);
#endif
  MB_CHK_SET_ERR(rval, "Failed to get obb for volume");
  return MB_SUCCESS;
}

ErrorCode DagMC::getobb(EntityHandle volume, double center[3], double axis1[3],
                        double axis2[3], double axis3[3]) {
#ifdef DOUBLE_DOWN
  ErrorCode rval = ray_tracer->get_obb(volume, center, axis1, axis2, axis3);
#else
  ErrorCode rval = GTT->get_obb(volume, center, axis1, axis2, axis3);
#endif
  MB_CHK_SET_ERR(rval, "Failed to get obb for volume");
  return MB_SUCCESS;
}

}  // namespace moab
