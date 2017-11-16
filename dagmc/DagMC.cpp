#include "DagMC.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <set>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#ifndef M_PI  /* windows */
# define M_PI 3.14159265358979323846
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define FACETING_TOL_TAG_NAME "FACETING_TOL"

#define SDF_PRECONDITIONER

#ifdef SDF_PRECONDITIONER
  #define SDF_BUILD
//  #define SDF_WRITE
  #define SDF_RF
//  #define SDF_DEBUG
  #define SDF_PIV
  #define SDF_CTL
#endif

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
DagMC::DagMC(Interface* mb_impl, double overlap_tolerance, double p_numerical_precision) {
  moab_instance_created = false;
  // if we arent handed a moab instance create one
  if (NULL == mb_impl) {
    mb_impl = new moab::Core();
    moab_instance_created = true;
  }

  // set the internal moab pointer
  MBI = mb_impl;

  // make new GeomTopoTool and GeomQueryTool
  GTT = new moab::GeomTopoTool(MBI, false);
  GQT = new moab::GeomQueryTool(GTT, overlap_tolerance, p_numerical_precision);

  // This is the correct place to uniquely define default values for the dagmc settings
  defaultFacetingTolerance = .001;
}

// Destructor
DagMC::~DagMC() {
  // delete the GeomTopoTool and GeomQueryTool
  delete GTT;
  delete GQT;

  // if we created the moab instance
  // clear it
  if (moab_instance_created) {
    MBI->delete_mesh();
    delete MBI;
  }
}

// get the float verision of dagmc version string
float DagMC::version(std::string* version_string) {
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

unsigned int DagMC::interface_revision() {
  unsigned int result = 0;
  const char* interface_string = DAGMC_INTERFACE_REVISION;
  if (strlen(interface_string) >= 5) {
    // start looking for the revision number after "$Rev: "
    result = strtol(interface_string + 5, NULL, 10);
  }
  return result;
}

/* SECTION I: Geometry Initialization and problem setup */

// the standard DAGMC load file method
ErrorCode DagMC::load_file(const char* cfile) {
  ErrorCode rval;
  std::cout << "Loading file " << cfile << std::endl;
  // load options
  char options[120] = {0};
  char file_ext[4] = "" ; // file extension

  // get the last 4 chars of file .i.e .h5m .sat etc
  memcpy(file_ext, &cfile[strlen(cfile) - 4], 4);

  EntityHandle file_set;
  rval = MBI->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != rval)
    return rval;

  rval = MBI->load_file(cfile, &file_set, options, NULL, 0, 0);

  if (MB_UNHANDLED_OPTION == rval) {
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end in '.h5m'
    std::string filename(cfile);
    if (filename.length() < 4 || filename.substr(filename.length() - 4) != ".h5m") {
      std::cerr << "DagMC warning: unhandled file loading options." << std::endl;
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
ErrorCode DagMC::load_existing_contents() {
  return finish_loading();
}

// setup the implicit compliment
ErrorCode DagMC::setup_impl_compl() {
  // If it doesn't already exist, create implicit complement
  // Create data structures for implicit complement
  ErrorCode rval = GTT->setup_implicit_complement();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle." << std::endl;
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
    std::cout << "Building OBB Tree..." << std::endl;
    rval = GTT->construct_obb_trees();
    MB_CHK_SET_ERR(rval, "Failed to build obb trees");
  }
  return MB_SUCCESS;
}

// setups of the indices for the problem, builds a list of
ErrorCode DagMC::setup_indices() {
  Range surfs, vols;
  ErrorCode rval = setup_geometry(surfs, vols);

  // build the various index vectors used for efficiency
  rval = build_indices(surfs, vols);
  MB_CHK_SET_ERR(rval, "Failed to build surface/volume indices");
  return MB_SUCCESS;
}

// initialise the obb tree
ErrorCode DagMC::init_OBBTree() {
  ErrorCode rval;

  // find all geometry sets
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "GeomTopoTool could not find the geometry sets");

  // implicit compliment
  // EntityHandle implicit_complement;
  //  rval = GTT->get_implicit_complement(implicit_complement, true);
  rval = setup_impl_compl();
  MB_CHK_SET_ERR(rval, "Failed to setup the implicit compliment");

  // build obbs
  rval = setup_obbs();
  MB_CHK_SET_ERR(rval, "Failed to setup the OBBs");

  // setup indices
  rval = setup_indices();
  MB_CHK_SET_ERR(rval, "Failed to setup problem indices");

#ifdef SDF_BUILD
  rval = build_preconditioner();
  MB_CHK_SET_ERR(rval, "Failed to setup preconditioner datastructure.");
#endif
  
  return MB_SUCCESS;
}

// helper function to finish setting up required tags.
ErrorCode DagMC::finish_loading() {
  ErrorCode rval;

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, NULL, false);

  facetingTolTag = get_tag(FACETING_TOL_TAG_NAME, 1, MB_TAG_SPARSE, MB_TYPE_DOUBLE);

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
    rval = MBI->tag_get_data(facetingTolTag, &(*tagged_sets.begin()), 1, &facet_tol_tagvalue);
    if (MB_SUCCESS != rval)
      return rval;
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

  // initialize GQT
  std::cout << "Initializing the GeomQueryTool..." << std::endl;
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "Failed to find the geometry sets");

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;

  return MB_SUCCESS;
}


/* SECTION II: Fundamental Geometry Operations/Queries */

ErrorCode DagMC::ray_fire(const EntityHandle volume, const double point[3],
                          const double dir[3], EntityHandle& next_surf,
                          double& next_surf_dist,
                          RayHistory* history,
                          double user_dist_limit, int ray_orientation,
                          OrientedBoxTreeTool::TrvStats* stats) {
  ErrorCode rval;
  
#ifdef SDF_RF
  bool precond_success;
  rval = precondition_ray_fire(volume, point, dir, user_dist_limit, next_surf, next_surf_dist, precond_success);
  MB_CHK_SET_ERR(rval, "Failed to precondition ray");
  if (precond_success) { return rval; }
  else { user_dist_limit = 0.0; }
#endif
  
  rval = GQT->ray_fire(volume, point, dir, next_surf, next_surf_dist,
                                 history, user_dist_limit, ray_orientation,
                                 stats);
  return rval;
}

ErrorCode DagMC::point_in_volume(const EntityHandle volume, const double xyz[3],
                                 int& result, const double* uvw,
                                 const RayHistory* history) {
  ErrorCode rval;
  
#ifdef SDF_PIV
  if(get_signed_distance_field(volume)) { 
    bool precond_success;
    rval = precondition_point_in_volume(volume, xyz, result, precond_success);
    MB_CHK_SET_ERR(rval, "Failed to precondition point in volume call for volume " << get_entity_id(volume));
    if (precond_success) return rval;
  }
#endif
  
  rval = GQT->point_in_volume(volume, xyz, result, uvw, history);
  return rval;
}

ErrorCode DagMC::test_volume_boundary(const EntityHandle volume,
                                      const EntityHandle surface,
                                      const double xyz[3], const double uvw[3],
                                      int& result,
                                      const RayHistory* history) {
  ErrorCode rval = GQT->test_volume_boundary(volume, surface, xyz, uvw, result,
                                             history);
  return rval;
}

// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMC::point_in_volume_slow(EntityHandle volume, const double xyz[3],
                                      int& result) {
  ErrorCode rval = GQT->point_in_volume_slow(volume, xyz, result);
  return rval;
}

// detemine distance to nearest surface
ErrorCode DagMC::closest_to_location(EntityHandle volume,
                                     const double coords[3], double& result,
                                     EntityHandle* surface) {
  ErrorCode rval;
  
#ifdef SDF_CTL
  bool precond_success;
  rval = precondition_closest_to_location(volume, coords, result, precond_success);
  MB_CHK_SET_ERR(rval, "Failed in preconditioning closest to location call");
  if (precond_success) return rval;
#endif
  
  rval = GQT->closest_to_location(volume, coords, result, surface);
  return rval;
}

// calculate volume of polyhedron
ErrorCode DagMC::measure_volume(EntityHandle volume, double& result) {
  ErrorCode rval = GQT->measure_volume(volume, result);
  return rval;
}

// sum area of elements in surface
ErrorCode DagMC::measure_area(EntityHandle surface, double& result) {
  ErrorCode rval = GQT->measure_area(surface, result);
  return rval;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense(EntityHandle volume, int num_surfaces,
                               const EntityHandle* surfaces, int* senses_out) {
  ErrorCode rval = GTT->get_surface_senses(volume, num_surfaces, surfaces,
                                           senses_out);
  return rval;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense(EntityHandle volume, EntityHandle surface,
                               int& sense_out) {
  ErrorCode rval = GTT->get_sense(surface, volume, sense_out);
  return rval;
}

ErrorCode DagMC::get_angle(EntityHandle surf, const double in_pt[3],
                           double angle[3],
                           const RayHistory* history) {
  ErrorCode rval = GQT->get_normal(surf, in_pt, angle, history);
  return rval;
}

ErrorCode DagMC::next_vol(EntityHandle surface, EntityHandle old_volume,
                          EntityHandle& new_volume) {
  ErrorCode rval = GTT->next_vol(surface, old_volume, new_volume);
  return rval;
}

/* SECTION III */

EntityHandle DagMC::entity_by_id(int dimension, int id) {
  return GTT->entity_by_id(dimension, id);
}

int DagMC::id_by_index(int dimension, int index) {
  EntityHandle h = entity_by_index(dimension, index);
  if (!h)
    return 0;

  int result = 0;
  MBI->tag_get_data(GTT->get_gid_tag(), &h, 1, &result);
  return result;
}

int DagMC::get_entity_id(EntityHandle this_ent) {
  return GTT->global_id(this_ent);
}

ErrorCode DagMC::build_indices(Range& surfs, Range& vols) {
  ErrorCode rval = MB_SUCCESS;

  // surf/vol offsets are just first handles
  setOffset = std::min(*surfs.begin(), *vols.begin());

  // max
  EntityHandle tmp_offset = std::max(surfs.back(), vols.back());

  // set size
  entIndices.resize(tmp_offset - setOffset + 1);

  // store surf/vol handles lists (surf/vol by index) and
  // index by handle lists
  surf_handles().resize(surfs.size() + 1);
  std::vector<EntityHandle>::iterator iter = surf_handles().begin();
  *(iter++) = 0;
  std::copy(surfs.begin(), surfs.end(), iter);
  int idx = 1;
  for (Range::iterator rit = surfs.begin(); rit != surfs.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  vol_handles().resize(vols.size() + 1);
  iter = vol_handles().begin();
  *(iter++) = 0;
  std::copy(vols.begin(), vols.end(), iter);
  idx = 1;
  for (Range::iterator rit = vols.begin(); rit != vols.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  // get group handles
  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                             MB_TAG_SPARSE, MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category + CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};
  Range groups;
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &category_tag,
                                           group_val, 1, groups);
  if (MB_SUCCESS != rval)
    return rval;
  group_handles().resize(groups.size() + 1);
  group_handles()[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles()[1]);

  return MB_SUCCESS;
}



/* SECTION IV */

void DagMC::set_overlap_thickness(double new_thickness) {
  GQT->set_overlap_thickness(new_thickness);
}

void DagMC::set_numerical_precision(double new_precision) {
  GQT->set_numerical_precision(new_precision);
}

ErrorCode DagMC::write_mesh(const char* ffile,
                            const int flen) {
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
  if (MB_SUCCESS != rval)
    return rval;
  name = static_cast<const char*>(v);
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_group_name(EntityHandle group_set, prop_map& result, const char* delimiters) {
  ErrorCode rval;
  std::string group_name;
  rval = get_group_name(group_set, group_name);
  if (rval != MB_SUCCESS)
    return rval;

  std::vector< std::string > group_tokens;
  tokenize(group_name, group_tokens, delimiters);

  // iterate over all the keyword positions
  // keywords are even indices, their values (optional) are odd indices
  for (unsigned int i = 0; i < group_tokens.size(); i += 2) {
    std::string groupkey = group_tokens[i];
    std::string groupval;
    if (i < group_tokens.size() - 1)
      groupval = group_tokens[i + 1];
    result[groupkey] = groupval;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::detect_available_props(std::vector<std::string>& keywords_list, const char* delimiters) {
  ErrorCode rval;
  std::set< std::string > keywords;
  for (std::vector<EntityHandle>::const_iterator grp = group_handles().begin();
       grp != group_handles().end(); ++grp) {
    std::map< std::string, std::string > properties;
    rval = parse_group_name(*grp, properties, delimiters);
    if (rval == MB_TAG_NOT_FOUND)
      continue;
    else if (rval != MB_SUCCESS)
      return rval;

    for (prop_map::iterator i = properties.begin();
         i != properties.end(); ++i) {
      keywords.insert((*i).first);
    }
  }
  keywords_list.assign(keywords.begin(), keywords.end());
  return MB_SUCCESS;
}

ErrorCode DagMC::append_packed_string(Tag tag, EntityHandle eh,
                                      std::string& new_string) {
  // When properties have multiple values, the values are tagged in a single character array
  // with the different values separated by null characters
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
  unsigned int tail_len = new_string.length() + 1;
  char* new_packed_string = new char[ len + tail_len ];
  memcpy(new_packed_string, str, len);
  memcpy(new_packed_string + len, new_string.c_str(), tail_len);

  int new_len = len + tail_len;
  p = new_packed_string;
  rval = MBI->tag_set_by_ptr(tag, &eh, 1, &p, &new_len);
  delete[] new_packed_string;
  return rval;
}

ErrorCode DagMC::unpack_packed_string(Tag tag, EntityHandle eh,
                                      std::vector< std::string >& values) {
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = MBI->tag_get_by_ptr(tag, &eh, 1, &p, &len);
  if (rval != MB_SUCCESS)
    return rval;
  str = static_cast<const char*>(p);
  int idx = 0;
  while (idx < len) {
    std::string item(str + idx);
    values.push_back(item);
    idx += item.length() + 1;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::parse_properties(const std::vector<std::string>& keywords,
                                  const std::map<std::string, std::string>& keyword_synonyms,
                                  const char* delimiters) {
  ErrorCode rval;

  // master keyword map, mapping user-set words in cubit to canonical property names
  std::map< std::string, std::string > keyword_map(keyword_synonyms);

  for (std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i) {
    keyword_map[*i] = *i;
  }

  // the set of all canonical property names
  std::set< std::string > prop_names;
  for (prop_map::iterator i = keyword_map.begin();
       i != keyword_map.end(); ++i) {
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
    if (MB_SUCCESS != rval)
      return rval;
    property_tagmap[(*i)] = new_tag;
  }

  // now that the keywords and tags are ready, iterate over all the actual geometry groups
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
    if (MB_SUCCESS != rval)
      return rval;
    if (grp_sets.size() == 0)
      continue;

    for (prop_map::iterator i = properties.begin();
         i != properties.end(); ++i) {
      std::string groupkey = (*i).first;
      std::string groupval = (*i).second;

      if (property_tagmap.find(groupkey) != property_tagmap.end()) {
        Tag proptag = property_tagmap[groupkey];
        const unsigned int groupsize = grp_sets.size();
        for (unsigned int j = 0; j < groupsize; ++j) {
          rval = append_packed_string(proptag, grp_sets[j], groupval);
        }
      }
    }
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_value(EntityHandle eh, const std::string& prop, std::string& value) {
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr(proptag, &eh, 1, &data, &ignored);
  if (rval != MB_SUCCESS)
    return rval;
  value = static_cast<const char*>(data);
  return MB_SUCCESS;
}

ErrorCode DagMC::prop_values(EntityHandle eh, const std::string& prop,
                             std::vector< std::string >& values) {

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


ErrorCode DagMC::get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &proptag, NULL, 1, all_ents);
  if (MB_SUCCESS != rval)
    return rval;

  std::set<std::string> unique_values;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (MB_SUCCESS != rval)
      return rval;
    unique_values.insert(values.begin(), values.end());
  }

  return_list.assign(unique_values.begin(), unique_values.end());
  return MB_SUCCESS;
}

ErrorCode DagMC::entities_by_property(const std::string& prop, std::vector<EntityHandle>& return_list,
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
  void* vals[2] = {NULL, (dimension != 0) ? &dimension : NULL };
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, tags, vals, 2, all_ents);
  if (MB_SUCCESS != rval)
    return rval;

  std::set<EntityHandle> handles;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (MB_SUCCESS != rval)
      return rval;
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

void DagMC::tokenize(const std::string& str,
                     std::vector<std::string>& tokens,
                     const char* delimiters) const {
  std::string::size_type last = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos  = str.find_first_of(delimiters, last);
  if (std::string::npos == pos)
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back(str.substr(last, pos - last));
      last = str.find_first_not_of(delimiters, pos);
      pos  = str.find_first_of(delimiters, last);
      if (std::string::npos == pos)
        pos = str.size();
    }
}

Tag DagMC::get_tag(const char* name, int size, TagType store,
                   DataType type, const void* def_value,
                   bool create_if_missing) {
  Tag retval = 0;
  unsigned flags = store | MB_TAG_CREAT;
  // NOTE: this function seems to be broken in that create_if_missing has
  // the opposite meaning from what its name implies.  However, changing the
  // behavior causes tests to fail, so I'm leaving the existing behavior
  // in place.  -- j.kraftcheck.
  if (!create_if_missing)
    flags |= MB_TAG_EXCL;
  ErrorCode result = MBI->tag_get_handle(name, size, type, retval, flags, def_value);
  if (create_if_missing && MB_SUCCESS != result)
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;

  return retval;
}

ErrorCode DagMC::build_preconditioner()
{
  ErrorCode rval;
  
  std::vector<std::string> keywords;
  std::string sdf = "sdf";
  keywords.push_back(sdf);
  rval = parse_properties(keywords);
  MB_CHK_SET_ERR(rval, "Failed to parse the signed distance field property");
  
  std::vector<EntityHandle> vols;
  rval = entities_by_property(sdf, vols, 3);
  MB_CHK_SET_ERR(rval, "Failed to find signed distance field volumes");

  // Range vols;
  // const int three = 3;
  // const void* const three_val[] = {&three};
  // rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag,
  //                                           three_val, 1, vols );
  // MB_CHK_SET_ERR(rval, "Could not get the model volumes.");

  sdfs.resize(num_entities(3));
  
  for(unsigned int i = 0; i < vols.size(); i++) {
    
    EntityHandle vol = vols[i];

    std::string val;
    rval = prop_value(vol, sdf, val);
    MB_CHK_SET_ERR(rval, "Could not get preconditioner step size from property tag");
    double sdf_step_size = std::atof(val.c_str());
    
    // never need this for the implicit compliment
    if(GTT->is_implicit_complement(vol)) continue;

    int vol_id = get_entity_id(vol);

    double x_min,x_max,y_min,y_max,z_min,z_max;
    CartVect lower_corner, upper_corner;
    rval = getobb(vol, lower_corner.array(), upper_corner.array());
    MB_CHK_SET_ERR(rval, "Could not get volume bounds.");

    x_min = lower_corner[0]; y_min = lower_corner[1]; z_min = lower_corner[2];
    x_max = upper_corner[0]; y_max = upper_corner[1]; z_max = upper_corner[2];
    
    int num_x_steps = (int)(x_max-x_min)/sdf_step_size;
    if ( (x_max-x_min)/sdf_step_size - (double)num_x_steps > 0.0 ) num_x_steps+=1/sdf_step_size;
    int num_y_steps = (int)(y_max-y_min)/sdf_step_size;
    if ( (y_max-y_min)/sdf_step_size - (double)num_y_steps > 0.0 ) num_y_steps+=1/sdf_step_size;
    int num_z_steps = (int)(z_max-z_min)/sdf_step_size;
    if ( (z_max-z_min)/sdf_step_size - (double)num_z_steps > 0.0 ) num_z_steps+=1/sdf_step_size;    
    
    std::cout << " x min/max " << x_min <<  "/" << x_max << std::endl;
    std::cout << " y min/max " << y_min <<  "/" << y_max << std::endl;
    std::cout << " z min/max " << z_min <<  "/" << z_max << std::endl;
    std::cout << " x steps " << num_x_steps << std::endl;
    std::cout << " y steps " << num_y_steps << std::endl;
    std::cout << " z steps " << num_z_steps << std::endl;
    
    // if any of these are zero, move on
    if( 0 == num_x_steps || 0 == num_y_steps || 0 == num_z_steps ) {
      continue;
    }
    
    int num_x_verts = num_x_steps + 2;
    int num_y_verts = num_y_steps + 2;
    int num_z_verts = num_z_steps + 2;
    
    num_x_verts++;
    num_y_verts++;
    num_z_verts++;
    
    x_min -= sdf_step_size;
    y_min -= sdf_step_size;
    z_min -= sdf_step_size;
    
    SignedDistanceField* sdf = new SignedDistanceField(x_min, y_min, z_min,
						       sdf_step_size,
						       num_x_verts, num_y_verts, num_z_verts);
    
    MB_CHK_SET_ERR(rval, "Could not create preconditioner ScdBox.");
    
    std::cout << "Populating preconditioner..." << std::endl;
    rval = populate_preconditioner_for_volume(vol, sdf);
    MB_CHK_SET_ERR(rval, "Could not populate preconditioner.");
    std::cout << "done." << std::endl;
}
  
  return MB_SUCCESS;
}

ErrorCode DagMC::populate_preconditioner_for_volume(EntityHandle &vol, SignedDistanceField* sdf) {
  ErrorCode rval;

  EntityHandle tree_root;
  rval = get_root(vol, tree_root);
  MB_CHK_SET_ERR(rval, "Could not retrieve the volume's OBB tree root.");
  
  int xpnts, ypnts, zpnts;  
  std::vector<double> sdvs;
  
  sdf->get_dims(xpnts,ypnts,zpnts);
  
  for(unsigned int k = 0 ; k < zpnts; k++){
    for(unsigned int j = 0 ; j < ypnts; j++){  
      for(unsigned int i = 0 ; i < xpnts; i++){
	CartVect vert_coords = sdf->get_coords(i,j,k);
	
	CartVect closest_loc;
	EntityHandle facet, surf;
	rval = GTT->obb_tree()->closest_to_location( vert_coords.array(), tree_root, closest_loc.array(), facet, &surf);	
	MB_CHK_SET_ERR(rval, "Could not get the closest location.");

	CartVect vec = closest_loc - vert_coords;

	//get facet connectivity
	const EntityHandle * con = NULL;
	int len;
	rval = MBI->get_connectivity(facet, con, len);
	MB_CHK_SET_ERR(rval, "Could not get connectivity of nearest facet.");
	CartVect coords[3];
	rval = MBI->get_coords( con, len, coords[0].array() );      
	MB_CHK_SET_ERR(rval, "Coult not get coordinates of nearest facet verts.");
	// get normal of triangle
	CartVect v0 = coords[1] - coords[0];
	CartVect v1 = coords[2] - coords[0];
	CartVect norm = (v0*v1);
	EntityHandle vols[2];
	//adjust normal for sense wrt the volume datastruct being populated
	rval = MBI->tag_get_data( sense_tag(), &surf, 1, vols );
	MB_CHK_SET_ERR(rval, "Could not get nearest surface sense wrt volume.");
	int sense = 1;
	if (vol != vols[0]) sense = -1;
	norm *= sense;

	// if direction to the facet opposes the adjusted facet normal this distance
	// value should be negative because it is on the outside
	//	double dot_prod = norm%vec;
	double sdv = fabs(vec.length());
	//	if( fabs(dot_prod) < 1e-8 ) {
	int inout = 0;
	rval = point_in_volume(vol, vert_coords.array(), inout);
	MB_CHK_SET_ERR(rval,"Point in volume failed");
	if(0 == inout) sdv*=-1;
	// }
	// else {
	//   if( dot_prod < 0 ) {
	//     sdv *= -1;
	//   }
	// }
	sdvs.push_back(sdv);
	MB_CHK_SET_ERR(rval, "Could not set SDF tag data.");
      }
    }
  }

  //set the data for the field
  sdf->set_data(sdvs);
  sdfs[index_by_handle(vol)] = sdf;

#ifdef SDF_WRITE
  ScdBox* b;
  rval = sdf->create_scdBox(b, moab_instance());
  MB_CHK_SET_ERR(rval, "Failed to create ScdBox from SDF.");
#endif
  
  return MB_SUCCESS;
}

ErrorCode DagMC::precondition_point_in_volume(EntityHandle volume, const double xyz[3], int& result, bool& preconditioned) {
  preconditioned = false;
  ErrorCode rval;
  
  // retrieve the signed distance value from the field
  double sdv, sdv_err;
  rval = find_sdv(volume, xyz, sdv, sdv_err);
  MB_CHK_SET_ERR(rval, "Failed to retrieve signed distance value for volume " << get_entity_id(volume));
    
  // if the signed distance value is large enough,
  // use its sign to determine the point containment
  if(fabs(sdv) > sdv_err) {
    preconditioned = true;
    result = sdv > 0.0;
  }
  
  return MB_SUCCESS;
  
}

  
ErrorCode DagMC::precondition_closest_to_location(EntityHandle volume, const double coords[3], double& result, bool& preconditioned) {
  ErrorCode rval;
  preconditioned = false;
  
  double sdv_err;
  
  rval = find_sdv(volume,coords,result, sdv_err);
  MB_CHK_SET_ERR(rval, "Failed to find signed distance value for volume " << get_entity_id(volume));
  
  // check that the nearest intersection can be considered valid
  // (is larger than interpolation error evaluation)
  if ( fabs(result) > sdv_err && result > 0.0) {
    preconditioned = true;
    result = result - sdv_err;
  }
  
  return MB_SUCCESS;
}

/** precondition ray using physical distance limit */
ErrorCode DagMC::precondition_ray_fire(const EntityHandle volume,
				       const double ray_start[3],
				       const double ray_dir[3],
				       const double ray_len,
				       EntityHandle& next_surf,
				       double& next_surf_dist,
				       bool& preconditioned) {

  if (ray_len == 0) {return MB_FAILURE;}
  
  CartVect ray_end = CartVect(ray_dir[0], ray_dir[1], ray_dir[2]);
  ray_end.normalize();
  ray_end *= ray_len;
  ray_end += CartVect(ray_start);

  return precondition_ray_fire(volume, ray_start, ray_end.array(), next_surf, next_surf_dist, preconditioned);
  
}
		     
   
ErrorCode DagMC::precondition_ray_fire(const EntityHandle volume,
				  const double ray_start[3],
				  const double ray_end[3],
				  EntityHandle& next_surf,
				  double& next_surf_dist,
				  bool& preconditioned) {
  preconditioned = false;
  // calculate ray length and try to account for all space in between
  double ray_len = (CartVect(ray_start)-CartVect(ray_end)).length();

  if(!get_signed_distance_field(volume)) return MB_SUCCESS;
  
  ErrorCode rval;
  
  double ssdv, ssdv_err, esdv, esdv_err;
  
  rval = find_sdv(volume,ray_start,ssdv, ssdv_err);
  MB_CHK_SET_ERR(rval,"Could not find start point sdv");

  rval = find_sdv(volume,ray_end,esdv, esdv_err);
  MB_CHK_SET_ERR(rval,"Could not find end point sdv");

  // if the starting point is outside of the volume, we have a problem, fire ray
  if(ssdv < 0) { preconditioned = false; return MB_SUCCESS; }
  // end point is outside the volume, fire ray to get intersection distance and surf
  if(esdv < 0) { preconditioned = false; return MB_SUCCESS; }
  
  if( (ray_len-ssdv-esdv) < -(ssdv_err + esdv_err) ) {
    preconditioned = true;
    Range surfs;
    rval = MBI->get_child_meshsets(volume, surfs);
    MB_CHK_SET_ERR(rval, "Failed to get child surfaces of volume " << get_entity_id(volume));
    // spoof next surface and distance
    next_surf_dist = 2.0 * ray_len;
    next_surf = surfs[0];
  }

  return MB_SUCCESS;
}

ErrorCode DagMC::find_sdv(const EntityHandle volume,
			  const double pnt[3],
			  double &sdv,
			  double &err) {
  ErrorCode rval;
  SignedDistanceField* sdf = get_signed_distance_field(volume);
  if( sdf != NULL) {
    sdv = sdf->find_sdv(pnt);
  }
  else {
    return MB_ENTITY_NOT_FOUND;
  }

  sdf->get_err(err);
  
  return MB_SUCCESS;
};
  

} // namespace moab
