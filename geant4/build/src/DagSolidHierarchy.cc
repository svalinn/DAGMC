#include <iostream>
#include "DagSolidHierarchy.hh"

/*
int main(int argc, char* argv[]) {

  moab::Core *mbi = new moab::Core();

  moab::EntityHandle file_set;
  moab::ErrorCode ec = mbi->create_meshset(moab::MESHSET_SET,file_set);
  ec = mbi->load_file("test.h5m",&file_set);

  DetermineHierarchy *dh = new DetermineHierarchy(mbi);

  ec = dh->DetermineTheHierarchy();

  return 0;
}
*/

// Constructor
DetermineHierarchy::DetermineHierarchy(moab::Interface *moab, moab::GeomTopoTool *geomtool) {
   mbi = moab;
   
   if(geomtool != NULL)
     gtt = geomtool;
   else
     gtt = new moab::GeomTopoTool(mbi);
   
   network = new Network(); // network to store relationships
}

// Destructor
DetermineHierarchy::~DetermineHierarchy(){};

// Determine the hierarchy in the problem
moab::ErrorCode DetermineHierarchy::DetermineTheHierarchy(bool loaded) {
  moab::ErrorCode ec; // for error checking

  // if data is already loaded 
  if(!loaded) {
    ec = gtt->setup_implicit_complement();
    MB_CHK_SET_ERR(ec,"Failed to setup the implicit complement");
  }
  
  moab::Range vols;
  ec = gtt->get_gsets_by_dimension(3,vols);
  MB_CHK_SET_ERR(ec,"Failed to get 3d entities from MOAB");

  range_it vol_it;

  std::map<moab::EntityHandle,moab::Range> vol_surfs;
  // loop over vols and get children
  for ( vol_it = vols.begin() ; vol_it != vols.end() ; vol_it++ ) {
    ec = DetermineVolume(*vol_it); // determine the children volume and
    // add to network
    MB_CHK_SET_ERR(ec,"Failed to get 3d entities from MOAB");
  }

  moab::EntityHandle implicit_complement;
  ec = gtt->get_implicit_complement(implicit_complement);
  MB_CHK_SET_ERR(ec,"Failed to get the implicit complement");

  // Remove all upward links, assuming implicit complement is the top vol
  network->Directionalise(implicit_complement);
  //
  return ec;
}

// Given a volume, determine which volumes are inside of it
moab::ErrorCode DetermineHierarchy::DetermineVolume(moab::EntityHandle volume) {
  // due to how geometry is set up we have a given volume, which has n surfaces
  // as children, those surfaces also belong to other volumes. If some of the n 
  // surfaces belong to the 'other' volume then it is a child of the searched volume

  moab::Range final_set;

  moab::Range child_surfaces; // children_surface sets
  moab::ErrorCode ec = mbi->get_child_meshsets(volume,child_surfaces); // get childen

  range_it surf_it;
  // loop over the children sets and determine which other volumes share the 
  moab::Range child_volumes;
  for (surf_it = child_surfaces.begin() ; surf_it != child_surfaces.end() ; ++surf_it) {
    ec = GetParentSets(volume,*surf_it,child_volumes);
    // 
  }
  // the parent volume range now includes all volumes that shares at least surface with the child surfaces
  // of the current volume

  // step through the list and remove any volumes that do not share all their child surfaces with the 
  // current volume 
  moab::Range parent_surfaces;
  ec = mbi->get_child_meshsets(volume,parent_surfaces);
  MB_CHK_SET_ERR(ec,"Failed to get child meshsets");

  range_it vol_it; 
  for ( vol_it = child_volumes.begin() ; vol_it != child_volumes.end() ; ++vol_it ) {

    moab::Range child_surfaces;  
    // get the children of the current volume
    ec = mbi->get_child_meshsets(*vol_it,child_surfaces);
    MB_CHK_SET_ERR(ec,"Failed to get child meshsets");
  
    // Note - if we do use all surfaces of the child volume we assume
    // volume not being inside of the parent.
   
    // if the intersect the ranges
    moab::Range overlap_set = intersect(parent_surfaces,child_surfaces);
    // if range not empty, we are not a child
    if( subtract(overlap_set,child_surfaces).empty() ) final_set.insert(*vol_it);
  }
  for ( vol_it = final_set.begin() ; vol_it != final_set.end() ; ++vol_it) 
    network->AddLink(volume,*vol_it);
  
  return moab::MB_SUCCESS;
}

// Given a surface set and volume set, return back other volumes that also share this surface
moab::ErrorCode DetermineHierarchy::GetParentSets(moab::EntityHandle volume, moab::EntityHandle surface, 
                                                  moab::Range &parent_volumes) {
  // get the parents of the surface set
  moab::ErrorCode ec = mbi->get_parent_meshsets(surface,parent_volumes);

  range_it surf_it; //oop over the surface sets
  for ( surf_it = parent_volumes.begin() ; surf_it != parent_volumes.end() ; ++surf_it ) {
    if(*surf_it == volume ) parent_volumes.erase(*surf_it); // remove the current volume from the range
  }

  return ec;
}

// Given the parent EH return its children
moab::ErrorCode DetermineHierarchy::GetChildren(moab::EntityHandle parent,
			    std::vector<moab::EntityHandle> &children) {

  std::vector<moab::EntityHandle> child_volumes;
  child_volumes = network->FindNeighbours(parent);
  children = child_volumes;
  if(children.size() == 0)
    return moab::MB_FAILURE;
  else
    return moab::MB_SUCCESS;
}

// Given a parent id, return vector of children ids
moab::ErrorCode DetermineHierarchy::GetChildren(int parent,
			    std::vector<int> &children) {
  moab::ErrorCode ec;
  moab::EntityHandle parent_eh = gtt->entity_by_id(3,parent);
  std::vector<int> children_id;
  std::vector<moab::EntityHandle> children_eh;
  ec = GetChildren(parent_eh,children_eh);

  // turn eh's into ids
  std::vector<moab::EntityHandle>::iterator it;
  for ( it = children_eh.begin() ; it != children_eh.end() ; ++it ) {
    children_id.push_back(gtt->global_id(*it));
  }
  children = children_id;
  return ec;
}

moab::ErrorCode DetermineHierarchy::FindParent(moab::EntityHandle child, moab::EntityHandle &parent) {
  parent = network->FindParent(child);
  if (parent != 0)
    return moab::MB_SUCCESS;
  else
    return moab::MB_FAILURE;
}

moab::ErrorCode DetermineHierarchy::FindParent(int child, int &parent) {
  moab::EntityHandle child_eh = gtt->entity_by_id(3,child);
  moab::EntityHandle parent_eh;
  moab::ErrorCode ec = FindParent(child,parent_eh);
  parent = gtt->global_id(parent_eh);
  return ec;
}

void DetermineHierarchy::PrintRange(moab::Range item) {
  range_it it;
  for ( it = item.begin() ; it != item.end() ; ++it ) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
}

void DetermineHierarchy::PrintTree(moab::EntityHandle vol, moab::Range item) {
  // int id = gtt->global_id(vol);
  range_it it;
  for ( it = item.begin() ; it != item.end() ; ++it ) {
    std::cout << gtt->global_id(vol) << "->" << gtt->global_id(*it) << std::endl;
  }
}
