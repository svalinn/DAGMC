#include "ExN01UserScoreWriter.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"
#include "G4SystemOfUnits.hh"
#include <map>

#include "iMesh.h"
#include "iMesh_extensions.h"

ExN01UserScoreWriter::ExN01UserScoreWriter()
  : G4VScoreWriter()
{;}

ExN01UserScoreWriter::~ExN01UserScoreWriter()
{;}

/* Takes the data contained within the mesh and dumps all scores into a
   file using the option provided, in this instance we always output as
   a MOAB mesh */
void ExN01UserScoreWriter::DumpAllQuantitiesToFile(const G4String& fileName,
						   const G4String& option)
{
  // get the mesh
  G4cout << "Dumping mesh " << fScoringMesh->GetWorldName() << " to file" << G4endl;

  // retrieve the map                                                                                                              
  MeshScoreMap scMap = fScoringMesh->GetScoreMap();

  // get the number of bins
  int num_bins[3];
  fScoringMesh->GetNumberOfSegments(num_bins);
  int extents[6] = {0,0,0,
		    num_bins[0],num_bins[1],num_bins[2]};
  
  G4ThreeVector size;
  size = fScoringMesh->GetSize();


  G4cout << size << " " << num_bins[0] << " " << num_bins[1] << " " << num_bins[2] << G4endl;

  std::vector<double> xvals_v = generate_bin_bounds(num_bins[0],size.getX()/cm);
  double* xvals = &xvals_v[0];
  
  std::vector<double> yvals_v = generate_bin_bounds(num_bins[1],size.getY()/cm);
  double* yvals = &yvals_v[0];

  std::vector<double> zvals_v = generate_bin_bounds(num_bins[2],size.getZ()/cm);
  double* zvals = &zvals_v[0];

  G4cout << size << " " << num_bins[0] << " " << num_bins[1] << " " << num_bins[2] << G4endl;

  iMesh_Instance imesh;
  iBase_EntitySetHandle root;
  iBase_EntitySetHandle mesh_handle;

  char* options = NULL;
  int options_len = 0;
  int ierr;

  // new imesh instance                            
  iMesh_newMesh(options, &imesh, &ierr, options_len);

  // create root set                              
  iMesh_getRootSet(imesh, &root, &ierr);

  // create the structure mesh
  iMesh_createStructuredMesh(imesh,extents,NULL,xvals,yvals,zvals,
			     0,-1,-1,-1,0,0,1,NULL,&ierr);

  // loop over the scores on this mesh and dump to moab
  MeshScoreMap::iterator it;
  for ( it = scMap.begin() ; it != scMap.end() ; ++it) {

    if ( it == scMap.end()) return;

    G4String score_name = it->first;
    std::map<G4int, G4double*> *score = it->second->GetMap();
    iBase_TagHandle tag;
    iBase_TagHandle iter_tag;

    G4cout << "Writing score for " << score_name << G4endl;
	  

    // create a tag
    iMesh_createTag(imesh,score_name.c_str(),1,iBase_DOUBLE,
		    &tag,&ierr, score_name.length());

    /*
    // create a tag
    iMesh_createTag(imesh,"iteration_order",1,iBase_INTEGER,
		    &iter_tag,&ierr, 16);
    */
    
    iBase_EntityHandle *hexes;
    int num_hexes;
    int num_hexes_alloc;

    // get the hex entities
    iMesh_getEntities(imesh,root,iBase_ALL_TYPES,iMesh_HEXAHEDRON,&hexes,
		      &num_hexes_alloc,&num_hexes,&ierr);

    G4cout << num_hexes_alloc << " " << num_hexes << G4endl;

    G4int idx; // mesh index
    double result = 0.0; // the result from the mesh
    // loop over each score 
    for(int x = 0; x < num_bins[0]; x++) {
      for(int y = 0; y < num_bins[1]; y++) {
	for(int z = 0; z < num_bins[2]; z++) {
	
	  idx = GetIndex(x,y,z);
	  //	  value = *((score->find(idx))->second);
	  std::map<G4int, G4double*>::iterator value = score->find(idx);
	  if(value != score->end()) result = *(value->second);
	  
	  iMesh_setArrData(imesh,&hexes[idx],1,tag,
	  		   (char*)(&result),sizeof(double),&ierr);

	  //	  double idx_d = static_cast<double>(idx);
	  //iMesh_setArrData(imesh,&hexes[idx],1,tag,
	  //		   (char*)(&idx_d),sizeof(double),&ierr);

	  //	  iMesh_setArrData(imesh,&hexes[idx],1,iter_tag,
	  //		  (char*)(&idx), sizeof(G4int), &ierr);

	  G4cout << x << " " << y << " " << z << " " << idx << " " << result << G4endl;
	  

	} // z                                                                                                                       
      } // y                                                                                                                         
    } // x   
  }

  // save the file
  iMesh_save(imesh,root,(fileName+".h5m").c_str(),options,&ierr,8,0);

  return;
}

// function to generate the bin bounds from -end_coord to end_coord with
// num_bound bins
std::vector<double> ExN01UserScoreWriter::generate_bin_bounds(int num_bounds, double end_coord)
{
  std::vector<double> bin_bounds;

  // width of each bin
  double bin_width = 2.0*end_coord / static_cast<double>(num_bounds);

  double boundary; // tmp calculated boundary
  // loop over the bins
  for ( int i = 0 ; i < num_bounds ; i++ )
  {
    boundary = -1.0*end_coord + static_cast<double>(i)*bin_width;
    bin_bounds.push_back(boundary);
  }
  // last bin is the upper extrema to avoid roundoff
  bin_bounds.push_back(end_coord);
  
  return bin_bounds;
}
