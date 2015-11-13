#include "ExN01UserScoreWriter.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"
#include "G4SystemOfUnits.hh"
#include <map>

ExN01UserScoreWriter::ExN01UserScoreWriter()
  : G4VScoreWriter()
{
  ;
}

ExN01UserScoreWriter::~ExN01UserScoreWriter()
{
  ;
}

/* Takes the data contained within the mesh and dumps all scores into a
   file using the option provided, in this instance we always output as
   a MOAB mesh */
void ExN01UserScoreWriter::DumpAllQuantitiesToFile(const G4String& fileName,
    const G4String& option)
{
  std::string file_name(fileName);
  // check to make sure the filename ends with .h5m
  if( fileName.substr(fileName.length()-4) != ".h5m")
    file_name += ".h5m";

  // get the mesh
  G4cout << "Dumping mesh " << fScoringMesh->GetWorldName() << " to file" << G4endl;

  // retrieve the map
  MeshScoreMap scMap = fScoringMesh->GetScoreMap();

  // get the number of bins
  int num_bins[3];
  fScoringMesh->GetNumberOfSegments(num_bins);
  int extents[6] = {0,0,0,
                    num_bins[0],num_bins[1],num_bins[2]
                   };

  G4ThreeVector size;
  size = fScoringMesh->GetSize();

  std::vector<double> x_vals = generate_bin_bounds(num_bins[0],size.getX()/cm);
  std::vector<double> y_vals = generate_bin_bounds(num_bins[1],size.getY()/cm);
  std::vector<double> z_vals = generate_bin_bounds(num_bins[2],size.getZ()/cm);

  // load the input file
  moab::ErrorCode rval;

  std::vector<moab::EntityHandle> mesh_elements;
  // create the mesh
  // mesh elements is vector of hex entity handles, ordered
  rval = generate_moab_mesh(x_vals, y_vals, z_vals, mesh_elements);

  // loop over the scores on this mesh and dump to moab
  MeshScoreMap::iterator it;
  for ( it = scMap.begin() ; it != scMap.end() ; ++it) {
    if ( it == scMap.end()) return;

    G4String score_name = it->first;
    std::map<G4int, G4double*> *score = it->second->GetMap();

    std::cout << "Writing score for " << score_name << std::endl;

    // create a tag
    moab::Tag tag_handle;
    rval = MBI()->tag_get_handle(score_name.c_str(), 1,
                                 moab::MB_TYPE_DOUBLE, tag_handle,
                                 moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);

    G4int idx; // mesh index
    double result = 0.0; // the result from the mesh
    // loop over each score

    for(int x = 0; x < num_bins[0]; x++) {
      for(int y = 0; y < num_bins[1]; y++) {
        for(int z = 0; z < num_bins[2]; z++) {

          idx = GetIndex(x,y,z);
          std::map<G4int, G4double*>::iterator value = score->find(idx);

          if(value != score->end())
            result = *(value->second);
          else
            result = 0.0;

          // set the tag data
          rval = MBI()->tag_set_data(tag_handle,&(mesh_elements[idx]), 1, &result);
        } // z
      } // y
    } // x
  }

// save the file
  rval = MBI()->write_mesh((fileName).c_str());

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
  for ( int i = 0 ; i < num_bounds ; i++ ) {
    boundary = -1.0*end_coord + static_cast<double>(i)*bin_width;
    bin_bounds.push_back(boundary);
  }
  // last bin is the upper extrema to avoid roundoff
  bin_bounds.push_back(end_coord);

  return bin_bounds;
}

// generates a structured moab mesh
moab::ErrorCode ExN01UserScoreWriter::generate_moab_mesh(std::vector<double> x_bins,
    std::vector<double> y_bins,
    std::vector<double> z_bins,
    std::vector<moab::EntityHandle> &mesh_elements)
{
  moab::ErrorCode rval;

  moab::EntityHandle input_set;
  // new mesh instance
  rval = MBI()->create_meshset( moab::MESHSET_SET, input_set );

  int num_x = x_bins.size();
  int num_y = y_bins.size();
  int num_z = z_bins.size();

  moab::EntityHandle eh;
  std::vector<moab::EntityHandle> verts;
  double coord[3] = {0.,0.,0.};
  // loop over the coordinates
  for ( int i = 0 ; i < num_x ; i++ ) {
    for ( int j = 0 ; j < num_y ; j++ ) {
      for ( int k = 0 ; k < num_z ; k++ ) {
        coord[0] = x_bins[i];
        coord[1] = y_bins[j];
        coord[2] = z_bins[k];
        rval = MBI()->create_vertex(coord,eh);
        rval = MBI()->add_entities(input_set,&(eh),1);
        verts.push_back(eh);
      }
    }
  }

  // create hexes
  moab::EntityHandle connect[8];
  moab::EntityHandle hex_eh;

  for ( int i = 0 ; i < num_x-1 ; i++ ) {
    for ( int j = 0 ; j < num_y-1 ; j++ ) {
      for ( int k = 0 ; k < num_z-1 ; k++ ) {
        // change this order at your peril
        // I know where you live
        connect[0] = verts[k + (j*num_z) + (i*num_z*num_y)];
        connect[1] = verts[k + (j*num_z) + (i*num_z*num_y) + num_z];
        connect[2] = verts[k + (j*num_z) + (i*num_z*num_y) + (num_z*num_y) + num_z];
        connect[3] = verts[k + (j*num_z) + (i*num_z*num_y) + (num_z*num_y)];
        connect[4] = verts[k + (j*num_z) + (i*num_z*num_y) + 1];
        connect[5] = verts[k + (j*num_z) + (i*num_z*num_y) + num_z + 1];
        connect[6] = verts[k + (j*num_z) + (i*num_z*num_y) + (num_z*num_y) + num_z + 1];
        connect[7] = verts[k + (j*num_z) + (i*num_z*num_y) + (num_z*num_y) + 1];

        rval = MBI()->create_element(moab::MBHEX,connect,8,hex_eh);
        if(rval != moab::MB_SUCCESS)
          std::cout << "Failed to create hex" << std::endl;
        rval = MBI()->add_entities(input_set,&(hex_eh),1);
        mesh_elements.push_back(hex_eh);
      }
    }
  }
  rval = MBI()->write_mesh("halfway.h5m");

  return rval;
}

// MOAB interface
moab::Interface *MBI()
{
  static moab::Core instance;
  return &instance;
}
