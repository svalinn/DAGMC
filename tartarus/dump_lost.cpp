#include <iostream>

void DumpParticles(std::vector<std::vector<Location> > AllLostParticles )
{
  unsigned int i;
  // routine variables
  unsigned int LostParticleCount = 0; 
  std::vector<std::vector<Location> >::const_iterator LostParticleNo;
  std::vector<Location>::const_iterator PositionHistory;

  std::stringstream intermediate;
  std::string FileString = "lostparticle_";
  std::string FilePost   = ".vtk";
  std::string Filename; 

  std::ofstream OutputFile; // Output file pointer

  std::cout << "Dumping lost particle information to file..." << std::endl;

  for ( LostParticleNo = AllLostParticles->begin() ; LostParticleNo != AllLostParticles->end() ; ++LostParticleNo )
    {
      ++LostParticleCount;
      intermediate << FileString << LostParticleCount << FilePost;
      Filename = intermediate.str();
     
      OutputFile.open(Filename.c_str(), ios::out | ios::trunc ); //open the file
      Filename = std::string(); // clear the filename string for reuse.
      intermediate.str(std::string()); // clear the strings for reuse

      OutputFile << "# vtk DataFile Version 2.0 " << std::endl;
      OutputFile << "Lost Particle Information" << std::endl;
      OutputFile << "ASCII" << std::endl;
      OutputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

      // Print out the points
      OutputFile << "POINTS " << AllLostParticles[LostParticleCount-1].size() << " FLOAT" << std::endl;
      for ( PositionHistory = LostParticleNo->begin() ; PositionHistory != LostParticleNo->end() ; ++PositionHistory )
	{
	   OutputFile << (PositionHistory->x) << " "
		      << (PositionHistory->y) << " "
		      << (PositionHistory->z) << std::endl;	  
	}

      OutputFile << std::endl;
      OutputFile << "CELLS " << (AllLostParticles[LostParticleCount-1].size())-1 << " " 
		             << (((AllLostParticles[LostParticleCount-1].size())-1)*3) << std::endl;
      for ( i = 0 ; i <= AllLostParticles[LostParticleCount-1].size()-2 ; i++ )
	{
	  OutputFile << "2 " << i << " " << i+1 << std::endl;
	}
      OutputFile << "CELL_TYPES " << AllLostParticles[LostParticleCount-1].size()-1 << std::endl;
      for ( i = 0 ; i <= AllLostParticles[LostParticleCount-1].size()-2 ; i++ )
	{
	  OutputFile << "3" << std::endl;
	}

      OutputFile.close(); // close the file
    }

  return;
}
