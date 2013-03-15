 /*
  Tartarus - cpp program to use DAG geometry and allows debugging of 
  geometry, runs in parallel
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <iterator>

#include "DagMC.hpp"
#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "mpi.h"

using moab::DagMC;
using namespace std;

#define DAG DagMC::instance()

  struct Location
  {
    double x;
    double y;
    double z;
  };

/* function prototype */
void PrintUsage(void); // prototype for print usage
void TransportCycle(int, int, int , double [3] ); // prototype for TransportCycle
int FindStartVolume(double[3],MBEntityHandle &); // prototype for FindStartVolume
void GetIsoDirection(double[3], int); // get a random isotropically distributed uvw
double RandomNumber(int); // random number generator
int CheckGraveyard(MBEntityHandle); // check if the current volume is the graveyard
int CheckForGraveyard(void);  // check to ensure graveyard exists anywhere in the model
void PrintLostParticles(std::vector< std::vector< Location > > ); // Dump the lost particle information 
void DumpParticles(std::vector< std::vector< Location > > ); // Dump the lost particle information to vtk
void EndMPI( void );
void GetNpsRange(int,int,int,int&,int&);
void CollectLostParticleData(int, std::vector<double>, std::vector<double>, std::vector<double>,
			      std::vector<double>&, std::vector<double>&, std::vector<double>&);


int main(int argc ,char *argv[])
{
  // generic variables
  int i; 
  // 
  double Position[3]; // particle position
  int NumPar = 0 ; // number of particles
  int ErrorCode ; //Error code returns from DAG functions
  std::string DagInputFile; // DAG filename

  MPI::Init (argc,argv); //start MPI

  int CpuId = MPI::COMM_WORLD.Get_rank( ); // CPU Id #
  int WorldSize = MPI::COMM_WORLD.Get_size( ); // Size of the world
 
  if( argc <= 1 ) // no arguments provided
    {
      if ( CpuId == 0 )
	{
	   PrintUsage(); // print how to use code	   
	}
      EndMPI();
    }
  else
    {
      for ( i = 1 ; i <= argc-1 ; i++)
	{
	  if ( argv[i] == std::string("--dag") )  // if dag keywrod
	    {
	      DagInputFile = argv[i+1];
	      std::cout << "Opening DAG input " << DagInputFile << std::endl;
	      ErrorCode = DAG -> load_file(DagInputFile.c_str()); // open the Dag file
	      if ( ErrorCode != MB_SUCCESS )
		{
		  if ( CpuId == 0 )
		    {
		      std::cout << "Could not load dag file" << std::endl;
		    }
		  EndMPI();
  		  exit(5);
		}
	    }
	  if ( argv[i] == std::string("--nps") ) // if nps keyword
	    {
	      NumPar = atoi(argv[i+1]);
	    }
	  if ( argv[i] == std::string("--pos") ) // if pos keyword
	    {
	      // Collect the position locations
	      Position[0] = atof(argv[i+1]);
	      Position[1] = atof(argv[i+2]);
	      Position[2] = atof(argv[i+3]);
	    }

	}

      if( NumPar == 0 )
	{
	  if ( CpuId == 0 )		     
	    {
	      std::cout << "Nps = 0!" << std::endl;
	    }
	  EndMPI();
	  exit(5);
	}



      /* otherwise everything ok, init dag */
      ErrorCode = DAG -> init_OBBTree(); // initialise the geometry

      /* ensure graveyard exists */
      ErrorCode = CheckForGraveyard();
      if ( CpuId == 0 )
	{
	  if ( ErrorCode == 0 )
	    {
	      std::cout << "No Graveyard Detected!" << std::endl;
	    }
	}

      if ( ErrorCode == 0 )
	{
	  EndMPI();
	  return 1;
	}

	      
      // determine which particle id's the sub tasks are running
      int NpStart = 0,NpFinish = 0 ;
      GetNpsRange(CpuId,WorldSize,NumPar,NpStart,NpFinish);
      // Run transport
      TransportCycle(123456789,NpStart,NpFinish,Position); // using arbitrary seed

      /* all done */
      EndMPI();
    }

  return 0;
}

void PrintUsage(void)
{
  std::cout << "tartarus - a program to determine gaps in geometry and thus collect  " <<  std::endl;
  std::cout << "lost particle information" << std::endl;
  std::cout << std::endl;
  std::cout << "--dag <filename> = determines which dag file the code will use" << std::endl;
  std::cout << std::endl;

  return;
}

/* 
Loop through the required number of particles
 */

// NumPar is the number of particles, and position is the start location (x,y,z)
void TransportCycle(int Seed, int StartPar, int EndPar, double Position[3])
{


  // Temporary Variables
  int i; //,j;
  int rval; 

  // mpi vars
  int CpuId,SizeWorld; // cpu id number and size of the mpi world

  // Routine Variables
  MBErrorCode ErrorCode;
  MBEntityHandle Surface = 0; // surface number
  MBEntityHandle Volume = 0; // volume number

  double Direction[3]; // vector containg the particle direction
  int stride = 5000000; // rn stride

  double HitDist = 10.0;
  double OriginalHistory[3]; // original location of start history

  int grave; // in graveyard flag 1 true 0 false
  //  int count = 0;

  // Lost particles
  int LostParticle = 0; // number of lost particles;
  Location CurrentPosition; //Struct of a position 

  // Collection of the lost particle data
  std::vector<Location>  ParticleHistory; // vector of particle locations
  std::vector<Location>::const_iterator position_iterator;
  // Lost Particle History
  std::vector<std::vector<Location>> LostParticleHistory;
  std::vector<std::vector<Location>> AllTheLostParticles;
  std::vector<double> x,y,z;
  std::vector<double> x_coord,y_coord,z_coord;


  static DagMC::RayHistory history; //keep the facets
 
  CpuId = MPI::COMM_WORLD.Get_rank( ); // CPU Id #
  SizeWorld = MPI::COMM_WORLD.Get_size( ); // Size of the world

  // Begining of Routine
  OriginalHistory[0] = Position[0], OriginalHistory[1] = Position[1], OriginalHistory[2] = Position[2];


  for ( i = StartPar ; i <= EndPar ; i++ )
    {
      //      if( i%(NumPar/10) == 0 )
      //	std::cout << "nps = " << i << std::endl;

      Position[0]=OriginalHistory[0], Position[1]=OriginalHistory[1], Position[2]=OriginalHistory[2] ; // reset location
      CurrentPosition.x=Position[0], CurrentPosition.y=Position[1], CurrentPosition.z=Position[2]; // record the start location
      ParticleHistory.push_back(CurrentPosition); // push it to the storage array

      //x.push_back(Position[0]); // initial array allocation
      //y.push_back(Position[1]); //
      //z.push_back(Position[2]); //
    
      rval = FindStartVolume(Position,Volume); //determine our start Volume
      if( rval == 1)
	{
	  GetIsoDirection(Direction,Seed+(i*stride));       // pick a random isotropic direction
	}
      else
	{
	  std::cout << "Particle not in any volume" << std::endl;
	  exit(0);
	}
 
      do 
      {
	ErrorCode = DAG->ray_fire(Volume,Position,Direction,Surface,HitDist,&history); // find next collision
	//	std::cout << history.prev_facets[1] << std::endl;
	
	if (MB_SUCCESS != ErrorCode )
	  {
	    std::cout << "ray_fire error" << std::endl;
	  }
	
	HitDist = HitDist + 1.0e-5; // bump the particle

	Position[0] = (Position[0]+(HitDist*Direction[0]));
	Position[1] = (Position[1]+(HitDist*Direction[1]));
	Position[2] = (Position[2]+(HitDist*Direction[2])); // update position

	CurrentPosition.x=Position[0], CurrentPosition.y=Position[1], CurrentPosition.z=Position[2];
	ParticleHistory.push_back(CurrentPosition); /* append to "ray history" */

	rval = FindStartVolume(Position,Volume); // get the current volume of position 

	if ( rval == 0) // particle is lost
	  {
	    LostParticle++; // increment the counter by one
	    // std::cout << "Lost particle " << LostParticle << std::endl;
	    LostParticleHistory.push_back(ParticleHistory); // add this particle to bank
	    x.push_back(Position[0]),y.push_back(Position[1]),z.push_back(Position[2]);
	    ParticleHistory.clear(); // clear the particle history
	    break; // new history
	  }  

	grave = CheckGraveyard(Volume);  //  see if its the graveyard
	if (grave == 1)
	  {
	    ParticleHistory.clear(); // wipe out the history
	    history.reset();
	    break; // we are in graveyard end of normal history
	  }

       } while( HitDist >= 0.0 );
 
   } // end of transport loop


  // transport is done 
  //Get the total number of lost particles
  int TotalNumLost = 0;
  MPI::COMM_WORLD.Barrier( ); 
  MPI::COMM_WORLD.Reduce(&LostParticle,&TotalNumLost,1,MPI::INT,MPI::SUM,0); // add up the local lost particle counts
  std::cout << "tot lost = " << TotalNumLost << " " << LostParticle << std::endl;
 
  // get all the lost particle data
  CollectLostParticleData(LostParticle,x,y,z,x_coord,y_coord,z_coord);

  if ( CpuId == 0 )
    {
      std::cout << "There were " << TotalNumLost << " lost particles " << std::endl;
      return;
      if ( LostParticleHistory.size() > 0 ) 
	{
	  DumpParticles(LostParticleHistory);
	  return;
	}
    }
  else
    {
      return;
    }
}

/* 
   for a given position determine which volume we start in
*/
int FindStartVolume(double Position[3], MBEntityHandle &Volume)
{
  int i;
  int NumVol; // the number of volumes in the problem
  MBErrorCode ErrorCode; 
  int result; // result of the point_in_volume query.
  
  NumVol = DAG->num_entities(3); // get the number of volumes
  
  for ( i = 1 ; i <= NumVol-1 ; i++) // loop over all volumes
    {
      Volume = DAG->entity_by_index(3,i); // convert index to MBHandle
      ErrorCode = DAG->point_in_volume(Volume,Position,result); // Find which volume we are in
      //std::cout << Volume << " " << Position[0] << " " << Position[1] << " " << Position[2] << " " << result << " " << ErrorCode << std::endl;
      if( ErrorCode != 0 )
	{
	  Volume = 0;
	  return 0;
	}
      if ( result == 1 )
	{
	  return 1;
	}
    }

  // if we are here then in no volumes
  Volume = 0;
  return 0;
}

/* 
Generates a random direction isotropic on the unit sphere.
*/

void GetIsoDirection(double Direction[3],int Seed)
{
  //  double u,v,w; // direction vectors
  double theta,phi; // components

  theta = RandomNumber(Seed)*2*3.14149;
  phi = RandomNumber(Seed+1)*3.14149;

  Direction[0] = cos(theta)*sin(phi);
  Direction[1] = sin(theta)*sin(phi);
  Direction[2] = cos(phi);

  return;
}

/*
  returns a double precision random number between 0 and 1
*/
double RandomNumber(int Seed)
{
  std::mt19937 gen(Seed);
  std::uniform_real_distribution<> uni(0,1);
  return uni(gen);
}

int CheckGraveyard(MBEntityHandle Volume)
{
  MBErrorCode rval;
  bool ReturnVal;
  std::vector<std::string> props;

  props.push_back("graveyard");

  rval = DAG->parse_properties(props);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  ReturnVal = DAG->has_prop(Volume,"graveyard");
 
  if ( ReturnVal )
    return 1;
  else
    return 0;
}

/* 
   Loop through all volumes to ensure we have a graveyard in the list
   This program relies on the fact that there is an all encompassing volume
   that particles must enter when traversing the geometry in question
*/
int CheckForGraveyard()
{
  int i;
  int rval;
  int NumVol;
  //  bool ReturnVal;
  std::vector<std::string> props;
  MBEntityHandle Volume;

  props.push_back("graveyard");
  rval = DAG->parse_properties(props);

  if (MB_SUCCESS != rval) 
    {
      std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
      exit(EXIT_FAILURE);
    }

  NumVol =  DAG->num_entities(3); // get the number of volumes
  
  for ( i = 1 ; i <= NumVol-1 ; i++ )
    {
      Volume = DAG->entity_by_index(3,i); // convert index to MBHandle
      rval = CheckGraveyard(Volume);
      if(rval == 1)
	{
	  return 1;
	}
    }

  //std::cout << "No Graveyard detected" << std::endl;
  return 0;
}

void PrintLostParticles(std::vector<std::vector<Location> > AllLostParticles )
{
  // routine variables
  int LostParticleCount = 0; 
  std::vector<std::vector<Location> >::const_iterator LostParticleNo;
  std::vector<Location>::const_iterator PositionHistory;

  for ( LostParticleNo = AllLostParticles.begin() ; LostParticleNo != AllLostParticles.end() ; ++LostParticleNo )
    {
      ++LostParticleCount;
      std::cout << LostParticleCount << std::endl;
      for ( PositionHistory = LostParticleNo->begin() ; PositionHistory != LostParticleNo->end() ; ++PositionHistory )
	{
	   std::cout << (PositionHistory->x) << " "
		     << (PositionHistory->y) << " "
		     << (PositionHistory->z) << std::endl;	  
	}
    }

  return;
}


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

  for ( LostParticleNo = AllLostParticles.begin() ; LostParticleNo != AllLostParticles.end() ; ++LostParticleNo )
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

/* terminates all MPI communications nicely */
void EndMPI( )
{
   MPI::COMM_WORLD.Barrier( ); 
   MPI::Finalize( );
   return;
}

/* get nps range */
// Determines the split of which particle ranges to simulate on which cpu
void GetNpsRange( int CpuId, int WorldSize, int NumPar, int &NumStart, int &NumFinish)
{
  int Nps,Remainder;
  int start,finish;

  Remainder =  NumPar%WorldSize;
  Nps = NumPar/WorldSize;
  //  std::cout << Remainder << std::endl;
  //  std::cout << Nps << std::endl;

  if ( CpuId == 0 )
    {
      start = 1;
      finish = ((CpuId+1)*Nps)+Remainder; // always give remainder to master
    }
  else // slave tasks
    {
      start = (CpuId*Nps)+2;
      finish = ((CpuId+1)*Nps)+1; 
    }

  std::cout << "CPU " << CpuId << " starts = " << start << " finshes = " << finish << std::endl;
  NumStart = start;
  NumFinish = finish;

  return;
}

void CollectLostParticleData(int NumPar, std::vector<double> x_local, std::vector<double> y_local, std::vector<double> z_local,
			     std::vector<double> &x_pos, std::vector<double> &y_pos, std::vector<double> &z_pos)
{
  int CpuId, NumProc;
  int NumParMaster;
    
  CpuId = MPI::COMM_WORLD.Get_rank( ); // CPU Id #
  NumProc = MPI::COMM_WORLD.Get_size( ); // Size of the world

  int i = 0;

  MPI::COMM_WORLD.Barrier( ); 

  /* Get size of the particle arrays */
  if ( CpuId == 0 )
    {

      int* num_lost = new int [NumProc];

      num_lost[0] = NumPar;
      //  std::cout << num_lost[i] << std::endl;
      for ( i = 1 ; i <= NumProc-1 ; i++)
	{
	  MPI::COMM_WORLD.Recv (&num_lost[i],1,MPI_INT,i,i);
	}

      // masters x data
      x_pos.resize(x_local.size());
      std::copy (x_local.begin(),x_local.end(), x_pos.begin() );
      // masters y data
      y_pos.resize(y_local.size());
      std::copy (y_local.begin(),y_local.end(), y_pos.begin() );
      // masters z data
      z_pos.resize(z_local.size());
      std::copy (z_local.begin(),z_local.end(), z_pos.begin() );

      // get the x data
      for ( i = 1 ; i <= NumProc-1 ; i++ )
	{
	  if ( num_lost[i] != 0 )
	    {
	      int tag;
	      std::cout << "Creating double array " << std::endl;
	      double* tmp_store = new double [num_lost[i]] ; // tmp array storage for x pos data
	      std::vector<double> tmp_vector;                // tmp vector for x pos data
	      MPI::COMM_WORLD.Irecv ( &tmp_store,num_lost[tag],MPI_DOUBLE,MPI_ANY_TAG,tag); // get the messages back
	      tmp_vector.resize( num_lost[i] );
	      std::copy ( tmp_store, tmp_store+num_lost[i],tmp_vector.begin() );
	      x_pos.insert(x_pos.end(),tmp_vector.begin(),tmp_vector.end() );
	      MPI::COMM_WORLD.Irecv ( &tmp_store,num_lost[tag],MPI_DOUBLE,MPI_ANY_TAG,tag); // get the messages back
	      tmp_vector.resize( num_lost[i] );
	      std::copy ( tmp_store, tmp_store+num_lost[i],tmp_vector.begin() );
	      y_pos.insert(y_pos.end(),tmp_vector.begin(),tmp_vector.end() );
	      MPI::COMM_WORLD.Irecv ( &tmp_store,num_lost[tag],MPI_DOUBLE,MPI_ANY_TAG,tag); // get the messages back
	      tmp_vector.resize( num_lost[i] );
	      std::copy ( tmp_store, tmp_store+num_lost[i],tmp_vector.begin() );
	      z_pos.insert(z_pos.end(),tmp_vector.begin(),tmp_vector.end() );

	      delete[] tmp_store; // delete the particle positions
	      tmp_vector.clear();
	      tmp_vector.shrink_to_fit();
	    }

	}

      delete[] num_lost; // delete the lost particle counts from each cpu


    }
  else
    {
      MPI::COMM_WORLD.Send (&NumPar,1,MPI_INT,0,CpuId); // send the number of lost particles (from each core)
      if ( NumPar != 0 )
	{
	  MPI::COMM_WORLD.Isend (&x_local,NumPar,MPI_DOUBLE,0,CpuId); // send the data per core 
	  MPI::COMM_WORLD.Isend (&y_local,NumPar,MPI_DOUBLE,0,CpuId); // send the data per core 
	  MPI::COMM_WORLD.Isend (&z_local,NumPar,MPI_DOUBLE,0,CpuId); // send the data per core 
	}
    }



  return;
}
