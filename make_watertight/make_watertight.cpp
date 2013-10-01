// ********************************************************************
// Brandon Smith
// August, 2009

/* _curve_to_be_tested_for_watertightness_
      vert1 X X vert1
            | |
      vert2 X |
  surf1     | |    surf2
            | |
      vert3 X X vert2
            | |
      vert4 X X vert3                   */

// input:  h5m filename, tolerance
// output: watertight h5m

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

// modified by Andrew Davis 2012
// Updated deprecated MOAB calls

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




MBInterface *MOAB();

int main(int argc, char **argv) 
  {

 

     MBErrorCode result;
     result= mw_func::make_model_watertight( argc, argv);
     if(gen::error(MB_SUCCESS!=result, "could not make model watertight")) return result;
       return 0;  
  }

MBInterface *MBI() 
{
    static MBCore instance;
    return &instance;
}
