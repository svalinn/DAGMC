#include <iostream>
#include <iomanip> // for setprecision                                  
#include <limits> // for min/max values  
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/GeomTopoTool.hpp"                               
#include "MBSkinner.hpp"

#include "gen.hpp"

namespace gen 
{

  bool error( const bool error_has_occured, const std::string message ) 
  {
    if(error_has_occured) 
      {
	if("" == message) 
	  {
	    std::cout << "Error at " << __FILE__ << ":" << __LINE__
		      << std::endl;
	  } 
	else 
	  {
	    std::cout << message << std::endl;
	  }
	return true;
      } 
    else 
      {
	return false;
      }
  }
}

