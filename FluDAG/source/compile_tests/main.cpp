#include <iostream>
#include "cpp/fluka_funcs.h"

int main()
{
  std::cout << "In main.cpp";
    
  char* filename;
  int clen;
  char* ftol;
  int ftlen;
  int parallel_read;
  double dagmc_version;
  int moab_version,max_pbl;

  dagmcinit_(filename,&clen,ftol,&ftlen,&parallel_read,&dagmc_version,&moab_version,&max_pbl);

  return 1;
}
