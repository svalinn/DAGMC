#include "fluka_funcs.h"

extern "C" void flukam(const int &flag);

int main(int argc, char* argv[]) 
{

	char *cfile = "myfile";
	if (argc > 2)
	{
		cfile = argv[1];
	}
	int clen = 1;  // geom
	char ftol;     // faceting tolerance
	int ftlen = 1; // faceting tolerance
	int parallel_file_mode = 0;  // parallel read mode
	

        dagmcinit_(cfile, &clen
                char *ftol,  int *ftlen, // faceting tolerance
                &parallel_file_mode, null, null, null);
                // double* dagmc_version, int* moab_version, int* max_pbl )

  return 0;
}




