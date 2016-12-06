#include "uwuw_preprocessor.hpp"
#include "moab/ProgOptions.hpp"

int main(int argc, char* argv[])
{

  ProgOptions po("uwuw_preproc: a tool for preprocessing DAGMC files to incorporate UWUW workflow information");

  bool verbose = false;
  bool fatal_errors = false;
  bool simulate = false;

  std::string lib_file;
  std::string dag_file;
  std::string out_file;

  po.addOpt<void>( "verbose,v", "Verbose output", &verbose);
  po.addOpt<void>( "fatal,f", "Fatal errors cause exit", &fatal_errors);
  po.addOpt<void>( "simulate,s", "Fatal errors cause exit", &simulate);

  po.addRequiredArg<std::string>("dag_file", "Path to DAGMC file to proccess", &dag_file);
  po.addOpt<std::string>("lib_file,l","Path to PyNE Material Library file to proccess", &lib_file);
  po.addOpt<std::string>("output,o", "Specify the output filename (default "")", &out_file);

  po.addOptionHelpHeading("Options for loading files");

  po.parseCommandLine(argc, argv);

  if(lib_file=="") {
    std::cout << "need to set the library" << std::endl;
    exit(1);
  }


  if(out_file== "" )
    out_file = dag_file;

  // make new preprocessor
  uwuw_preprocessor *uwuw_preproc = new uwuw_preprocessor(lib_file,dag_file,out_file,verbose,fatal_errors);

  // process the materials
  uwuw_preproc->process_materials();

  // process the tallies
  uwuw_preproc->process_tallies();

  if(simulate) {
    uwuw_preproc->print_summary();
    return 0;
  }

  // write the material data
  uwuw_preproc->write_uwuw_materials();

  // write the tally data
  uwuw_preproc->write_uwuw_tallies();

  // thats all
  return 0;
}
