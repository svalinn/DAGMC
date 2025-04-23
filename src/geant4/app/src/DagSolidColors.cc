#include "DagSolidColors.hh"

// class to manage colour generation
UniformColourGenerator::UniformColorGenerator(int num_colours) {
  num_colours = num_colours;
}

// destructor
UniformColourGenerator::~UniformColorGenerator(){
  delete colours;
}

// generator some colours
void UniformColourGenerator::Generate() {
  // angles
  std::vector<double> angles = linspace(0.,360.,num_colours);
  for ( int i = 0 ; i < num_colours ; i++ ) {
    RGB rgb = hslToRgb(angles[i],0.75,0.75);
    colours.append(rgb);
  }
}

// return the colours
std::vector<RGB> UniformColorGenerator::GetColours(std::map<int,RGB> colours){
  return colours;
}
