#include "DagSolidColors.hh"
#include <iostream>
// Generate a linearly spaced array of values
std::vector<double> linspace(double start, double stop, int num) {
  std::vector<double> result;
  if (num <= 0) {
    return result;
  }
    
  if (num == 1) {
    result.push_back(start);
    return result;
  }

  double step = (stop - start) / (num - 1);
  for (int i = 0; i < num; ++i) {
    result.push_back(start + step * i);
  }
  return result;
}

RGB hslToRgb(double h, double s, double l) {
  double c = (1. - std::abs(2.*l - 1.0))*s;
  double x = c * (1. - std::abs(fmod(h/60., 2.) - 1.));  
  double m = l - c/2.;
  
  if ( h < 60.) {
    return { c+m, x+m, m };
  } else if (h < 120.) {
    return { x+m, c+m, m };
  } else if (h < 180.) {
    return { m, c+m, x+m };
  } else if (h < 240.) {
    return { m, x+m, c+m };
  } else if (h < 300.) {
    return { x+m, m, c+m };
  } else {
    return { c+m, m, x+m };
  }

}

// class to manage colour generation
UniformColorGenerator::UniformColorGenerator(int num_colors) {
  this->num_colors = num_colors;
}

// destructor
UniformColorGenerator::~UniformColorGenerator(){
}

// generator some colours
void UniformColorGenerator::Generate() {
  // angles
  std::vector<double> angles = linspace(0.,360.,num_colors);
  for ( int i = 0 ; i < num_colors ; i++ ) {
    RGB rgb = hslToRgb(angles[i],0.75,0.75);
    colors.push_back(rgb);
  }
}

// return the colours
std::vector<RGB> UniformColorGenerator::GetColors(){
  return colors;
}
