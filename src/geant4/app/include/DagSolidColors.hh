#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#ifndef UNIFORM_COLOR_GEN_HPP
#define UNIFORM_COLOR_GEN_HPP

// struct for RGB
struct RGB {
    double r, g, b;
};

class UniformColorGenerator {
  public:
  UniformColorGenerator(int num_colors);
  ~UniformColorGenerator();
  void Generate();
  std::vector<RGB> GetColors();

  private:
  int num_colors;
  std::vector<RGB> colors;
};

#endif //UNIFORM_COLOR_GEN_HPP
