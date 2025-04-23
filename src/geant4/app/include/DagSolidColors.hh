#include <cmath>
#include <tuple>
#include <vector>

#ifndef UNIFORM_COLOUR_GEN_HPP
#define UNIFORM_COLOUR_GEN_HPP

// struct for RGB
struct RGB {
    int r, g, b;
};

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

double hue2rgb(double p, double q, double t) {
    if (t < 0.0) t += 1.0;
    if (t > 1.0) t -= 1.0;
    if (t < 1.0/6.0) return p + (q - p) * 6.0 * t;
    if (t < 1.0/2.0) return q;
    if (t < 2.0/3.0) return p + (q - p) * (2.0/3.0 - t) * 6.0;
    return p;
}

RGB hslToRgb(double h, double s, double l) {
    double r, g, b;

    h = fmod(h, 360.0) / 360.0;  // Normalize to [0, 1]

    if (s == 0.0) {
        r = g = b = l;  // achromatic
    } else {
        double q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        double p = 2 * l - q;
        r = hue2rgb(p, q, h + 1.0/3.0);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1.0/3.0);
    }

    return {
        static_cast<int>(std::round(r * 255)),
        static_cast<int>(std::round(g * 255)),
        static_cast<int>(std::round(b * 255))
    };
}

class UniformColorGenerator {
  public:
  UniformColorGenerator(int num_colours);
  ~UniformColorGenerator();
  void Generate();
  std::vector<RGB> GetColours();
  private:
  int num_colours;
  std::vector<RGB> colours;
};

#endif //UNIFORM_COLOUR_GEN_HPP
