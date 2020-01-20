
#ifndef DAGMC_PROGRESSBAR_H
#define DAGMC_PROGRESSBAR_H

#include <string>

class ProgressBar {

 public:
  // constructor
  ProgressBar();

  // destructor
  ~ProgressBar();

  void set_value(double val);

  static bool is_terminal();

 private:
  int current {0};
};

#endif // HEADER GUARD
