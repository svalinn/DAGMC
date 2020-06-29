
#ifndef DAGMC_PROGRESSBAR_H
#define DAGMC_PROGRESSBAR_H

class ProgressBar {

 public:
  // constructor
  ProgressBar() {
    // initialize bar
    set_value(0.0);
  };

  // destructor
  ~ProgressBar();

  void set_value(double val);

  static bool is_terminal();

 private:
  int current {0};
  bool need_final_newline {true};
};

#endif // HEADER GUARD
