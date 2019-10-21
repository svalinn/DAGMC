

#include "ProgressBar.hpp"

#include <sstream>
#include <iostream>
#include <iomanip>

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#endif

#define BAR_WIDTH 72

ProgressBar::ProgressBar() {
  // initialize bar
  set_value(0.0);
}

ProgressBar::~ProgressBar() {
  std::cout << "\n";
}

bool ProgressBar::is_terminal() {
#ifdef _POSIX_VERSION
  return isatty(STDOUT_FILENO) != 0;
#else
  return false;
#endif
}


void ProgressBar::set_value(double val) {

  if ((int)val == current) {
    return;
  } else {
    current = (int)val;
  }

  std::string bar;

  if (!is_terminal())
    return;

  // set the bar percentage
  if (val >= 100.0) {
    bar.append("100");
  } else if (val <= 0.0) {
    bar.append("  0");
  } else {
    std::stringstream ss;
    ss << std::setfill(' ') << std::setw(3) << (int)val;
    bar.append(ss.str());
  }

  bar.append("% |");
  // remaining width of the bar
  int remaining_width = BAR_WIDTH - 8;

  // set the bar width
  if (val >= 100.0) {
    bar.append(remaining_width, '=');
  } else if (val < 0.0) {
    bar.append(remaining_width, ' ');
  } else {
    int width = (int)((double)remaining_width * val / 100);
    bar.append(width, '=');
    bar.append(1, '>');
    bar.append(remaining_width - width - 1, ' ');
  }

  bar.append("|+");

  // write the bar
  std::cout << '\r' << bar << std::flush;
  if (val >= 100.0) {
    std::cout << "\n";
  }

  // reset the bar value
  bar = "";
}
