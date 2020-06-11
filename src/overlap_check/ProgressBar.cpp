

#include "ProgressBar.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#endif

#define BAR_WIDTH 72

void ProgressBar::~ProgressBar() {
  if (need_final_newline) {
    std::cout << std::endl;
  }
}

bool ProgressBar::is_terminal() {
#ifdef _POSIX_VERSION
  return isatty(STDOUT_FILENO) != 0;
#else
  return false;
#endif
}


void ProgressBar::set_value(double val) {

  if (!is_terminal()) {
    return;
  }

  val = std::max(std::min(100.0, val), 0.0);

  if ((int)val == current) {
    return;
  } else {
    current = (int)val;
  }

  // current value and opening bracket
  std::stringstream bar;
  bar <<  std::setfill(' ') << std::setw(3) << (int)val;
  bar << "% |";

  // remaining width of the bar, leaving room for the closing characters
  // and the arrowhead
  int remaining_width = BAR_WIDTH - 9;
  int width = (int)((double)remaining_width * val / 100);
  bar << std::string(width, '=');
  bar << std::string(1, '>');
  bar << std::string(remaining_width - width, ' ');

  // closing bracket
  bar << "|+";

  // write the bar to screen
  std::cout << '\r' << bar.str() << std::flush;
  need_final_newline = true;
  if (val >= 100.0) {
    std::cout << "\n";
    need_final_newline = false;
  }
}
