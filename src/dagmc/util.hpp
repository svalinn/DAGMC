#ifndef _DAGMC_UTIL
#define _DAGMC_UTIL

#include <string>

namespace dagmc_util {

inline void lowercase_str(std::string& input) {
  std::transform(input.begin(), input.end(), input.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

}  // namespace dagmc_util
#endif