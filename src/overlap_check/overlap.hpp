
#include <map>
#include <set>
#include "moab/Core.hpp"

using namespace moab;

ErrorCode
check_file_for_overlaps(std::string filename,
                        std::map<std::set<int>,std::array<double,3>>& overlap_map);
