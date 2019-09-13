
#include <map>
#include <memory>
#include <set>
#include "moab/Core.hpp"

using namespace moab;

using OverlapMap = std::map<std::set<int>,std::array<double,3>>;

ErrorCode
check_file_for_overlaps(std::shared_ptr<Interface> MBI,
                        OverlapMap& overlap_map);
