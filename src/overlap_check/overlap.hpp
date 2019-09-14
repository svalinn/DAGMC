
#include <map>
#include <memory>
#include <set>
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"

using namespace moab;

using OverlapMap = std::map<std::set<int>, CartVect>;

ErrorCode
check_file_for_overlaps(std::shared_ptr<Interface> MBI,
                        OverlapMap& overlap_map);

void report_overlaps(const OverlapMap& overlap_map);
