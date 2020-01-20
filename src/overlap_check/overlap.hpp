
#include <map>
#include <memory>
#include <set>
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"

using namespace moab;

using OverlapMap = std::map<std::set<int>, CartVect>;

// checks MOAB instance with a loaded file for overlaps
ErrorCode
check_instance_for_overlaps(std::shared_ptr<Interface> MBI,
                            OverlapMap& overlap_map,
                            int pnts_per_edge = 0);

// a convenience function for reporting overlaps in the model
void report_overlaps(const OverlapMap& overlap_map);
