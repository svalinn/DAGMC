#include <iostream>
#include <string>
#include <vector>
#include <map>

class HistogramManager
{
 public:
  // Constructor
  HistogramManager();
  // Destructor
  ~HistogramManager();
  // Add histogram collection
  void add_histogram(int volume_id, std::string particle_name);
  void add_histogram(int volume_id, int pdc_number);
  int get_histogram_id(int volume_id, std::string particle_name);
  int get_histogram_id(int volume_id, int pdc_number);
  std::vector<int> get_senstitive_particles(int volume_id);
  void print_histogram_collection();
 private:
  std::map<int,std::map<int, int> >         histogram_collection;
  int histogram_id;
};
