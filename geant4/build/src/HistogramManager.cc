#include "HistogramManager.hh"

HistogramManager::HistogramManager()
{
  histogram_id = 1;
  return;
}

HistogramManager::~HistogramManager()
{
  return;
}

void HistogramManager::add_histogram(int volume_id, std::string particle_name)
{
  // convert particle name to pdc
  // note, by using pdc, we cannot index by the heavy ion type
  // need to map to zzaaa or something
  int pdc_number;
  add_histogram(volume_id,pdc_number);
  return;
}

void HistogramManager::add_histogram(int volume_id, int pdc_number)
{
  // convert particle name to pdc
  // note, by using pdc, we cannot index by the heavy ion type
  // need to map to zzaaa or something
  std::map<int,int> histogram_pdcnum_map;
  histogram_pdcnum_map = histogram_collection[volume_id];
  histogram_pdcnum_map[pdc_number] = histogram_id;
  histogram_collection[volume_id] = histogram_pdcnum_map;
  histogram_id++;
  return;
}

std::vector<int> HistogramManager::get_senstitive_particles(int volume_id)
{
  std::map<int,int> histogram_pdcnum_map = histogram_collection[volume_id];
  std::vector<int> pdcs;
  std::map<int,int>::iterator it;
  for ( it = histogram_pdcnum_map.begin() ; it != histogram_pdcnum_map.end() ; ++it ) {
    pdcs.push_back(it->first);
  }
  return pdcs;
}

int HistogramManager::get_histogram_id(int volume_id, std::string particle_name)
{
  // convert particle name to pdc
  int pdc_number;
  int id_num = get_histogram_id(volume_id,pdc_number);
  return id_num;
}

int HistogramManager::get_histogram_id(int volume_id, int pdc_number)
{
  std::map<int,int> histogram_pdcnum_map;
  histogram_pdcnum_map = histogram_collection[volume_id];
  int id_num  = histogram_pdcnum_map[pdc_number];
  return id_num;
}

void HistogramManager::print_histogram_collection()
{
  // map of histrogram
  std::map<int,std::map<int,int> >::iterator it;
  for ( it = histogram_collection.begin() ; it != histogram_collection.end() ;
        ++it ) {
    std::cout << "Volume id " << it->first << " has tally collections " << std::endl;
    std::map<int,int>::iterator iter;
    for ( iter = it->second.begin() ; iter != it->second.end() ; ++iter) {
      std::cout << "  pdc number = " << iter->first << " histo idx  = ";
      std::cout << iter->second << std::endl;
    }
  }
}
