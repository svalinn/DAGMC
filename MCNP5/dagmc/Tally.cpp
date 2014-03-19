// MCNP5/dagmc/Tally.cpp

#include <cassert>
#include <iostream>
#include <cmath>

#include "Tally.hpp"
#include "TrackLengthMeshTally.hpp"
#include "KDEMeshTally.hpp"
#include "CellTally.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Tally::Tally(const TallyInput& input) 
    : input_data(input), data(NULL)
{
    bool total_energy_bin = true;
    unsigned int num_energy_bins = 0;

    int num_boundaries = input_data.energy_bin_bounds.size();

    // turn off total energy bin if only one bin exists
    if (num_boundaries == 2)
    {
        total_energy_bin = false;
        num_energy_bins = 1;
    }
    else
    {
        // determine total number of energy bins requested
        assert(num_energy_bins > 2);
        num_energy_bins = num_boundaries;
    }

    data = new TallyData(num_energy_bins, total_energy_bin); 
}
//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
Tally::~Tally()
{
   delete data;
}
//---------------------------------------------------------------------------//
// FACTORY METHOD
//
// 3D Flux Tally Types:
//    Other      | Estimator type | Generic Type     Name         Status
//                                 Mesh, Cell, Surf 
//   ------        --------------   -------------   ----------    -----------
//  Unstructured | Track Length   | Mesh Tally   || unstr_track   implemented
//  KDE          | Integral Track | Mesh Tally   || kde_track     KD's Thesis
//  KDE          | SubTrack       | Mesh Tally   || kde_subtrack  implemented
//  KDE          | Collision      | Mesh Tally   || kde_coll      implemented
//               | Track Length   | Cell         || cell_track    implemented 
//               | Collision      | Cell         || cell_coll     implemented 
//---------------------------------------------------------------------------//
Tally *Tally::create_tally(const TallyInput& input)
{
    Tally *newTally = NULL;
          
    if (input.tally_type == "unstr_track")
    {
        newTally = new moab::TrackLengthMeshTally(input);
    }
    else if (input.tally_type == "kde_track")
    {
        KDEMeshTally::Estimator estimator = KDEMeshTally::INTEGRAL_TRACK;
        newTally = new KDEMeshTally(input, estimator); 
    }
    else if (input.tally_type == "kde_subtrack")
    {
        KDEMeshTally::Estimator estimator = KDEMeshTally::SUB_TRACK;
        newTally = new KDEMeshTally(input, estimator); 
    }
    else if (input.tally_type == "kde_coll")
    {
        // This line is not necessary because COLLISION is default estimator.
        KDEMeshTally::Estimator estimator = KDEMeshTally::COLLISION;
        newTally = new KDEMeshTally(input, estimator); 
    }
    else if (input.tally_type == "cell_track")
    {
        newTally = new CellTally(input, TallyEvent::TRACK);
    }
    else if (input.tally_type == "cell_coll")
    {
        newTally = new CellTally(input, TallyEvent::COLLISION);
    }
    else 
    {
        std::cout << "Warning: " << input.tally_type
                  << " is not a valid tally type." << std::endl;
    }
         
    return newTally;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void Tally::end_history()
{
    data->end_history();
}

const TallyData& Tally::getTallyData()
{
      return *data;
}
//---------------------------------------------------------------------------//
// PROTECTED INTERFACE
//---------------------------------------------------------------------------//
bool Tally::get_energy_bin(double energy, unsigned int& ebin)
{
    bool bin_found = false;

    if (energy_in_bounds(energy))
    {
        if (data->get_num_energy_bins() == 1)
        {
            ebin = 0;
            bin_found = true;
        }
        else  // in bounds and more than one energy bin
        {
            // Case where we are close to the highest energy
            double tol = 1e-6;
            unsigned int maxbin = input_data.energy_bin_bounds.size() - 1;

            if (fabs(energy - input_data.energy_bin_bounds.at(maxbin)) < tol)
            {
                ebin =  maxbin;
                bin_found = true;
            }

            unsigned int i = 0;

            while (!bin_found)
            {
                if (input_data.energy_bin_bounds.at(i) <= energy &&
                    energy < input_data.energy_bin_bounds.at(i+1))
                {
                    ebin = i;
                    bin_found = true;
                }
                else
                {
                    ++i;
                }
            }  // end while
        }  // end else in bounds and >1 energy bin
    }  // end if in bounds

    return bin_found;
}
//---------------------------------------------------------------------------//
bool Tally::energy_in_bounds(double energy)
{
    unsigned int maxbin = input_data.energy_bin_bounds.size() - 1;

    return !(energy < input_data.energy_bin_bounds.at(0) ||
             energy > input_data.energy_bin_bounds.at(maxbin));
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/Tally.cpp
