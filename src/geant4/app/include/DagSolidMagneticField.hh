#ifndef DagSolidMagneticField_HH
#define DagSolidMagneticField_HH 1

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"

#include <array>
#include <vector>
#include <string>

class G4MagneticField;

// enum for the sampling methods
enum Sampling { NEAREST, BILINEAR_INTERPOLATION };

double BilinearInterpolation(const std::array<double,2> x_stencil,
			     const std::array<double,2> y_stencil,
			     const std::array<double,4> values,
			     const double x,
			     const double y);

class MagneticField : public G4MagneticField {
    public:
        // constructor
        MagneticField(const Sampling mode = BILINEAR_INTERPOLATION);
        // destructor
        ~MagneticField() override;
        // main get field value for Geant4
        void GetFieldValue(const G4double Point[4], double *Bfield) const;
        // load the field data
        void LoadFile(const std::string filename);
    private:
        // process the field ready to use
        void ProcessField();

        // main access method
        void InterpolateField(const double x, const double y, const double z,
            double &field_x, double &field_y, double &field_z);
        
        // field lookup functions
        void FieldByNearest(const double r, const double z,
            double &final_r, double &field_z, double &field_t);

        void FieldByInterpolation(const double r, const double z,
				  double &field_x, double &field_y, double &field_z); 

        // data access functions
  std::array<std::array<double,4>, 3> GetFourFieldValuesByIndex(const int r_index,
            const int z_index);
  std::array<std::array<double,4>, 3> GetFourFieldValuesByValue(const double r,
            const double z);

        // indexing functions
        int FindRadialIndex(const double r);
        int FindVerticalIndex(const double z);
        void FindRadialandVerticalIndices(const double r, const double z,
            int &r_ind, int &z_ind);
        int GlobalIndexByRadialandVerticalIndex(const int r_index, 
                const int z_index);
        int GlobalIndexByRadialandVerticalPosition(const double, 
                const double z);  
        std::array<double,2> LookupRPtsByIndex(const int r_index);
        std::array<double,2> LookupZPtsByIndex(const int z_index);
        std::array<double,2> LookupRPtsByValue(const double r);
        std::array<double,2> LookupZPtsByValue(const double z);  

    private:
        Sampling mode; // the mode of sampling
        double scale; // the scale factor
        int n_r; // number of r points
        int n_z; // number of z points
        std::vector<double> r_point; // the unique grid points in r
        std::vector<double> z_point; // the unique grid points in z

        std::vector<double> radial_pos; /// temp
        std::vector<double> vertical_pos; /// temp
        std::vector<std::array<double,3> > b_field; // the field data

};
#endif
