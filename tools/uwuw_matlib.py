#!/usr/bin/env python
################################################################################
# This script read an MCNP5 input file and creates a material library for use
# with the UWUW workflow. The naming scheme in the material library is 
# compatible with the mcnp2cad naming scheme, when the mcnp2cad flag "-U" is 
# used. Note that there will be one entry per material, even when materials have
# multiple densities.
#
# This script is dependent on PyNE (http://pyne.io).
################################################################################
import argparse

from pyne.material import Material, MultiMaterial, MaterialLibrary
from pyne.mcnp import mats_from_inp

def uwuw_matlib(inp, out):
    mats = mats_from_inp(inp)
    uwuw_mats = {}
    for mat in mats.values():
        if isinstance(mat, MultiMaterial):
            m = mat._mats.keys()[0]
        else:
            m = mat
        mat_name = "m{0}".format(m.metadata["mat_number"])
        m.metadata["name"] = mat_name
        uwuw_mats[mat_name] = m
           
    ml = MaterialLibrary(uwuw_mats) 
    ml.write_hdf5(out, datapath='/material_library/materials', 
                       nucpath='/material_library/nucid')

def main():
    description = ("This program read MCNP5 input files and creates UWUW material libraries"
                   " with naming schemes that match mcnp2cad with the -U flag.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("inp", action='store', help="Name of the MCNP input file.")
    parser.add_argument("-o", "--output", action='store', dest='out',
                       default="matlib.h5m", help="The name of the .h5m file to produce.")
    args = parser.parse_args()
    uwuw_matlib(args.inp, args.out)

if __name__ == '__main__':
    main()
