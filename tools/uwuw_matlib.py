#!/usr/bin/env python
#
# add option to make more compatible with uwuw_prepoc
#     -use default nucpath,datapath
#     -use material names like 11 instead of m11
#
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

def uwuw_matlib(inp, out, enhanced_uwuw):
    mats = mats_from_inp(inp)
    uwuw_mats = {}
    for mat in mats.values():
        if isinstance(mat, MultiMaterial):
            m = mat._mats.keys()[0]
        else:
            m = mat
        if enhanced_uwuw:
            print "Using enhanced uwuw_preproc compatibility mat_name...(no leading m)"
            mat_name = str(m.metadata["mat_number"])
        else:
            mat_name = "m{0}".format(m.metadata["mat_number"]) # adds "m" to mat_name
        m.metadata["name"] = mat_name
        uwuw_mats[mat_name] = m
           
    ml = MaterialLibrary(uwuw_mats)
    if enhanced_uwuw:
        print "\n Setting enhanced uwuw_preproc compatibility for writing material library..."
        ml.write_hdf5(out) # don't set datapath,nucpath...will be pyne default values
    else:
        ml.write_hdf5(out, datapath='/material_library/materials', 
                       nucpath='/material_library/nucid')

def main():
    description = ("This program read MCNP5 input files and creates UWUW material libraries"
                   " with naming schemes that match mcnp2cad with the -U flag.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-u","--uwuw_preproc", help="produce library more compatible with uwuw_preproc",action="store_true")
    parser.add_argument("inp", action='store', help="Name of the MCNP input file.")
    parser.add_argument("-o", "--output", action='store', dest='out',
                       default="matlib.h5m", help="The name of the .h5m file to produce.")
    args = parser.parse_args()
    uwuw_matlib(args.inp, args.out, args.uwuw_preproc)

if __name__ == '__main__':
    main()
