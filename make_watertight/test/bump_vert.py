#! /usr/bin/env python

import sys
from itaps import iMesh
import argparse as ap

parser = ap.ArgumentParser(description = "Short script for bumping a single vertex out of position in a mesh file.")

parser.add_argument('--filename', type = str, dest = 'fn', help = "Mesh file to modify.")

args = parser.parse_args(sys.argv[1:])
mesh = iMesh.Mesh()
mesh.load(args.fn)

print mesh.getNumOfTopo(0)
print mesh.getNumOfTopo(1)
print mesh.getNumOfTopo(2)
print mesh.getNumOfTopo(3)

vertices=mesh.getEntities(0)

for i in vertices:
    pass

x, y, z = mesh.getVtxCoords(i)
print "%f, %f, %f" % (x, y, z)
coords=[x,y,z+1.e-2]

mesh.setVtxCoords(i,coords)
x, y, z = mesh.getVtxCoords(i)
print "%f, %f, %f" % (x, y, z)

basename = args.fn.split('.')[0]
out_fn = basename + "_mod.h5m"
mesh.save(out_fn)

#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


