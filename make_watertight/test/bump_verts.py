#! /usr/bin/env python

import sys
from itaps import iMesh, iBase
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
j=0
bump_start=(mesh.getNumOfTopo(0)-201)

basename = args.fn.split('.')[0]
out_fn = basename + "_mod.h5m"

for i in vertices:
	if j > bump_start :
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		coords=[x,y,z+1e-01]

		mesh.setVtxCoords(i,coords)
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		mesh.save(out_fn)
		j+=1
	else:
		pass
		j+=1




#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


