#! /usr/bin/env python

from itaps import iMesh

datafile="cyl.h5m"

mesh = iMesh.Mesh()
mesh.load(datafile)

print mesh.getNumOfTopo(0)
print mesh.getNumOfTopo(1)
print mesh.getNumOfTopo(2)
print mesh.getNumOfTopo(3)

vertices=mesh.getEntities(0)

for i in vertices:
    pass

x, y, z = mesh.getVtxCoords(i)
print "%f, %f, %f" % (x, y, z)
coords=[x,y,z+1.e-1]

mesh.setVtxCoords(i,coords)
x, y, z = mesh.getVtxCoords(i)
print "%f, %f, %f" % (x, y, z)

newfile="cyl_mod.h5m"
mesh.save(newfile)

#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


