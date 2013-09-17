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
j=0
for i in vertices:
	if j == 0 :
		mesh.remove(i)
	else:   
		pass


newfile="cyl_mod.h5m"
mesh.save(newfile)

#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


