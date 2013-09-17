#! /usr/bin/env python

from itaps import iMesh, iBase

datafile="cyl.h5m"

mesh = iMesh.Mesh()
mesh.load(datafile)

print mesh.getNumOfTopo(0)
print mesh.getNumOfTopo(1)
print mesh.getNumOfTopo(2)
print mesh.getNumOfTopo(3)

vertices=mesh.getEntities(0)
j=0
bump_start=(mesh.getNumOfTopo(0)-10)
for i in vertices:
	if j > bump_start and j < (bump_start+4):
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		coords=[x,y,z]

		mesh.setVtxCoords(i,coords)
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		j+=1

        elif j == (bump_start+4):
                x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		coords=[x-4e-02,y+9e-02,z+1e-02]

		mesh.setVtxCoords(i,coords)
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		j+=1
        elif j == (bump_start+5):
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		coords=[x,y,z+1e-02]

		mesh.setVtxCoords(i,coords)
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		j+=1
		
	else:
		pass
		j+=1


newfile="cyl_mod.h5m"
mesh.save(newfile)
#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


