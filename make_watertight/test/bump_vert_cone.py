#! /usr/bin/env python

from itaps import iMesh

datafile="cone.h5m"

mesh = iMesh.Mesh()
mesh.load(datafile)

print mesh.getNumOfTopo(0)
print mesh.getNumOfTopo(1)
print mesh.getNumOfTopo(2)
print mesh.getNumOfTopo(3)

vertices=mesh.getEntities(0)
j=0
for i in vertices:
   
	
	x, y, z = mesh.getVtxCoords(i)
	print "%f, %f, %f" % (x, y, z)
	if (x==0 and y==0 and j==0):
		coords=[x,y,z+1.e-1]
		mesh.setVtxCoords(i,coords)
		x, y, z = mesh.getVtxCoords(i)
		print "%f, %f, %f" % (x, y, z)
		j+=1
	else:
		pass
             
	
		
	
newfile="cone_mod.h5m"
mesh.save(newfile)

#for i in vertices:
 #   for j in i:
#    x, y, z = mesh.getVtxCoords(i)
#    print "%f, %f, %f" % (x, y, z)

 


