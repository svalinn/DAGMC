import argparse
import subprocess
# Only used for get_rand_dirs
import btrn
import numpy as np
# import tempfile
import os
import sys

try:
    from scipy.spatial import ConvexHull
except ImportError:
    raise ImportError("scipy.spatial.ConvexHull could not be imported.")

try:
    from itaps import iMesh, iBase
except ImportError:
    raise ImportError("Could not import iMesh and/or iBase.")


# Number of columns devoted to direction and depth
num_meta = 6


############################################################3
# Functions used for testing and development
"""
Order the set of points as a triangular mesh 
xyz	 nx3 np array with points on sphere
return   triangular facets as sets of points
"""
def get_verts_on_sphere(xyz):
    hull = ConvexHull(xyz)
    indices = hull.simplices
    return xyz[indices]

"""
facets an nx3x3 array of triangulated points on a sphere
"""
def create_meshed_sphere(facets, scale):
    msph = iMesh.Mesh()
 
    for facet in facets:
	verts = msph.createVtx(facet*scale)
	tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph

"""
    xyz	 set of points on a sphere: nx3x3
"""
def mesh_sphere(xyz, scale):
    facets = get_verts_on_sphere(xyz)
    return create_meshed_sphere(facets, scale)

"""
    Create and write out a mesh constructed from random points 
    on a sphere
"""
def mesh_write_rand(filename, numOfPts):
    xyz = btrn.get_rand_dirs(numOfPts)
    mesh = mesh_sphere(xyz, 1.0)
    mesh.save(filename)
    return 
# end optional functions
############################################################3

"""
    data   data indexed the same as the facets
    facets nx3x3 array of triangulated points on a sphere
"""
"""
def create_tagged_meshed(data, facets, scale):
    msph = iMesh.Mesh()
 
    for facet in facets:
	# ToDo:  is it possible to createVtx with vtxs, i.e.
	#        an array of facets?
	verts = msph.createVtx(facet*scale)

	tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph
"""

"""
    Read the lines of a  database file, stripping 
    off the header and statistics lines
"""
def get_row_data(infile):

    rows = np.array([])
    with open(infile) as fp:
        for index, line in enumerate(fp):
	    i = index + 1
	    if i == 1:
	        header = line.split()
	    else:
	        try: 
	            row = map(float, line.split())
	            if len(rows) == 0:
	                rows = row
	            else:
	                rows = np.vstack((rows, row))
	        except ValueError:
		    pass
	           
    # Determine how many columns of the header are non-numeric
    max_non_numeric_col = -1    
    for field in header:
        try:
	    num = float(field)
	    break
	except ValueError:
	    max_non_numeric_col = max_non_numeric_col + 1
    return header, max_non_numeric_col, rows
    
"""
    Read the indicated columns from the infile.
    Using the direction vectors on each line of the file, and data
    from the selected columns, create a tagged iMesh.Mesh that can be
    directly writen to a *.vtk file for viewing with VisIt.
    The direction vectors are expected to be unit length and will
    be scaled by scale.
"""
def tag_mesh(header, data_rows, columns, scale):
    "Create a scaled, tagged mesh from data in a file."

    vtxd = np.array(data_rows)
   
    # Separate into vertexes (1st 3) and data
    vtcs = vtxd[:,:3]

    start = columns[0]
    end = 0
    if len(columns) == 1:
	end = start + 1
    else:
        end = columns[1] + 1

    # Triangulate 
    hull       = ConvexHull(vtcs)
    indices    = hull.simplices
    facets     = vtcs[indices]
    num_facets = facets.shape[0]

    # Create the mesh, and a data tag on it
    msph = iMesh.Mesh()

    tag_map = {}

    for col in range(start, end): 
        # get the desired column of data
	data = vtxd[:,col]
	# Reorder the data column in terms of vertices
        di = data[indices]
        # creat a header-field-based name for the current column
        data_name = header[col]+"_{}".format(col)
	# Associate the hashable Tag array object with the ascii data_name
        tag_map[data_name] = msph.createTag(data_name, 1, float)

        # Use indexing to get matching data at each vertex
        for i in range(num_facets):
            facet        = facets[i]
	    data_at_vtcs = di[i]
	    # Create an entity handle for each point in the facet
	    verts = msph.createVtx(facet*scale)
	    # Tag each entity handle (lhs key) with 
	    # corresponding data value (rhs)
	    tag_map[data_name][verts[0]] = data_at_vtcs[0]
	    tag_map[data_name][verts[1]] = data_at_vtcs[1]
	    tag_map[data_name][verts[2]] = data_at_vtcs[2]

	    # Tell the mesh about this triangular facet
	    tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph

"""
Argument parsing
returns : args: -d for the run directory
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='infile',
	help='The name of the .dat file with one line per direction.')

    parser.add_argument(
        '-c', action='store', dest='data_column', 
	help='The column, column range, or "a" for all columns.  Must be {} or greater'.format(num_meta))

    parser.add_argument(
        '-s', action='store', dest='scale', 
	help='Scale to use when creating the mesh.', type=float)
    parser.set_defaults(scale=1.0)

    parser.add_argument(
        '-o', action='store', dest='outfile',
	help='The name of the .vtk file, with extension.')

    args = parser.parse_args()
    if not args.infile:
        raise Exception('Input file not specified. [-f] not set.')

    if not args.outfile:
        raise Exception('Output file not specified.  [-o] not set.')
    
    return args

def main():

    args = parsing()
    if not args.data_column:
        args.data_column = 'header'

    ########################
    # Get the header, the highest-numbered column with a non-numeric
    # header, and the floating point row data
    header, num_last_dose_col, row_data = get_row_data(args.infile)

    ########################
    # get columns
    ########################
    if args.data_column == 'h' or args.data_column == 'header':
        columns = [num_meta, num_last_dose_col]
	print "Header columns being used:", columns
    elif args.data_column == 'a' or args.data_column == 'all':
        columns = [num_meta, row_data.shape[1]]
	print "All columns being used:", columns
    elif len(args.data_column.split(',')) > 1:
        sys.exit("Please enter only one '-' separated field range")
    else:
        columns = map(int, args.data_column.split('-'))
	if len(columns) > 2:
	    sys.exit("For a range of columns please enter two integers separated by a '-'.")
        if min(columns) < num_meta:
	    sys.exit("The minimum column must be {} or greater".format(num_meta))
	# Reverse the order if needed
	if len(columns) == 2 and columns[1] < columns[0]:
	    columns = columns[::-1]
        print "Column range:", columns

    mesh = tag_mesh(header, row_data, columns, args.scale)
    mesh.save(args.outfile)

if __name__ == '__main__':
    main()
