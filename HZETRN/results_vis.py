import argparse
import logging
import subprocess
# Only used for get_rand_dirs
import tag_utils as tu
import numpy as np
import tempfile
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


num_non_data = 4
config_log  = 'config.log'


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
facets nx3x3 array of triangulated points on a sphere
"""
def create_meshed_sphere(facets, scale):
    msph = iMesh.Mesh()
 
    for facet in facets:
	verts = msph.createVtx(facet*scale)
	tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph

def create_tagged_meshed(data, facets, scale):
    """
    data   data indexed the same as the facets
    facets nx3x3 array of triangulated points on a sphere
    """
    msph = iMesh.Mesh()
 
    for facet in facets:
	# ToDo:  is it possible to createVtx with vtxs, i.e.
	#        an array of facets?
	verts = msph.createVtx(facet*scale)

	tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph

def mesh_sphere(xyz, scale):
    """
    xyz	 set of points on a sphere: nx3x3
    """
    facets = get_verts_on_sphere(xyz)
    return create_meshed_sphere(facets, scale)


"""
    Read a databasefile and strip of the header and 
    statistics lines
"""
def get_row_data(infile):
    global num_non_data 

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
	            print "Have reached an ascii line at i =", i

    # Determine how many columns of the header are non-numeric
    max_non_numeric_col = -1    
    for field in header:
        try:
	    num = float(field)
	    break
	except ValueError:
	    max_non_numeric_col = max_non_numeric_col + 1
	    # print field, max_non_numeric_col	    
    return header, max_non_numeric_col, rows
    
"""
    Read certain columns from the infile.
    Using direction vectors on each line of the file, and data
    from dat_column, create a tagged iMesh.Mesh that can be
    directly writen to a *.vtk file for viewing with VisIt.
    The direction vectors are expected to be unit length and will
    be scaled by scale.
"""
def tag_mesh(header, data_rows, columns, scale):
    "Create a scaled, tagged mesh from data in a file."

    vtxd = np.array(data_rows)
    print "rows shape", vtxd.shape
    print "rows", data_rows[0, 0:10]
   
    # Separate into vertexes (1st 3) and data
    vtcs = vtxd[:,:3]
    print "vtcs, row 0", vtcs[0,:]

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
    
def mesh_write_rand(filename, numOfPts):
    xyz = tu.get_rand_dirs(numOfPts)
    mesh = mesh_sphere(xyz, 1.0)
    mesh.save(filename)
    return 

def dmesh_write_from_db(data_rows, data_column, scale, outfile):
    mesh = tag_mesh(data_rows, data_column, scale)
    mesh.save(outfile)
    return

"""
Argument parsing
returns : args: -d for the run directory
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    global num_non_data
    # Change this to 4 if depth column is removed
    min_col = num_non_data
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='infile',
	help='The name of the .dat file with one line per direction.')

    parser.add_argument(
        '-c', action='store', dest='data_column', 
	help='The column, column range, or "a" for columns to use.  Must be {} or greater'.format(num_non_data))

    parser.add_argument(
        '-s', action='store', dest='scale', 
	help='Scale to use when creating the mesh.', type=float)
    parser.set_defaults(scale=1.0)

    parser.add_argument(
        '-n', action='store', dest='num_species',
	help='Number of species: 6 for spe, 59 for gcr')
    parser.set_defaults(num_species=6)

    parser.add_argument(
        '-o', action='store', dest='outfile',
	help='The name of the .vtk file, with extension.')

    args = parser.parse_args()
    if not args.infile:
        raise Exception('Input file not specified. [-f] not set.')
    """
    if not args.data_column:
        last_header_col = 2*args.num_species + num_non_data
	columns = str(num_non_data) + "-{}".format(last_header_col)
        parser.set_defaults(data_column=columns)
    """    
    if not args.outfile:
        raise Exception('Output file not specified.  [-o] not set.')
    
    return args

def main():
    global config_log
    global num_non_data

    args = parsing()

    ########################
    # logging
    ########################
    logging.basicConfig(filename=config_log, level=logging.INFO, \
                        format='%(asctime)s %(message)s')
    message = 'python results_vis.py -f ' + args.infile + \
              ' -c ' + str(args.data_column) + \
              ' -s ' + "{0:.2f}".format(args.scale) + \
              ' -o ' + args.outfile
    logging.info(message)
    ########################
    header, num_last_dose_col, row_data = get_row_data(args.infile)

    ########################
    # get columns
    ########################
    if not args.data_column \
        or args.data_column == 'h' \
	or args.data_column == 'header':
        columns = [num_non_data, num_last_dose_col]
    elif args.data_column == 'a' or args.data_column == 'all':
        columns = [num_non_data, row_data.shape[1]]
    elif len(args.data_column.split(',')) > 1:
        sys.exit("Please enter only one '-' separated field range")
    else:
        columns = map(int, args.data_column.split('-'))
	if len(columns) > 2:
	    sys.exit("For a range of columns please enter two integers separated by a '-'.")
        if min(columns) < num_non_data:
	    sys.exit("The minimum column must be {} or greater".format(num_non_data))
	if len(columns) == 2 and columns[1] < columns[0]:
	    columns = columns[::-1]
    print "columns", columns
    mesh = tag_mesh(header, row_data, columns, args.scale)
    mesh.save(args.outfile)

if __name__ == '__main__':
    main()
