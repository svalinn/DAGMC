import argparse
import logging
import subprocess
import tag_utils as tu
import numpy as np
try:
    from scipy.spatial import ConvexHull
except ImportError:
    raise ImportError("scipy.spatial.ConvexHull could not be imported.")

try:
    from itaps import iMesh, iBase
except ImportError:
    raise ImportError("Could not import iMesh and/or iBase.")


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


def get_cols_from_file(infile, data_column):
    "Create a scaled, tagged mesh from data in a file."
    # extract direction and one data column
    # us columns rather than fields to count.
    field_width = 13
    col_start = str((data_column-1)*field_width)
    col_end   = str(data_column*field_width)
    # Columns 1-39 contain the three direction values
    column_list = "1-39,{}-{}".format(col_start, col_end)
    print "cut column values", column_list
    # Call a linux command to get desired columns; capture output
    args = ["cut", "-c", column_list, infile]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    all_lines = output.splitlines()
    for line in all_lines:
        print line
    # Skip first (header) and last two (statistics) lines
    lines = all_lines[1:len(all_lines)-2]
    return lines

"""
    Read certain columns from the infile.
    Using direction vectors on each line of the file, and data
    from dat_column, create a tagged iMesh.Mesh that can be
    directly writen to a *.vtk file for viewing with VisIt.
    The direction vectors are expected to be unit length and will
    be scaled by scale.
"""
def tag_mesh(infile, data_column, scale):
    "Create a scaled, tagged mesh from data in a file."

    lines = get_cols_from_file(infile, data_column)
    vals = []
    for line in lines:
        vals.append(map(float, line.split()))
    vtxd = np.array(vals)
    
    # Separate into vertexes (1st 3) and data
    vtcs = vtxd[:,:3]
    data = vtxd[:,3]

    # Triangulate 
    hull       = ConvexHull(vtcs)
    indices    = hull.simplices
    facets     = vtcs[indices]
    di         = data[indices]
    num_facets = facets.shape[0]

    # Create the mesh, and a data tag on it
    msph = iMesh.Mesh()
    point_data = msph.createTag("point_data", 1, float)

    # Use indexing to get matching data at each vertex
    for i in range(num_facets):
        facet        = facets[i]
	data_at_vtcs = di[i]
	# Create an entity handle for each point in the facet
	verts = msph.createVtx(facet*scale)
	# Tag each entity handle (lhs key) with 
	# corresponding data value (rhs)
	point_data[verts[0]] = data_at_vtcs[0]
	point_data[verts[1]] = data_at_vtcs[1]
	point_data[verts[2]] = data_at_vtcs[2]

	# Tell the mesh about this triangular facet
	tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
    return msph
    
def mesh_write_rand(filename, numOfPts):
    xyz = tu.get_rand_dirs(numOfPts)
    mesh = mesh_sphere(xyz, 1.0)
    mesh.save(filename)
    return 

def dmesh_write_from_db(infile, data_column, scale, outfile):
    mesh = tag_mesh(infile, data_column, scale)
    mesh.save(outfile)
    return

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
	help='The column of data to use.  Must be greater than 3', type=int)

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
    if not args.data_column:
        raise Exception('Data column not specified. [-c] not set.')
    elif args.data_column < 4:
        raise Exception('Data column specified must be greater than 3.')
    if not args.outfile:
        raise Exception('Output file not specified.  [-o] not set.')
    
    return args

def main():
    global config_log

    args = parsing()

    logging.basicConfig(filename=config_log, level=logging.INFO, format='%(asctime)s %(message)s')
    message = 'python results_vis.py -f ' + args.infile + \
              ' -c ' + str(args.data_column) + \
              ' -s ' + ('not given - default=1.0' if not args.scale else "{0:.2f}".format(args.scale)) + \
              ' -o ' + args.outfile
    logging.info(message)

    dmesh_write_from_db(args.infile, args.data_column, args.scale, args.outfile)

if __name__ == '__main__':
    main()
