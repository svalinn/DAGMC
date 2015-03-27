import subprocess
import numpy as np

try:
    from scipy.spatial import ConvexHull
except ImportError:
    raise ImportError("scipy.spatial.ConvexHull could not be imported.")

try:
    from itaps import iMesh, iBase

    HAVE_PYTAPS = True
except ImportError:

    HAVE_PYTAPS = False

try:
   from pyne import material
   from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependencies could not be imported.")    


"""
Produce length-one tuples randomly oriented on a sphere
"""
def get_rand_dirs(number):
    rays = []
    if number > 0:
        for i in range(number):
            # [0.0,1.0)
            rnum = np.random.uniform()
            # map rnum onto [0.0,1.0)
            z = 2*rnum - 1

	    # Call a new randome number
            theta = 2*np.pi*np.random.uniform()
            norm_fac = np.sqrt(1 - z*z)
            y = norm_fac*np.sin(theta)
            x = norm_fac*np.cos(theta)
	    rays.append([x, y, z])
    return np.array(rays)

def get_verts_on_sphere(xyz):
    """
    Order the set of points as a triangular mesh 
    xyz	     nx3 np array with points on sphere
    return   triangular facets as sets of points
    """
    hull = ConvexHull(xyz)
    indices = hull.simplices
    return xyz[indices]

def create_meshed_sphere(facets, scale):
    """
    facets nx3x3 array of triangulated points on a sphere
    """
    msph = iMesh.Mesh()
 
    for facet in facets:
	# ToDo:  is it possible to createVtx with vtxs, i.e.
	#        an array of facets?
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


def mesh_write(filename, numOfPts):
    xyz = get_rand_dirs(numOfPts)
    mesh = mesh_sphere(xyz, 1.0)
    mesh.save(filename)
    return 

"""
    Read certain columns from the infile.
    Using direction vectors on each line of the file, and data
    from the given column, create a tagged iMesh.Mesh that can be
    directly writen to a *.vtk file for viewing with VisIt.
    The direction vectors are expected to be unit length and will
    be scaled by scale.
"""
def tag_mesh(infile, data_column, scale):
    "Create a scaled, tagged mesh from data in a file."
    # extract direction and one data column
    # us columns rather than fields to count.
    field_width = 13
    col_start = str((data_column-1)*field_width)
    col_end   = str(data_column*field_width)
    # Columns 1-39 contain the three direction values
    column_list = "1-39,{}-{}".format(col_start, col_end)
    # Call a linux command to get desired columns; capture output
    args = ["cut", "-c", column_list, infile]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    all_lines = output.splitlines()
    # Skip first (header) and last two (statistics) lines
    lines = all_lines[1:len(all_lines)-2]
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
    

"""
function to transform the tags into strings
tag : string of the tag to add to tag_list
tag_list : vector of tags in the problem
returns tag_list
"""
def _stringify(tag):
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii
            a.append(str(unichr(part)))
            # join to end string
            string = ''.join(a)
            # the the string we are testing for is not in the list of found
            # tag values, add to the list of tag_values

    return string

"""
function that gets all the material tags on dagmc geometry
------------------------------
path : the dagmc filename
return vector of tag_values
"""
def get_mat_tag_values(path):
    
    vol_mat_dict = {}

    # create imesh instance  and load the file
    dag_geom = iMesh.Mesh()
    dag_geom.load(path)

    # create a  mesh set         
    mesh_sets = dag_geom.getEntSets()

    cat_tag  = dag_geom.getTagHandle('CATEGORY')
    id_tag   = dag_geom.getTagHandle('GLOBAL_ID')
    name_tag = dag_geom.getTagHandle('NAME')

    # get the list we need    
    for mesh_set in mesh_sets:
        tags = dag_geom.getAllTags(mesh_set)
	if name_tag in tags and cat_tag in tags and _stringify(cat_tag[mesh_set])=="Group":
	    name = _stringify(name_tag[mesh_set])
	    # Now get the volumes of this group
	    child_sets = mesh_set.getEntSets()
	    for child in child_sets:
	        vol_id = id_tag[child]
	        if 'mat:' in name:
	            vol_mat_dict[vol_id] = name

	# else this is a meshset that doesn't have volume

    return vol_mat_dict

def get_fnames_for_vol(path):
    vol_fname_dict = {}
    vol_mat_dict = get_mat_tag_values(path)
    # Cross-Section: load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(path)
    
    for key in mat_lib:
        material_obj = mat_lib[key]
	mat_name = material_obj.metadata['name']
	for vol in vol_mat_dict:
	    if mat_name == vol_mat_dict[vol]:
	        vol_fname_dict[vol] = material_obj.metadata['fluka_name']
	        

    # Since 'graveyard' is not a 'fluka_name', add it specifically to the dictionary
    for vol in vol_mat_dict:
        if 'graveyard' in vol_mat_dict[vol].lower():
	   vol_fname_dict[vol] = 'graveyard' 
	   print 'Setting name for vol', vol, 'graveyard'
    return vol_fname_dict 

