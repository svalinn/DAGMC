   
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

