   
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
def tag_to_script(tag):
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii
            a.append(str(unichr(part)))
            # join to end string
            test = ''.join(a)
            # the the string we are testing for is not in the list of found
            # tag values, add to the list of tag_values

    return test

"""
function that gets all the material tags on dagmc geometry
------------------------------
path : the dagmc filename
return vector of tag_values
"""
def get_mat_tag_values(path):
    # dag_properties = set()
    # material tags
    # dag_material_tags = []

    # material_vol 
    # mat_assigns=[]
    # mat_vol_dict = {}
    vol_mat_dict = {}

    # create imesh instance  and load the file
    dag_geom = iMesh.Mesh()
    dag_geom.load(path)

    # get all entities      
    ents = dag_geom.getEntities()
    # create a  mesh set         
    mesh_set = dag_geom.getEntSets()
    # list of volume ent handles and geom_dimension ents
    mat_list = []
    geom_list = []

    vol_tag  = dag_geom.getTagHandle('GEOM_DIMENSION')
    name_tag = dag_geom.getTagHandle('GLOBAL_ID')
    mat_tag  = dag_geom.getTagHandle('NAME')

    # get the list we need    
    for i in mesh_set:
        tags = dag_geom.getAllTags(i)
        for tag in tags:
            if tag == vol_tag:
                geom_list.append(i)
            if tag == mat_tag:
                mat_list.append(i)


    found_meshset_for_vol = False
    # loop over all the volumes 
    for entity in geom_list:
        # loop over the material sets
        for meshset in mat_list:
	    # reset
            found_meshset_for_vol = False
            # if volume in set       
            if meshset.contains(entity):
	        if found_meshset_for_vol:
		    # ToDo: throw exception:, should not be here again for same vol
	            print 'Error: a volume cannot belong to more than one material'
		# Expect a volume to be in only one mat_list set
		found_meshset_for_vol = True
                mat_name = mat_tag[meshset]
		# print ('in entity-meshset loop: ', mat_tag, meshset, mat_tag[meshset])
                volume_name = name_tag[entity]
                # dag_properties.add(tag_to_script(mat_name))

		# uwuw_preproc's version, get_tag_values(), uses 'tally:' here
		if 'mat:' in tag_to_script(mat_name):
		    # pair = (volume_name, tag_to_script(mat_name))
		    # mat_assigns.append(pair)
                    # mat_vol_dict[tag_to_script(mat_name)] = volume_name
		    vol_mat_dict[volume_name] = tag_to_script(mat_name)
		  

    # now we have dag properties
    # for tag in dag_properties:
    #     if 'mat:' in tag:
    #        dag_material_tags.append(tag)

    # a dictiounary will be more convenient
    # mat_vol_dict = {}
    # pair 0 is the vol, par 1 is the material tag name
    # for pair in mat_assigns:
    #    mat_vol_dict[pair[1]]=pair[0]

    # return mat_vol_dict
    return vol_mat_dict

def get_fnames_for_vol(path):
    vol_fname_dict = {}
    vol_mat_dict = get_mat_tag_values(path)
    print 'materials and volumes', vol_mat_dict
    # Cross-Section: load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(path)
    # cs_file_mats = {}
    # num_materials = len(mat_lib.keys())
    for key in mat_lib:
        material_obj = mat_lib[key]
	mat_name = material_obj.metadata['name']
	for vol in vol_mat_dict:
	    if mat_name == vol_mat_dict[vol]:
	        vol_fname_dict[vol] = material_obj.metadata['fluka_name']
	        

	# This is the whole point of this method:  join two tables
	# name       = material_obj.metadata['name']
	#vol = mat_vol_dict[material_obj.metadata['name']]
	
        # fluka_name = material_obj.metadata['fluka_name']
        # Associate the vol and the fluka name	
	#vol_fname_dict[vol] = material_obj.metadata['fluka_name']

    # Since 'graveyard' is not a 'fluka_name', add it specifically to the dictionary
    # for mat in mat_vol_dict:
    for vol in vol_mat_dict:
        if 'graveyard' in vol_mat_dict[vol].lower():
	   vol_fname_dict[vol] = 'graveyard' 
	   print 'Setting name for vol', vol, 'graveyard'
    return vol_fname_dict 

