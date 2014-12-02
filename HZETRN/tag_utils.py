   
try:
    from itaps import iMesh, iBase

    HAVE_PYTAPS = True
except ImportError:

    HAVE_PYTAPS = False

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
    dag_vol_names = [] # list of dag volume names (Cubit id)
    dag_properties = set()
    # material tags
    dag_material_tags = []

    # material_vol 
    mat_assigns=[]

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

    # for the 3d entities     
    for entity in geom_list:
        if vol_tag[entity] == 3:
             dag_vol_names.append(str(name_tag[entity]))

    # loop over all the volumes 
    for entity in geom_list:
        # loop over the material sets
        for meshset in mat_list:
            # if volume in set       
            if meshset.contains(entity):
                mat_name = mat_tag[meshset]
		# print ('in entity-meshset loop: ', mat_tag, meshset, mat_tag[meshset])
                volume_name = name_tag[entity]
                dag_properties.add(tag_to_script(mat_name))

		# uwuw_preproc's version, get_tag_values(), uses 'tally:' here
		if 'mat:' in tag_to_script(mat_name):
		    pair = (volume_name, tag_to_script(mat_name))
		    mat_assigns.append(pair)

    # now we have dag properties
    for tag in dag_properties:
        if 'mat:' in tag:
            dag_material_tags.append(tag)
    # a dictiounary will be more convenient
    vol_mat_dict = {}
    for pair in mat_assigns:
        vol_mat_dict[pair[0]]=pair[1]

    # return dag_material_tags, vol_mat_dict
    return vol_mat_dict

