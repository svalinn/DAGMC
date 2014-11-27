   
try:
    from itaps import iMesh, iBase

    HAVE_PYTAPS = True
except ImportError:

    HAVE_PYTAPS = False

try:
    from pyne.material import Material, MaterialLibrary
    HAVE_MATERIAL = True
except ImportError:
    HAVE_MATERIAL = False

try:
    from pyne.tally import Tally
    HAVE_TALLY = True
except ImportError:
    HAVE_TALLY = False

try:
    from pyne.particle import is_valid,name
    HAVE_PARTICLE = True
except ImportError:
    HAVE_PARTICLE = False


import string
import argparse
import os

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
function that gets all tags on dagmc geometry
------------------------------
path : the dagmc filename
return vector of tag_values
"""
def get_tag_values(path):
    ##########################################################
    dag_vol_names = [] # list of dag volume names (Cubit id)
    dag_properties = set()
    # material tags
    dag_material_tags = []
    # tally tags
    dag_tally_tags=[]
    # tally assignments
    tally_assigns=[]

    # material_vol 
    mat_assigns=[]

    # create imesh instance  
    dag_geom = iMesh.Mesh()
    # load the file
    dag_geom.load(path)

    # get all entities      
    ents = dag_geom.getEntities()
    # create a  mesh set         
    mesh_set = dag_geom.getEntSets()
    # list of volume ent handles 
    mat_list = []
    # get all geom_dimension ents
    geom_list = []
    cat_list  = []

    vol_tag  = dag_geom.getTagHandle('GEOM_DIMENSION')
    name_tag = dag_geom.getTagHandle('GLOBAL_ID')
    mat_tag  = dag_geom.getTagHandle('NAME')
    cat_tag  = dag_geom.getTagHandle('CATEGORY')

    # get the list we need    
    for i in mesh_set:
        tags = dag_geom.getAllTags(i)
        for tag in tags:
            if tag == vol_tag:
                geom_list.append(i)
            if tag == mat_tag:
                mat_list.append(i)
            if tag == cat_tag:
                cat_list.append(i)

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
                volume_name = name_tag[entity]
                # dag_materials[volume_name]="".join( chr( val ) for val in mat_name )
                # dag_properties.add("".join( chr( val ) for val in mat_name))
                dag_properties.add(tag_to_script(mat_name))

		if 'mat:' in tag_to_script(mat_name):
		    pair = (volume_name, tag_to_script(mat_name))
		    mat_assigns.append(pair)
		"""
                if 'tally:' in tag_to_script(mat_name):
                    pair = (volume_name,tag_to_script(mat_name))
                    tally_assigns.append(pair)
		"""

    # now we have dag properties, create one with materials and one with tallies
    for tag in dag_properties:
        if 'mat:' in tag:
            dag_material_tags.append(tag)
        if 'tally:' in tag:
            dag_tally_tags.append(tag)

    print ('#####################################################################')
    print dag_vol_names, dag_properties
    print ('#####################################################################')
    return dag_material_tags,mat_assigns
