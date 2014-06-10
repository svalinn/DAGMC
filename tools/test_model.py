#!/usr/bin/python

import os, subprocess
import nose
import dagmc_get_materials as gtag

#filename is the .h5m CAD file
filename = 'test.h5m'

"""
Function that loads the CAD h5m file and compares the group names with the obtained group names using get_tag_values on the dagmc_get_materials script 
"""
def test_model():
    # using subprocess to run mbsize on the h5m CAD file and assign the output to "out"
    proc = subprocess.Popen(["mbsize -ll %s | grep \"mat:\"" %filename], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    #the output is of <str> type and contains all group names.
    #split to seperate group names
    psuedo_list=out.split(' ')
    group_list=[]
    #loop over the pseudo_list to find group names and edit them to be in a proper format
    #group_list is the final group names list that is found on the CAD h5m file
    for k in psuedo_list:
        if "mat:" in k:
            group_list.append(k.rstrip("\n"))
    #tag_list is the list of group names obtained by running get_tag_values function
    tag_list=gtag.get_tag_values(filename)
    #comparing group names on both the tag_list "from dagmc_get_materials script" and mat_list "from CAM h5m file"
    for group in group_list:
        if group in tag_list:
            continue     
        else:
           raise ValueError('Mismatch between group names on the model and output of get_tag_values function!!')           
