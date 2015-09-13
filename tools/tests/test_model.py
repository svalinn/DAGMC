#!/usr/bin/python

import os
import subprocess
import nose
import uwuw_preproc as gtag


filename = 'model.h5m'


def test_model_material():
    """
    Function that loads the model h5m file and compares the group names 
    with the obtained group names using get_tag_values on uwuw_preproc script 
    """
    output_filename = 'output.h5m'
    # using subprocess to run mbsize on the h5m CAD file and assign the output
    # to "out"
    proc = subprocess.Popen(
        ["mbsize -ll %s | grep \"mat:\"" % filename], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    # the output is of <str> type and contains all group names.
    # split to seperate group names
    pseudo_list = out.split('NAME = ')
    group_list = []
    # loop over the pseudo_list to find group names and edit them to be in a proper format
    # group_list is the final list of group names found on the CAD h5m file
    for k in pseudo_list:
        if ("mat:" in k):
            index = k.index('\n')
            k = k[0:index]
            group_list.append(k.rstrip("\n"))
    # tag_list is the list of group names obtained by running get_tag_values
    # function
    tag_list = gtag.get_tag_values(filename, output_filename)
    # comparing group names on both the tag_list "from dagmc_get_materials
    # script" and mat_list "from CAM h5m file"
    for group in group_list:
        if group in tag_list:
            continue
        else:
            raise ValueError('Mismatch in %s' % group)


def test_model_tally():
    output_filename = 'output.h5m'
    # using subprocess to run mbsize on the h5m CAD file and assign the output
    # to "out"
    proc = subprocess.Popen(
        ["mbsize -ll %s | grep \"tally:\"" % filename], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    # the output is of <str> type and contains all group names.
    # split to seperate group names
    pseudo_list = out.split('NAME = ')
    group_list = []
    # loop over the pseudo_list to find group names and edit them to be in a proper format
    # group_list is the final group names list that is found on the CAD h5m
    # file
    for k in pseudo_list:
        if ("tally:" in k):
            index = k.index('\n')
            k = k[0:index]
            group_list.append(k.rstrip("\n"))
    # tag_list is the list of group names obtained by running get_tag_values
    # function
    tag_list = gtag.get_tag_values(filename, output_filename)
    # comparing group names on both the tag_list "from dagmc_get_materials
    # script" and mat_list "from CAM h5m file"
    for group in group_list:
        if group in tag_list:
            continue
        else:
            raise ValueError('Mismatch in %s' % group)
