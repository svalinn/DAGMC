import os, subprocess
import nose
import dagmc_get_materials as gtag

filename = 'test.h5m'
def test_model():
    proc = subprocess.Popen(["mbsize -ll %s | grep \"mat:\"" %filename], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    psuedo_list=out.split(' ')
    group_list=[]
    for k in psuedo_list:
        if "mat:" in k:
            group_list.append(k.rstrip("\n"))
   
    tag_list=gtag.get_tag_values(filename)
    for group in group_list:
        if group in tag_list:
            continue     
        else:
           raise ValueError('Mismatch between group names on the model and output of get_tag_values function!!')           
