import sys
import os

import pandas as pd
import numpy as np

from pymoab import core
from pymoab import types
from pymoab.tag import Tag



class h5mexplorer:

    def __init__(self, filename):

        self.mb = core.Core()
        self.h5m_file = filename
        self.mb.load_file(filename)

        self.root_set = self.mb.get_root_set()
        self.handles = self.mb.get_entities_by_handle(self.root_set)
        
        self.materials = BuildMaterialTable()
        self.volumes = BuildVolumeTable()


    def GetCATEGORYHandle(self, category):
        category_tg_name=self.mb.tag_get_handle(
            "CATEGORY", 32, types.MB_TYPE_OPAQUE, False)
        return self.mb.get_entities_by_type_and_tag(self.root_set, types.MBENTITYSET,
                category_tg_name, np.array([category]))


    def BuildMaterialTable(self):
        mat_handle = []
        mat_name = []
        mat_assignment = []
        mat_id = []

        group_handle = self.GetCATEGORYHandle("Group")
        name_tg_name = self.mb.tag_get_handle(
            "NAME", 32, types.MB_TYPE_OPAQUE, False)
        id_tg_name = self.mb.tag_get_handle(
            "GLOBAL_ID", 1, types.MB_TYPE_INTEGER, False)

        for handle in group_handle:
            try:
                name = self.mb.tag_get_data(name_tg_name, handle)[
                                            0][0].decode('utf-8')
                if name[0:3] == "mat":
                    mat_handle.append(handle)
                    mat_name.append(name)

                    assignent = []
                    for vol_hld in self.mb.get_entities_by_handle(handle):
                        assignent.append(vol_hld)
                    mat_assignment.append(assignent)

                    mat_id.append(self.mb.tag_get_data(id_tg_name, handle)[0][0])
            except RuntimeError:
                continue

        index = range(0, len(mat_handle))
        columns = ('ID', 'Handle', 'Name', 'Assignment')
        mat_pdf = pd.DataFrame( columns=columns)
        for i in range(len(mat_handle)):
            for assigment in mat_assignment[i]:
                row = [mat_id[i], mat_handle[i], mat_name[i], assigment]
                mat_pdf.loc[len(mat_pdf)] = row

        return mat_pdf

    def BuildVolumeTable(self):
        vol_id = []
        vol_mat_names = []
        vol_mat_handles = []
        vol_surface_handle = []
        
        vol_handle = self.GetCATEGORYHandle("Volume")
        
        tg_name = self.mb.tag_get_handle(
            "GLOBAL_ID", 1, types.MB_TYPE_INTEGER, False)
        
        for i in vol_handle:
            try:
                vol_id.append(self.mb.tag_get_data(tg_name, i)[0][0])
                
                vol_mat_handles.append(self.GetMaterialHandleInVolume(i))
                vol_mat_names.append(self.GetMaterialNameInVolume(i))
                
            except RuntimeError:
                continue
        
        index = range(0, len(vol_handle))
        columns = ('ID', 'Handle', 'MaterialName', 'MaterialHandle')
        vol_pdf = pd.DataFrame( columns=columns)
        for i in range(len(vol_handle)):
            for j in range(len(vol_mat_names[i])):
                row = [vol_id[i], vol_handle[i], vol_mat_names[i][j],
                        vol_mat_handles[i][j]]
                vol_pdf.loc[len(vol_pdf)] = row

        return vol_pdf


    def GetMaterialHandleInVolume(self, volume_handle):
        df_vol = self.materials[self.materials['Assignment'] == volume_handle]

        if len(df_vol) == 0:
            return [None]
        elif len(df_vol) == 1:
            return [df_vol['Handle']]
        else:
            print("!!WARNING!! Multiple material assignements for a single Volume")
            return df_vol['Handle'].tolist()
    
    
    def GetMaterialNameInVolume(self, volume_handle):
        df_vol = self.materials[self.materials['Assignment'] == volume_handle]

        if len(df_vol) == 0:
            return [None]
        elif len(df_vol) == 1:
            return [df_vol['Name']]
        else:
            print("!!WARNING!! Multiple material assignements for a single Volume")
            return df_vol['Name'].tolist()

