import sys
import os

import pandas as pd
import pymoab

from pymoab import core
from pymoab import types
from pymoab.tag import Tag

import numpy as npi


class h5mexplorer:

    def __init__(self, filename):

        self.mb = core.Core()
        self.h5m_file = filename
        self.mb.load_file(filename)

        self.root_set = self.mb.get_root_set()
        self.handles = self.mb.get_entities_by_handle(self.root_set)
        self.materials = BuildMaterial()

    def BuildMaterial(self):
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
        for x in range(len(mat_handle)):
            for assigment in mat_assignment[x]:
                row = [mat_id[x], mat_handle[x], mat_name[x], assigment]
                mat_pdf.loc[len(mat_pdf)] = row


        return mat_pdf



    def GetCATEGORYHandle(self, category):
        category_tg_name=self.mb.tag_get_handle(
            "CATEGORY", 32, types.MB_TYPE_OPAQUE, False)
        return self.mb.get_entities_by_type_and_tag(self.root_set, types.MBENTITYSET,
                category_tg_name, np.array([category]))


    def GetMaterialNameInVolume(self, volume_handle):
        # Get assignment and Handle
        mat_assigment=GetMaterialAssignment(self)
        mat_handle=GetMaterialHandle(self)

        # Search for matching mat assigment
        mat_in_vol=[]
        for i, asgmts in enumerate(mat_assignment):
            if vol_handle in asgmts:
                mat_in_vol.append(mat_handle[i])

        if len(mat_in_vol) == 0:
            return NULL
        elif len(mat_in_vol) == 1:
            return mat_in_vol[0]
        else:
            print("!!WARNING!! Multiple material assignements for a single Volume")
            return mat_in_vol


    def GetVolumeId(self):
        rtn=[]
        volume_identities=[]
        vol_hdl=self.GetCATEGORYHandle("Volume")
        tg_name=self.mb.tag_get_handle(
            "GLOBAL_ID", 1, types.MB_TYPE_INTEGER, False)
        for i in vol_hdl:
            try:
                volume_identities.append(
                    self.mb.tag_get_data(tg_name, i)[0][0])
            except RuntimeError:
                continue
        return volume_identities



    def GetVolumeMaterialAssignment_handle(self):
        mat_assigment=GetMaterialAssignment(self)
        volume_handle=GetVolumeHandle(self)
        mat_handle=GetMaterialHandle(self)

        mat_in_vol=[]
        for vol_hdl in volume_handle:
            for i, asgmts in enumerate(mat_assignment):
                if vol_hdl in asgmts:
                    mat_in_vol.append(mat_handle[i])
        return mat_in_vol

    def GetVolumeMaterialAssignment_name(self):
        mat_assigment=GetMaterialAssignment(self)
        volume_handle=GetVolumeHandle(self)
        mat_name=GetMaterialList(self)

        mat_in_vol=[]
        for vol_hdl in volume_handle:
            for i, asgmts in enumerate(mat_assignment):
                if vol_hdl in asgmts:
                    mat_in_vol.append(mat_name[i])
        return mat_in_vol
