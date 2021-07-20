import sys
import os

import pandas as pd
import numpy as np

from pymoab import core
from pymoab import types
from pymoab.tag import Tag


class h5mexplorer:

    def __init__(self, filename):
        """ Constructor """

        self.mb = core.Core()
        self.h5m_file = filename
        self.mb.load_file(filename)

        self.root_set = self.mb.get_root_set()
        self.handles = self.mb.get_entities_by_handle(self.root_set)

        self.materials = BuildMaterialTable()
        self.volumes = BuildVolumeTable()

    def GetCATEGORYHandles(self, category):
        """
        Get the list of handle with the CATEGORY TAG matching the user
        spacification

        Parameters
        ----------
        category : string
            Value of the category tag to select

        Returns
        -------
            MOAB Range of EntityHandles

        """
        category_tg_name = self.mb.tag_get_handle(
            "CATEGORY", 32, types.MB_TYPE_OPAQUE)
        return self.mb.get_entities_by_type_and_tag(self.root_set, types.MBENTITYSET,
                                                    category_tg_name, np.array([category]))

    def BuildMaterialTable(self):
        """
        Build the Material Table containing the list of material and one line
        per assignment.

        Returns
        -------
            Pandas DataFrame
        """
        mat_handle = []
        mat_name = []
        mat_assignment = []
        mat_id = []

        group_handle = self.GetCATEGORYHandles("Group")
        name_tg_name = self.mb.tag_get_handle(
            "NAME", 32, types.MB_TYPE_OPAQUE)
        id_tg_name = self.mb.tag_get_handle(
            "GLOBAL_ID", 1, types.MB_TYPE_INTEGER)

        for handle in group_handle:
            try:
                name = self.mb.tag_get_data(name_tg_name, handle, flat=True)[0].decode('utf-8')
                if name[0:3] == "mat":
                    mat_handle.append(handle)
                    mat_name.append(name)

                    assignent = []
                    for vol_hld in self.mb.get_entities_by_handle(handle):
                        assignent.append(vol_hld)
                    mat_assignment.append(assignent)

                    mat_id.append(self.mb.tag_get_data(
                        id_tg_name, handle, flat=True)[0])
            except RuntimeError:
                continue

        index = range(0, len(mat_handle))
        columns = ('ID', 'Handle', 'Name', 'Assignment')
        mat_pdf = pd.DataFrame(columns=columns)
        for i in range(len(mat_handle)):
            for assigment in mat_assignment[i]:
                row = [mat_id[i], mat_handle[i], mat_name[i], assigment]
                mat_pdf.loc[len(mat_pdf)] = row

        return mat_pdf

    def BuildVolumeTable(self):
        """
        Build the Volume Table containing the list of volume, their ID, handle,
        material name and handle.

        Returns
        -------
        Pandas DataFrame
        """
        vol_id = []
        vol_mat_names = []
        vol_mat_handles = []
        vol_surface_handle = []

        vol_handle = self.GetCATEGORYHandles("Volume")

        tg_name = self.mb.tag_get_handle(
            "GLOBAL_ID", 1, types.MB_TYPE_INTEGER)

        for i in vol_handle:
            try:
                vol_id.append(self.mb.tag_get_data(tg_name, i, flat=True)[0])

                vol_mat_handles.append(self.GetMaterialHandleInVolume(i))
                vol_mat_names.append(self.GetMaterialNameInVolume(i))

            except RuntimeError:
                continue

        index = range(0, len(vol_handle))
        columns = ('ID', 'Handle', 'MaterialName', 'MaterialHandle')
        vol_pdf = pd.DataFrame(columns=columns)
        for i in range(len(vol_handle)):
            for j in range(len(vol_mat_names[i])):
                row = [vol_id[i], vol_handle[i], vol_mat_names[i][j],
                       vol_mat_handles[i][j]]
                vol_pdf.loc[len(vol_pdf)] = row

        return vol_pdf

    def GetMaterialHandleInVolume(self, volume_handle):
        """
        Retrieve the Material Handle for a volume

        Parameters
        ----------
        volume_handle : MOAB EntityHandle

        Returns
        -------
        list of material handles
        """
        df_vol = self.materials[self.materials['Assignment'] == volume_handle]

        if len(df_vol) == 0:
            return [None]
        elif len(df_vol) > 1:
            print("!!WARNING!! Multiple material assignements for a single Volume")
        return df_vol['Handle'].tolist()

    def GetMaterialNameInVolume(self, volume_handle):
        """
        Retrieve the Material Name for a Volume

        Parameters
        ----------
        volume_handle : MOAB EntityHandle (long)

        Returns
        -------
        list of material names
        """
        df_vol = self.materials[self.materials['Assignment'] == volume_handle]

        if len(df_vol) == 0:
            return [None]
        elif len(df_vol) > 1:
            print("!!WARNING!! Multiple material assignements for a single Volume")
        return df_vol['Name'].tolist()
