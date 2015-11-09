import sys 
import os
from pyne.material import MaterialLibrary
from nose.tools import assert_equal, assert_almost_equal

sys.path.append("..")
from uwuw_matlib import uwuw_matlib

def test_uwuw_matlib():
    inp = "test.inp"
    out = "out.h5m"
    if os.path.exists(out):
        os.remove(out)
    uwuw_matlib(inp, out)
    ml = MaterialLibrary()
    ml.from_hdf5(out, datapath='/material_library/materials',
                      nucpath='/material_library/nucid')
    mats = ml.values()
    assert_equal(len(mats), 2)
    
    assert_equal(mats[0].metadata["mat_number"], "1")
    assert_equal(mats[0].metadata["name"], "m1")
    assert_almost_equal(mats[0].density, 9.0)
    assert_equal(mats[0].comp, {10010000: 1.0})

    assert_equal(mats[1].metadata["mat_number"], "2")
    assert_equal(mats[1].metadata["name"], "m2")
    assert_almost_equal(mats[1].density, 21.0)
    assert_equal(mats[1].comp, {10020000: 1.0})
    os.remove(out)

