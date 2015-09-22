#! /usr/bin/env python
import unittest
import sys
import os
import os.path
from tempfile import NamedTemporaryFile as NTF
import ConfigParser
import btrn
import hzetrn as one_d_tool
import nose
import subprocess

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
   assert_almost_equal, assert_true, assert_false, assert_in


class TestBTRN(unittest.TestCase):
    rundir = './tdir' 
    resp_rundir   = '3_Response_functions/'
    resp_input    = 'resp_input.dat'


    def setUp(self):
        pass

    def tearDown(self):
        pass

    def verify_response(self, target):
        infile = self.rundir + '/' + self.resp_rundir + self.resp_input

	with open(infile) as f:
	    for i, line in enumerate(f):
	       if i == 0:
	           self.assertEqual("common_static_data",line.split()[0])
	       elif i == 1:
	           self.assertEqual("cross_sections_out",line.split()[0])
	       elif i == 2:
	           self.assertEqual("transport_data_out",line.split()[0])
	       elif i == 3:
	           self.assertEqual("response_data_out",line.split()[0])
	       elif i == 7:
	           self.assertEqual("flux.dat",line.split()[0])
	       elif i == 11:
	           self.assertEqual(target,line.split()[0])
        return

    def test_target_parameter(self):

#	subprocess.call(['python', '../btrn.py', '-f', './tgeom/brick.h5m', '-d', self.rundir, '-n 1'])
#	self.verify_response('WATERLIQ')

#	subprocess.call(['python', '../btrn.py', '-f', 'tgeom/brick.h5m', '-d', self.rundir, '-n 1', '-target', 'WATERLIQ'])
#	self.verify_response('WATERLIQ')

	subprocess.call(['python', '../btrn.py', '-f',  './tgeom/brick.h5m', '-d', self.rundir, '-n 1', '-target', 'ALUMINA1'])
	self.verify_response('ALUMINA1')


