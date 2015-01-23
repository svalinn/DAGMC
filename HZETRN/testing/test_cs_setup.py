#! /usr/bin/env python
import unittest
import sys
import os
import os.path
from tempfile import NamedTemporaryFile as NTF
import ConfigParser
import cs
import hzetrn as one_d_tool
import nose

# from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
#   assert_almost_equal, assert_true, assert_false, assert_in
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
   assert_almost_equal, assert_true, assert_false, assert_in


class TestLoadConfigs(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_load_configs(self):
        """ Simulate a .cfg file and check that everything is read correctly
	    jcz note:  this is like a protected 'cfgNTF = NTF()'
	"""
	with NTF() as cfgNTF:
	    # Create placeholder for hze.cfg file
	    cfgNTF.write("[common]\n" \
                         "run_directory = tdir\n" \
                         "common_data = common_static_data\n" \
			 "[cs]\n" \
                         "cross_dir = cross\n" \
                         "cs_out = cross_sections_out\n" \
			 )
	    cfgNTF.seek(0) # Goes to beginning

	    config = ConfigParser.SafeConfigParser()
	    config.read(cfgNTF.name)

        rundir, cs_common, cross_dir, cs_outdir = cs.load_config_params(config)
        # Check for correctness
	self.assertEqual(config.get('common', 'run_directory'), 'tdir')
	self.assertEqual(cs_common, 'common_static_data')
	self.assertEqual(config.get('cs', 'cross_dir'), 'cross')
	self.assertEqual(cs_outdir, 'cross_sections_out')

    def test_cs(self):
        header_lines = cs.xs_create_header('test_common_data', 'test_output')
	expect = "test_common_data                    # Name of folder where static data is stored\ntest_output                         # Name of folder to put data in (folder must already exist)\n\n------------------------------------------------------------------------------------\n\n"
	assert_equal(header_lines, expect)

        test_dir = 'test_dir'
        cur_path = os.path.dirname(os.path.abspath(__file__)) + '/'
        cross_path = cur_path + test_dir + '/'
	
	# mat_lib = 
        # cs_file_dict = write_cs_input(header_lines, mat_lib, cross_path) 

    
