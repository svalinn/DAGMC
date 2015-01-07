#! /usr/bin/env python
import unittest
import sys
import os
import os.path
from tempfile import NamedTemporaryFile as NTF
import ConfigParser
import cs
import hzetrn as one_d_tool

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

