#!/usr/bin/env python3

"""
The standard run script for the calibration of SSOC ROSA and Zyla
datasets.

-------------------------------------------------------------------------

Usage
-----

	standardCalScript.py <instrument name> <configuration file>

	instrument name : any of the following: ROSA_3500, ROSA_4170,
		ROSA_CAK, ROSA_GBAND, ZYLA.
	
	configuration file : path to an instrument configuration file.
		Consult the ssosoft documents for more information about
		this file.

-------------------------------------------------------------------------

This script completes all the steps necessary for an end-to-end
calibration for ROSA and Zyla datasets including speckle reconstruction
with the Kiepenheuer-Institut Speckle Interferometry Package (KISIP).
During the run, the RosaZylaCal class will use an instance of the logger
class to write to a log file named '<obsTime>_<instrument name>.log'.
In that log, the user will find informative messages regarding the
progress of the calibration.

-------------------------------------------------------------------------

Author
------

	G.A. MacDonald, NMSU, gordonm@nmsu.edu

-------------------------------------------------------------------------

Website
-------

	https://github.com/SSOCsoft/SSOsoft

-------------------------------------------------------------------------
"""

import ssosoft
import sys

assert len(sys.argv)==3, "Usage: {0} <instrument> <config file>".format(sys.argv[0])

r=ssosoft.rosaZylaCal(*sys.argv[1:])
r.rosa_zyla_run_calibration()
k=ssosoft.kisipWrapper(r)
k.kisip_despeckle_all_batches()
r.rosa_zyla_save_despeckled_as_fits()

