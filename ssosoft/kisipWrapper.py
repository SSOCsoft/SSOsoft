"""
This is a wrapper library for running KISIP with
regular Rosa/Hardcam data reduction runs.

Initialize with an instance of RosaHardCamCal class.
"""
import configparser
import logging, logging.config
import os

class kisipWrapper:

	from . import ssosoftConfig

	def __init__(self, RHCC):
		self.burstFileForm=RHCC.burstFileForm
		self.burstNumber=RHCC.burstNumber
		self.configFile=RHCC.configFile
		self.imageShape=RHCC.imageShape
		self.instrument=RHCC.instrument.upper()
		self.kisipPreSpeckleStartInd=0
		self.kisipPreSpeckleEndInd=0
		self.logFile=RHCC.logFile
		self.noiseFile=RHCC.noiseFile
		self.obsDate=RHCC.obsDate
		self.obsTime=RHCC.obsTime
		self.preSpeckleBase=RHCC.preSpeckleBase
		self.workBase=RHCC.workBase
	
	def kisip_configure_run(self):
		config=configparser.ConfigParser()
		config.read(self.configFile)
	
		self.kisipArcsecPerPixX=config[self.instrument]['kisipArcsecPerPixX']
		self.kisipArcsecPerPixY=config[self.instrument]['kisipArcsecPerPixY']
		self.speckledFileForm=config[self.instrument]['speckledFileForm']
		self.wavelengthnm=config[self.instrument]['wavelengthnm']
	
		self.speckleBase=os.path.join(self.workBase, 'speckle')
	
		self.kisipMethodMethod=config['KISIP_METHOD']['kisipMethodMethod']
		self.kisipMethodSubfieldArcsec=config['KISIP_METHOD']['kisipMethodSubfieldArcsec']
		self.kisipMethodPhaseRecLimit=config['KISIP_METHOD']['kisipMethodPhaseRecLimit']
		self.kisipMethodUX=config['KISIP_METHOD']['kisipMethodUX']
		self.kisipMethodUV=config['KISIP_METHOD']['kisipMethodUV']
		self.kisipMethodMaxIter=config['KISIP_METHOD']['kisipMethodMaxIter']
		self.kisipMethodSNThresh=config['KISIP_METHOD']['kisipMethodSNThresh']
		self.kisipMethodWeightExp=config['KISIP_METHOD']['kisipMethodWeightExp']
		self.kisipMethodPhaseRecApod=config['KISIP_METHOD']['kisipMethodPhaseRecApod']
		self.kisipMethodNoiseFilter=config['KISIP_METHOD']['kisipMethodNoiseFilter']
		self.kisipPropsHeaderOff=config['KISIP_PROPS']['kisipPropsHeaderOff']
		self.kisipPropsTelescopeDiamm=config['KISIP_PROPS']['kisipPropsTelescopeDiamm']
		self.kisipPropsAoLockX=config['KISIP_PROPS']['kisipPropsAoLockX']
		self.kisipPropsAoLockY=config['KISIP_PROPS']['kisipPropsAoLockY']
		self.kisipPropsAoUsed=config['KISIP_PROPS']['kisipPropsAoUsed']
		self.kisipEnvLib=config['KISIP_ENV']['kisipEnvLib']
		self.kisipEnvBin=config['KISIP_ENV']['kisipEnvBin']
		self.kisipEnvMpiNproc=config['KISIP_ENV']['kisipMpiNproc']
	
		## Set-up logging.
		logging.config.fileConfig(self.configFile,
				defaults={'logfilename': self.logFile}
				)
		self.logger=logging.getLogger('kisipLog')
	
		self.logger.info("This is kisipWrapper, part of SSOsoft "
				"version {0}".format(self.ssosoftConfig.__version__)
				)
		self.logger.info("Contact {0} to report bugs, make suggestions, "
				"or contribute.".format(self.ssosoftConfig.__email__)
				)
		self.logger.info(
				"kisipWrapper is designed to be used in conjunction with "
				"the SSOsoft RosaHardCamCal class. Any other way of using "
				"this package is at the user's risk."
				)
		self.logger.info("Now configuring this KISIP run.")
	
		## Directories preSpeckleBase and speckleBase mut exist or be
		## created.
		if not os.path.isdir(self.preSpeckleBase):
			self.logger.info("os.mkdir: attempting to create directory: "
					"{0}".format(self.preSpeckleBase)
					)
			try:
				os.mkdir(self.preSpeckleBase)
			except Exception as err:
				self.logger.critical("CRITICAL: {0}".format(err))
				raise
		if not os.path.isdir(self.speckleBase):
			self.logger.info("os.mkdir: attempting to create directory: "
					"{0}".format(self.speckleBase)
					)
			try:
				os.mkdir(self.speckleBase)
			except Exception as err:
				self.logger.critical("CRITICAL: {0}".format(err))
				raise
	
	def kisip_set_environment(self):
		self.logger.info("Setting system environment variables for KISIP.")
		self.logger.info("Pre-appending to PATH: {0}".format(self.kisipEnvBin))
		os.environ['PATH']="{0}{1}{2}".format(
				self.kisipEnvBin, os.pathsep, os.environ['PATH'])
		self.logger.info("Pre-appending to LD_LIBRARY_PATH: {0}".format(self.kisipEnvLib))
		os.environ['LD_LIBRARY_PATH']="{0}{1}{2}".format(
				self.kisipEnvLib, os.pathsep, os.environ['PATH']
				)

	def kisip_spawn_kisip(self):
		pass

	def kisip_write_init_files(self):
		self.logger.info("Preparing to write KISIP init files.")
		self.logger.info("Writing KISIP config file: "
				"{0}".format(os.path.join(self.workBase, 'init_file.dat'))
				)
		try:
			with open(os.path.join(self.workBase, 'init_file.dat'), mode='wt') as f:
				f.write("{0}{1}".format(
					os.path.join(self.preSpeckleBase, self.burstFileForm),
					os.linesep
					)
					)
				f.write("{0}{1}".format(self.kisipPreSpeckleStartInd, os.linesep))
				f.write("{0}{1}".format(self.kisipPreSpeckleEndInd, os.linesep))
				f.write("{0}{1}".format(
					os.path.join(self.speckleBase, self.speckledFileForm),
					os.linesep
					)
					)
				f.write("{0}{1}".format(
					os.path.join(self.preSpeckleBase, self.noiseFile),
					os.linesep
					)
					)
		except Exception as err:
			self.logger.critical("CRITICAL: {0}".format(err))
			raise

		self.logger.info("Writing KISIP config file: "
				"{0}".format(os.path.join(self.workBase, 'init_method.dat'))
				)
		try:
			with open(os.path.join(self.workBase, 'init_method.dat'), mode='wt') as f:
				f.write("{0}{1}".format(self.kisipMethodMethod, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodSubfieldArcsec, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodPhaseRecLimit, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodUX, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodUV, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodMaxIter, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodSNThresh, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodWeightExp, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodPhaseRecApod, os.linesep))
				f.write("{0}{1}".format(self.kisipMethodNoiseFilter, os.linesep))
		except Exception as err:
			self.logger.critical("CRITICAL: {0}".format(err))
			raise

		self.logger.info("Writing KISIP config file: "
				"{0}".format(os.path.join(self.workBase, 'init_props.dat'))
				)
		try:
			with open(os.path.join(self.workBase, 'init_props.dat'), mode='wt') as f:
				f.write("{0}{1}".format(self.imageShape[1], os.linesep))
				f.write("{0}{1}".format(self.imageShape[0], os.linesep))
				f.write("{0}{1}".format(self.burstNumber, os.linesep))
				f.write("{0}{1}".format(self.kisipPropsHeaderOff, os.linesep))
				f.write("{0}{1}".format(self.kisipArcsecPerPixX, os.linesep))
				f.write("{0}{1}".format(self.kisipArcsecPerPixY, os.linesep))
				f.write("{0}{1}".format(self.kisipPropsTelescopeDiamm, os.linesep))
				f.write("{0}{1}".format(self.wavelengthnm, os.linesep))
				f.write("{0}{1}".format(self.kisipPropsAoLockX, os.linesep))
				f.write("{0}{1}".format(self.kisipPropsAoLockY, os.linesep))
				f.write("{0}{1}".format(self.kisipPropsAoUsed, os.linesep))
		except Exception as err:
			self.logger.critical("CRITICAL: {0}".format(err))
			raise

		self.logger.info("Successfully wrote KISIP init files.")

