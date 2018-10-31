"""
This is a wrapper library for running KISIP with
regular Rosa/Hardcam data reduction runs.

Initialize with an instance of RosaHardCamCal class.
"""
import configparser
import glob
import logging, logging.config
import os
import subprocess

class kisipWrapper:

	from . import ssosoftConfig

	def __init__(self, RHCC):
		self.batchList=RHCC.batchList
		self.burstFileForm=RHCC.burstFileForm
		self.burstNumber=RHCC.burstNumber
		self.configFile=RHCC.configFile
		self.imageShape=RHCC.imageShape
		self.instrument=RHCC.instrument.upper()
		self.kisipPreSpeckleBatch=0
		self.kisipPreSpeckleStartInd=0
		self.kisipPreSpeckleEndInd=0
		self.logFile=RHCC.logFile
		self.noiseFile=RHCC.noiseFile
		self.obsDate=RHCC.obsDate
		self.obsTime=RHCC.obsTime
		self.postSpeckleBase=RHCC.postSpeckleBase
		self.preSpeckleBase=RHCC.preSpeckleBase
		self.speckleBase=RHCC.speckleBase
		self.workBase=RHCC.workBase
	
	def kisip_configure_run(self):
		config=configparser.ConfigParser()
		config.read(self.configFile)
	
		self.kisipArcsecPerPixX=config[self.instrument]['kisipArcsecPerPixX']
		self.kisipArcsecPerPixY=config[self.instrument]['kisipArcsecPerPixY']
		self.speckledFileForm=config[self.instrument]['speckledFileForm']
		self.wavelengthnm=config[self.instrument]['wavelengthnm']
	
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
		self.kisipEnvMpiNproc=config['KISIP_ENV']['kisipEnvMpiNproc']
		self.kisipEnvMpirun=config['KISIP_ENV']['kisipEnvMpirun']
		self.kisipEnvKisipExe=config['KISIP_ENV']['kisipEnvKisipExe']
	
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
	
		## Directories preSpeckleBase and speckleBase must exist or be
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

	def kisip_despeckle_all_batches(self):
		self.logger.info("Preparing to run KISIP on batches: "
				"{0}".format(self.batchList)
				)
		for batch in self.batchList:
			self.kisip_set_batch_start_end_inds(batch)
			self.kisip_set_environment()
			self.kisip_write_init_files()
			self.kisip_spawn_kisip()
	
	def kisip_set_environment(self):
		self.logger.info("Setting system environment variables for KISIP.")
		self.logger.info("Pre-appending to PATH: {0}".format(self.kisipEnvBin))
		os.environ['PATH']="{0}{1}{2}".format(
				self.kisipEnvBin, os.pathsep, os.environ['PATH'])
		self.logger.info("Pre-appending to LD_LIBRARY_PATH: {0}".format(self.kisipEnvLib))
		os.environ['LD_LIBRARY_PATH']="{0}{1}{2}".format(
				self.kisipEnvLib, os.pathsep, os.environ['PATH']
				)

	def kisip_set_batch_start_end_inds(self, batch):
		self.logger.info("Setting batch number: {0}".format(batch))
		self.kisipPreSpeckleBatch=batch
		self.logger.info("Searching for files: {0}".format(
			os.path.join(
				self.preSpeckleBase,
				self.burstFileForm.format(
					self.obsDate,
					self.obsTime,
					self.kisipPreSpeckleBatch,
					0
					)[:-3]+'*' ## Remove last '000' add '*'.
				)
			)
			)
			
		fList=glob.glob(
				os.path.join(
					self.preSpeckleBase,
					self.burstFileForm.format(
						self.obsDate,
						self.obsTime,
						self.kisipPreSpeckleBatch,
						0
						)[:-3]+'*' ## Remove last '000' add '*'.
					)
				)
		nFile=len(fList)
		if nFile==0:
			self.logger.error("ERROR: batch {0} has {1} files.".format(
				self.kisipPreSpeckleBatch,
				nFile
				)
				)
			self.logger.warning("WARNING: KISIP might run, "
					"but will ultimately do nothing."
					)
		else:
			self.logger.info("Batch {0} has {1} files.".format(
				self.kisipPreSpeckleBatch,
				nFile
				)
				)
		## RISKY ASSUMPTION! kisipPreSpeckleStartInd=0.
		self.logger.info("Batch {0}: setting "
				"start index: {1} "
				"end index: {2}.".format(
					self.kisipPreSpeckleBatch,
					0,
					nFile-1
					)
				)
		self.kisipPreSpeckleStartInd=0	##Assumption.
		self.kisipPreSpeckleEndInd=nFile-1	## Following assumption.

	def kisip_spawn_kisip(self):
		kisipCommand="{0} {1} {2} {3}".format(
				os.path.join(self.kisipEnvBin, self.kisipEnvMpirun),
				'-np',
				self.kisipEnvMpiNproc,
				os.path.join(self.kisipEnvBin, self.kisipEnvKisipExe)
				)
		self.logger.info("KISIP command: {0}".format(kisipCommand))
		self.logger.info("KISIP log will be in directory: {0}".format(self.speckleBase))
		self.logger.info("Now running KISIP for batch: {0} on: "
				"{1} threads.".format(
					self.kisipPreSpeckleBatch,
					self.kisipEnvMpiNproc
					)
				)
		try:
			os.chdir(self.workBase)
			returnCode=subprocess.call([
						os.path.join(self.kisipEnvBin, self.kisipEnvMpirun),
						'-np',
						self.kisipEnvMpiNproc,
						os.path.join(self.kisipEnvBin, self.kisipEnvKisipExe)
						]
						)
		except Exception as err:
			self.logger.critical("CRITICAL: KISIP run failed: {0}".format(err))
			raise
		if returnCode:
			self.logger.error("Something went wrong with KISIP run. "
					"Check logfile. Code: {0}".format(returnCode)
					)
		else:
			self.logger.info(
					"KISIP batch: {0} exited with code: "
					"{1}".format(
						self.kisipPreSpeckleBatch,
						returnCode
						)
					)

	def kisip_write_init_files(self):
		self.logger.info("Preparing to write KISIP init files.")
		self.logger.info("Writing KISIP config file: "
				"{0}".format(os.path.join(self.workBase, 'init_file.dat'))
				)
		try:
			with open(os.path.join(self.workBase, 'init_file.dat'), mode='wt') as f:
				f.write("{0}{1}".format(
					os.path.join(
						self.preSpeckleBase,
						self.burstFileForm.format(
							self.obsDate,
							self.obsTime,
							self.kisipPreSpeckleBatch,
							0
							)[:-4] ## Remove from end '.000'
						),
					os.linesep
					)
					)
				f.write("{0}{1}".format(
					"{:03d}".format(self.kisipPreSpeckleStartInd),
					os.linesep
					)
					)
				f.write("{0}{1}".format(
					"{:03d}".format(self.kisipPreSpeckleEndInd),
					os.linesep
					)
					)
				f.write("{0}{1}".format(
					os.path.join(
						self.speckleBase,
						self.speckledFileForm.format(
							self.obsDate,
							self.obsTime,
							self.kisipPreSpeckleBatch,
							0
							)[:-4] ## Remove from end '.000'
						),
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

