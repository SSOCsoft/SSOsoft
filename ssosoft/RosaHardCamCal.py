import astropy.io.fits as fits
import configparser
import glob
import logging, logging.config
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import sys

class RosaHardCamCal:
	
	from . import ssosoftConfig
	
	def __init__(self, instrument, configFile):
		## Test to see if configFile exists.
		try:
			f=open(configFile, mode='r')
			f.close()
		except Exception as err:
			print("Exception: {0}".format(err))
			raise
		
		self.avgDark=None
		self.avgFlat=None
		self.batchList=[]
		self.burstNumber=0
		self.configFile=configFile
		self.darkBase=""
		self.darkList=[""]
		self.dataBase=""
		self.dataList=[""]
		self.dataShape=None
		self.filePattern=""
		self.flatBase=""
		self.flatList=[""]
		self.gain=None
		self.imageShape=None
		self.instrument=instrument.upper()
		self.logFile=""
		self.noise=None
		self.noiseFile=""
		self.obsDate=""
		self.obsTime=""
		self.preSpeckleBase=""
		self.workBase=""

	def rosa_hardcam_average_image_from_list(self, fileList):
		def rosa_hardcam_print_average_image_progress():
			if not fNum % 100:
				self.logger.info("Progress: "
						"{:0.1%}.".format(fNum/len(fileList))
						)

		self.logger.info("Computing average image from {0} files "
				"in directory: {1}".format(
					len(fileList), os.path.dirname(fileList[0])
					)
				)
		avgIm=np.zeros(self.imageShape, dtype=np.float32)
		fNum=0
		for file in fileList:
			## Read-in image files, add them up.
			avgIm=avgIm+self.rosa_hardcam_read_binary_image(file)
			rosa_hardcam_print_average_image_progress()
			fNum+=1
		avgIm=avgIm/len(fileList)

		self.logger.info("Average complete, directory: "
				"{0}".format(os.path.dirname(fileList[0]))
				)
		return avgIm

	def rosa_hardcam_check_dark_data_flat_shapes(self):
		pass

	def rosa_hardcam_compute_gain(self):
		self.logger.info("Computing gain table.")
		darkSubtracted=self.avgFlat-self.avgDark
		try:
			self.gain=np.median(darkSubtracted)/(darkSubtracted)
		except Exception as err:
			self.error("Error computing gain table: "
					"{0}".format(err)
					)

		self.logger.info("Gain table computed.")

	def rosa_hardcam_compute_noise_file(self):
		self.logger.info("Computing noise cube: shape: "
				"{0}".format(self.imageShape+(self.burstNumber,))
				)
		self.noise=np.zeros((self.burstNumber,)+self.imageShape,
				dtype=np.float32
				)
		i=0
		for flat in self.flatList[0:self.burstNumber]:
			self.noise[i, :, :]=(self.rosa_hardcam_read_binary_image(flat)
					-self.avgDark)
			i+=1
		self.noise=np.multiply(self.noise, self.gain)
		
		self.logger.info("Noise cube complete. "
				"Saving to noise file: "
				"{0}".format(
					os.path.join(self.preSpeckleBase, self.noiseFile)
					)
				)
		self.rosa_hardcam_save_binary_image_cube(
				self.noise,
				os.path.join(self.preSpeckleBase, self.noiseFile)
				)
		self.logger.info("Saved noise file: "
				"{0}".format(os.path.join(self.preSpeckleBase, self.noiseFile))
				)
	
	def rosa_hardcam_configure_run(self):
		def rosa_hardcam_assert_base_dirs(baseDir):
			assert(os.path.isdir(baseDir)), (
					"Directory does not exist: {0}".format(baseDir)
					)

		config=configparser.ConfigParser()
		config.read(self.configFile)
		
		self.burstNumber=np.int(config[self.instrument]['burstNumber'])
		self.burstFileForm=config[self.instrument]['burstFileForm']
		self.darkBase=config[self.instrument]['darkBase']
		self.dataBase=config[self.instrument]['dataBase']
		self.filePattern=config[self.instrument]['filePattern']
		self.flatBase=config[self.instrument]['flatBase']
		self.noiseFile=config[self.instrument]['noiseFile']
		self.obsDate=config[self.instrument]['obsDate']
		self.obsTime=config[self.instrument]['obsTime']
		self.speckledFileForm=config[self.instrument]['speckledFileForm']
		self.workBase=config[self.instrument]['workBase']

		self.preSpeckleBase=os.path.join(self.workBase, 'preSpeckle')
		self.speckleBase=os.path.join(self.workBase, 'speckle')
		self.postSpeckleBase=os.path.join(self.workBase, 'postSpeckle')
		self.darkFile=os.path.join(self.workBase, '{0}_dark.fits'.format(self.instrument))
		self.flatFile=os.path.join(self.workBase, '{0}_flat.fits'.format(self.instrument))
		self.gainFile=os.path.join(self.workBase, '{0}_gain.fits'.format(self.instrument))
		self.noiseFileFits=os.path.join(self.workBase, '{0}_noise'.format(self.instrument))

		## Directories preSpeckleBase, speckleBase, and postSpeckle
		## must exist or be created in order to continue.
		for dirBase in [self.preSpeckleBase, self.speckleBase, self.postSpeckleBase]:
			if not os.path.isdir(dirBase):
				print("{0}: os.mkdir: attempting to create directory:"
						"{1}".format(__name__, dirBase)
						)
				try:
					os.mkdir(dirBase)
				except Exception as err:
					print("An exception was raised: {0}".format(err))
					raise

		## Set-up logging.
		self.logFile='{0}{1}'.format(
				os.path.join(self.workBase,
					'{0}_{1}'.format(self.obsTime, self.instrument.lower())
					),
				'.log'
				)
		logging.config.fileConfig(self.configFile,
				defaults={'logfilename': self.logFile}
				)
		self.logger=logging.getLogger('{0}Log'.format(self.instrument.lower()))

		## Print an intro message.
		self.logger.info("This is SSOsoft version {0}".format(self.ssosoftConfig.__version__))
		self.logger.info("Contact {0} to report bugs, make suggestions, "
				"or contribute.".format(self.ssosoftConfig.__email__))
		self.logger.info("Now configuring this {0} data calibration run.".format(self.instrument))

		## darkBase, dataBase, and flatBase directories must exist.
		try:
			rosa_hardcam_assert_base_dirs(self.darkBase)
		except AssertionError as err:
			self.logger.critical("Fatal: {0}".format(err))
			raise
		else:
			self.logger.info("Using dark directory: {0}".format(self.darkBase))
		
		try:
			rosa_hardcam_assert_base_dirs(self.dataBase)
		except AssertionError as err:
			self.logger.critical("Fatal: {0}".format(err))
			raise
		else:
			self.logger.info("Using data directory: {0}".format(self.dataBase))
		
		try:
			rosa_hardcam_assert_base_dirs(self.flatBase)
		except AssertionError as err:
			self.logger.critical("Fatal: {0}".format(err))
			raise
		else:
			self.logger.info("Using flat directory: {0}".format(self.flatBase))

	def rosa_hardcam_detect_dims(self, imageData):
		## Detects data then usable image dimensions.
		## Assumes all overscan regions within the raw
		## image has a zero value. If no overscan,
		## this function will fail. Will also fail if
		## dead pixels are present.
		def rosa_hardcam_detect_overscan():
			## Borrowed from a Stackoverflow article
			## titled "Finding the consecutive zeros
			## in a numpy array." Author unknown.
			## Create an array that is 1 where imageData is 0,
			## and pad each end with an extra 0.
			## LATER: check for identical indices which would
			##	be indicative of non-contiguous dead
			##	pixels.
			iszero = np.concatenate(([0],
				np.equal(imageData, 0).view(np.int8),
				[0]))
			absdiff = np.abs(np.diff(iszero))
			# Runs start and end where absdiff is 1.
			ranges = np.where(absdiff == 1)[0]
			self.logger.info("Zeros boundary detected at: "
					"{0}".format(np.unique(ranges))
					)
			return ranges
	
		## Detects image columns by looking for the
		## last overscan pixel in the first row. Rows
		## are computed from the quotient of the
		## number of pixels and the number of columns.
		## Get data dimensions.
		## LATER: add function to check for overscan.
		##	Could be done by examining the number
		## 	of results in ovrScn.
		self.logger.info("Attempting to detect overscan and data shape.")
		ovrScn=rosa_hardcam_detect_overscan()
		datDim=(np.uint16(imageData.size/ovrScn[1]),
				ovrScn[1])
		## Detects usable image columns by using
		## the first overscan index. Finds first
		## overscan row by looking for the first
		## occurance where overscan indices are
		## not separated by dx1 or dx2.
		## Get image dimensions.
		self.logger.info("Attempting to detect image shape.")
		dx1=ovrScn[1]-ovrScn[0]
		dx2=ovrScn[2]-ovrScn[1]
		DeltaX=np.abs(np.diff(ovrScn))
		endRow=(np.where(
				np.logical_and(
					DeltaX != dx1,
					DeltaX != dx2
					)
				)[0])[0]/2+1
		imgDim=(np.uint16(endRow), ovrScn[0])
		self.dataShape=datDim
		self.imageShape=imgDim
		self.logger.info("Auto-detected data dimensions "
				"(rows, cols): {0}".format(self.dataShape))
		self.logger.info("Auto-detected image dimensions "
				"(rows, cols): {0}".format(self.imageShape))

	def rosa_hardcam_display_image(self,im):
		plt.imshow(im, origin='upper',
				interpolation='none',
				cmap='hot'
				)
		plt.show()

	def rosa_hardcam_get_cal_images(self):
		if os.path.exists(self.darkFile):
			self.logger.info("Average dark file found: {0}".format(self.darkFile))
			self.logger.info("Reading average dark.")
			with fits.open(self.darkFile) as hdu:
				self.avgDark=hdu[0].data
		else:
			self.avgDark=self.rosa_hardcam_average_image_from_list(
					self.darkList
					)
		if os.path.exists(self.flatFile):
			self.logger.info("Average flat file found: {0}".format(self.flatFile))
			self.logger.info("Reading average flat.")
			with fits.open(self.flatFile) as hdu:
				self.avgFlat=hdu[0].data
		else:
			self.avgFlat=self.rosa_hardcam_average_image_from_list(
					self.flatList
					)
		if os.path.exists(self.gainFile):
			self.logger.info("Gain file found: {0}".format(self.gainFile))
			self.logger.info("Reading gain file.")
			with fits.open(self.gainFile) as hdu:
				self.gain=hdu[0].data
		else:
			self.rosa_hardcam_compute_gain()


	def rosa_hardcam_get_file_lists(self):
		def rosa_hardcam_assert_file_list(fList):
			assert(len(fList)!=0), "List contains no matches."

		self.logger.info("Searching for darks, flats, and data "
				"ending with {0}".format(self.filePattern))
		self.logger.info("Searching for dark image files: {0}".format(self.darkBase))
		self.darkList=glob.glob(
			os.path.join(self.darkBase, self.filePattern)
			)
		try:
			rosa_hardcam_assert_file_list(self.darkList)
		except AssertionError as err:
			self.logger.critical("Error: darkList: {0}".format(err))
			raise
		else:
			self.logger.info("Files in darkList: {0}".format(len(self.darkList)))

		self.logger.info("Searching for data image files: {0}".format(self.dataBase))
		self.dataList=glob.glob(
				os.path.join(self.dataBase, self.filePattern)
				)
		try:
			rosa_hardcam_assert_file_list(self.dataList)
		except AssertionError as err:
			self.logger.critical("Error: dataList: {0}".format(err))
			raise
		else:
			self.logger.info("Files in dataList: {0}".format(len(self.dataList)))

		self.logger.info("Searching for flat image files: {0}".format(self.flatBase))
		self.flatList=glob.glob(
				os.path.join(self.flatBase, self.filePattern)
				)
		try:
			rosa_hardcam_assert_file_list(self.flatList)
		except AssertionError as err:
			self.logger.critical("Error: flatList: {0}".format(err))
			raise
		else:
			self.logger.info("Files in flatList: {0}".format(len(self.flatList)))

	def rosa_hardcam_get_data_image_shapes(self, file):
		self.logger.info("Detecting image and data dimensions in "
				"binary file: {0}".format(file)
				)
		try:
			with open(file, mode='rb') as imageFile:
				imageData=np.fromfile(imageFile,
						dtype=np.uint16
						)
		except Exception as err:
			self.logger.critical("Could not get image or data "
					"shapes: {0}".format(err)
					)
			raise
		self.rosa_hardcam_detect_dims(imageData)

	def rosa_hardcam_order_files(self):
		def rosa_hardcam_order_file_list(fList):
			orderList=['']*len(fList) ## List length same as fList
			ptrn='[0-9]+'	#Match any digit one or more times.
			p=re.compile(ptrn)
			for f in fList:
				head, tail=os.path.split(f)
				match=p.match(tail)
				if match:
					digitNew=match.group()[::-1]
					orderList[int(digitNew)]=f
				else:
					self.logger.error(
							"Unexpected filename format: "
							"{0}".format(f)
							)
			try:
				assert(all(orderList)), "List could not be ordered."
			except AssertionError as err:
				self.logger.critical("Error: {0}".format(err))
				raise
			else:
				return orderList

		self.logger.info("Sorting darkList.")
		self.darkList=rosa_hardcam_order_file_list(self.darkList)
		self.logger.info("Sorting flatList.")
		self.flatList=rosa_hardcam_order_file_list(self.flatList)
		self.logger.info("Sorting dataList.")
		self.dataList=rosa_hardcam_order_file_list(self.dataList)

	def rosa_hardcam_read_binary_image(self, file, dataShape=None, imageShape=None):
		if dataShape is None:
			dataShape=self.dataShape
		if imageShape is None:
			imageShape=self.imageShape

		try:
			with open(file, mode='rb') as imageFile:
				imageData=np.fromfile(imageFile,
						dtype=np.uint16
						)
		except Exception as err:
			self.logger.critical("Could not open/read binary image file: "
					"{0}".format(err)
					)
			raise

		im=imageData.reshape((dataShape))
		## Generate a tuple of slice objects, s[i]~ 0:imageShape[i]
		s=tuple()
		for t in imageShape:
			s=s+np.index_exp[0:t]
		im=im[s]
		return np.float32(im)

	def rosa_hardcam_run_calibration(self):
		self.rosa_hardcam_configure_run()
		self.logger.info("Starting standard HARDCAM calibration.")
		self.rosa_hardcam_get_file_lists()
		self.rosa_hardcam_order_files()
		self.rosa_hardcam_get_data_image_shapes(self.flatList[0])
		self.rosa_hardcam_get_cal_images()
		self.rosa_hardcam_save_cal_images()
		self.rosa_hardcam_save_bursts()
		self.logger.info("Finished standard HARDCAM calibration.")

	def rosa_hardcam_save_binary_image_cube(self, data, file):
		try:
			with open(file, mode='wb') as f:
				data.tofile(f)
		except Exception as err:
			self.logger.critical("Could not save binary file: {0}".format(err))
			raise

	def rosa_hardcam_save_bursts(self):
		def rosa_hardcam_flatfield_correction():
			data=self.rosa_hardcam_read_binary_image(file)
			corrected=self.gain*(data-self.avgDark)
			return corrected

		def rosa_hardcam_print_progress_save_bursts():
			self.logger.info("Progress: {:0.2%} "
					"with file: {:s}".format(
						burst/lastBurst,
						burstFile
						)
					)

		lastBurst=len(self.dataList)//self.burstNumber
		burstShape=(self.burstNumber,)+self.imageShape
		self.logger.info("Preparing burst files, saving in directory: "
				"{0}".format(self.preSpeckleBase)
				)
		self.logger.info("Number of files to be saved: "
				"{0}".format(len(self.dataList))
				)
		self.logger.info("Flat-fielding and saving data to burst files "
				"with burst number: {0}: shape: {1}".format(
					self.burstNumber, burstShape
					)
				)
		burstCube=np.zeros(burstShape, dtype=np.float32)
		burst=0
		batch=-1
		while burst < lastBurst:
			burstList=self.dataList[burst*self.burstNumber:(burst+1)*self.burstNumber]
			i=0
			for file in burstList:
				burstCube[i, :, :]=rosa_hardcam_flatfield_correction()
				i+=1
			## Construct filename, KISIP takes batches of 1000 or less.
			burstThsnds=burst//1000
			burstHndrds=burst%1000
			burstFile=os.path.join(
					self.preSpeckleBase,
					(self.burstFileForm).format(
						self.obsDate, self.obsTime, burstThsnds, burstHndrds
						)
					)
			self.rosa_hardcam_save_binary_image_cube(
					burstCube, 
					burstFile
					)
			rosa_hardcam_print_progress_save_bursts()
			## Add burstThsnds to batchList if needed.
			if burstThsnds != batch:
				batch=burstThsnds
				(self.batchList).append(batch)
			burst+=1

		self.logger.info("Burst files complete: {0}".format(self.preSpeckleBase))

	def rosa_hardcam_save_cal_images(self):
		if os.path.exists(self.darkFile):
			self.logger.info("Dark file already exists: {}".format(self.darkFile))
		else:
			self.logger.info("Saving average dark: "
					"{0}".format(
						os.path.join(
							self.workBase,
							self.darkFile
							)
						)
					)
			self.rosa_hardcam_save_fits_image(self.avgDark,
					os.path.join(
						self.workBase,
						self.darkFile
						)
					)
		if os.path.exists(self.flatFile):
			self.logger.info("Flat file already exists: {0}".format(self.flatFile))
		else:
			self.logger.info("Saving average flat: "
					"{0}".format(
						os.path.join(
							self.workBase,
							self.flatFile
							)
						)
					)
			self.rosa_hardcam_save_fits_image(self.avgFlat,
					os.path.join(
						self.workBase,
						self.flatFile
						)
					)
		if os.path.exists(self.gainFile):
			self.logger.info("Gain file already exists: {0}".format(self.gainFile))
		else:
			self.logger.info("Saving gain: "
					"{0}".format(
						os.path.join(
							self.workBase,
							self.gainFile
							)
						)
					)
			self.rosa_hardcam_save_fits_image(self.gain,
					os.path.join(
						self.workBase,
						self.gainFile
						)
					)
		if os.path.exists(self.noiseFileFits):
			self.logger.info("Noise FITS file already exists: {0}".format(self.noiseFileFits))
		else:
			self.logger.info("Saving noise: "
					"{0}".format(
						os.path.join(
							self.workBase,
							self.noiseFileFits
							)
						)
					)
			self.rosa_hardcam_save_fits_image(self.noise,
					os.path.join(
						self.workBase,
						self.noiseFileFits
						)
					)

	def rosa_hardcam_save_despeckled_as_fits(self):
		self.logger.info("Saving despeckled binary image files to FITS.")
		self.logger.info("Searching for files: "
				"{0}".format(
					os.path.join(self.speckleBase,
						self.speckledFileForm.format(
							self.obsDate, self.obsTime, 0, 0
							)[:-7]+'*.final'
						)
					)
				)
		fList=glob.glob(os.path.join(self.speckleBase, 
			self.speckledFileForm.format(
				self.obsDate, self.obsTime, 0, 0
				)[:-7]+'*.final'
			)
			)
		try:
			assert(len(fList) != 0)
		except Exception as err:
			self.logger.critical("CRITICAL: no files found: {0}".format(err))
			raise
		else:
			self.logger.info("Found {0} files.".format(len(fList)))
		for file in fList:
			im=self.rosa_hardcam_read_binary_image(file,
					imageShape=self.imageShape[::-1]+(2,),
					dataShape=self.imageShape[::-1]+(2,)
					) ## KISIP saves with Fortran ordering.
			fName=os.path.basename(file)
			self.rosa_hardcam_save_fits_image(im[:, :, 1],
					os.path.join(
						self.postSpeckleBase,
						fName+'.fits'
						)
					)
		self.logger.info("Finished saving despeckled images as FITS"
				"in directory: {0}".format(self.postSpeckleBase))

	def rosa_hardcam_save_fits_image(self, image, file, clobber=True):
		hdu=fits.PrimaryHDU(image)
		hdul=fits.HDUList([hdu])
		try:
			hdul.writeto(file, clobber=clobber)
		except Exception as err:
			self.logger.warning("Could not write FITS file: "
					"{0}".format(file)
					)
			self.logger.warning("FITS write warning: continuing, but "
					"this could cause problems later."
					)
			

