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

	def __init__(self, configFile): ## Later on, input a named instrument.
		## Test to see if configFile exists.
		try:
			f=open(configFile, mode='r')
			f.close()
		except FileNotFoundError as err:
			print("File not found error: {0}".format(err))
			raise
		
		self.avgDark=None
		self.avgFlat=None
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
		self.noise=None
		self.noiseFile=""
		self.obsDate=""
		self.obsTime=""
		self.preSpeckleBase=""

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
		
		self.burstNumber=np.int(config['HARDCAM']['burstNumber'])
		self.darkBase=config['HARDCAM']['darkBase']
		self.dataBase=config['HARDCAM']['dataBase']
		self.filePattern=config['HARDCAM']['filePattern']
		self.flatBase=config['HARDCAM']['flatBase']
		self.noiseFile=config['HARDCAM']['noiseFile']
		self.obsDate=config['HARDCAM']['obsDate']
		self.obsTime=config['HARDCAM']['obsTime']
		self.preSpeckleBase=config['HARDCAM']['preSpeckleBase']

		## Directory preSpeckleBase must exist or be created in order
		## to continue.
		if not os.path.isdir(self.preSpeckleBase):
			print("{0}: os.mkdir: attempting to create directory:"
					"{1}".format(__name__, self.preSpeckleBase)
					)
			try:
				os.mkdir(self.preSpeckleBase)
			except FileExistsError as err:
				print("File exists error: {0}".format(err))
				raise
			except PermissionError as err:
				print("Permission error: {0}".format(err))
				raise

		## Set-up logging.
		logging.config.fileConfig(self.configFile,
				defaults={'logfilename':
					'{0}/hardcam.log'.format(self.preSpeckleBase)
					}
				)
		self.logger=logging.getLogger('HardcamLog')

		## Print an intro message.
		self.logger.info("This is SSOsoft version {0}".format(self.ssosoftConfig.__version__))
		self.logger.info("Contact {0} to report bugs, make suggestions, "
				"or contribute.".format(self.ssosoftConfig.__email__))
		self.logger.info("Now configuring this {0} data calibration run.".format('HARDCAM'))

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
		self.logger.info("Starting standard HARDCAM run.")
		self.rosa_hardcam_configure_run()
		self.rosa_hardcam_get_file_lists()
		self.rosa_hardcam_order_files()
		self.rosa_hardcam_get_data_image_shapes(self.flatList[0])
		self.avgDark=self.rosa_hardcam_average_image_from_list(
				self.darkList
				)
		self.avgFlat=self.rosa_hardcam_average_image_from_list(
				self.flatList
				)
		self.rosa_hardcam_compute_gain()
		self.rosa_hardcam_save_cal_images()
		self.rosa_hardcam_save_bursts()
		self.logger.info("Finished standard HARDCAM run.")

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

		lastBurst=np.int(len(self.dataList)/self.burstNumber)
		burstShape=(self.burstNumber,)+self.imageShape
		self.logger.info("Preparing burst files, saving in directory: "
				"{0}".format(self.preSpeckleBase)
				)
		self.logger.info("Flat-fielding and saving data to burst files "
				"with burst number: {0}: shape: {1}".format(
					self.burstNumber, burstShape
					)
				)
		burstCube=np.zeros(burstShape, dtype=np.float32)
		burst=0
		while burst < lastBurst:
			burstList=self.dataList[burst*self.burstNumber:(burst+1)*self.burstNumber]
			i=0
			for file in burstList:
				burstCube[i, :, :]=rosa_hardcam_flatfield_correction()
				i+=1
			## Construct filename.
			burstThsnds=burst//1000
			burstHndrds=burst%1000
			burstFile=os.path.join(self.preSpeckleBase,
					"%s_%s_kisip.raw.batch.%02d.%03d" %
					(self.obsDate, self.obsTime, burstThsnds, burstHndrds)
					)
			### Save burstCube
			#try:
			#	with open(burstFile, mode='wb') as f:
			#		## Reverse axis order to save fortran-style.
			#		burstCube.swapaxes(0,2).tofile(f)
			#except Exception as err:
			#	self.logger.critical("Could not save burst file: {0}".format(err))
			#	raise
			#else:
			#	rosa_hardcam_print_progress_save_bursts()
			#	burst+=1
			self.rosa_hardcam_save_binary_image_cube(
					burstCube,#.swapaxes(0,2), 
					burstFile
					)
			rosa_hardcam_print_progress_save_bursts()
			burst+=1

		self.logger.info("Burst files complete: {0}".format(self.preSpeckleBase))

	def rosa_hardcam_save_cal_images(self):
		self.logger.info("Saving average dark: "
				"{0}".format(
					os.path.join(self.preSpeckleBase, 'Halpha_dark.fits')
					)
				)
		self.rosa_hardcam_save_fits_image(self.avgDark,
				os.path.join(self.preSpeckleBase, 'Halpha_dark.fits')
				)

		self.logger.info("Saving average flat: "
				"{0}".format(
					os.path.join(self.preSpeckleBase, 'Halpha_flat.fits')
					)
				)
		self.rosa_hardcam_save_fits_image(self.avgFlat,
				os.path.join(self.preSpeckleBase, 'Halpha_flat.fits')
				)

		self.logger.info("Saving gain: "
				"{0}".format(
					os.path.join(self.preSpeckleBase, 'Halpha_gain.fits')
					)
				)
		self.rosa_hardcam_save_fits_image(self.gain,
				os.path.join(self.preSpeckleBase, 'Halpha_gain.fits')
				)

		self.logger.info("Saving noise: "
				"{0}".format(
					os.path.join(self.preSpeckleBase, 'Halpha_noise.fits')
					)
				)
		self.rosa_hardcam_save_fits_image(self.noise,
				os.path.join(self.preSpeckleBase, 'Halpha_noise.fits')
				)


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
			

