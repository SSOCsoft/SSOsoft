import astropy.io.fits as fits
import configparser
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import re

class RosaHardCamCal:

	def __init__(self, configFile): ## Later on, input a named instrument.
		## Test to see if configFile exists.
		try:
			f=open(configFile, mode='r')
			f.close()
		except FileNotFoundError as err:
			print("File not found error: {0}".format(err))
		
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
		self.gain=[]
		self.imageShape=None
		self.obsDate=""
		self.obsTime=""
		self.preSpeckleBase=""

		
	def rosa_hardcam_average_image_from_list(self, fileList):
		avgIm=np.zeros(self.imageShape, dtype=np.float32)
		for file in fileList:
			## Read-in image files, add them up.
			avgIm=avgIm+self.rosa_hardcam_read_binary_image(file)
		avgIm=avgIm/len(fileList)
		return avgIm

	def rosa_hardcam_check_dark_data_flat_shapes(self):
		pass

	def rosa_hardcam_configure_run(self):
		config=configparser.ConfigParser()
		config.read(self.configFile)
		self.burstNumber=np.int(config['HARDCAM']['burstNumber'])
		self.darkBase=config['HARDCAM']['darkBase']
		self.dataBase=config['HARDCAM']['dataBase']
		self.filePattern=config['HARDCAM']['filePattern']
		self.flatBase=config['HARDCAM']['flatBase']
		self.obsDate=config['HARDCAM']['obsDate']
		self.obsTime=config['HARDCAM']['obsTime']
		self.preSpeckleBase=config['HARDCAM']['preSpeckleBase']

		## LATER: check to be sure all directories and files exist,
		## if not, create them.

	def rosa_hardcam_compute_gain(self):
		darkSubtracted=self.avgFlat-self.avgDark
		self.gain=np.median(darkSubtracted)/(darkSubtracted)
	
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
			return ranges
	
		## Detects image columns by looking for the
		## last overscan pixel in the first row. Rows
		## are computed from the quotient of the
		## number of pixels and the number of columns.
		## Get data dimensions.
		## LATER: add function to check for overscan.
		##	Could be done by examining the number
		## 	of results in ovrScn.
		ovrScn=rosa_hardcam_detect_overscan()
		datDim=(np.uint16(imageData.size/ovrScn[1]+1),
				ovrScn[1])
		## Detects usable image columns by using
		## the first overscan index. Finds first
		## overscan row by looking for the first
		## occurance where overscan indices are
		## not separated by dx1 or dx2.
		## Get image dimensions.
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
		print('Auto-detected data dimensions: ', self.dataShape)
		print('Auto-detected image dimensions: ', self.imageShape)

	
	def rosa_hardcam_display_image(self,im):
		plt.imshow(im, origin='upper',
				interpolation='none',
				cmap='hot'
				)
		plt.show()

	def rosa_hardcam_get_file_lists(stlf):
		self.darkList=sorted(glob.glob(self.darkBase+self.filePattern))
		self.dataList=sorted(glob.glob(self.dataBase+self.filePattern))
		self.flatList=sorted(glob.glob(self.flatBase+self.filePattern))

	def rosa_hardcam_get_data_image_shapes(self, file):
		imageFile=open(file, mode='rb')
		imageData=np.fromfile(imageFile, dtype=np.uint16)
		self.rosa_hardcam_detect_dims(imageData)
		imageFile.close()

	def rosa_hardcam_read_binary_image(self, file):
		imageFile=open(file, mode='rb')
		imageData=np.fromfile(imageFile, dtype=np.uint16)
		im=np.empty(self.dataShape, dtype=np.uint16)
		for i in np.arange(self.dataShape[0]-1):
			im[i, :]=imageData[i*self.dataShape[1]:(i+1)*self.dataShape[1]]
		im=im[0:self.imageShape[0], 0:self.imageShape[1]]
		imageFile.close()
		return np.float32(im)

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
					print("Invalid filename format for dark.")
			return orderList

		self.darkList=rosa_hardcam_order_file_list(self.darkList)
		self.flatList=rosa_hardcam_order_file_list(self.flatList)
		self.dataList=rosa_hardcam_order_file_list(self.dataList)

	def rosa_hardcam_run_calibration(self):
		self.rosa_hardcam_configure_run()
		self.rosa_hardcam_get_file_lists()
		self.rosa_hardcam_order_files()
		#self.rosa_hardcam_get_file_lists()
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

	def rosa_hardcam_save_cal_images(self):
		self.rosa_hardcam_save_fits_image(self.avgDark,
				os.path.join(self.preSpeckleBase, 'Halpha_dark.fits')
				)
		self.rosa_hardcam_save_fits_image(self.avgFlat,
				os.path.join(self.preSpeckleBase, 'Halpha_flat.fits')
				)
		self.rosa_hardcam_save_fits_image(self.gain,
				os.path.join(self.preSpeckleBase, 'Halpha_gain.fits')
				)

	def rosa_hardcam_save_bursts(self):
		def rosa_hardcam_flatfield_correction():
			data=self.rosa_hardcam_read_binary_image(file)
			corrected=self.gain*(data-self.avgDark)
			return corrected

		lastBurst=np.int(len(self.dataList)/self.burstNumber)
		burstShape=self.imageShape+(self.burstNumber,)
		burstCube=np.zeros(burstShape, dtype=np.float32)
		burst=0
		while burst < lastBurst:
			burstList=self.dataList[burst*self.burstNumber:(burst+1)*self.burstNumber]
			i=0
			for file in burstList:
				burstCube[:, :, i]=rosa_hardcam_flatfield_correction()
				i+=1
			## Construct filename.
			burstFile=os.path.join(self.preSpeckleBase,
					"%s_%s_kisip.raw.batch.%04d" %
					(self.obsDate, self.obsTime, burst)
					)
			## Save burstCube
			f=open(burstFile, mode='wb')
			burstCube.tofile(f)	## Saving the file using Numpy.
			f.close()
			burst+=1

	def rosa_hardcam_save_fits_image(self, image, file, clobber=True):
		hdu=fits.PrimaryHDU(image)
		hdul=fits.HDUList([hdu])
		hdul.writeto(file, clobber=clobber)

