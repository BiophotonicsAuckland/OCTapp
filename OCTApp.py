from PyQt5 import QtCore, QtGui, QtWidgets
from UI import Ui_MainWindow
import sys

from multiprocessingFunc import threeD_parallelCalculation
from multiprocessing import freeze_support
freeze_support()

import numpy as np
from matplotlib import pyplot as plt

import multiprocessing
from functools import partial

from scipy.interpolate import interp1d
import os

import locale
locale.setlocale(locale.LC_NUMERIC, '') #Set up decimal mark conversion in case of european ',' decimal mark, etc.

## Multiprocessing can cause some problems, especially with IDEs. Set to 'false' to revert to a backup version.
multithreadingEnabled = True

## For the 3D tab, the calculation can be multithreaded. Set up the number of processes to be created here
processCount = multiprocessing.cpu_count() - 1

## This ought to be a class variable, but refactoring would take time.
imageOCT = []

#Misc Global Variables
transformationHistory = []

# Settings
CURRENT_SETTINGS_VERSION = 1.3 # Last changed 27 June 2018
oldSettingsFile = False
settings = {}

def initialiseSettingsFile():
	global settings
	settings = {'spectra_line':'', 'linearphase_line':'', 'nonlinearphase_line':'', 'referencespectrum_line':'', 'dispersionvectors_line':'', 'linearisation_check':False, 'dispersioncompensation_check':False, 'backgroundsubtraction_check':True, 'referencefromfile_radio':False, 'transformationHistory_transform':[], 'rotateState_state':0, 'upperthreshold_dspin':1, 'lowerthreshold_dspin':0, 'version':CURRENT_SETTINGS_VERSION}

#If settings file doesn't exist, create it with default values.
if (not os.path.isfile( os.getcwd() + '\\'+'settings.npy')): 
	initialiseSettingsFile()
	np.save('settings', settings)
else: # If the settings file does exist, load it
	settings = np.load('settings.npy').item()

#Check to see if settings file is outdated
try: 
	if (settings['version'] != CURRENT_SETTINGS_VERSION):
		oldSettingsFile = True
except KeyError:
	oldSettingsFile = True

# If settings file is old, reset to defaults
if oldSettingsFile: 
	initialiseSettingsFile()
	np.save('settings', settings)

class Prog(QtWidgets.QMainWindow):
	def __init__(self):
		global settings
		global transformationHistory
		global oldSettingsFile
		super().__init__();
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		
		self.onClickPos = [] # Stores mouse coordinates when click-drag zooming
		self.onReleasePos = []
		self.lastFFTValue = None
		self.firstCalculation = True
		self.extent = None #Image dimensions
		self.FFTResChanged = False # If the FFT resolution changed since the last calculation
		self.FFTRatio = None #Ratio between the previously used FFT resolution and the current resolution
		self.lim = None # Currently viewed area
		self.colourmap = 'Greys_r' # Currently selected colourmap
		
		self.currentProcessingImage = [] # The image currently displayed in the image processing tab
		self.preFFTImage = [] # Image prior to FFT, used in subimages tab
		self.subimage1 = [] # Bottom left image on the subimages tab
		self.subimage2 = [] # Bottom right image on subimages tab
		self.rawSpectra = []			# Black on subimages grah
		self.spectrum_sub1_disp = []	# Blue
		self.spectrum_sub2_disp = []	# Red
		self.ui.scrollArea_2.setObjectName("ScrollArea_2") # For some reason, the scroll background becomes grey without these lines
		self.ui.scrollArea_2.setStyleSheet("QScrollArea#ScrollArea_2 { background-color: rgb(255, 255,255)}")
		self.ui.scrollArea_2.show()
		self.aspect = "auto" #aspect ratio used for image
		self.processingFlag = False # True if the image displayed in image processing tab has been enhanced, if false then it's only been copied from standard tab
		
		self.processingExtent = None # dimensions, used only for axis recalculation in processing tab
		self.processingAspect = None # used to determine image ratio for axis recalculation
		
		self.threeD_currentImage = []
		self.lastCalculationWasOnlyOneFace = True # Used to determine if the last calculation done on the 3D tab was one face or all cross-sections
		self.threeD_extent = None # Image dimensions of 3D image
		self.threeD_lim = None # Currently viewed area on 3D tab
		self.threeD_aspect = "auto"
		
		## Used for referencing tab indices
		self.standardImageTab = 0
		self.imageProcessingTab = 1
		self.subImagesTab = 2
		self.dispersionTab = 3
		self.segmentationTab = 4
		self.threeDTab = 5
		
		for i in [self.subImagesTab, self.dispersionTab, self.segmentationTab]:
			self.ui.tabwidget.setTabEnabled(i, False)
		
		########## Open File Dialogs ##########
		self.ui.spectra_open.clicked.connect(self.getFile) # on button click, call getFile function
		self.ui.referencespectrum_open.clicked.connect(self.getFile)
		self.ui.linearphase_open.clicked.connect(self.getFile)
		self.ui.nonlinearphase_open.clicked.connect(self.getFile)
		self.ui.dispersioncompensation_open.clicked.connect(self.getFile)
		self.ui.threeD_datafolder_open.clicked.connect(self.getFolder)
		
		########## Open Save Dialogs ##########
		self.ui.standardimage_savePNG_button.clicked.connect(self.saveAsPNG) # on button click call saveAs PNG function
		self.ui.threeD_savePNG_button.clicked.connect(self.saveAsPNG)
		
		# self.ui.dispersion_cornea_opticalthickness_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_cornea_geometricalthickness_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_cornea_walkoff_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_cornea_beta2_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_aqueoushumor_opticalthickness_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_aqueoushumor_geometricalthickness_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_aqueoushumor_walkoff_savePNG_button.clicked.connect(self.saveAsPNG)
		# self.ui.dispersion_aqueoushumor_beta2_savePNG_button.clicked.connect(self.saveAsPNG)
		
		self.ui.standardimage_saveTXT_button.clicked.connect(self.saveAsTXT) # etc...
		self.ui.threeD_saveTXT_button.clicked.connect(self.saveAsTXT)
		
		# self.ui.dispersion_cornea_opticalthickness_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_cornea_geometricalthickness_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_cornea_walkoff_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_cornea_beta2_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_aqueoushumor_opticalthickness_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_aqueoushumor_geometricalthickness_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_aqueoushumor_walkoff_saveTXT_button.clicked.connect(self.saveAsTXT)
		# self.ui.dispersion_aqueoushumor_beta2_saveTXT_button.clicked.connect(self.saveAsTXT)
		
		########## Standard Image Setup ##########
		self.ui.calculate_button.clicked.connect(self.calculateImage)
		
		self.ui.hendpix_spin.valueChanged.connect(self.updateStartZoomMax) # Update spinbox maximum when start pos is changed
		self.ui.vendpix_spin.valueChanged.connect(self.updateStartZoomMax)
		
		self.ui.hstartpix_spin.valueChanged.connect(self.updateEndZoomMin) # Update spinbox minimum when end is changed
		self.ui.vstartpix_spin.valueChanged.connect(self.updateEndZoomMin)
		
		self.ui.resetzoom_button.clicked.connect(self.resetZoom)
		self.ui.rotate_button.clicked.connect(self.rotateImage)
		self.ui.mirror_button.clicked.connect(self.mirrorImage)
		self.ui.zoom_button.clicked.connect(self.zoom)
		
		self.ui.upperthreshold_dspin.valueChanged.connect(self.redrawIntensity) # make sure you can't set lower intensity higher than upper intensity
		self.ui.upperthreshold_slider.valueChanged.connect(self.redrawIntensity)
		self.ui.lowerthreshold_dspin.valueChanged.connect(self.redrawIntensity) # make sure you can't set upper intensity lower than lower intensity
		self.ui.lowerthreshold_slider.valueChanged.connect(self.redrawIntensity)
		
		self.ui.vstartpix_spin.valueChanged.connect(self.zoom)
		self.ui.vendpix_spin.valueChanged.connect(self.zoom)
		self.ui.hstartpix_spin.valueChanged.connect(self.zoom)
		self.ui.hendpix_spin.valueChanged.connect(self.zoom)
		
		self.ui.colourmap_combo.currentIndexChanged.connect(self.changeColourmap)
		
		self.ui.imagewidget_mpl.canvas.fig.delaxes(self.ui.imagewidget_mpl.canvas.ax) # reduce widget margins
		self.ui.imagewidget_mpl.canvas.ax = self.ui.imagewidget_mpl.canvas.fig.add_axes([0.1,0.08,0.86,0.9]) # new margins as a fraction of widget size: [left_start, bottom_start, horiz_length, vert_length]
		self.ui.axisrecalculation_check.stateChanged.connect(self.axisRecalculation)
		
		self.ui.imagewidget_mpl.canvas.mpl_connect('draw_event', self.synchroniseZoom)
		
		########## Image Processing Tab Setup ##########
		self.ui.processing_colourmap_combo.currentIndexChanged.connect(self.changeColourmap)
		self.ui.processing_removelastpixels.valueChanged.connect(lambda: self.processing_enhanceDeeperLayers(resetThresholds = False))
		self.ui.processing_removefirstpixels.valueChanged.connect(lambda: self.processing_enhanceDeeperLayers(resetThresholds = False))
		self.ui.processing_processimage_button.clicked.connect(self.processing_enhanceDeeperLayers)
		
		self.ui.processing_upperthreshold_dspin.valueChanged.connect(self.processing_redrawIntensity)
		self.ui.processing_upperthreshold_slider.valueChanged.connect(self.processing_redrawIntensity)
		self.ui.processing_lowerthreshold_dspin.valueChanged.connect(self.processing_redrawIntensity)
		self.ui.processing_lowerthreshold_slider.valueChanged.connect(self.processing_redrawIntensity)
		
		self.ui.processing_imagewidget_mpl.canvas.fig.delaxes(self.ui.processing_imagewidget_mpl.canvas.ax) # reduce widget margins
		self.ui.processing_imagewidget_mpl.canvas.ax = self.ui.processing_imagewidget_mpl.canvas.fig.add_axes([0.1,0.08,0.86,0.9]) # new margins as a fraction of widget size: [left_start, bottom_start, horiz_length, vert_length]
		
		self.ui.processing_saveTXT_button.clicked.connect(self.saveAsTXT)
		self.ui.processing_savePNG_button.clicked.connect(self.saveAsPNG)
		self.ui.processing_copyimagefromstandard_button.clicked.connect(self.processing_copyImage)
		self.processing_updateCid = self.ui.processing_imagewidget_mpl.canvas.mpl_connect('draw_event', self.processing_cropImageFromZoom)
		self.ui.processing_enhancedeeperlayers_check.stateChanged.connect(self.processing_toggleProcessImageButton)
		self.ui.processing_progressbar.hide()
		self.ui.processing_axisrecalculation_check.stateChanged.connect(self.processing_axisRecalculation)
		
		########## Sub-Images Tab Setup ##########
		self.ui.tabwidget.currentChanged.connect(lambda: self.subimages_calculateImages(calculateTopGraphOnly = True))
		self.ui.subimages_calculate_button.clicked.connect(self.subimages_calculateImages)
		self.ui.subimages_subimage1upperthreshold_dspin.valueChanged.connect(self.subimages_updateGraph)
		self.ui.subimages_subimage1lowerthreshold_dspin.valueChanged.connect(self.subimages_updateGraph)
		self.ui.subimages_subimage2upperthreshold_dspin.valueChanged.connect(self.subimages_updateGraph)
		self.ui.subimages_subimage2lowerthreshold_dspin.valueChanged.connect(self.subimages_updateGraph)
		self.ui.subimages_spectrumnum_spin.valueChanged.connect(self.subimages_updateSpectrumNum)
		
		self.ui.subimages_main_mpl.canvas.fig.delaxes(self.ui.subimages_main_mpl.canvas.ax)
		self.ui.subimages_main_mpl.canvas.ax = self.ui.subimages_main_mpl.canvas.fig.add_axes([0.12,0.14,0.84,0.84])
		
		########## 3D Tab Setup ##########
		self.ui.threeD_one_button.clicked.connect(self.threeD_calculateOneImage)
		self.ui.threeD_all_button.clicked.connect(self.threeD_calculateAllImages)
		
		self.ui.threeD_all_radio.toggled.connect(self.threeD_toggleStorage)
		
		self.ui.threeD_all_spin.valueChanged.connect(self.threeD_changeViewedFace)
		self.ui.threeD_hendpix_spin.valueChanged.connect(self.threeD_updateStartZoomMax)
		self.ui.threeD_vendpix_spin.valueChanged.connect(self.threeD_updateStartZoomMax)
		
		self.ui.threeD_hstartpix_spin.valueChanged.connect(self.threeD_updateEndZoomMin)
		self.ui.threeD_vstartpix_spin.valueChanged.connect(self.threeD_updateEndZoomMin)
		
		self.ui.threeD_vstartpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_vendpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_hstartpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_hendpix_spin.valueChanged.connect(self.threeD_zoom)
		
		self.ui.threeD_upperthreshold_dspin.valueChanged.connect(self.threeD_redrawIntensity) # make sure you can't set lower intensity higher than upper intensity
		self.ui.threeD_upperthreshold_slider.valueChanged.connect(self.threeD_redrawIntensity)
		self.ui.threeD_lowerthreshold_dspin.valueChanged.connect(self.threeD_redrawIntensity) # make sure you can't set upper intensity lower than lower intensity
		self.ui.threeD_lowerthreshold_slider.valueChanged.connect(self.threeD_redrawIntensity)
		
		self.ui.threeD_colourmap_combo.currentIndexChanged.connect(self.changeColourmap)
		self.ui.threeD_imagewidget_mpl.canvas.mpl_connect('draw_event', self.threeD_synchroniseZoom)
		self.ui.threeD_rotate_button.clicked.connect(self.threeD_rotateImage)
		self.ui.threeD_mirror_button.clicked.connect(self.threeD_mirrorImage)
		
		self.ui.threeD_progressbar.hide()
		
		########## Miscellaneous Setup ##########
		self.updateLegend()
		self.ui.backgroundsubtraction_check.stateChanged.connect(self.toggleRadio)
		self.ui.dispersioncompensation_check.stateChanged.connect(self.toggleRadio)
		self.ui.std_progressbar.hide()
		
		########## Load previous Settings ##########
		if (oldSettingsFile == True):
			OldSettingsFileBox = QtWidgets.QMessageBox()
			OldSettingsFileBox.setText("Outdated Settings File Detected. New Settings File Created. Settings reverted to Defaults.")
			OldSettingsFileBox.setIcon(QtWidgets.QMessageBox.Warning)
			OldSettingsFileBox.setWindowTitle("Old Settings File Detected")
			OldSettingsFileBox.exec()
		
		for i in settings:
			if (settings[i] != ''):
				if (i[-5:] == "_line"):
					eval('self.ui.' + i + '.setText(settings[i])') # This does open up arbitrary code execution as the settings file is partly human-readable. I doubt it would be a problem though.
				elif ((i[-6:] == "_check") or (i[-6:] == "_radio")) and (eval('self.ui.' + i + '.isChecked()') != settings[i]):
					eval('self.ui.' + i + '.setChecked(settings[i])')
				elif (i[-10:] == '_transform'):
					transformationHistory = settings[i]
				elif (i[-6:] == '_dspin'):
					eval('self.ui.' + i + '.setValue(float(settings[i]))')
		
		########## Link Parameter Updates to Settings File ##########
		self.ui.upperthreshold_dspin.valueChanged.connect(self.updateSettingsFileDspin)
		self.ui.lowerthreshold_dspin.valueChanged.connect(self.updateSettingsFileDspin)
		
		self.ui.upperthreshold_dspin.valueChanged.connect(self.redrawIntensity)
		self.ui.lowerthreshold_dspin.valueChanged.connect(self.redrawIntensity)
		self.ui.upperthreshold_slider.valueChanged.connect(self.redrawIntensity)
		self.ui.lowerthreshold_slider.valueChanged.connect(self.redrawIntensity)
		
		self.ui.spectra_line.textChanged.connect(self.resetTransformationHistory) #Reset transformation history on file change
		self.ui.linearphase_line.textChanged.connect(self.resetTransformationHistory)
		self.ui.nonlinearphase_line.textChanged.connect(self.resetTransformationHistory)
		self.ui.referencespectrum_line.textChanged.connect(self.resetTransformationHistory)
		self.ui.dispersioncompensation_line.textChanged.connect(self.resetTransformationHistory)
		
		self.ui.spectra_line.textChanged.connect(self.updateSettingsFileFilename) #Update the settings file if the filename is manually typed in/cleared
		self.ui.linearphase_line.textChanged.connect(self.updateSettingsFileFilename)
		self.ui.nonlinearphase_line.textChanged.connect(self.updateSettingsFileFilename)
		self.ui.referencespectrum_line.textChanged.connect(self.updateSettingsFileFilename)
		self.ui.dispersioncompensation_line.textChanged.connect(self.updateSettingsFileFilename)
		
		self.ui.linearphase_line.textChanged.connect(self.autoUpdateParameters) #Toggle the relevant parameter if a filepath is entered
		self.ui.nonlinearphase_line.textChanged.connect(self.autoUpdateParameters)
		self.ui.referencespectrum_line.textChanged.connect(self.autoUpdateParameters)
		self.ui.dispersioncompensation_line.textChanged.connect(self.autoUpdateParameters)
		
		self.ui.linearisation_check.toggled.connect(self.updateSettingsFileCheckbox)
		self.ui.dispersioncompensation_check.toggled.connect(self.updateSettingsFileCheckbox)
		self.ui.backgroundsubtraction_check.toggled.connect(self.updateSettingsFileCheckbox)
		self.ui.referencefromfile_radio.toggled.connect(self.updateSettingsFileCheckbox)
	
	#Disconnect the spinboxes from their slots, in case we want to update their values without promping a redraw
	def disconnectZoomUpdates(self): 
		self.ui.vstartpix_spin.disconnect()
		self.ui.vendpix_spin.disconnect()
		self.ui.hstartpix_spin.disconnect()
		self.ui.hendpix_spin.disconnect()
	
	#Reconnect the spinboxes once we've done what we want
	def reconnectZoomUpdates(self): 
		self.ui.vstartpix_spin.valueChanged.connect(self.zoom)
		self.ui.vendpix_spin.valueChanged.connect(self.zoom)
		self.ui.hstartpix_spin.valueChanged.connect(self.zoom)
		self.ui.hendpix_spin.valueChanged.connect(self.zoom)
		
		self.ui.hendpix_spin.valueChanged.connect(self.updateStartZoomMax)
		self.ui.vendpix_spin.valueChanged.connect(self.updateStartZoomMax)
		self.ui.hstartpix_spin.valueChanged.connect(self.updateEndZoomMin)
		self.ui.vstartpix_spin.valueChanged.connect(self.updateEndZoomMin)
	
	#Make sure that if the user pans/zooms the image with the toolbar, these changes are applied to the start/end pixel spinboxes too
	def synchroniseZoom(self, event): 
		if (self.ui.imagewidget_mpl.mpl_toolbar._active != None): #If the current toolbar mode is either PAN or ZOOM, then update the coordinates of the spinboxes with the coordinates of the on-screen image
			vstart = min(self.ui.imagewidget_mpl.canvas.ax.get_ylim())
			vend = max(self.ui.imagewidget_mpl.canvas.ax.get_ylim())
			hstart = min(self.ui.imagewidget_mpl.canvas.ax.get_xlim())
			hend = max(self.ui.imagewidget_mpl.canvas.ax.get_xlim())
			
			self.disconnectZoomUpdates()
			self.ui.vstartpix_spin.setMaximum(vend - 1)
			self.ui.vstartpix_spin.setValue(vstart)
			
			self.ui.vendpix_spin.setMinimum(vstart + 1)
			self.ui.vendpix_spin.setValue(vend)
			
			self.ui.hstartpix_spin.setMaximum(hend - 1)
			self.ui.hstartpix_spin.setValue(hstart)
			
			self.ui.hendpix_spin.setMinimum(hstart + 1)
			self.ui.hendpix_spin.setValue(hend)
			
			self.reconnectZoomUpdates()
	
	#If we are reading a string and want an int... Make sure decimal mark is '.', etc. (just in case)
	def cleanInt(self, string): 
		return locale.atoi(string)
	
	#If we are reading a string and want a float... Make sure decimal mark is '.'
	def cleanFloat(self, string): 
		return locale.atof(string)
	
	#If input is already a float or int, return same type
	def cleanNum(self, num): 
		if (type(num) is int):
			return int(locale.format('%d', num))
		elif (type(num) is float):
			return float(locale.format('%g', num))
		else:
			raise TypeError("Warning: Unrecognised variable type being cleaned: " + type(num))
	
	# Update the image widget on the standard image tab. This needs to be called whenever things like intensity thresholds or zoomed areas change.
	def updateGraph(self, image):
		if image == []:
			return
		
		global transformationHistory
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		
		if (not self.firstCalculation):
			if (self.lim is not None): # then stay zoomed on the same spot
				[hstart, hend, vstart, vend] = self.lim
			else: # set zoom to image size
				[hstart, hend, vstart, vend] = self.ui.imagewidget_mpl.canvas.ax.axis()
		else:
			[hstart, hend, vstart, vend] = self.extent # Otherwise just set the zoom on the entire image
		
		self.ui.imagewidget_mpl.canvas.ax.clear() # Clears the canvas between draw commands.
		
		#vmin, vmax set the intensity thresholds
		#extent dictates the image size (hence the tick values on the axes)
		#aspect controls the canvas scaling: square, 1:1, etc.
		#origin decides whether (0,0) is bottom left, like graphs, or upper left, like arrays
		self.ui.imagewidget_mpl.canvas.ax.imshow(image, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = self.extent, aspect = self.aspect, origin = "upper")
		
		# If the fft size was changed, we need to correct, since the length of one axis is determined by the fft size
		if (self.FFTResChanged):
			#If an even number of rotations were done, we need to correct the xaxis, otherwise it's the yaxis
			total = 0
			for transform in transformationHistory:
				if (transform > 0): # Only sum the rotations, not the mirrors
					total = (total + transform) % 2
			if (total == 0): 
				self.ui.imagewidget_mpl.canvas.ax.set_xlim(hstart * self.FFTRatio, hend * self.FFTRatio)
				self.ui.imagewidget_mpl.canvas.ax.set_ylim(vstart, vend)
			else:
				self.ui.imagewidget_mpl.canvas.ax.set_xlim(hstart, hend)
				self.ui.imagewidget_mpl.canvas.ax.set_ylim(vstart * self.FFTRatio, vend * self.FFTRatio)
		else: # If the FFT size wasn't changed, just update as usual
			self.ui.imagewidget_mpl.canvas.ax.set_xlim(hstart, hend)
			self.ui.imagewidget_mpl.canvas.ax.set_ylim(vstart, vend)
		
		# Redraw the canvas to apply changes
		self.ui.imagewidget_mpl.canvas.draw()
		
		# Update the colourbar to match changes to image
		self.updateLegend()
	
	# Updates the colourbar beneath the image
	def updateLegend(self): # I know this is slipshod, I couldn't get colormap() to work.
		vmin, vmax = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		
		colorbar = np.linspace(vmin, vmax, num = 200)
		colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
		
		self.ui.imagelegend_mpl.canvas.ax.cla() #clear the canvas between draws
		self.ui.imagelegend_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [vmin, vmax, vmin, vmax], aspect = 0.02)
		self.ui.imagelegend_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.imagelegend_mpl.canvas.draw()
		self.ui.calculate_button.setEnabled(True)
	
	#Shows a blank canvas
	def showBlankImage(self):
		self.ui.imagelegend_mpl.canvas.ax.clear()
		self.ui.imagewidget_mpl.canvas.ax.imshow([[1,1],[1,1]], 'Greys', extent = [0,1,0,1], aspect = 'auto', origin = 'upper')
		self.ui.imagewidget_mpl.canvas.ax.set_xlim(0,1)
		self.ui.imagewidget_mpl.canvas.ax.set_ylim(0,1)
		self.ui.imagewidget_mpl.canvas.draw()
		self.firstCalculation = True
		self.ui.axisrecalculation_check.setChecked(False)
		self.ui.axisrecalculation_check.setDisabled(True)
	
	# Clear transforms in memory.
	# Since mirror/rotations should be maintained between recalculations, we need to store what transforms are applied to the current image so we can apply them again if the image is recalculated
	def resetTransformationHistory(self):
		global transformationHistory
		global settings
		global imageOCT
		imageOCT = []
		
		settings['transformationHistory_transform'] = []
		transformationHistory = []
		settings['upperthreshold_dspin'] = 1
		self.ui.upperthreshold_dspin.setValue(1)
		settings['lowerthreshold_dspin'] = 0
		self.ui.lowerthreshold_dspin.setValue(0)
		
		np.save('settings', settings)
		
		self.showBlankImage()
	
	# Get the location to cut the spectra for linearisation or dispersion compensation.
	def getCuts(self, filename): 
		start_start_loc = filename.find("_st") + 3
		start_end_loc = filename.find("_len")
		length_start_loc = filename.find("_len") + 4
		length_end_loc = filename.find(".txt")
		start = self.cleanInt(filename[start_start_loc:start_end_loc])
		length = self.cleanInt(filename[length_start_loc:length_end_loc])
		return (start, length)
	
	# Every time we do a new image calculation, we go through all the previous transformations so the image is in the same orientiation
	# This makes keeping track of axes positions and such much easier.
	def applyPreviousTransformations(self, catchup = False): 
		global transformationHistory 
		
		if (len(transformationHistory) == 0):
			return
		else:
			for transform in transformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						self.rotateImage(catchup = True)
				elif (np.sign(transform) < 0):
						self.mirrorImage(catchup = True) 
	
	# Calculate the image from the files provided
	def calculateImage(self):
		global imageOCT
		
		# Error box preparation
		FileNotFoundErrorBox = QtWidgets.QMessageBox()
		FileNotFoundErrorBox.setText("One or more files were not found. The relevant files may have been moved or might not exist.")
		FileNotFoundErrorBox.setWindowTitle("File Not Found Error")
		PermissionErrorBox = QtWidgets.QMessageBox()
		PermissionErrorBox.setText("Not all file paths were given (or read permissions was denied). Please enter all relevant file paths or deselect the relevant parameter(s).")
		PermissionErrorBox.setWindowTitle("Permission Error")
		
		# Progress bar setup
		self.ui.calculate_button.setEnabled(False)
		self.ui.std_progressbar.reset()
		self.ui.std_progressbar.show()
		self.ui.std_progressbar.setValue(1)
		
		try: #Data validation try block
			try: #File location try block
				spectra_raw = np.loadtxt(self.ui.spectra_line.text())
				
				reference = None
				phase_lin = None
				phase_nonlin = None
				dispersion_compensation = None
				
				if (self.ui.backgroundsubtraction_check.isChecked()):
					if (self.ui.referencefromfile_radio.isChecked()):
						reference = np.loadtxt(self.ui.referencespectrum_line.text())
					elif (self.ui.referencesynthesize_radio.isChecked()):
						reference = spectra_raw.mean(axis=0)
				
				if (self.ui.linearisation_check.isChecked()):
					phase_lin = np.loadtxt(self.ui.linearphase_line.text())
					phase_nonlin = np.loadtxt(self.ui.nonlinearphase_line.text())
				
				if (self.ui.dispersioncompensation_check.isChecked()):
					dispersion_compensation = np.loadtxt(self.ui.dispersioncompensation_line.text())
				
			except PermissionError: # This error occurs if the path is blank, or attempts to access a file without proper permissions
				PermissionErrorBox.exec()
				self.ui.std_progressbar.reset()
				self.ui.std_progressbar.hide()
				self.ui.calculate_button.setEnabled(True)
				return
			except FileNotFoundError: # This error occurs when a file cannot be found.
				FileNotFoundErrorBox.exec()
				self.ui.std_progressbar.reset()
				self.ui.std_progressbar.hide()
				self.ui.calculate_button.setEnabled(True)
				return
			
			spectra = spectra_raw
			self.rawSpectra = spectra # This will be useful in the subimages tab
			
			########## Calculate/Load reference ##########
			if (self.ui.backgroundsubtraction_check.isChecked()):
				spectra = spectra_raw - reference
			
			self.ui.std_progressbar.setValue(10)
			
			########## Load phases and linearise ##########
			if (self.ui.linearisation_check.isChecked()):
				linear_cut_start, linear_cut_length = self.getCuts(self.ui.linearphase_line.text())
				spectra = spectra[::, linear_cut_start : (linear_cut_start+linear_cut_length)]
				
				self.rawSpectra = self.rawSpectra[::, linear_cut_start : (linear_cut_start+linear_cut_length)] # Cut raw spectra down to size for subimages tab
				
				f = interp1d(phase_nonlin, spectra, fill_value = 'extrapolate')
				spectra = f(phase_lin)
			
			self.ui.std_progressbar.setValue(20)
			
			########## Dispersion Compensation ##########
			if (self.ui.dispersioncompensation_check.isChecked()):
				mirrorHalfCoeff = -1 if self.ui.mirrorHalf_check.isChecked() else 1
				dispersion_cut_start, dispersion_cut_length = self.getCuts(self.ui.dispersioncompensation_line.text())
				spectra = spectra[::, dispersion_cut_start : (dispersion_cut_start+dispersion_cut_length)] * np.exp(1j*mirrorHalfCoeff*dispersion_compensation)
				
				self.rawSpectra = self.rawSpectra[::, dispersion_cut_start : (dispersion_cut_start+dispersion_cut_length)] # Cut raw spectra down to size for subimages tab
			
			self.preFFTImage = spectra
			
			########## Calculate Fourier Transform ##########
			fftres = self.cleanInt(self.ui.fftres_spin.text())
			fft_data_one = np.fft.fft(spectra, n=fftres)
			imageOCT = np.absolute(fft_data_one[::, 0:int(fftres/2)])
			del fft_data_one
			if (self.ui.linearisation_check.isChecked()): 
				del phase_nonlin, phase_lin
			
			self.ui.std_progressbar.setValue(30)
			
			########## Save intermediate file for image processing ###########
			self.unnormalisedImage = imageOCT
			
			########## Normalisation ##########
			# Normalise values to a minimum of 0 and a maximum of 1
			imageOCT = imageOCT - np.amin(imageOCT)
			imageOCT = imageOCT / np.amax(imageOCT)
			
			self.ui.std_progressbar.setValue(70)
			
			########## Image Cleaning ##########
			imageOCT = imageOCT[::, 0:int(fftres/2)] # prune some pixels that we don't need
			self.ui.std_progressbar.setValue(90)
		except (ArithmeticError, LookupError, NameError, ValueError, TypeError, EOFError) as error: #Generic error in case something goes wrong in the calculation steps
			CalculationErrorBox = QtWidgets.QMessageBox()
			CalculationErrorBox.setText("There was an error when calculating the image. Please check that all input files are in the correct format.\n \nError Info:\n" + type(error).__name__ + ": " + str(error))
			CalculationErrorBox.setWindowTitle("Calculation Error")
			CalculationErrorBox.setIcon(QtWidgets.QMessageBox.Critical)
			self.ui.calculate_button.setEnabled(True)
			self.ui.std_progressbar.setValue(100)
			self.ui.std_progressbar.reset()
			self.ui.std_progressbar.hide()
			CalculationErrorBox.exec()
			return
		
		########## Miscellania ##########
		if ((fftres != self.lastFFTValue) and (self.lastFFTValue is not None) ): # Check to see if the fft size has changed since the last calculation
			self.FFTResChanged = True
			self.FFTRatio = fftres / self.lastFFTValue
		else:
			self.FFTResChanged = False
			self.FFTRatio = 1
		
		# Every time we click 'calculate image', we reset the view to default and run through all the transforms in transformationHistory
		# It makes things a lot easier in the long run, and the transformations are fairly inexpensive to compute anyway
		self.extent = [0, imageOCT.shape[1], 0, imageOCT.shape[0]]
		
		########## Display Image ##########
		#If this isn't the first calculation, then save the current zoom location to a variable that can be read after applyPreviousTransformations
		if (not self.firstCalculation): 
			self.lim = self.ui.imagewidget_mpl.canvas.ax.axis()
		
		if (self.ui.automaticthreshold_check.isChecked() == True): # Set upper threshold to a sensible number
			self.ui.upperthreshold_dspin.setValue(np.sqrt(imageOCT.mean())) #Note, must be done after self.extent is first defined
			
			
			self.ui.upperthreshold_slider.disconnect()
			self.ui.upperthreshold_slider.setValue(10**(2-2*np.sqrt(imageOCT.mean()))) # If previously stored intensity value is the same as the new calculated one, the synchronising script between the spinbox and the slider isn't called, so synchronise manually here
			self.ui.lowerthreshold_dspin.setMaximum(min(np.sqrt(imageOCT.mean()), 0.99))
			self.ui.lowerthreshold_slider.setMinimum(min(10**(2-2*np.sqrt(imageOCT.mean())), 99))
			
			self.ui.upperthreshold_slider.valueChanged.connect(self.redrawIntensity)
		
		self.applyPreviousTransformations(catchup = self.firstCalculation)
		
		self.ui.vendpix_spin.setMaximum(np.shape(imageOCT)[0])
		self.ui.hendpix_spin.setMaximum(np.shape(imageOCT)[1])
		
		# If no prior calculations have been done, then set the zoom to the whole image
		if (self.firstCalculation): 
			self.ui.hendpix_spin.setValue(np.shape(imageOCT)[1])
			self.ui.vendpix_spin.setValue(np.shape(imageOCT)[0])
			self.ui.hstartpix_spin.setValue(0)
			self.ui.vstartpix_spin.setValue(0)
		else: # just read the zoom values from the zoom fields
			self.ui.hendpix_spin.setValue(max(self.lim[0], self.lim[1]))
			self.ui.vendpix_spin.setValue(max(self.lim[2], self.lim[3]))
			
			self.ui.hstartpix_spin.setValue(min(self.lim[0], self.lim[1]))
			self.ui.vstartpix_spin.setValue(min(self.lim[2], self.lim[3]))
		
		self.updateStartZoomMax(fftres/2) # Reset zoom limts
		
		self.lastFFTValue = fftres
		self.updateGraph(imageOCT)
		self.firstCalculation = False
		self.ui.axisrecalculation_check.setDisabled(False)
		self.ui.calculate_button.setEnabled(True)
		
		self.ui.std_progressbar.setValue(100)
		self.ui.std_progressbar.reset()
		self.ui.std_progressbar.hide()
	
	# Update spinbox start maximums so you can't have nonsense values
	def updateStartZoomMax(self, val):
		if (self.sender() == self.ui.calculate_button): # If you're coming from calculation, set max of both startspinboxes
			self.ui.hstartpix_spin.setMaximum(self.ui.hendpix_spin.value() - 1)
			self.ui.vstartpix_spin.setMaximum(self.ui.vendpix_spin.value() - 1)
		
		# Else if one of the spinboxes sent you, just change it's neighbour
		if (self.sender() == self.ui.hendpix_spin): 
			self.ui.hstartpix_spin.setMaximum(max(val-1,0))
		elif ((self.sender() == self.ui.vendpix_spin) and (self.cleanNum(self.ui.vendpix_spin.value()) > 0)):
			self.ui.vstartpix_spin.setMaximum(max(val-1,0))
	
	# Similar to previous, but for setting a minimum value on the end position
	def updateEndZoomMin(self, val):
		global imageOCT
		
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		if ((self.sender() == self.ui.vstartpix_spin) and (self.cleanNum(self.ui.vstartpix_spin.value()) < imageOCT.shape[0])):
			self.ui.vendpix_spin.setMinimum(val+1)
		elif ((self.sender() == self.ui.hstartpix_spin) and (self.cleanNum(self.ui.hstartpix_spin.value()) < imageOCT.shape[1])):
			self.ui.hendpix_spin.setMinimum(val+1)
	
	def resetZoom(self): # Show the entire image
		## Because of how the imagewidgets work, every instance of an imagewidget will have it's reset zoom button go here.
		## We can use which tab is active to figure out where we should send it
		if (self.ui.tabwidget.currentIndex() == self.threeDTab):
			self.threeD_resetZoom()
			return
		
		global imageOCT
		
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		## IF Axis Recalculation is Enabled, disable it for now and re-enable it afterwards
		resetRecalc = False
		if (self.ui.axisrecalculation_check.isChecked()):
			self.axisRecalculation(False)
			resetRecalc = True
		
		# Disconnect zoom updates so we can change values in peace without prompting heaps of redraws
		self.disconnectZoomUpdates()
		
		# Reset all values to sensible ones
		self.ui.hendpix_spin.setMinimum(1)
		
		self.ui.hendpix_spin.setMaximum(imageOCT.shape[1])
		self.ui.hendpix_spin.setValue(imageOCT.shape[1])
		
		self.ui.vendpix_spin.setMinimum(1)
		self.ui.vendpix_spin.setMaximum(imageOCT.shape[0])
		self.ui.vendpix_spin.setValue(imageOCT.shape[0])
		
		self.ui.hstartpix_spin.setMaximum(self.ui.hendpix_spin.value()-1)
		self.ui.hstartpix_spin.setMinimum(0)
		self.ui.hstartpix_spin.setValue(0)
		
		self.ui.vstartpix_spin.setMaximum(self.ui.vendpix_spin.value()-1)
		self.ui.vstartpix_spin.setMinimum(0)
		self.ui.vstartpix_spin.setValue(0)
		
		self.reconnectZoomUpdates()
		
		# If axes are inverted, our limits need to go from end -> start, not start -> end
		if ( self.ui.imagewidget_mpl.canvas.ax.xaxis_inverted() ): 
			self.ui.imagewidget_mpl.canvas.ax.set_xlim(imageOCT.shape[1], 0)
		else:
			self.ui.imagewidget_mpl.canvas.ax.set_xlim(0, imageOCT.shape[1])
		
		if ( self.ui.imagewidget_mpl.canvas.ax.yaxis_inverted() ):
			self.ui.imagewidget_mpl.canvas.ax.set_ylim(imageOCT.shape[0], 0)
		else:
			self.ui.imagewidget_mpl.canvas.ax.set_ylim(0, imageOCT.shape[0])
		
		## Re-enable axis recalc if it was toggled off before
		if resetRecalc:
			self.axisRecalculation(True)
		
		# Prompt a redraw so changes take effect
		self.ui.imagewidget_mpl.canvas.draw()
	
	# Whenever a transform is applied to the image, we need to save it here so a recalculation can be set to the same orientation
	def updateTransformationHistory(self, type):
		global transformationHistory
		
		# Mirrors are saved as -1, rotations are 1,2,3 where the magnitude is the number of 90 degree clockwise rotations
		# i.e [-1, 2, -1] is a mirror, a 180 degree rotation, and then a mirror again.
		if (type == 'mirror'):
			sgn = -1
		elif (type == 'rotate'):
			sgn = 1
		else:
			print ('Error - invalid transform type parameter. Should be "mirror" or "rotate" only')
			return
		
		# Try to combine transforms together i.e [1, 1, 1] = [3], [3, 1] = [4] = [0] = []
		if (len(transformationHistory) > 0) and (np.sign(transformationHistory[-1]) == sgn):
			transformationHistory[-1] += sgn
			if (sgn == 1) and (transformationHistory[-1] % 4 == 0):
				transformationHistory = transformationHistory[:-1]
			elif (sgn == -1) and (transformationHistory[-1] %2 == 0):
				transformationHistory = transformationHistory[:-1]
		else: # just append if we can't simplify in any way
			transformationHistory.append(sgn)
		
		# Save transformation history to settings file
		settings['transformationHistory_transform'] = transformationHistory
		np.save('settings', settings)
	
	#The rotation code is quite complex. The image, axes, and start/end pixel posiitons need to be constantly synchronised.
	def rotateImage(self, catchup = False): 
		global imageOCT
		
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		# Store currently viewed area
		limits = self.ui.imagewidget_mpl.canvas.ax.axis() 
		
		# Reset zoom and clear image
		self.resetZoom()
		self.ui.imagewidget_mpl.canvas.ax.clear()
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		
		# Rotate image clockwise 
		imageOCT = np.rot90(imageOCT, -1)
		
		#Rotate axes with image
		self.extent = [self.extent[2], self.extent[3], self.extent[1], self.extent[0]] 
		
		self.ui.imagewidget_mpl.canvas.ax.imshow(imageOCT, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.extent, aspect = self.aspect, origin = "upper")
		
		if (not catchup):
			self.updateTransformationHistory('rotate')
		
		# Not sure why, but we need to reset the zoom again so it displays properly
		self.resetZoom() 
		
		########## Re-zoom to previous spot ##########
		limits = [limits[2], limits[3], limits[1], limits[0]] #Rotate viewing area with image
		self.ui.imagewidget_mpl.canvas.ax.set_xlim(limits[0], limits[1]) # Reapply limits to viewing area
		self.ui.imagewidget_mpl.canvas.ax.set_ylim(limits[2], limits[3])
		
		self.updatePixValues(limits)
		
		self.ui.imagewidget_mpl.canvas.draw()
	
	def updatePixValues(self, limits):
		
		self.disconnectZoomUpdates()
		
		self.ui.vendpix_spin.setMinimum(min(limits[2], limits[3]))
		self.ui.vendpix_spin.setMaximum(max(self.extent[2], self.extent[3]))
		self.ui.vendpix_spin.setValue(max(limits[2], limits[3]))
		
		self.ui.hendpix_spin.setMinimum(min(limits[0], limits[1]))
		self.ui.hendpix_spin.setMaximum(max(self.extent[0], self.extent[1]))
		self.ui.hendpix_spin.setValue(max(limits[0], limits[1]))
		
		self.ui.vstartpix_spin.setMaximum(max(limits[2], limits[3]))
		self.ui.vstartpix_spin.setValue(min(limits[2], limits[3]))
		
		self.ui.hstartpix_spin.setMaximum(max(limits[0], limits[1]))
		self.ui.hstartpix_spin.setValue(min(limits[0], limits[1]))
		
		
		self.reconnectZoomUpdates()
	
	def mirrorImage(self, catchup = False):
		global imageOCT
		
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		# Store currently viewed area
		limits = self.ui.imagewidget_mpl.canvas.ax.axis()
		
		# Reset zoom and clear image
		self.resetZoom()
		self.ui.imagewidget_mpl.canvas.ax.clear()
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		
		# Flip image left/right
		imageOCT = np.fliplr(imageOCT)
		
		#Flip axes with image
		self.extent = [self.extent[1], self.extent[0], self.extent[2], self.extent[3]]
		
		self.ui.imagewidget_mpl.canvas.ax.imshow(imageOCT, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.extent, aspect = self.aspect, origin = "upper")
		if (not catchup):
			self.updateTransformationHistory('mirror')
		
		########## Re-zoom to previous spot ##########
		limits = [limits[1], limits[0], limits[2], limits[3]] #Rotate viewing area with image
		self.ui.imagewidget_mpl.canvas.ax.set_xlim(limits[0], limits[1])
		self.ui.imagewidget_mpl.canvas.ax.set_ylim(limits[2], limits[3])
		
		self.updatePixValues(limits)
		
		self.ui.imagewidget_mpl.canvas.draw()
	
	def zoom(self):
		global imageOCT
		
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		# Get limits from spinboxes
		vstart, vend = self.cleanInt(self.ui.vstartpix_spin.text()), self.cleanInt(self.ui.vendpix_spin.text()) 
		hstart, hend = self.cleanInt(self.ui.hstartpix_spin.text()), self.cleanInt(self.ui.hendpix_spin.text())
		
		# Check for inverted axes, as this inverts the view limits too.
		# Zoom image accordingly
		if (self.ui.imagewidget_mpl.canvas.ax.xaxis_inverted()): 
			self.ui.imagewidget_mpl.canvas.ax.set_xlim(hend, hstart)
		else:
			self.ui.imagewidget_mpl.canvas.ax.set_xlim(hstart, hend) 
		
		if (self.ui.imagewidget_mpl.canvas.ax.yaxis_inverted()):
			self.ui.imagewidget_mpl.canvas.ax.set_ylim(vend, vstart)
		else:
			self.ui.imagewidget_mpl.canvas.ax.set_ylim(vstart, vend)
		
		#Redraw so changes take effect
		self.ui.imagewidget_mpl.canvas.draw() 
	
	def redrawIntensity(self, val):
		global imageOCT
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		# Synchronise slider and spinbox
		# This is a real mess, made worse by the fact that the sliders are logarithmic and the spinboxes are linear
		if (self.sender() == self.ui.upperthreshold_dspin):
			# Set minimum value of the other pair to correct value
			self.ui.lowerthreshold_dspin.setMaximum(min(val, 0.99))
			self.ui.lowerthreshold_slider.setMinimum(min(10**(2-2*val), 99))
			
			# Try to set correct value on matching slider
			self.ui.upperthreshold_slider.disconnect()
			try:
				self.ui.upperthreshold_slider.setValue(10**(2-2*self.ui.upperthreshold_dspin.value()))
			except ValueError: # obviously this was important at some point, but I can't see why this is here.
				pass
			self.ui.upperthreshold_slider.valueChanged.connect(self.redrawIntensity)
			
		elif (self.sender() == self.ui.upperthreshold_slider):
			self.ui.lowerthreshold_dspin.setMaximum(min(1-np.log10(val)/2, 0.99))
			self.ui.lowerthreshold_slider.setMinimum(min(val, 99))
			
			self.ui.upperthreshold_dspin.disconnect()
			try:
				self.ui.upperthreshold_dspin.setValue(1-np.log10(self.ui.upperthreshold_slider.value())/2)
			except ValueError:
				pass
			self.ui.upperthreshold_dspin.valueChanged.connect(self.redrawIntensity)
			
		elif (self.sender() == self.ui.lowerthreshold_dspin):
			self.ui.upperthreshold_dspin.setMinimum(max(val, 0.01))
			self.ui.upperthreshold_slider.setMaximum(max(10**(2-2*val), 1))
			
			self.ui.lowerthreshold_slider.disconnect()
			try:
				self.ui.lowerthreshold_slider.setValue(10**(2-2*self.ui.lowerthreshold_dspin.value()))
			except ValueError:
				pass
			self.ui.lowerthreshold_slider.valueChanged.connect(self.redrawIntensity)
			
		elif (self.sender() == self.ui.lowerthreshold_slider):
			self.ui.upperthreshold_dspin.setMinimum(max(1-np.log10(val)/2, 0.01))
			self.ui.upperthreshold_slider.setMaximum(max(val, 1))
			
			self.ui.lowerthreshold_dspin.disconnect()
			try:
				self.ui.lowerthreshold_dspin.setValue(1-np.log10(self.ui.lowerthreshold_slider.value())/2)
			except ValueError:
				pass
			self.ui.lowerthreshold_dspin.valueChanged.connect(self.redrawIntensity)
		self.updateGraph(imageOCT)
	
	# Reset intensity to [0, 1]
	def resetIntensity(self):
		global imageOCT
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		self.ui.upperthreshold_dspin.setValue(1)
		self.ui.lowerthreshold_dspin.setValue(0)
	
	def updateIntensityMax(self, val):
		pass #old function, shouldn't be used anywhere
	
	def updateIntensityMin(self, val):
		pass #old function, shouldn't be used anywhere
	
	# When parent checkbox is changed, enable/disable child objects
	def toggleRadio(self):
		if self.sender() == self.ui.backgroundsubtraction_check:
			if (self.ui.backgroundsubtraction_check.isChecked()):
				self.ui.referencefromfile_radio.setDisabled(False)
				self.ui.referencesynthesize_radio.setDisabled(False)
			elif (not self.ui.backgroundsubtraction_check.isChecked()):
				self.ui.referencefromfile_radio.setDisabled(True)
				self.ui.referencesynthesize_radio.setDisabled(True)
		elif self.sender() == self.ui.dispersioncompensation_check:
			if (self.ui.dispersioncompensation_check.isChecked()):
				self.ui.mirrorHalf_check.setDisabled(False)
			elif (not self.ui.dispersioncompensation_check.isChecked()):
				self.ui.mirrorHalf_check.setDisabled(True)
	
	# Store filenames in settings file
	def updateSettingsFileFilename(self, filename):
		global settings
		origin = str(self.sender().objectName())[0:-5] # This gives us a string like 'dispersioncompensation'
		settings[origin + '_line'] = filename
		np.save('settings', settings)
	
	# If a filename is changed, turn on the related parameter automatically
	def autoUpdateParameters(self, filename):
		if (self.sender() == self.ui.dispersioncompensation_line):
			if (filename == ''):
				self.ui.dispersioncompensation_check.setChecked(False)
				return
			self.ui.dispersioncompensation_check.setChecked(True)
		
		elif (self.sender() == self.ui.referencespectrum_line):
			if (filename == ''):
				self.ui.referencesynthesize_radio.setChecked(True)
				return
			self.ui.backgroundsubtraction_check.setChecked(True)
			self.ui.referencefromfile_radio.setChecked(True)
		
		elif (self.sender() == self.ui.linearphase_line or self.sender() == self.ui.nonlinearphase_line):
			if (filename == ''):
				self.ui.linearisation_check.setChecked(False)
				return
			self.ui.linearisation_check.setChecked(True)
	
	# Open file dialog
	def getFile(self):
		options = QtWidgets.QFileDialog.Options()
		options |= QtWidgets.QFileDialog.DontUseNativeDialog
		filename, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Open File", "","Text Files (*.txt)", options=options)
		if (filename == ""):
			return
		
		origin = str(self.sender().objectName())[0:-5] 
		eval('self.ui.' + origin + '_line.setText(filename)') # This does open up arbitrary code execution via fname. I doubt it would be a problem though.
		self.updateSettingsFileFilename(filename);
	
	# Open folder select dialog
	def getFolder(self):
		filename = QtWidgets.QFileDialog.getExistingDirectory(self,"Select Folder", "", QtWidgets.QFileDialog.ShowDirsOnly)
		if (filename == ""):
			return
		
		origin = str(self.sender().objectName())[0:-5] 
		eval('self.ui.' + origin + '_line.setText(filename)') # This does open up arbitrary code execution via fname. I doubt it would be a problem though.
		self.updateSettingsFileFilename(filename);
	
	# Update settings file when a checkbox is changed
	def updateSettingsFileCheckbox(self):
		global settings
		origin = str(self.sender().objectName())
		settings[origin] = self.sender().isChecked()
		np.save('settings', settings)
	
	# Update settings file when a spinbox is changed
	def updateSettingsFileDspin(self):
		global settings
		origin = str(self.sender().objectName())
		settings[origin] = self.sender().value()
		np.save('settings', settings)
	
	# Saves current image to a png
	def saveAsPNG(self): 
		SaveErrorBox = QtWidgets.QMessageBox()
		SaveErrorBox.setText("I forgot to code this button properly. Nothing was saved to file.")
		SaveErrorBox.setIcon(QtWidgets.QMessageBox.Critical)
		SaveErrorBox.setWindowTitle("Error")
		
		fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', "", "PNG File  (*.png)")
		if (fname == ''):
			return
		
		
		#Figure out which save button was pressed to figure out what image to save. I'll make this better eventually. A similar setup to getFile could be used instead of a mountain of else if's
		sending_button = self.sender().objectName()
		
		if (sending_button.find("standardimage") >= 0): 
			self.ui.imagewidget_mpl.canvas.figure.savefig(fname, bbox_inches = 'tight')
		elif(sending_button.find("processing") >= 0):
			self.ui.processing_imagewidget_mpl.canvas.figure.savefig(fname, bbox_inches = 'tight')
		elif (sending_button.find("threeD") >= 0):
			self.ui.threeD_imagewidget_mpl.canvas.figure.savefig(fname, bbox_inches = 'tight')
		# elif (sending_button.find("cornea_opticalthickness") >= 0): 
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_geometricalthickness") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_walkoff") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_beta2") >= 0):
			# #Do a thing
			# return
			
		# elif (sending_button.find("aqueoushumor_opticalthickness") >= 0): 
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_geometricalthickness") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_walkoff") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_beta2") >= 0):
			# #Do a thing
			# return
		else:
			SaveErrorBox.exec()
	# Saves current image to a text file
	def saveAsTXT(self):
		SaveErrorBox = QtWidgets.QMessageBox()
		SaveErrorBox.setText("I forgot to code this button properly. Nothing was saved to file.")
		SaveErrorBox.setIcon(QtWidgets.QMessageBox.Critical)
		SaveErrorBox.setWindowTitle("Error")
		
		fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', "", "Text File  (*.txt)")
		if (fname == ''):
			return
		
		#Figure out which save button was pressed to figure out what image to save I'll make htis better eventually
		sending_button = self.sender().objectName()
		
		if (sending_button.find("standardimage") >= 0): 
			np.savetxt(fname, imageOCT, delimiter = '\t')
		elif(sending_button.find("processing") >= 0):
			np.savetxt(fname, self.currentProcessingImage, delimiter = '\t')
		elif (sending_button.find("threeD") >= 0):
			np.savetxt(fname, self.threeD_currentImage, delimiter = '\t')
		# elif (sending_button.find("cornea_opticalthickness") >= 0): 
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_geometricalthickness") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_walkoff") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("cornea_beta2") >= 0):
			# #Do a thing
			# return
			
		# elif (sending_button.find("aqueoushumor_opticalthickness") >= 0): 
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_geometricalthickness") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_walkoff") >= 0):
			# #Do a thing
			# return
		# elif (sending_button.find("aqueoushumor_beta2") >= 0):
			# #Do a thing
			# return
		else:
			SaveErrorBox.exec()
	
	# Change Colourmap on images
	def changeColourmap(self, index):
		global imageOCT
		
		oldColourmap = self.colourmap
		
		if (index == 0):
			self.colourmap = 'Greys_r'
		elif (index == 1):
			self.colourmap = 'viridis'
		elif (index == 2):
			self.colourmap = 'plasma'
		
		if (self.sender() == self.ui.processing_colourmap_combo):
			self.ui.colourmap_combo.setCurrentIndex(index)
			self.ui.threeD_colourmap_combo.setCurrentIndex(index)
		elif (self.sender() == self.ui.colourmap_combo):
			self.ui.processing_colourmap_combo.setCurrentIndex(index)
			self.ui.threeD_colourmap_combo.setCurrentIndex(index)
		elif (self.sender() == self.ui.threeD_colourmap_combo):
			self.ui.processing_colourmap_combo.setCurrentIndex(index)
		
		try:
			if (self.colourmap == oldColourmap):
				return
		except NameError:
			return
		self.lim = self.ui.imagewidget_mpl.canvas.ax.axis()
		
		self.updateGraph(imageOCT)
		self.processing_updateGraph(self.currentProcessingImage)
		self.threeD_updateGraph()
		
		try:
			if (self.subimage1 != [] and self.subimage2 != []):
				self.subimages_updateGraph(updateBoth = True)
		except NameError:
			pass
	
	# Convert axes from pixels into physical dimensions
	def axisRecalculation(self, state):
		global imageOCT
		global transformationHistory
		
		try:
			if (self.extent is None or self.extent == []):
				self.ui.axisrecalculation_check.setChecked(False)
				return
		except Exception:
			self.ui.axisrecalculation_check.setChecked(False)
			return
		
		if (state == 0): # state signal comes through as an int, convert it to bool for convenience
			state = False
		else:
			state = True
		
		#Disallow changing refr_ind and vpp while currently applied because it's a nightmare
		self.ui.refrindex_dspin.setDisabled(state)
		self.ui.vpp_dspin.setDisabled(state)
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		refrIndex = self.ui.refrindex_dspin.value()
		vpp = self.ui.vpp_dspin.value()
		tempExtent = self.extent
		
		# choose whether we're doing or undoing recalculation
		if (state):
			longMult = 0.0025125/refrIndex # forward transform
			shortMult = (0.053*vpp)/2
			self.aspect = "equal"
		else:
			longMult = refrIndex/0.0025125 # backward transform
			shortMult = (2*vpp)/0.053
			self.aspect = "auto"
		
		# If not zoomed in, take whole image as the zoomed area
		#if(self.lim is None):
		#	tempLim = self.extent[:]
		#else:
		tempLim = list(self.ui.imagewidget_mpl.canvas.ax.axis()[:])
		
		######### Apply axis scaling
		if ( (len(transformationHistory) == 0) ): # undo rotations so we know which is the long side, etc
			pass
		else:
			for transform in transformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						tempExtent = [tempExtent[3], tempExtent[2], tempExtent[0], tempExtent[1]]
						tempLim = [tempLim[3], tempLim[2], tempLim[0], tempLim[1]]
				elif (np.sign(transform) < 0):
						tempExtent = [tempExtent[1], tempExtent[0], tempExtent[2], tempExtent[3]]
						tempLim = [tempLim[1], tempLim[0], tempLim[2], tempLim[3]]
		
		
		[tempExtent[0], tempExtent[1]] = [tempExtent[0] * longMult, tempExtent[1] * longMult]
		[tempLim[0], tempLim[1]] = [tempLim[0] * longMult, tempLim[1] * longMult]
		
		[tempExtent[2], tempExtent[3]] = [tempExtent[2] * shortMult, tempExtent[3] * shortMult]
		[tempLim[2], tempLim[3]] = [tempLim[2] * shortMult, tempLim[3] * shortMult]
		
		
		if ( (len(transformationHistory) == 0) ): # redo rotations
			pass
		else:
			for transform in transformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						tempExtent = [tempExtent[2], tempExtent[3], tempExtent[1], tempExtent[0]]
						tempLim = [tempLim[2], tempLim[3], tempLim[1], tempLim[0]]
				elif (np.sign(transform) < 0):
						tempExtent = [tempExtent[1], tempExtent[0], tempExtent[2], tempExtent[3]]
						tempLim = [tempLim[1], tempLim[0], tempLim[2], tempLim[3]]
		
		self.extent = tempExtent
		self.lim = tempLim
		
		self.ui.imagewidget_mpl.canvas.ax.clear()
		self.ui.imagewidget_mpl.canvas.ax.imshow(imageOCT, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.extent, aspect = self.aspect, origin = "upper")
		self.ui.imagewidget_mpl.canvas.ax.set_xlim(self.lim[0], self.lim[1]) 
		self.ui.imagewidget_mpl.canvas.ax.set_ylim(self.lim[2], self.lim[3])
		
		self.updatePixValues(self.lim)
		
		self.ui.imagewidget_mpl.canvas.draw()
	
	def processing_updateIntensityMin(self):
		pass
	
	def processing_updateIntensityMax(self):
		pass
	
	# Similar to redrawing intensity in standard image tab. See comments there
	def processing_redrawIntensity(self, val):
		### Redraw image and synchronise spinbox and slider
		if (self.sender() == self.ui.processing_upperthreshold_dspin):
			self.ui.processing_lowerthreshold_dspin.setMaximum(min(val, 0.999))
			self.ui.processing_lowerthreshold_slider.setMinimum(min(10**(3-3*val), 999))
			
			self.ui.processing_upperthreshold_slider.disconnect()
			try:
				self.ui.processing_upperthreshold_slider.setValue(10**(3-3*self.ui.processing_upperthreshold_dspin.value()))
			except ValueError:
				pass
			self.ui.processing_upperthreshold_slider.valueChanged.connect(self.processing_redrawIntensity)
			
		elif (self.sender() == self.ui.processing_upperthreshold_slider):
			self.ui.processing_lowerthreshold_dspin.setMaximum(min(1-np.log10(val)/3, 0.999))
			self.ui.processing_lowerthreshold_slider.setMinimum(min(val, 999))
			
			self.ui.processing_upperthreshold_dspin.disconnect()
			try:
				self.ui.processing_upperthreshold_dspin.setValue(1-np.log10(self.ui.processing_upperthreshold_slider.value())/3)
			except ValueError:
				pass
			self.ui.processing_upperthreshold_dspin.valueChanged.connect(self.processing_redrawIntensity)
			
		elif (self.sender() == self.ui.processing_lowerthreshold_dspin):
			self.ui.processing_upperthreshold_dspin.setMinimum(max(val, 0.001))
			self.ui.processing_upperthreshold_slider.setMaximum(max(10**(3-3*val), 1))
			
			self.ui.processing_lowerthreshold_slider.disconnect()
			try:
				self.ui.processing_lowerthreshold_slider.setValue(10**(3-3*self.ui.processing_lowerthreshold_dspin.value()))
			except ValueError:
				pass
			self.ui.processing_lowerthreshold_slider.valueChanged.connect(self.processing_redrawIntensity)
			
		elif (self.sender() == self.ui.processing_lowerthreshold_slider):
			self.ui.processing_upperthreshold_dspin.setMinimum(max(1-np.log10(val)/3, 0.001))
			self.ui.processing_upperthreshold_slider.setMaximum(max(val, 1))
			
			self.ui.processing_lowerthreshold_dspin.disconnect()
			try:
				self.ui.processing_lowerthreshold_dspin.setValue(1-np.log10(self.ui.processing_lowerthreshold_slider.value())/3)
			except ValueError:
				pass
			self.ui.processing_lowerthreshold_dspin.valueChanged.connect(self.processing_redrawIntensity)
		self.processing_updateGraph(self.currentProcessingImage)
	
	# Copy the image currently displayed in the standard image tab into the processing tab
	def processing_copyImage(self):
		global imageOCT
		global transformationHistory
		try:
			if (imageOCT == [] or self.unnormalisedImage == []):
				return
		except (NameError, AttributeError):
			return
		
		#### Synchronise axis recalculation between standard image and image processing
		self.processingTransformationHistory = transformationHistory # store rotation/mirror state of current image so we can do axis recalculation later
		if (self.ui.axisrecalculation_check.isChecked()):
			axisState = True
		else:
			axisState = False
		
		self.processingExtent = self.extent
		self.ui.processing_refrindex_dspin.setValue(self.ui.refrindex_dspin.value())
		self.ui.processing_vpp_dspin.setValue(self.ui.vpp_dspin.value())
		self.processingAspect = self.aspect
		
		self.ui.processing_axisrecalculation_check.disconnect()
		self.ui.processing_axisrecalculation_check.setChecked(axisState)
		self.ui.processing_axisrecalculation_check.stateChanged.connect(self.processing_axisRecalculation)
		self.ui.processing_axisrecalculation_check.setDisabled(False)
		
		### Copy values from standard image tab
		limits = self.ui.imagewidget_mpl.canvas.ax.axis()
		self.copiedLimits = limits # Used to restore limits after removing last x pixels
		self.processing_transformationHistory = transformationHistory
		xMin = self.ui.hstartpix_spin.value()
		xMax = self.ui.hendpix_spin.value()
		yMin = self.ui.vstartpix_spin.value()
		yMax = self.ui.vendpix_spin.value()
		
		currentImage = self.unnormalisedImage
		
		if (len(transformationHistory) == 0):
			pass
		else:
			for transform in transformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						currentImage = np.rot90(currentImage, -1)
				elif (np.sign(transform) < 0):
						currentImage = np.fliplr(currentImage)
		
		currentImage = currentImage[yMin:yMax, xMin:xMax]
		
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.lowerthreshold_dspin.text()), self.cleanFloat(self.ui.upperthreshold_dspin.text())
		
		self.ui.processing_lowerthreshold_dspin.setValue(lowerThreshold)
		self.ui.processing_upperthreshold_dspin.setValue(upperThreshold)
		
		self.processing_imageToBeEnhanced = currentImage
		
		displayImage = currentImage - np.amin(currentImage)
		displayImage = displayImage / np.amax(displayImage)
		
		self.currentProcessingImage = displayImage
		
		## Disconnect and update image widget
		self.ui.processing_imagewidget_mpl.canvas.mpl_disconnect(self.processing_updateCid)
		
		self.ui.processing_imagewidget_mpl.canvas.ax.clear()
		self.ui.processing_imagewidget_mpl.canvas.ax.imshow(displayImage, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = limits, aspect = self.processingAspect, origin = "upper")
		self.ui.processing_imagewidget_mpl.canvas.draw()
		
		self.processing_updateCid = self.ui.processing_imagewidget_mpl.canvas.mpl_connect('draw_event', self.processing_cropImageFromZoom)
		
		## Colour bar stuff ##
		colorbar = np.linspace(lowerThreshold, upperThreshold, num = 200)
		colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
		
		self.ui.processing_imagelegend_mpl.canvas.ax.clear()
		self.ui.processing_imagelegend_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [lowerThreshold, upperThreshold, lowerThreshold, upperThreshold], aspect = 0.02)
		self.ui.processing_imagelegend_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.processing_imagelegend_mpl.canvas.draw()
		
		# Clear flag so we know that currently displayed image has only been copied.
		self.processingFlag = False
	
	def processing_enhanceDeeperLayers(self, resetThresholds = True,*args):
		global imageOCT
		global transformationHistory
		
		if (not self.ui.processing_enhancedeeperlayers_check.isChecked()):
			return
		
		try:
			if (imageOCT == [] or self.processing_imageToBeEnhanced == []):
				return
		except (NameError, AttributeError):
			return
		
		## If axis recalculation is on, turn off for the mean time
		axisRecalculated = self.ui.processing_axisrecalculation_check.isChecked()
		if ( axisRecalculated == True):
			self.ui.processing_axisrecalculation_check.setChecked(False)
		
		self.ui.processing_progressbar.setValue(0)
		self.ui.processing_progressbar.show()
		
		# Since the function is called before the mouseRelease event is handled on the spinbox, the spinbox thinks that the spinbox arrow is still being held down while the function is being run and calls valueChanged again, so we disable it while the function is doing it's thing
		self.ui.processing_removefirstpixels.setDisabled(True) 
		self.ui.processing_removelastpixels.setDisabled(True)
		
		n = 2
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.processing_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.processing_upperthreshold_dspin.text())
		
		limits = list(self.copiedLimits)
		
		image = self.processing_imageToBeEnhanced
		
		### Undo rotations before applying algorithm
		if ( (len(transformationHistory) == 0) ): 
			pass
		else:
			for transform in transformationHistory[::-1]:
				if (np.sign(transform) > 0):
					for i in range(transform):
						image = np.rot90(image, 1)
						limits = [limits[3], limits[2], limits[0], limits[1]]
				elif (np.sign(transform) < 0):
						image = np.fliplr(image)
						limits = [limits[1], limits[0], limits[2], limits[3]]
		self.ui.processing_progressbar.setValue(20)
		
		### Crop image before algorithm
		cropStart = self.ui.processing_removefirstpixels.value()
		cropEnd = self.ui.processing_removelastpixels.value()
		if ((cropStart != 0) or (cropEnd != 0)):
			limits[0] += cropStart
			limits[1] -= cropEnd
			if cropEnd != 0:
				image = image[:, cropStart:-cropEnd]
			else:
				image = image[:, cropStart:]
		self.ui.processing_progressbar.setValue(30)
		
		### Apply algorithm
		image_enhanced = []
		for row in image:
			k = 0
			row2 = []
			for pixel in row:
				end = int(len(row))
				
				pix2 = ( pixel / ((2 * (end - k) )**n) )
				row2.append(pix2)
				k = k + 1
				
			image_enhanced.append(row2)
		
		### Crop Image After algorithm
		image_enhanced = np.array(image_enhanced)
		image_enhanced = image_enhanced[:, 0:-600]
		self.ui.processing_progressbar.setValue(70)
		
		### Redo rotations after applying algorithm
		if ( (len(transformationHistory) == 0) ):
			pass
		else:
			for transform in transformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						image_enhanced = np.rot90(image_enhanced, -1)
						limits = [limits[2], limits[3], limits[1], limits[0]]
				elif (np.sign(transform) < 0):
						image_enhanced = np.fliplr(image_enhanced)
						limits = [limits[1], limits[0], limits[2], limits[3]]
		
		if ((self.ui.processing_automaticthreshold_check.isChecked() == True) and (resetThresholds == True or self.processingFlag == False)):
			upperThreshold = upperThreshold / 10
			self.ui.processing_upperthreshold_dspin.setValue(upperThreshold)
		
		## Update image widget
		self.ui.processing_imagewidget_mpl.canvas.mpl_disconnect(self.processing_updateCid)
		
		self.ui.processing_imagewidget_mpl.canvas.ax.clear()
		self.ui.processing_imagewidget_mpl.canvas.ax.imshow(image_enhanced, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = limits, aspect = "auto", origin = "upper")
		self.ui.processing_imagewidget_mpl.canvas.draw()
		
		self.processing_updateCid = self.ui.processing_imagewidget_mpl.canvas.mpl_connect('draw_event', self.processing_cropImageFromZoom)
		
		## Tidy-up
		self.currentProcessingImage = image_enhanced
		
		self.ui.processing_progressbar.setValue(90)
		self.ui.processing_progressbar.hide()
		
		self.ui.processing_removefirstpixels.setDisabled(False)
		self.ui.processing_removelastpixels.setDisabled(False)
		self.processingExtent = limits
		
		## Re-enable axis recalculation if it was on beforehand
		if ( axisRecalculated == True):
			self.ui.processing_axisrecalculation_check.setChecked(True)
		
		self.ui.processing_imagewidget_mpl.mpl_toolbar._active = None
		
		#Set flag so we know the currently displayed image has been processed.
		self.processingFlag = True
	
	# Similar to updateGraph, see comments there
	def processing_updateGraph(self, image):
		try:
			if (self.currentProcessingImage == []):
				return
		except (NameError, AttributeError):
			return
		
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.processing_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.processing_upperthreshold_dspin.text())
		limits = self.ui.processing_imagewidget_mpl.canvas.ax.axis()
		
		self.ui.processing_imagewidget_mpl.canvas.mpl_disconnect(self.processing_updateCid)
		
		self.ui.processing_imagewidget_mpl.canvas.ax.clear() # Clears the canvas between draw commands.
		
		#vmin, vmax set the intensity thresholds
		#extent dictates the image size
		#aspect controls the canvas scaling, square, 1:1, etc.
		#origin decides whether (0,0) is bottom left, like graphs, or upper left, like arrays
		self.ui.processing_imagewidget_mpl.canvas.ax.imshow(self.currentProcessingImage, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = limits, aspect = "auto", origin = "upper")
		self.ui.processing_imagewidget_mpl.canvas.draw()
		
		self.processing_updateCid = self.ui.processing_imagewidget_mpl.canvas.mpl_connect('draw_event', self.synchroniseZoom)
		
		##### Colour bar #####
		colorbar = np.linspace(lowerThreshold, upperThreshold, num = 200)
		colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
		
		self.ui.processing_imagelegend_mpl.canvas.ax.clear()
		self.ui.processing_imagelegend_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [lowerThreshold, upperThreshold, lowerThreshold, upperThreshold], aspect = 0.02)
		self.ui.processing_imagelegend_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.processing_imagelegend_mpl.canvas.draw()
	
	# Similar to cropImageFromZoom, see comments there.
	def processing_cropImageFromZoom(self, event):
		try:
			if (self.currentProcessingImage == []):
				return
		except (NameError, AttributeError):
			return
		
		if (self.ui.processing_imagewidget_mpl.mpl_toolbar._active != 'ZOOM'):
			return
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.processing_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.processing_upperthreshold_dspin.text())
		limits = list(self.ui.processing_imagewidget_mpl.canvas.ax.axis())
		originalLimits = limits
		
		if (len(self.processing_transformationHistory) > 0): 
			rotateState = 0
			for transform in self.processing_transformationHistory:
				if (transform > 0):
					rotateState += transform
			if (rotateState % 2 == 1): # If the x and y axes have been swapped from image rotations, unswap them
				[limits[0], limits[1], limits[2], limits[3]] = [limits[2], limits[3], limits[0], limits[1]]
		
		if (self.ui.processing_imagewidget_mpl.canvas.ax.yaxis_inverted()):
			[limits[0], limits[1], limits[2], limits[3]] = [limits[1], limits[0], limits[2], limits[3]]
		if (self.ui.processing_imagewidget_mpl.canvas.ax.xaxis_inverted()):
			[limits[0], limits[1], limits[2], limits[3]] = [limits[0], limits[1], limits[3], limits[2]]
		
		self.currentProcessingImage = self.currentProcessingImage[limits[0]:limits[1], limits[2]:limits[3]]
		self.processing_imageToBeEnhanced = self.processing_imageToBeEnhanced[limits[0]:limits[1], limits[2]:limits[3]]
		
		self.ui.processing_imagewidget_mpl.canvas.mpl_disconnect(self.processing_updateCid)
		
		self.ui.processing_imagewidget_mpl.canvas.ax.clear()
		self.ui.processing_imagewidget_mpl.canvas.ax.imshow(self.currentProcessingImage, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = originalLimits, aspect = "auto", origin = "upper")
		self.ui.processing_imagewidget_mpl.canvas.draw()
		
		self.processing_updateCid = self.ui.processing_imagewidget_mpl.canvas.mpl_connect('draw_event', self.processing_cropImageFromZoom)
	
	def processing_toggleProcessImageButton(self, state):
		if (state):
			self.ui.processing_processimage_button.setDisabled(False)
		else:
			self.ui.processing_processimage_button.setDisabled(True)
	
	def processing_axisRecalculation(self, state):
		global transformationHistory
		try:
			if (self.extent is None or self.extent == []):
				self.ui.processing_axisrecalculation_check.setChecked(False)
				return
		except Exception:
			self.ui.processing_axisrecalculation_check.setChecked(False)
			return
		if (state == 0): # state signal comes through as an int, convert it to bool for convenience
			state = False
		else:
			state = True
		
		#Disallow changing refr_ind and vpp while currently applied because it's a nightmare
		self.ui.processing_refrindex_dspin.setDisabled(state)
		self.ui.processing_vpp_dspin.setDisabled(state)
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.processing_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.processing_upperthreshold_dspin.text())
		refrIndex = self.ui.processing_refrindex_dspin.value()
		vpp = self.ui.processing_vpp_dspin.value()
		tempExtent = self.processingExtent
		
		# choose whether we're doing or undoing recalculation
		if (state):
			longMult = 0.0025125/refrIndex # forward transform
			shortMult = (0.053*vpp)/2
			self.processingAspect = "equal"
		else:
			longMult = refrIndex/0.0025125 # backward transform
			shortMult = (2*vpp)/0.053
			self.processingAspect = "auto"
		
		# If not zoomed in, take whole image as the zoomed area
		if(self.lim is None):
			tempLim = self.processingExtent[:]
		else:
			tempLim = self.copiedLimits[:]
		
		######### Apply axis scaling
		if ( (len(self.processingTransformationHistory) == 0) ): # undo rotations so we know which is the long side, etc
			pass
		else:
			for transform in self.processingTransformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						tempExtent = [tempExtent[3], tempExtent[2], tempExtent[0], tempExtent[1]]
						tempLim = [tempLim[3], tempLim[2], tempLim[0], tempLim[1]]
				elif (np.sign(transform) < 0):
						tempExtent = [tempExtent[1], tempExtent[0], tempExtent[2], tempExtent[3]]
						tempLim = [tempLim[1], tempLim[0], tempLim[2], tempLim[3]]
		
		
		[tempExtent[0], tempExtent[1]] = [tempExtent[0] * longMult, tempExtent[1] * longMult]
		[tempLim[0], tempLim[1]] = [tempLim[0] * longMult, tempLim[1] * longMult]
		
		[tempExtent[2], tempExtent[3]] = [tempExtent[2] * shortMult, tempExtent[3] * shortMult]
		[tempLim[2], tempLim[3]] = [tempLim[2] * shortMult, tempLim[3] * shortMult]
		
		
		if ( (len(self.processingTransformationHistory) == 0) ): # redo rotations
			pass
		else:
			for transform in self.processingTransformationHistory:
				if (np.sign(transform) > 0):
					for i in range(transform):
						tempExtent = [tempExtent[2], tempExtent[3], tempExtent[1], tempExtent[0]]
						tempLim = [tempLim[2], tempLim[3], tempLim[1], tempLim[0]]
				elif (np.sign(transform) < 0):
						tempExtent = [tempExtent[1], tempExtent[0], tempExtent[2], tempExtent[3]]
						tempLim = [tempLim[1], tempLim[0], tempLim[2], tempLim[3]]
		
		self.processingExtent = tempExtent
		#self.copiedLimits = tempLim
		self.ui.processing_imagewidget_mpl.canvas.ax.clear()
		self.ui.processing_imagewidget_mpl.canvas.ax.imshow(self.currentProcessingImage, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.processingExtent, aspect = self.processingAspect, origin = "upper")
		#self.ui.imagewidget_mpl.canvas.ax.set_xlim(self.lim[0], self.lim[1]) # Weird things happen when we maintain zoom, probably something to do with inverted axes
		#self.ui.imagewidget_mpl.canvas.ax.set_ylim(self.lim[2], self.lim[3])
		self.ui.processing_imagewidget_mpl.canvas.draw()
		self.ui.processing_imagewidget_mpl.mpl_toolbar._active = None
	
	def subimages_gauss(self, x, a, x0, sigma):
		return a*np.exp(-(x-x0)**2/(2*sigma**2))
	
	def subimages_calculateImages(self, calculateTopGraphOnly = False):
		global transformationHistory
		if (calculateTopGraphOnly and (self.ui.tabwidget.currentIndex() != 2)):
			return
		
		self.calculateImage()
		
		centre = self.cleanNum(self.ui.subimages_centre_spin.value()) # "Centre of split"
		spacing = self.cleanNum(self.ui.subimages_spacing_spin.value()) # "Distance between peaks"
		width = self.cleanNum(self.ui.subimages_width_spin.value()) # "Width of peak"
		ss = self.cleanNum(self.ui.subimages_spectrumnum_spin.value()) # spectrum number
		
		# create Gaussian filters
		xax = np.array(range(len(self.preFFTImage[0])))
		height = 1
		fit_gauss1 = self.subimages_gauss(xax, height, centre - np.floor(spacing/2), width)
		fit_gauss2 = self.subimages_gauss(xax, height, centre + np.floor(spacing/2), width)
		
		# Filtering - for further processing
		spectrum_sub1 = fit_gauss1 * self.preFFTImage
		spectrum_sub2 = fit_gauss2 * self.preFFTImage
		
		# Filtering - to show the user where and how the spectrum is cut
		self.spectrum_sub1_disp = fit_gauss1 * self.rawSpectra
		self.spectrum_sub2_disp = fit_gauss2 * self.rawSpectra
		
		self.ui.subimages_main_mpl.canvas.ax.clear()
		self.ui.subimages_main_mpl.canvas.ax.plot(self.rawSpectra[ss], color='black')
		self.ui.subimages_main_mpl.canvas.ax.plot(self.spectrum_sub1_disp[ss], color='blue')
		self.ui.subimages_main_mpl.canvas.ax.plot(self.spectrum_sub2_disp[ss], color='red')
		self.ui.subimages_main_mpl.canvas.draw()
		
		if (calculateTopGraphOnly):
			return
		
		''' Fourier transformation '''
		mm = 4
		N = int(2048*mm)
		
		fft_data_whole = np.fft.fft(self.preFFTImage, n=N)
		fft_data_1 = np.fft.fft(spectrum_sub1, n=N)
		fft_data_2 = np.fft.fft(spectrum_sub2, n=N)
		
		abs_FFT_whole = np.absolute(fft_data_whole)[:,0:N/2]
		abs_FFT_1 = np.absolute(fft_data_1)[:,0:N/2]
		abs_FFT_2 = np.absolute(fft_data_2)[:,0:N/2]
		
		''' Normalisation '''
		# image is horizontal heading up
		imageOCT_whole = np.rot90(abs_FFT_whole, -1)
		imageOCT_1 = np.rot90(abs_FFT_1, -1)
		imageOCT_2 = np.rot90(abs_FFT_2, -1)
		
		min_val_whole = np.amin(imageOCT_whole)
		max_val_whole = np.amax(imageOCT_whole)
		imageOCT_whole = imageOCT_whole - min_val_whole
		imageOCT_whole = imageOCT_whole / max_val_whole
		
		vmin1 = self.cleanNum(self.ui.subimages_subimage1lowerthreshold_dspin.value())
		vmax1 = self.cleanNum(self.ui.subimages_subimage1upperthreshold_dspin.value())
		
		min_val_1 = np.amin(imageOCT_1) # Normalise
		max_val_1 = np.amax(imageOCT_1)
		imageOCT_1 = imageOCT_1 - min_val_1
		imageOCT_1 = imageOCT_1 / max_val_1
		
		self.subimage1 = imageOCT_1 # Display image
		self.ui.subimages_subimage1_mpl.canvas.ax.clear()
		self.ui.subimages_subimage1_mpl.canvas.ax.imshow(imageOCT_1, self.colourmap, aspect='auto', vmin = vmin1, vmax = vmax1)
		self.ui.subimages_subimage1_mpl.canvas.draw()
		
		
		vmin2 = self.cleanNum(self.ui.subimages_subimage2lowerthreshold_dspin.value())
		vmax2 = self.cleanNum(self.ui.subimages_subimage2upperthreshold_dspin.value())
		
		min_val_2 = np.amin(imageOCT_2) # Normalise
		max_val_2 = np.amax(imageOCT_2)
		imageOCT_2 = imageOCT_2 - min_val_2
		imageOCT_2 = imageOCT_2 / max_val_2
		
		self.subimage2 = imageOCT_2 # Display image
		self.ui.subimages_subimage2_mpl.canvas.ax.clear()
		self.ui.subimages_subimage2_mpl.canvas.ax.imshow(imageOCT_2, self.colourmap, aspect='auto', vmin = vmin2, vmax = vmax2)
		self.ui.subimages_subimage2_mpl.canvas.draw()
		
		colorbar1 = np.linspace(vmin1, vmax1, num = 200)
		colorbar1 = np.rot90(np.transpose(np.tile(colorbar1, (4,1))))
		colorbar2 = np.linspace(vmin2, vmax2, num = 200)
		colorbar2 = np.rot90(np.transpose(np.tile(colorbar2, (4,1))))
		
		self.ui.subimages_subimagelegend1_mpl.canvas.ax.clear() #clear the canvas between draws
		self.ui.subimages_subimagelegend1_mpl.canvas.ax.imshow(colorbar1, self.colourmap, extent = [vmin1, vmax1, vmin1, vmax1], aspect = 0.04)
		self.ui.subimages_subimagelegend1_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.subimages_subimagelegend1_mpl.canvas.draw()
		self.ui.subimages_subimagelegend2_mpl.canvas.ax.clear() #clear the canvas between draws
		self.ui.subimages_subimagelegend2_mpl.canvas.ax.imshow(colorbar2, self.colourmap, extent = [vmin2, vmax2, vmin2, vmax2], aspect = 0.04)
		self.ui.subimages_subimagelegend2_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.subimages_subimagelegend2_mpl.canvas.draw()
	
	def subimages_updateGraph(self, updateBoth = False):
		if (self.subimage1 == [] or self.subimage2 == []):
			return
		if (updateBoth):
			vmin1, vmax1 = self.cleanFloat(self.ui.subimages_subimage1lowerthreshold_dspin.text()), self.cleanFloat(self.ui.subimages_subimage1upperthreshold_dspin.text())
			vmin2, vmax2 = self.cleanFloat(self.ui.subimages_subimage2lowerthreshold_dspin.text()), self.cleanFloat(self.ui.subimages_subimage2upperthreshold_dspin.text())
			
			colorbar1 = np.linspace(vmin1, vmax1, num = 200)
			colorbar1 = np.rot90(np.transpose(np.tile(colorbar1, (4,1))))
			self.ui.subimages_subimagelegend1_mpl.canvas.ax.clear() #clear the canvas between draws
			self.ui.subimages_subimagelegend1_mpl.canvas.ax.imshow(colorbar1, self.colourmap, extent = [vmin1, vmax1, vmin1, vmax1], aspect = 0.04)
			self.ui.subimages_subimagelegend1_mpl.canvas.ax.get_yaxis().set_visible(False)
			self.ui.subimages_subimagelegend1_mpl.canvas.draw()
			self.ui.subimages_subimage1_mpl.canvas.ax.clear()
			self.ui.subimages_subimage1_mpl.canvas.ax.imshow(self.subimage1, self.colourmap, aspect='auto', vmin = vmin1, vmax = vmax1)
			self.ui.subimages_subimage1_mpl.canvas.draw()
			
			colorbar2 = np.linspace(vmin2, vmax2, num = 200)
			colorbar2 = np.rot90(np.transpose(np.tile(colorbar2, (4,1))))
			self.ui.subimages_subimagelegend2_mpl.canvas.ax.clear() #clear the canvas between draws
			self.ui.subimages_subimagelegend2_mpl.canvas.ax.imshow(colorbar2, self.colourmap, extent = [vmin2, vmax2, vmin2, vmax2], aspect = 0.04)
			self.ui.subimages_subimagelegend2_mpl.canvas.ax.get_yaxis().set_visible(False)
			self.ui.subimages_subimagelegend2_mpl.canvas.draw()
			self.ui.subimages_subimage2_mpl.canvas.ax.clear()
			self.ui.subimages_subimage2_mpl.canvas.ax.imshow(self.subimage2, self.colourmap, aspect='auto', vmin = vmin2, vmax = vmax2)
			self.ui.subimages_subimage2_mpl.canvas.draw()
		else:
			if ( (self.sender() == self.ui.subimages_subimage1lowerthreshold_dspin) or (self.sender() == self.ui.subimages_subimage1upperthreshold_dspin) ):
				vmin, vmax = self.cleanFloat(self.ui.subimages_subimage1lowerthreshold_dspin.text()), self.cleanFloat(self.ui.subimages_subimage1upperthreshold_dspin.text())
			elif ( (self.sender() == self.ui.subimages_subimage2lowerthreshold_dspin) or (self.sender() == self.ui.subimages_subimage2upperthreshold_dspin)):
				vmin, vmax = self.cleanFloat(self.ui.subimages_subimage2lowerthreshold_dspin.text()), self.cleanFloat(self.ui.subimages_subimage2upperthreshold_dspin.text())
			
			colorbar = np.linspace(vmin, vmax, num = 200)
			colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
			if ( (self.sender() == self.ui.subimages_subimage1lowerthreshold_dspin) or (self.sender() == self.ui.subimages_subimage1upperthreshold_dspin) ):
				self.ui.subimages_subimagelegend1_mpl.canvas.ax.clear() #clear the canvas between draws
				self.ui.subimages_subimagelegend1_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [vmin, vmax, vmin, vmax], aspect = 0.04)
				self.ui.subimages_subimagelegend1_mpl.canvas.ax.get_yaxis().set_visible(False)
				self.ui.subimages_subimagelegend1_mpl.canvas.draw()
				
				self.ui.subimages_subimage1_mpl.canvas.ax.clear()
				self.ui.subimages_subimage1_mpl.canvas.ax.imshow(self.subimage1, self.colourmap, aspect='auto', vmin = vmin, vmax = vmax)
				self.ui.subimages_subimage1_mpl.canvas.draw()
			
			elif ( (self.sender() == self.ui.subimages_subimage2lowerthreshold_dspin) or (self.sender() == self.ui.subimages_subimage2upperthreshold_dspin) ):
				self.ui.subimages_subimagelegend2_mpl.canvas.ax.clear() #clear the canvas between draws
				self.ui.subimages_subimagelegend2_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [vmin, vmax, vmin, vmax], aspect = 0.04)
				self.ui.subimages_subimagelegend2_mpl.canvas.ax.get_yaxis().set_visible(False)
				self.ui.subimages_subimagelegend2_mpl.canvas.draw()
				
				self.ui.subimages_subimage2_mpl.canvas.ax.clear()
				self.ui.subimages_subimage2_mpl.canvas.ax.imshow(self.subimage2, self.colourmap, aspect='auto', vmin = vmin, vmax = vmax)
				self.ui.subimages_subimage2_mpl.canvas.draw()
	
	def subimages_updateSpectrumNum(self, num):
		if (self.rawSpectra == [] or self.spectrum_sub1_disp == [] or self.spectrum_sub2_disp == []):
			return
		else:
			self.ui.subimages_main_mpl.canvas.ax.clear()
			self.ui.subimages_main_mpl.canvas.ax.plot(self.rawSpectra[num], color='black')
			self.ui.subimages_main_mpl.canvas.ax.plot(self.spectrum_sub1_disp[num], color='blue')
			self.ui.subimages_main_mpl.canvas.ax.plot(self.spectrum_sub2_disp[num], color='red')
			self.ui.subimages_main_mpl.canvas.draw()
	
	# Choose to store all cross sections or just one. Disable other option.
	def threeD_toggleStorage(self, value):
		self.ui.threeD_one_spin.setEnabled(not value)
		self.ui.threeD_one_button.setEnabled(not value)
		self.ui.threeD_one_spinlabel1.setEnabled(not value)
		
		self.ui.threeD_all_spinlabel1.setEnabled(value)
		self.ui.threeD_all_spin.setEnabled(value)
		self.ui.threeD_all_button.setEnabled(value)
	
	# Calculates the image cross sections.
	# IS multithreaded (multiprocessed) to speed up calculation
	def threeD_calculateImage(self):
		self.threeD_3DImage = []
		self.threeD_currentImage = []
		
		## Folder path of data files
		folderPath = self.ui.threeD_datafolder_line.text()
		
		## Create and sort list of files in directory
		onlyfiles = [os.path.join(folderPath, f) for f in os.listdir(folderPath) if os.path.isfile(os.path.join(folderPath, f))]
		onlyfiles = sorted(onlyfiles, key=lambda x: int(os.path.splitext(x)[0].split("_")[1].split("b")[1]))
		
		## Using calibration files from regular tab for now
		try: #to read files
			FileNotFoundErrorBox = QtWidgets.QMessageBox()
			FileNotFoundErrorBox.setText("One or more files were not found. The relevant files may have been moved or might not exist.")
			FileNotFoundErrorBox.setWindowTitle("File Not Found Error")
			PermissionErrorBox = QtWidgets.QMessageBox()
			PermissionErrorBox.setText("Not all file paths were given (or read permissions was denied). Please enter all relevant file paths or deselect the relevant parameter(s).")
			PermissionErrorBox.setWindowTitle("Permission Error")
			
			
			spectra_raw = np.loadtxt(self.ui.spectra_line.text())
			onlyOne = self.ui.threeD_one_radio.isChecked() # Store all or just one cross section?
			
			reference = None
			phase_lin = None
			phase_nonlin = None
			dispersion_compensation = None
			
			if (onlyOne):
				pix = self.cleanNum(self.ui.threeD_one_spin.value())
			else:
				pix = self.cleanNum(self.ui.threeD_all_spin.value())
			
			if (self.ui.backgroundsubtraction_check.isChecked()):
				if (self.ui.referencefromfile_radio.isChecked()):
					reference = np.loadtxt(self.ui.referencespectrum_line.text())
				elif (self.ui.referencesynthesize_radio.isChecked()):
					reference = spectra_raw.mean(axis=0)
			
			if (self.ui.linearisation_check.isChecked()):
				phase_lin = np.loadtxt(self.ui.linearphase_line.text())
				phase_nonlin = np.loadtxt(self.ui.nonlinearphase_line.text())
			
			if (self.ui.dispersioncompensation_check.isChecked()):
				dispersion_compensation = np.loadtxt(self.ui.dispersioncompensation_line.text())
			
		except PermissionError: # This error occurs if the path is blank, or attempts to access a file without proper permissions
			PermissionErrorBox.exec()
			self.ui.std_progressbar.reset()
			self.ui.std_progressbar.hide()
			self.ui.calculate_button.setEnabled(True)
			return
		except FileNotFoundError: # This error occurs when a file cannot be found.
			FileNotFoundErrorBox.exec()
			self.ui.std_progressbar.reset()
			self.ui.std_progressbar.hide()
			self.ui.calculate_button.setEnabled(True)
			return
		
		self.ui.threeD_progressbar.setValue(1)
		numFiles = len(onlyfiles)
		
		i=0
		
		## This is a dict that will be passed into each process. It has everything that should be needed for independent processing
		settings = {'onlyOne':onlyOne, 'pix':pix, 'fftres':self.cleanInt(self.ui.fftres_spin.text())}
		
		backSubSettings = {}
		if (self.ui.backgroundsubtraction_check.isChecked()):
			backSubSettings.update({'enabled':True})
			backSubSettings.update({'file':reference})
		else:
			backSubSettings.update({'enabled':False})
		settings.update({'backgroundSubtraction':backSubSettings})
		
		linSettings = {}
		if (self.ui.linearisation_check.isChecked()):
			linSettings.update({'enabled':True})
			linSettings.update({'cuts':self.getCuts(self.ui.linearphase_line.text())})
			linSettings.update({'files':{'lin':phase_lin, 'nonlin':phase_nonlin}})
		else:
			linSettings.update({'enabled':False})
		settings.update({'linearisation':linSettings})
		
		dispSettings = {}
		if (self.ui.dispersioncompensation_check.isChecked()):
			dispSettings.update({'enabled':True})
			dispSettings.update({'cuts':self.getCuts(self.ui.dispersioncompensation_line.text())})
			dispSettings.update({'file':dispersion_compensation})
		else:
			dispSettings.update({'enabled':False})
		settings.update({'dispersionCompensation':dispSettings})
		
		if multithreadingEnabled:
			
			global processCount
			
			# Don't make less than one process
			poolSize = max(1, processCount)
			
			# Create pool object which delegates jobs to the processes
			pool = multiprocessing.Pool(processes=poolSize)
			
			# Supply parameters that don't change. Every process uses the dict created above
			fixedSettings = partial(threeD_parallelCalculation, inputDict = settings)
			
			# Give the pool a function (and it's dynamic parameters) to assign to processes
			jobs = pool.starmap_async(fixedSettings, enumerate(onlyfiles), 1)
			
			# Signal to pool that no more jobs will be assigned
			pool.close()
			
			# While pool still has jobs remaining
			while True:
				if jobs.ready():
					break
				else:
					# Report to OS and update progress bar
					QtWidgets.QApplication.instance().processEvents()
					self.ui.threeD_progressbar.setValue(1+84*(float(numFiles - jobs._number_left)/numFiles))
			
			# Get data from processes
			image = jobs.get()
			
			# The multiprocessing objects used to eat up memory because of how I first implemented it, but I fixed that. This can stay anyway.
			del pool, jobs, reference, phase_lin, phase_nonlin, dispersion_compensation
		else:
			image = []
			for index, item in enumerate(onlyfiles):
				image.append(threeD_parallelCalculation(index, item, inputDict = settings))
				QtWidgets.QApplication.instance().processEvents()
				self.ui.threeD_progressbar.setValue(1+84*(float(index)/numFiles))
		
		# numpy-ify array and normalise
		image = np.array(image)
		
		image = image - np.amin(image)
		self.ui.threeD_progressbar.setValue(90)
		QtWidgets.QApplication.instance().processEvents()
		
		image = image / np.amax(image)
		self.ui.threeD_progressbar.setValue(95)
		QtWidgets.QApplication.instance().processEvents()
		
		return image
	
	# Called when storing only one face is checked and calulate is pressed
	def threeD_calculateOneImage(self):
		self.ui.threeD_progressbar.show()
		self.threeD_currentImage = self.threeD_calculateImage()
		
		## Sets zoom area to image size
		self.threeD_extent = [0, self.threeD_currentImage.shape[1], 0, self.threeD_currentImage.shape[0]]
		self.ui.threeD_vendpix_spin.setMaximum(self.threeD_currentImage.shape[0])
		self.ui.threeD_hendpix_spin.setMaximum(self.threeD_currentImage.shape[1])
		self.ui.threeD_vendpix_spin.setValue(self.threeD_currentImage.shape[0])
		self.ui.threeD_hendpix_spin.setValue(self.threeD_currentImage.shape[1])
		
		self.lastCalculationWasOnlyOneFace = True
		self.threeD_updateGraph()
		
		## for some reason, the image disappears the first time you zoom, etc. but works fine after.
		## Force a refresh now to get it out of the way
		self.threeD_zoom()
		self.threeD_updateGraph()
		
		self.ui.threeD_progressbar.hide()
	
	# Called when storing all faces is checked and calulate is pressed
	def threeD_calculateAllImages(self):
		self.ui.threeD_progressbar.show()
		pix = self.cleanNum(self.ui.threeD_all_spin.value())
		
		self.threeD_3DImage = self.threeD_calculateImage()
		self.threeD_currentImage = self.threeD_3DImage[:,:,pix]
		
		##Set zoom area to image size
		self.threeD_extent = [0, self.threeD_currentImage.shape[1], 0, self.threeD_currentImage.shape[0]]
		self.ui.threeD_vendpix_spin.setMaximum(self.threeD_currentImage.shape[0])
		self.ui.threeD_hendpix_spin.setMaximum(self.threeD_currentImage.shape[1])
		self.ui.threeD_vendpix_spin.setValue(self.threeD_currentImage.shape[0])
		self.ui.threeD_hendpix_spin.setValue(self.threeD_currentImage.shape[1])
		
		self.lastCalculationWasOnlyOneFace = False
		self.threeD_updateGraph()
		
		## for some reason, the image disappears the first time you zoom, etc. but works fine after.
		## Force a refresh now to get it out of the way
		self.threeD_zoom()
		self.threeD_updateGraph()
		
		if self.ui.threeD_automaticthreshold_check.isChecked():
			self.ui.threeD_upperthreshold_dspin.setValue(np.sqrt(self.threeD_currentImage.mean()))
		
		self.ui.threeD_progressbar.hide()
	
	# Called when the En Face Number spinbox is changed in the 'store all cross ections' setting
	def threeD_changeViewedFace(self):
		if (self.lastCalculationWasOnlyOneFace):
			return
		else:
			pix = self.cleanNum(self.ui.threeD_all_spin.value())
			self.threeD_currentImage = self.threeD_3DImage[:,:,pix]
			self.threeD_updateGraph()
	
	## Most functions below this are essentially copies of functions found in the standard tab, since they share a lot of functionality.
	## Ideally I would have subclassed the imagewidget class and implemented most of these functions there to reduce duplicate code, 
	## but this tab was added much later in development and I hadn't considered what might happen if I needed to add a new tab with similiar functionality.
	
	# Make sure that if the user pans/zooms the image with the toolbar, these changes are applied to the start/end pixel spinboxes too
	def threeD_synchroniseZoom(self, event): 
		if (self.ui.threeD_imagewidget_mpl.mpl_toolbar._active != None): #If the current toolbar mode is either PAN or ZOOM, then update the coordinates of the spinboxes with the coordinates of the on-screen image
			vstart = min(self.ui.threeD_imagewidget_mpl.canvas.ax.get_ylim())
			vend = max(self.ui.threeD_imagewidget_mpl.canvas.ax.get_ylim())
			hstart = min(self.ui.threeD_imagewidget_mpl.canvas.ax.get_xlim())
			hend = max(self.ui.threeD_imagewidget_mpl.canvas.ax.get_xlim())
			
			self.threeD_disconnectZoomUpdates()
			self.ui.threeD_vstartpix_spin.setMaximum(vend - 1)
			self.ui.threeD_vstartpix_spin.setValue(vstart)
			
			self.ui.threeD_vendpix_spin.setMinimum(vstart + 1)
			self.ui.threeD_vendpix_spin.setValue(vend)
			
			self.ui.threeD_hstartpix_spin.setMaximum(hend - 1)
			self.ui.threeD_hstartpix_spin.setValue(hstart)
			
			self.ui.threeD_hendpix_spin.setMinimum(hstart + 1)
			self.ui.threeD_hendpix_spin.setValue(hend)
			
			self.threeD_reconnectZoomUpdates()
	
	def threeD_disconnectZoomUpdates(self): #Disconnect the spinboxes from their slots, in case we want to update their values without promping a redraw
		self.ui.threeD_vstartpix_spin.disconnect()
		self.ui.threeD_vendpix_spin.disconnect()
		self.ui.threeD_hstartpix_spin.disconnect()
		self.ui.threeD_hendpix_spin.disconnect()
	
	def threeD_reconnectZoomUpdates(self): #Reconnect the spinboxes once we've done what we want
		self.ui.threeD_vstartpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_vendpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_hstartpix_spin.valueChanged.connect(self.threeD_zoom)
		self.ui.threeD_hendpix_spin.valueChanged.connect(self.threeD_zoom)
		
		self.ui.threeD_hendpix_spin.valueChanged.connect(self.threeD_updateStartZoomMax)
		self.ui.threeD_vendpix_spin.valueChanged.connect(self.threeD_updateStartZoomMax)
		self.ui.threeD_hstartpix_spin.valueChanged.connect(self.threeD_updateEndZoomMin)
		self.ui.threeD_vstartpix_spin.valueChanged.connect(self.threeD_updateEndZoomMin)
	
	def threeD_zoom(self):
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		vstart, vend = self.cleanInt(self.ui.threeD_vstartpix_spin.text()), self.cleanInt(self.ui.threeD_vendpix_spin.text()) #get limits from spinboxes
		hstart, hend = self.cleanInt(self.ui.threeD_hstartpix_spin.text()), self.cleanInt(self.ui.threeD_hendpix_spin.text())
		
		if (self.ui.threeD_imagewidget_mpl.canvas.ax.xaxis_inverted()): # Zoom image accordingly
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(hend, hstart)
		else:
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(hstart, hend) 
		
		if (self.ui.threeD_imagewidget_mpl.canvas.ax.yaxis_inverted()):
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(vend, vstart)
		else:
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(vstart, vend)
		
		self.ui.threeD_imagewidget_mpl.canvas.draw() #Redraw so changes take effect
	
	def threeD_updateStartZoomMax(self, val):
		if (self.sender() == self.ui.threeD_one_button or self.sender() == self.ui.threeD_all_button):
			self.ui.threeD_hstartpix_spin.setMaximum(self.ui.threeD_hendpix_spin.value() -1)
			self.ui.threeD_vstartpix_spin.setMaximum(self.ui.threeD_vendpix_spin.value() -1)
		if (self.sender() == self.ui.threeD_hendpix_spin):
			self.ui.threeD_hstartpix_spin.setMaximum(max(val-1,0))
		elif ((self.sender() == self.ui.threeD_vendpix_spin) and (self.cleanNum(self.ui.threeD_vendpix_spin.value()) > 0)):
			self.ui.threeD_vstartpix_spin.setMaximum(max(val-1,0))
	
	def threeD_updateEndZoomMin(self, val):
		try:
			if (imageOCT == []):
				return
		except NameError:
			return
		
		if ((self.sender() == self.ui.threeD_vstartpix_spin) and (self.cleanNum(self.ui.threeD_vstartpix_spin.value()) < imageOCT.shape[0])):
			self.ui.threeD_vendpix_spin.setMinimum(val+1)
		elif ((self.sender() == self.ui.threeD_hstartpix_spin) and (self.cleanNum(self.ui.threeD_hstartpix_spin.value()) < imageOCT.shape[1])):
			self.ui.threeD_hendpix_spin.setMinimum(val+1)
	
	def threeD_resetZoom(self):
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		self.threeD_disconnectZoomUpdates()
		
		self.ui.threeD_hendpix_spin.setMinimum(1)
		self.ui.threeD_hendpix_spin.setMaximum(self.threeD_currentImage.shape[1])
		self.ui.threeD_hendpix_spin.setValue(self.threeD_currentImage.shape[1])
		
		self.ui.threeD_vendpix_spin.setMinimum(1)
		self.ui.threeD_vendpix_spin.setMaximum(self.threeD_currentImage.shape[0])
		self.ui.threeD_vendpix_spin.setValue(self.threeD_currentImage.shape[0])
		
		self.ui.threeD_hstartpix_spin.setMaximum(self.ui.hendpix_spin.value()-1)
		self.ui.threeD_hstartpix_spin.setMinimum(0)
		self.ui.threeD_hstartpix_spin.setValue(0)
		
		self.ui.threeD_vstartpix_spin.setMaximum(self.ui.vendpix_spin.value()-1)
		self.ui.threeD_vstartpix_spin.setMinimum(0)
		self.ui.threeD_vstartpix_spin.setValue(0)
		
		self.threeD_reconnectZoomUpdates()
		
		if (self.ui.threeD_imagewidget_mpl.canvas.ax.xaxis_inverted() ): #Check for inverted axes and set the limits accordingly
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(self.threeD_currentImage.shape[1], 0)
		else:
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(0, self.threeD_currentImage.shape[1])
		
		if (self.ui.threeD_imagewidget_mpl.canvas.ax.yaxis_inverted() ):
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(self.threeD_currentImage.shape[0], 0)
		else:
			self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(0, self.threeD_currentImage.shape[0])
		
		self.ui.threeD_imagewidget_mpl.canvas.draw()
	
	def threeD_redrawIntensity(self, val):
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		### Synchronise slider and spinbox
		if (self.sender() == self.ui.threeD_upperthreshold_dspin):
			self.ui.threeD_lowerthreshold_dspin.setMaximum(min(val, 0.99))
			self.ui.threeD_lowerthreshold_slider.setMinimum(min(10**(2-2*val), 99))
			
			self.ui.threeD_upperthreshold_slider.disconnect()
			try:
				self.ui.threeD_upperthreshold_slider.setValue(10**(2-2*self.ui.threeD_upperthreshold_dspin.value()))
			except ValueError:
				pass
			self.ui.threeD_upperthreshold_slider.valueChanged.connect(self.threeD_redrawIntensity)
			
		elif (self.sender() == self.ui.threeD_upperthreshold_slider):
			self.ui.threeD_lowerthreshold_dspin.setMaximum(min(1-np.log10(val)/2, 0.99))
			self.ui.threeD_lowerthreshold_slider.setMinimum(min(val, 99))
			
			self.ui.threeD_upperthreshold_dspin.disconnect()
			try:
				self.ui.threeD_upperthreshold_dspin.setValue(1-np.log10(self.ui.threeD_upperthreshold_slider.value())/2)
			except ValueError:
				pass
			self.ui.threeD_upperthreshold_dspin.valueChanged.connect(self.threeD_redrawIntensity)
			
		elif (self.sender() == self.ui.threeD_lowerthreshold_dspin):
			self.ui.threeD_upperthreshold_dspin.setMinimum(max(val, 0.01))
			self.ui.threeD_upperthreshold_slider.setMaximum(max(10**(2-2*val), 1))
			
			self.ui.threeD_lowerthreshold_slider.disconnect()
			try:
				self.ui.threeD_lowerthreshold_slider.setValue(10**(2-2*self.ui.threeD_lowerthreshold_dspin.value()))
			except ValueError:
				pass
			self.ui.threeD_lowerthreshold_slider.valueChanged.connect(self.threeD_redrawIntensity)
			
		elif (self.sender() == self.ui.threeD_lowerthreshold_slider):
			self.ui.threeD_upperthreshold_dspin.setMinimum(max(1-np.log10(val)/2, 0.01))
			self.ui.threeD_upperthreshold_slider.setMaximum(max(val, 1))
			
			self.ui.threeD_lowerthreshold_dspin.disconnect()
			try:
				self.ui.threeD_lowerthreshold_dspin.setValue(1-np.log10(self.ui.threeD_lowerthreshold_slider.value())/2)
			except ValueError:
				pass
			self.ui.threeD_lowerthreshold_dspin.valueChanged.connect(self.threeD_redrawIntensity)
		self.threeD_updateGraph()
	
	def threeD_resetIntensity(self):
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		self.ui.upperthreshold_dspin.setValue(1)
		self.ui.lowerthreshold_dspin.setValue(0)
	
	def threeD_updateGraph(self):
		
		try:
			if (self.threeD_currentImage == []):
				return
		except (NameError, AttributeError):
			return
		
		lowerThreshold, upperThreshold = self.cleanFloat(self.ui.threeD_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.threeD_upperthreshold_dspin.text())
		limits = self.ui.threeD_imagewidget_mpl.canvas.ax.axis()
		
		self.ui.threeD_imagewidget_mpl.canvas.ax.clear() # Clears the canvas between draw commands.
		
		#vmin, vmax set the intensity thresholds
		#extent dictates the image size
		#aspect controls the canvas scaling, square, 1:1, etc.
		#origin decides whether (0,0) is bottom left, like graphs, or upper left, like arrays
		self.ui.threeD_imagewidget_mpl.canvas.ax.imshow(self.threeD_currentImage, self.colourmap, vmin = lowerThreshold, vmax = upperThreshold, extent = self.threeD_extent, aspect = "auto", origin = "upper")
		
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(limits[0], limits[1]) # Reapply limits to viewing area
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(limits[2], limits[3])
		self.ui.threeD_imagewidget_mpl.canvas.draw()
		
		##### Colour bar #####
		colorbar = np.linspace(lowerThreshold, upperThreshold, num = 200)
		colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
		
		self.ui.threeD_imagelegend_mpl.canvas.ax.clear()
		self.ui.threeD_imagelegend_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [lowerThreshold, upperThreshold, lowerThreshold, upperThreshold], aspect = 0.02)
		self.ui.threeD_imagelegend_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.threeD_imagelegend_mpl.canvas.draw()
		
		self.threeD_updateLegend()
	
	def threeD_updateLegend(self):
		vmin, vmax = self.cleanFloat(self.ui.threeD_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.threeD_upperthreshold_dspin.text())
		
		colorbar = np.linspace(vmin, vmax, num = 200)
		colorbar = np.rot90(np.transpose(np.tile(colorbar, (4,1))))
		
		self.ui.threeD_imagelegend_mpl.canvas.ax.cla() #clear the canvas between draws
		self.ui.threeD_imagelegend_mpl.canvas.ax.imshow(colorbar, self.colourmap, extent = [vmin, vmax, vmin, vmax], aspect = 0.02)
		self.ui.threeD_imagelegend_mpl.canvas.ax.get_yaxis().set_visible(False)
		self.ui.threeD_imagelegend_mpl.canvas.draw()
		
		if (self.ui.threeD_one_radio.isChecked()):
			self.ui.threeD_one_button.setEnabled(True)
		else:
			self.ui.threeD_all_button.setEnabled(True)
	
	def threeD_rotateImage(self):#, catchup = False):
		
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		limits = self.ui.threeD_imagewidget_mpl.canvas.ax.axis() # Store currently viewed area
		
		self.threeD_resetZoom()
		self.ui.threeD_imagewidget_mpl.canvas.ax.clear()
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.threeD_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.threeD_upperthreshold_dspin.text())
		
		if self.ui.threeD_all_radio.isChecked():
			self.threeD_3DImage = np.rot90(self.threeD_3DImage, -1) # Flip image left/right
			pix = self.cleanNum(self.ui.threeD_all_spin.value())
			self.threeD_currentImage = self.threeD_3DImage[:,:,pix]
		else:
			self.threeD_currentImage = np.rot90(self.threeD_currentImage, -1)
		
		self.threeD_extent = [self.threeD_extent[2], self.threeD_extent[3], self.threeD_extent[1], self.threeD_extent[0]] #Rotate axes with image
		
		self.ui.threeD_imagewidget_mpl.canvas.ax.imshow(self.threeD_currentImage, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.threeD_extent, aspect = self.threeD_aspect, origin = "upper")
		
		#if (not catchup):
		#	self.updateTransformationHistory('rotate')
		
		
		self.threeD_resetZoom() # Not sure why, but we need to reset the zoom again so it displays properly
		
		########## Re-zoom to previous spot ##########
		limits = [limits[2], limits[3], limits[1], limits[0]] #Rotate viewing area with image
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(limits[0], limits[1]) # Reapply limits to viewing area
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(limits[2], limits[3])
		self.ui.threeD_imagewidget_mpl.canvas.draw()
	
	def threeD_mirrorImage(self):#, catchup = False):
		try:
			if (self.threeD_currentImage == []):
				return
		except NameError:
			return
		
		limits = self.ui.threeD_imagewidget_mpl.canvas.ax.axis()
		
		self.threeD_resetZoom()
		self.ui.threeD_imagewidget_mpl.canvas.ax.clear()
		
		thresholdMin, thresholdMax = self.cleanFloat(self.ui.threeD_lowerthreshold_dspin.text()), self.cleanFloat(self.ui.threeD_upperthreshold_dspin.text())
		if self.ui.threeD_all_radio.isChecked():
			self.threeD_3DImage = np.fliplr(self.threeD_3DImage) # Flip image left/right
			pix = self.cleanNum(self.ui.threeD_all_spin.value())
			self.threeD_currentImage = self.threeD_3DImage[:,:,pix]
		else:
			self.threeD_currentImage = np.fliplr(self.threeD_currentImage) # Flip image left/right
		
		self.threeD_extent = [self.threeD_extent[1], self.threeD_extent[0], self.threeD_extent[2], self.threeD_extent[3]] #Rotate axes with image
		
		self.ui.threeD_imagewidget_mpl.canvas.ax.imshow(self.threeD_currentImage, self.colourmap, vmin = thresholdMin, vmax = thresholdMax, extent = self.threeD_extent, aspect = self.threeD_aspect, origin = "upper")
		
		#if (not catchup):
		#	self.updateTransformationHistory('mirror')
		
		self.threeD_resetZoom()
		
		########## Re-zoom to previous spot ##########
		limits = [limits[1], limits[0], limits[2], limits[3]] #Mirror viewing area with image
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_xlim(limits[0], limits[1])
		self.ui.threeD_imagewidget_mpl.canvas.ax.set_ylim(limits[2], limits[3])
		self.ui.threeD_imagewidget_mpl.canvas.draw()

if __name__=='__main__':
	Program =  QtWidgets.QApplication(sys.argv)
	MyProg = Prog()
	MyProg.show()
	sys.exit(Program.exec_())

