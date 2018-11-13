# Imports
from PyQt5 import QtWidgets
import matplotlib
# Ensure using PyQt5 backend
matplotlib.use('qt5agg')

from matplotlib.figure import Figure
from backend_qt5aggDrawRect import FigureCanvasQTAggDrawRect as Canvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


# Matplotlib canvas class to create figure
class MplCanvas(Canvas):
	def __init__(self, tight = True, tiny = False):
		
		if (tight):
			self.fig = Figure(tight_layout = True)
		else:
			self.fig = Figure()
		
		if (tiny):
			self.ax = self.fig.add_axes([0.12,0.5,0.80,0.38])
		else:
			self.ax = self.fig.add_axes([0.12,0.1,0.84,0.88])
		
		Canvas.__init__(self)
		Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
		Canvas.updateGeometry(self)

class NavigationToolbarBasic(NavigationToolbar):
	# only display the buttons we need
	toolitems = [t for t in NavigationToolbar.toolitems if
				 t[0] in ('Home', 'Pan', 'Zoom')]
	def home(self, *args):#Reset zoom
		
		#self._views.home()
		#self._positions.home()
		#self.set_history_buttons()
		#self._update_view()	
		self.window().resetZoom()

# Matplotlib widget
class MplWidgetNavbar(QtWidgets.QWidget):
	def __init__(self, parent=None):
		QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
		self.canvas = MplCanvas()				  # Create canvas object
		self.mpl_toolbar = NavigationToolbarBasic(self.canvas, self) # Create toolbar object
		self.vbl = QtWidgets.QVBoxLayout()		 # Set box for plotting
		self.vbl.addWidget(self.mpl_toolbar)
		self.vbl.addWidget(self.canvas)
		
		self.setLayout(self.vbl)

class MplWidget(QtWidgets.QWidget): # Image widget without toolbar
	def __init__(self, parent=None):
		QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
		self.canvas = MplCanvas(tight = False)				  # Create canvas object
		self.vbl = QtWidgets.QVBoxLayout()		 # Set box for plotting
		self.vbl.addWidget(self.canvas)
		
		self.setLayout(self.vbl)

class MplWidgetLegend(QtWidgets.QWidget): # Image widget without toolbar
	def __init__(self, parent=None):
		QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
		self.canvas = MplCanvas(tight = True, tiny = True)				  # Create canvas object
		self.vbl = QtWidgets.QVBoxLayout()		 # Set box for plotting
		self.vbl.addWidget(self.canvas)
		
		self.setLayout(self.vbl)