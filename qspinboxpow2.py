from PyQt5 import QtWidgets

class QSpinBoxPow2(QtWidgets.QSpinBox):
	def __init__(self, parent=None):
		QtWidgets.QSpinBox.__init__(self, parent)   # Inherit from QSpinBox
	
	def stepBy(self, steps):
		i = 1
		pow2 = 0
		while ( i <= self.value() ):
			if (i == self.value()):
				pow2=1
			i *= 2
		
		if ((steps > 0) and (i < 999999999)):
			self.setValue(i)
			return
		elif (steps < 0):
			if (i==1):
				self.setValue(1)
				return
			if( pow2==1 ):
				self.setValue(i/4)
				return
			self.setValue (i/2)