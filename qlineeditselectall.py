from PyQt5 import QtWidgets

class QLineEditSelectAll(QtWidgets.QLineEdit):
	def __init__(self, parent=None):
		QtWidgets.QLineEdit.__init__(self, parent)   # Inherit from QLinEdit
	def mousePressEvent(self, event):
		self.selectAll()