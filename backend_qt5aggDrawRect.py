import six

import ctypes
import sys

from PyQt5 import QtCore, QtGui

from PyQt5.QtWidgets import QSizePolicy, QWidget, QVBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

QT_API = 'PyQt5'
DEBUG = False

_decref = ctypes.pythonapi.Py_DecRef
_decref.argtypes = [ctypes.py_object]
_decref.restype = None

class FigureCanvasQTAggDrawRect(FigureCanvas):
    def __init__(self):
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


    #Change de color of rectangle zoom toolbar rewriting painEvent
    #the original code is in the backend_qt5agg.py file inside
    #matplotlib/backends directory
    def paintEvent(self, e):
        """
        Copy the image from the Agg canvas to the qt.drawable.
        In Qt, all drawing should be done inside of here when a widget is
        shown onscreen.
        """
        # if the canvas does not have a renderer, then give up and wait for
        # FigureCanvasAgg.draw(self) to be called
        if not hasattr(self, 'renderer'):
            return

        if DEBUG:
            print('FigureCanvasQtAgg.paintEvent: ', self,
                  self.get_width_height())

        if len(self.blitbox) == 0:
            # matplotlib is in rgba byte order.  QImage wants to put the bytes
            # into argb format and is in a 4 byte unsigned int.  Little endian
            # system is LSB first and expects the bytes in reverse order
            # (bgra).
            if QtCore.QSysInfo.ByteOrder == QtCore.QSysInfo.LittleEndian:
                stringBuffer = self.renderer._renderer.tostring_bgra()
            else:
                stringBuffer = self.renderer._renderer.tostring_argb()

            refcnt = sys.getrefcount(stringBuffer)

            # convert the Agg rendered image -> qImage
            qImage = QtGui.QImage(stringBuffer, self.renderer.width,
                                  self.renderer.height,
                                  QtGui.QImage.Format_ARGB32)
            if hasattr(qImage, 'setDevicePixelRatio'):
                # Not available on Qt4 or some older Qt5.
                qImage.setDevicePixelRatio(self._dpi_ratio)
            # get the rectangle for the image
            rect = qImage.rect()
            p = QtGui.QPainter(self)
            # reset the image area of the canvas to be the back-ground color
            p.eraseRect(rect)
            # draw the rendered image on to the canvas
            p.drawPixmap(QtCore.QPoint(0, 0), QtGui.QPixmap.fromImage(qImage))

            # draw the zoom rectangle to the QPainter
            ########################################################
            ## Edit: Draw a solid black rectangle, then a white dashed line
            ########################################################
            if self._drawRect is not None:
                pen = QtGui.QPen(QtCore.Qt.black, 1 / self._dpi_ratio,
                                 QtCore.Qt.SolidLine)
                p.setPen(pen)
                x, y, w, h = self._drawRect
                p.drawRect(x, y, w, h)
                
                pen.setColor(QtCore.Qt.white)
                pen.setStyle(QtCore.Qt.DashLine)
                p.setPen(pen)
                x, y, w, h = self._drawRect
                p.drawRect(x, y, w, h)
            p.end()

            # This works around a bug in PySide 1.1.2 on Python 3.x,
            # where the reference count of stringBuffer is incremented
            # but never decremented by QImage.
            # TODO: revert PR #1323 once the issue is fixed in PySide.
            del qImage
            if refcnt != sys.getrefcount(stringBuffer):
                _decref(stringBuffer)
        else:
            p = QtGui.QPainter(self)

            while len(self.blitbox):
                bbox = self.blitbox.pop()
                l, b, r, t = bbox.extents
                w = int(r) - int(l)
                h = int(t) - int(b)
                t = int(b) + h
                reg = self.copy_from_bbox(bbox)
                stringBuffer = reg.to_string_argb()
                qImage = QtGui.QImage(stringBuffer, w, h,
                                      QtGui.QImage.Format_ARGB32)
                if hasattr(qImage, 'setDevicePixelRatio'):
                    # Not available on Qt4 or some older Qt5.
                    qImage.setDevicePixelRatio(self._dpi_ratio)
                # Adjust the stringBuffer reference count to work
                # around a memory leak bug in QImage() under PySide on
                # Python 3.x
                if QT_API == 'PySide' and six.PY3:
                    ctypes.c_long.from_address(id(stringBuffer)).value = 1

                origin = QtCore.QPoint(l, self.renderer.height - t)
                pixmap = QtGui.QPixmap.fromImage(qImage)
                p.drawPixmap(origin / self._dpi_ratio, pixmap)

            # draw the zoom rectangle to the QPainter
            if self._drawRect is not None:
                pen = QtGui.QPen(QtCore.Qt.black, 1 / self._dpi_ratio,
                                 QtCore.Qt.DotLine)
                p.setPen(pen)
                x, y, w, h = self._drawRect
                p.drawRect(x, y, w, h)

            p.end()