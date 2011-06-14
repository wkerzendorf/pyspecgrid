from PyQt4 import QtGui, QtCore
import pylab
import copy
import minuit
import numpy as np

class converter(object):
    def __init__(self, limit, sliderResolution=1000):
        self.sliderResolution = sliderResolution
        self.limit = limit
        self.delta = limit[1] - limit[0]
        self.curVal = limit[0]
        
            
            
    def convertFromSlider(self, value):
        return self.limit[0] + \
           value*self.delta/float(self.sliderResolution)
    def convertToSlider(self, value):
        return int((value-self.limit[0])*self.sliderResolution/self.delta)
        
    def valueUpdate(val):
        spec = specGrid.interpGrid(tuple([slider.val for slider in sliders]))
        specPlot.set_ydata(spec)
        if autoScale:
            plotAxis.relim()
            plotAxis.autoscale_view()
        pylab.draw()
    
    def textboxChanged(self, value):
        try:
            value = float(value)
        except ValueError:
            self.textbox.setText(str(self.curValue))
            return textboxChanged(self.curValue)
        self.slider.setValue(self.convertToSlider(value))
        return
        
    def sliderChanged(self, value):
        value = self.convertFromSlider(value)
        self.curValue = value
        self.textbox.setText(str(value))
        return

class FitterGUI(QtGui.QWidget):
    def __init__(self, specGrid):
        QtGui.QWidget.__init__(self)
        
        self.specGrid = specGrid
        self.figure = pylab.figure(1)
        self.figure.clf()
        self.sampleSpec = None

        self.fitChiSq = None
        self.ax = self.figure.add_subplot(111)
        spec = self.specGrid.getSpec(*(item[0] for item in self.specGrid.limits))
        self.specPlot, = self.ax.plot(spec.wave, spec.flux)
        self.sampleSpecPlot, = self.ax.plot([])
        layout = QtGui.QGridLayout()
        self.sliders = []
        self.textBoxes = []
        self.converters = []
        self.autoScaleView = True
        self.curValue = [item[0] for item in self.specGrid.limits]
        for i, (label, limit) in enumerate(zip(self.specGrid.params, self.specGrid.limits)):    
            slider = QtGui.QSlider(QtCore.Qt.Horizontal)
            slider.setMinimum(0)
            slider.setMaximum(1000)
            self.sliders.append(slider)
            
            tbox = QtGui.QLineEdit(str(limit[0]))
            self.textBoxes.append(tbox)
            
            self.converters.append(converter(limit=limit))
            self.connect(slider,
            QtCore.SIGNAL("sliderMoved(int)"),
            self._sliderChanged)
            
            self.connect(tbox,
            QtCore.SIGNAL("editingFinished()"),
            self._textBoxChanged)
            
            qlabel = QtGui.QLabel(label)
        
            layout.addWidget(qlabel, i, 1, 1, 1)
            layout.addWidget(slider, i, 2, 1, 10)
            layout.addWidget(tbox, i, 12, 1, 1)
            
        self.buttonChiSq = QtGui.QPushButton('Chi Square')
        self.buttonChiSq.setDisabled(True)
        self.connect(self.buttonChiSq,
            QtCore.SIGNAL("clicked()"),
            self._fitChiSq)
        
        
        self.buttonSimplex = QtGui.QPushButton('Simplex')
        self.buttonSimplex.setDisabled(True)
        self.connect(self.buttonSimplex,
            QtCore.SIGNAL("clicked()"),
            self._fitSimplex)
        
        self.buttonMigrad = QtGui.QPushButton('Migrad')
        self.buttonMigrad.setDisabled(True)
        self.connect(self.buttonMigrad,
            QtCore.SIGNAL("clicked()"),
            self._fitMigrad)
        self.fvalLabel = QtGui.QLabel('fval = None')
        layout.addWidget(self.buttonChiSq, i+1, 1)
        layout.addWidget(self.buttonSimplex, i+1, 2)
        layout.addWidget(self.buttonMigrad, i+1, 3)
        layout.addWidget(self.fvalLabel, i+1, 4)
        self.setLayout(layout)
        self.show()
        self.figure.show()
    def _sliderChanged(self, value):

        for i, (slider, textBox, converter) in enumerate(zip(self.sliders, self.textBoxes, self.converters)):
            self.curValue[i] = converter.convertFromSlider(slider.value())
            textBox.setText(str(self.curValue[i]))
        if self.sampleSpec != None:
            self.fvalLabel.setText('fval=%s' % self.fminuit.fcn(*self.curValue))
        self._updateSpecPlot()
        
    def _updateSpecPlot(self):
        spec = self.specGrid.interpGrid(tuple(self.curValue))
        self.specPlot.set_ydata(spec)
        if self.autoScaleView:
            self.ax.relim()
            self.ax.autoscale_view()
        pylab.draw()
    
    def _textBoxChanged(self):
        for i, (slider, textBox, converter) in enumerate(zip(self.sliders, self.textBoxes, self.converters)):
            try:
                value = float(textBox.text())
            except ValueError:
                value = self.curValue[i]
                textBox.setText(str(value))
            if value < self.specGrid.limits[i][0] or value > self.specGrid.limits[i][1]:
                value = self.curValue[i]
                textBox.setText(str(value))
            
            self.curValue[i] = value
            slider.setValue(converter.convertToSlider(value))
        if self.sampleSpec != None:
            self.fvalLabel.setText('fval=%s' % self.fminuit.fcn(*self.curValue))
        self._updateSpecPlot()
    
    def _updateFromCurValue(self):
        for slider, textbox, converter, value in zip(self.sliders, self.textBoxes, self.converters, self.curValue):
            slider.setValue(converter.convertToSlider(value))
            textbox.setText(str(value))
        if self.sampleSpec != None:
            self.fvalLabel.setText('fval=%s' % self.fminuit.fcn(*self.curValue))
        self._updateSpecPlot()
        
    def _fitChiSq(self):
        if self.fitChiSq ==None:
            fitParams = self.specGrid.fitChiSq(self.sampleSpec)
            self.fitChiSq = fitParams
        
        self.curValue = copy.copy(self.fitChiSq)
        self._updateFromCurValue()
        
    def _fitSimplex(self):
        for param, value in zip(self.specGrid.params, self.curValue):
            self.fminuit.values[param] = value
        try:
            self.fminuit.simplex()
        except minuit.MinuitError as me:
            pass
        
        for i, (param, limit) in enumerate(zip(self.specGrid.params, self.specGrid.limits)):
            value = self.fminuit.values[param]
            self.curValue[i] = np.max((np.min((value, limit[1])), limit[0]))
        
        self._updateFromCurValue()
        
        
    def _fitMigrad(self):
        for param, value in zip(self.specGrid.params, self.curValue):
            self.fminuit.values[param] = value
        try:
            self.fminuit.migrad()
        except minuit.MinuitError as me:
            pass
        
        for i, (param, limit) in enumerate(zip(self.specGrid.params, self.specGrid.limits)):
            value = self.fminuit.values[param]
            self.curValue[i] = np.max((np.min((value, limit[1])), limit[0]))
        
        self._updateFromCurValue()
        
    def setSampleSpec(self, sampleSpec):
        self.sampleSpec = sampleSpec
        self.sampleSpecPlot.set_data(self.sampleSpec.wave, self.sampleSpec.flux)
        self.fitChiSq = None
        self.fminuit = self.specGrid.getMinuit(self.sampleSpec)
        self.buttonChiSq.setEnabled(True)
        self.buttonMigrad.setEnabled(True)
        self.buttonSimplex.setEnabled(True)
        pylab.draw()
        
        
        