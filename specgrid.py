import pyspec
import os
from scipy import interpolate, ndimage
import numpy as np
import sqlite3
from glob import glob
import pysynphot
import pydib
import ConfigParser

from matplotlib.widgets import Slider
specGridDir = '/Users/wkerzend/projects/grids'

def getSimpleSpecs(gridName, fnames, config, **kwargs):

    #reading config
    waveFName = config.get('wave', 'wave')
    specDType = config.get('structure', 'datatype')
    specSize = config.getint('structure', 'specsize')
    
    #Reading wave solution
    wave = np.fromfile(os.path.join(specGridDir, gridName, waveFName))
    
    gridSize = len(fnames) * specSize / 1024**2
    print "Processing %.3f MB in to Memory." % gridSize
    
    specs = []
    if kwargs.has_key('normRange'):
        if kwargs.has_key('normMode'):
            normMode = kwargs['normMode']
        else:
            normMode = 'simple'
        
        if normMode == 'simple':
            normSlice = slice(*[wave.searchsorted(item) for item in kwargs['normRange']])
            
    for specFName in fnames:
        spec = np.fromfile(os.path.join(specGridDir, gridName, 'data', specFName))
        if kwargs.has_key('normRange'):
            normFac = np.mean(spec[normSlice])
            spec /= normFac
        if kwargs.has_key('newWave'):
            newWave = kwargs['newWave']
            f = interpolate.interp1d(wave, spec)
            spec = f(kwargs['newWave'])
        if kwargs.has_key('smoothRes'):
            if kwargs.has_key('newWave'):
                deltaLambda = ((newWave.max()+newWave.min())/2.)/kwargs['smoothRes']
            else:
                deltaLambda = ((wave.max()+wave.min())/2.)/kwargs['smoothRes']
            sigma = deltaLambda / 2.3548
            spec = ndimage.gaussian_filter1d(spec, sigma)
        specs.append(spec)
    if kwargs.has_key('newWave'):
        return newWave, np.array(specs)
    else:
        return wave, np.array(specs)
        
def getMunariSpecs(gridName, fnames, config, **kwargs):

    #reading config
    waveFName = config.get('wave', 'wave')
    specDType = config.get('structure', 'datatype')
    specSize = config.getint('structure', 'specsize')
    
    #Reading wave solution
    wave = np.fromfile(os.path.join(specGridDir, gridName, waveFName))
    
    gridSize = len(fnames) * specSize / 1024**2
    print "Processing %.3f MB in to Memory." % gridSize
    
    specs = []
    if kwargs.has_key('normRange'):
        if kwargs.has_key('normMode'):
            normMode = kwargs['normMode']
        else:
            normMode = 'simple'
        
        if normMode == 'simple':
            normSlice = slice(*[wave.searchsorted(item) for item in kwargs['normRange']])
            
    for specFName in fnames:
        spec = np.fromfile(os.path.join(specGridDir, gridName, 'data', specFName))
        if kwargs.has_key('normRange'):
            normFac = np.mean(spec[normSlice])
            spec /= normFac
        if kwargs.has_key('newWave'):
            newWave = kwargs['newWave']
            f = interpolate.interp1d(wave, spec)
            spec = f(kwargs['newWave'])
        if kwargs.has_key('smoothRes'):
            if kwargs.has_key('newWave'):
                deltaLambda = ((newWave.max()+newWave.min())/2.)/kwargs['smoothRes']
            else:
                deltaLambda = ((wave.max()+wave.min())/2.)/kwargs['smoothRes']
            sigma = deltaLambda / 2.3548
            spec = ndimage.gaussian_filter1d(spec, sigma)
        specs.append(spec)
    if kwargs.has_key('newWave'):
        return newWave, np.array(specs)
    else:
        return wave, np.array(specs)

def getMarcsSpecs(gridName, fnames, config, **kwargs):

    #reading config
    waveFName = config.get('wave', 'wave')
    specDType = config.get('structure', 'datatype')
    specSize = config.getint('structure', 'specsize')
    
    #Reading wave solution
    wave = np.fromfile(os.path.join(specGridDir, gridName, waveFName))
    
    gridSize = len(fnames) * specSize / 1024**2
    print "Processing %.3f MB in to Memory." % gridSize
    
    specs = []
    if kwargs.has_key('normRange'):
        if kwargs.has_key('normMode'):
            normMode = kwargs['normMode']
        else:
            normMode = 'simple'
        
        if normMode == 'simple':
            normSlice = slice(*[wave.searchsorted(item) for item in kwargs['normRange']])
            
    for specFName in fnames:
        spec = np.fromfile(os.path.join(specGridDir, gridName, 'data', specFName))
        if kwargs.has_key('normRange'):
            normFac = np.mean(spec[normSlice])
            spec /= normFac
        if kwargs.has_key('newWave'):
            newWave = kwargs['newWave']
            f = interpolate.interp1d(wave, spec)
            spec = f(kwargs['newWave'])
        if kwargs.has_key('smoothRes'):
            if kwargs.has_key('newWave'):
                deltaLambda = ((newWave.max()+newWave.min())/2.)/kwargs['smoothRes']
            else:
                deltaLambda = ((wave.max()+wave.min())/2.)/kwargs['smoothRes']
            sigma = deltaLambda / 2.3548
            spec = ndimage.gaussian_filter1d(spec, sigma)
        specs.append(spec)
    if kwargs.has_key('newWave'):
        return newWave, np.array(specs)
    else:
        return wave, np.array(specs)



specGridFunc = {
        'munari': getMunariSpecs,
        'marcs' : getMarcsSpecs,
        'specgrid1' : getSimpleSpecs
}

def getGridNames():
    return [os.path.basename(item.strip('.db3')) for item in glob(os.path.join(specGridDir, '*.db3'))]
    
    
def getGridDBConnection(gridName):
#    if gridName not in getGridNames():
#        raise ValueError("%s does not exist" % gridName)
#    else:
        return sqlite3.connect(os.path.join(specGridDir, gridName, 'index.db3'), detect_types=sqlite3.PARSE_DECLTYPES)

#

                
        
class specGrid(object):
    def __init__(self, gridName, params, paramLimits, mode='spec', **kwargs):
        """
            Function specgrid
            
            modify kwargs
            normRange
            normMode
            newWave
            smoothRes
        """
        
        config = ConfigParser.ConfigParser()
        config.read(os.path.join(specGridDir, gridName, 'config.ini'))
        
        #get spectralGrid
        # modelGrid
        # params specified ((param1,(lowerLimit, upperLimit)), (param2,
        sqlStmt = 'select %s ' % (', '.join(list(params) + ['fname']),)
        sqlStmt += ' from %s ' % gridName
        sqlCond = 'where '
        self.params = params
        self.paramLimits = paramLimits
        for param, paramLimit in zip(params, paramLimits):
            if paramLimit[0] != None:
                sqlCond += "%s>=%s and " % (param, paramLimit[0])
            if paramLimit[1] != None:
                sqlCond += "%s<=%s and " % (param, paramLimit[1])
            initParams = config.items('param-default')
            
        #converting the initial values to floats (if possible)
        for i in xrange(len(initParams)):
            try:
                initParams[i] = (initParams[i][0], float(initParams[i][1]))
            except ValueError:
                pass
            
        sqlCond += ' and '.join([ '%s=%s' % (initParam, initValue) for initParam, initValue in initParams if initParam not in params])
        sqlCond = sqlCond.strip().strip('and')
        sqlCond += " order by %s" % (', '.join(params),) 
        
        sqlStmt += sqlCond
        print "Reading from Database:\n%s\n\n" % sqlStmt
        conn = getGridDBConnection(gridName)
        print "Reading Points"
        DBData = conn.execute(sqlStmt).fetchall()
        fnames = zip(*DBData)[-1]
        self.points = np.array(zip(*zip(*DBData)[:-1]))
        print "Reading Values"
        self.wave, self.values = specGridFunc[gridName](gridName, fnames, config)
        self.interpGrid = interpolate.LinearNDInterpolator(self.points, self.values, fill_value=-1.)
        return None
    
    
    def __call__(self, *args):
        return self.interpGrid(args)
        
class reddenGrid(object):
    def __init__(self, ebvRange, wave, enableDIB=False, enableFlux=True, extinctionLaw = 'gal3'):
        self.points = np.array(ebvRange)
        values = []
        for ebv in ebvRange:
            extinctFlux = np.ones(wave.shape)
            if enableFlux:
                extinct = pysynphot.Extinction(ebv, extinctionLaw)    
                extinctThroughPut = extinct.GetThroughput()[::-1]
                f = interpolate.interp1d(extinct.wave[::-1], extinctThroughPut)
                extinctFlux *= f(wave)
            
            if enableDIB:
                extinctFlux *= pydib.makeSpectrum(wave, ebv).flux
            
            values.append(extinctFlux)
        self.values = np.array(values)
        self.interpGrid = interpolate.interp1d(self.points, self.values.transpose(), fill_value=1, bounds_error=False)
    
    def __call__(self, *args):
        return self.interpGrid(args)[:,0]
        
    
    
    
def showGrid(specGrid, sampleSpec=None, autoScale=True, fig=None):
    import pylab
    if fig==None:
        fig = pylab.figure(1)
    def sliderUpdate(val):
        spec = specGrid.interpGrid(tuple([slider.val for slider in sliders]))
        specPlot.set_ydata(spec)
        if autoScale:
            plotAxis.relim()
            plotAxis.autoscale_view()
        pylab.draw()
        
    
    #adding sliders
    sliders=[]
    paramNames = specGrid.params
    #paramInit = [5, 0]
    axes = [np.unique(item) for item in specGrid.points.transpose()]
    for i,axis in enumerate(axes):
        ax = fig.add_axes([0.1,i*0.05+0.05,0.7,0.03])
        sliders.append(Slider(ax, paramNames[i], min(axis), max(axis)))
        
        
    for slider in sliders:
        slider.on_changed(sliderUpdate)
        
        
    
    plotAxis = fig.add_axes([0.1, 0.1+0.15, 0.8, (1-0.1+0.1)*0.7])
#    if sampleStar != None:
#        plotAxis.plot(sampleStar.x, sampleStar.y, lw=3, color='black')
    spec = specGrid.interpGrid(tuple([slider.val for slider in sliders]))
    
    if sampleSpec != None:
        plotAxis.plot(sampleSpec.wave, sampleSpec.flux, color='black', lw=2)
    
    specPlot, = plotAxis.plot(specGrid.wave, spec)
    
    plotAxis.set_xlabel('Wavelength')
    plotAxis.set_ylabel('Intensity')
    pylab.show()    
    

sqlite3.register_converter("npmap", np.fromstring)