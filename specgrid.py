import pyspec
from pyspec import oned
import os
from scipy import interpolate, ndimage, optimize
from scipy.interpolate import LinearNDInterpolator
import numpy as np
import sqlite3
from glob import glob
import pysynphot
#import pydib
import ConfigParser
import pdb
import minuit
import time
from matplotlib.widgets import Slider
from iwfinterpolate import Invdisttree
try:
    import pydib
except ImportError:
    print "PyDIB not available"
    pass
defaultSpecGridDir = os.path.expanduser('~/.specgrids')


#specGridDir = '/Users/wkerzend/projects/grids'

def getSpecs(specDataDir, fnames, config, **kwargs):

    #reading config
    waveFName = config.get('wave', 'wave')
    specDType = config.get('structure', 'datatype')
    specSize = config.getint('structure', 'specsize')
    specsize = config.getint('wave', 'islog')
    R = config.getfloat('params', 'r')
    
    if kwargs.has_key('smoothres'):
        if kwargs['smoothres']>R:
            raise ValueError('requested resolution (R=%f) higher than'
                             'models intrinsic resolution (R=%f)' % (kwargs['smoothres'], R))
    #Reading wave solution
    wave = np.fromfile(os.path.join(specDataDir, waveFName))
    
    if kwargs.has_key('wave'):
        gridSize = len(fnames) * len(kwargs['wave'].tostring()) / 1024**2
    else:
        gridSize = len(fnames) * specSize / 1024**2
    print "Processing %d spectra" % len(fnames)
    
    print "Processing %.3f MB in to Memory." % gridSize
    
    specs = []
    if kwargs.has_key('normrange'):
        if kwargs.has_key('normmode'):
            normMode = kwargs['normmode']
        else:
            normMode = 'simple'
        
        if normMode == 'simple':
            normSlice = slice(*[wave.searchsorted(item) for item in kwargs['normrange']])
        else:
            raise NotImplementedError('Normalisation mode %s has not been implemented yet' % normMode)
    startTime = time.time()
    for i, specFName in enumerate(fnames):
        flux = np.fromfile(os.path.join(specDataDir, 'data', specFName))
        spec = oned.onedspec(wave, flux, mode='waveflux')
        if i%100 == 0:
            print "@%d took %.2f s" % (i, time.time() - startTime)
            startTime = time.time()
        
            
        if kwargs.has_key('normrange'):
            normFac = np.mean(spec[normSlice].flux)
            spec /= normFac
            
            
        if kwargs.has_key('smoothres') or kwargs.has_key('smoothrot'):
            if kwargs.has_key('wave'):
                tmpSpec = spec[float(kwargs['wave'].min()):float(kwargs['wave'].max())]
                logDelta, logSpec = tmpSpec.interpolate_log()
            else:
                logDelta, logSpec = spec.interpolate_log()
            if kwargs.has_key('smoothres'):
                logSpec = logSpec.convolve_profile(kwargs['smoothres'], smallDelta=logDelta)
            if kwargs.has_key('smoothrot'):
                logSpec = logSpec.convolve_rotation(kwargs['smoothrot'], smallDelta=logDelta)
            spec = logSpec
            
        if kwargs.has_key('wave'):
            spec = spec.interpolate(kwargs['wave'])

            
        specs.append(spec.flux)
    if kwargs.has_key('wave'):
        return kwargs['wave'], np.array(specs)
    else:
        return wave, np.array(specs)
        
    

def getGridNames():
    return [os.path.basename(item.strip('.db3')) for item in glob(os.path.join(specGridDir, '*.db3'))]
    
    
def getGridDBConnection(gridName, specGridDir=None):
    if specGridDir == None:
        specGridDir = defaultSpecGridDir
#    if gridName not in getGridNames():
#        raise ValueError("%s does not exist" % gridName)
#    else:
    return sqlite3.connect(os.path.join(specGridDir, gridName, 'index.db3'), detect_types=sqlite3.PARSE_DECLTYPES)

#

class minuitFunction2(object):
    def __init__(self, specGrid, sampleSpec):
        self.specGrid = specGrid
        self.sampleSpec = sampleSpec
    class func_code:
        co_varnames = []
        co_argcount = 0

    def varnames(self, *args):
        self.func_code.co_varnames = args
        self.func_code.co_argcount = len(args)

    def __call__(self, *args):
        if len(args) != self.func_code.co_argcount:
            raise TypeError, "wrong number of arguments"
        
        sampleSpec = self.sampleSpec
        
        
        if 'vrot' in self.func_code.co_varnames:
            vrot = args[-1]
            vrot = np.max((2., np.abs(vrot)))
            gridSpec = self.specGrid.getSpec(*args[:-1])
            gridSpec = gridSpec.convolve_rotation(abs(vrot))
            gridSpec = gridSpec.interpolate(sampleSpec.wave)
        else:
            gridSpec = self.specGrid.getSpec(*args)
            gridSpec = gridSpec.interpolate(sampleSpec.wave)
            
        
        
        
        if sampleSpec.var != None:
            var = sampleSpec.var
        else:
            var = 1.

        if  sampleSpec.dq != None:
            dqMask = sampleSpec.dq
        else:
            dqMask = np.ones(grid.shape[1]).astype(bool)
        
        
            
        chiSq = ((gridSpec.flux[dqMask]-sampleSpec.flux[dqMask])/var[dqMask])**2
        nu = gridSpec.wave.shape[0] - self.func_code.co_argcount - 1
        redChiSq = np.sum(chiSq) / nu
        
        return redChiSq
    
    
class minuitFunction(object):
    def __init__(self, specGrid, sampleSpec, priors=None):
        self.specGrid = specGrid
        self.sampleSpec = sampleSpec
        self.priors = priors
    class func_code:
        co_varnames = []
        co_argcount = 0

    def varnames(self, *args):
        self.func_code.co_varnames = args
        self.func_code.co_argcount = len(args)

    def __call__(self, *args):
        if len(args) == 1:
            print "warning attempting to unwrap args"
            args = args[0]
        if len(args) != self.func_code.co_argcount:
            raise TypeError, "wrong number of arguments"
        
        sampleSpec = self.sampleSpec
        priors = self.priors
        
        if 'vrot' in self.func_code.co_varnames:
            vrot = args[-1]
            vrot = np.max((2., np.abs(vrot)))
            gridSpec = self.specGrid.getSpec(*args[:-1])
            if vrot != 0.:
                gridSpec = gridSpec.convolve_rotation(abs(vrot))
            gridSpec = gridSpec.interpolate(sampleSpec.wave)
        else:
            #pdb.set_trace()
            gridSpec = self.specGrid.getSpec(*args)
            gridSpec = gridSpec.interpolate(sampleSpec.wave)
            
        
        
        
        if sampleSpec.var != None:
            var = sampleSpec.var
        else:
            var = 1.

        if  sampleSpec.dq != None:
            dqMask = sampleSpec.dq
        else:
            dqMask = np.ones(grid.shape[1]).astype(bool)
        
        
            
        chiSq = ((gridSpec.flux[dqMask]-sampleSpec.flux[dqMask])/var[dqMask])**2
        nu = gridSpec.wave.shape[0] - self.func_code.co_argcount - 1
        
        chi_prior = 0.
        if priors != None:
            for i, prior in enumerate(priors):
                if prior == None: continue
                chi_prior += ((args[i] - prior[0]) /  prior[1])**2
                nu += 1
                
        
        redChiSq = (np.sum(chiSq) + chi_prior) / nu
        
        return redChiSq

        
class specGrid(object):
    def __init__(self, gridName, params, paramLimits,
                 mode='spec', specGridDir=None, interpolator=Invdisttree,
                 **kwargs):
        """
            Function specgrid
            
            modify kwargs
            normRange
            normMode
            newWave
            smoothRes
        """
        if specGridDir==None:
            specGridDir = defaultSpecGridDir
        
        
        config = ConfigParser.ConfigParser()
        config.read(os.path.join(specGridDir, gridName, 'config.ini'))
        gridTableName = config.get('structure', 'table')
        #get spectralGrid
        # modelGrid
        # params specified ((param1,(lowerLimit, upperLimit)), (param2,
        self.params = params
        if 'ebv' in params:
            if params.index('ebv') != len(params) - 1:
                raise IndexError('ebv param must be at the end')
            if len(paramLimits[-1]) != 3:
                raise ValueError('paramLimits for ebv must consist of start end and steps')
            params = self.params[:-1]
        
        #SQL Statement assembly
        sqlStmt = 'select %s ' % (', '.join(list(params) + ['fname']),)
        sqlStmt += ' from %s ' % gridTableName
        sqlCond = ''
        self.paramLimits = paramLimits
        for param, paramLimit in zip(params, paramLimits):
            if paramLimit[0] != None:
                if sqlCond == '': sqlCond = 'where '
                sqlCond += "%s>=%s and " % (param, paramLimit[0])
            if paramLimit[1] != None:
                if sqlCond == '': sqlCond = 'where '
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
        conn = getGridDBConnection(gridName, specGridDir)
        print "Reading Points"
        DBData = conn.execute(sqlStmt).fetchall()
        fnames = zip(*DBData)[-1]
        self.points = np.array(zip(*zip(*DBData)[:-1]))
        self.limits = [(np.min(item), np.max(item)) for item in self.points.transpose()] 
        print "Reading Values"

        specGridFunction = getSpecs
        self.wave, self.values = specGridFunction(os.path.join(specGridDir, gridName), fnames, config, **kwargs)
        self.interpGrid = interpolator(self.points, self.values)
        
        if 'ebv' in self.params:
            self.redGrid = reddenGrid(np.linspace(*paramLimits[-1]), self.wave, enableFlux=True, enableDIB=True)            
            
            
            
        return None
    
    
    def __call__(self, *args):
        if 'ebv' in self.params:
            curParams = args[:-1]
            redFactor = self.redGrid(args[-1])
            return redFactor*self.interpGrid(curParams)
        else:    
            return self.interpGrid(args)
    def getSpec(self, *args):
        return oned.onedspec(self.wave, self(args), mode='waveflux')
    
    def fitChiSq(self, sampleSpec, returnChiSq=False, priors = None):

        if np.all(sampleSpec.wave == self.wave):
            newSampleSpec = sampleSpec
            grid = self.values
        else:
            newSampleSpec = sampleSpec.interpolate(self.wave)          
            grid = self.values

        
        if newSampleSpec.var != None:
            var = newSampleSpec.var
        else:
            var = 1.

        if  newSampleSpec.dq != None:
            dqMask = newSampleSpec.dq
        else:
            dqMask = np.ones(grid.shape[1]).astype(bool)
        
            
        chiSq = np.sum(((grid[:,dqMask]-newSampleSpec.flux[dqMask])/var[dqMask])**2, axis=1)
        nu = (np.ones(grid.shape[0])*grid.shape[1]) - len(self.params) - 1
        chi_prior = np.zeros(self.points.shape[0])
        if priors != None:
            for i, prior in enumerate(priors):
                if prior == None: continue
                chi_prior += ((self.points[:,i] - prior[0]) /  prior[1])**2
                nu += 1
                
        
        redChiSq = (chiSq + chi_prior) / nu
        
        
            

        if returnChiSq:
            return np.min(redChiSq), self.points[np.argmin(redChiSq)]
        else:
            return self.points[np.argmin(redChiSq)]
            
    def fitBFGS(self, sampleSpec, initValues=None, vrot=None, priors=None):
        
        f = minuitFunction(self, sampleSpec, priors=priors)
        
        if vrot != None:
            varnames = list(self.params) + ['vrot']
            f.varnames(*varnames)
            return optimize.fmin_bfgs(f, initValues)
        else:
            f.varnames(*self.params)
            return optimize.fmin_bfgs(f, initValues)
    
    def getMinuitFunction(self, sampleSpec):
        pass
    def getMinuit(self, sampleSpec, init=None, vrot=None, priors=None):
        
        f = minuitFunction(self, sampleSpec, priors=priors)
        
        if vrot != None:
            varnames = list(self.params) + ['vrot']
            f.varnames(*varnames)
            fmin = minuit.Minuit(f)
            fmin.values['vrot']=vrot
        else:
            f.varnames(*self.params)
            fmin = minuit.Minuit(f)
        
        if init!=None:
            fmin.values.update(dict(zip(self.params, init)))
        
        return fmin
    
        
        
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