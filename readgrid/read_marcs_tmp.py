import sqlite3
import gzip
import re
import os
from glob import glob

def makeMarcsDB(dbName):
    os.system('rm %s' % dbName)
    conn = sqlite3.connect(dbName)
    
    sqlStmt = ("create table %s( id integer primary key,"
            "teff float,"
            "logg float,"
            "feh float,"
            "alpha float,"
            "c float,"
            "n float,"
            "o float,"
            "r float,"
            "s float,"
            "mass float,"
            "micro float,"
            "geometry text,"
            "spec npmap);"
            "create table wave (wave npmap);" % dbName.strip('.db3'))
    conn.executescript(sqlStmt)
    conn.commit()
    conn.close()
    



def getSpecList2(specPath,):
    specPattern = re.compile('([sp])(\d+)g([pm])(\d+)([pm])(\d+)a(\d+)t(\d+).flx')
    convSign = {'m':-1, 'p':1}
    params = []
    specFiles = []
    for fname in glob(os.path.join(specPath, '*.flx')):
        specMatch = specPattern.match(os.path.basename(fname))
        if specMatch == None: continue
        geom, teff, loggSign, logg, fehSign, feh, alpha, micro = specMatch.groups()
        teff = float(teff)
        logg = convSign[loggSign] * (float(logg)/10)
        feh = convSign[fehSign] * float(feh)/1000.
        alpha = float(alpha) / 10
        micro = float(micro)
        specFiles.append(fname)
        params.append((teff, logg, feh, alpha, c, n, o, r, s, mass, micro, geom))
    return specFiles, params

def getSpecList (specPath, separator='_'):
    specPattern = re.compile('([sp])(\d+){sep}g([+-]?\d+\.\d+){sep}m(\d+\.\d+){sep}t(\d+){sep}st{sep}'
                             'z([+-]?\d+\.\d+){sep}'
                             'a([+-]?\d+\.\d+){sep}'
                             'c([+-]?\d+\.\d+){sep}'
                             'n([+-]?\d+\.\d+){sep}'
                             'o([+-]?\d+\.\d+){sep}'
                             'r([+-]?\d+\.\d+){sep}'
                             's([+-]?\d+\.\d+).flx'.format(sep=separator))
    
    params = []
    specFiles = []
    for fname in glob(os.path.join(specPath, '*.flx')):
        specMatch = specPattern.match(os.path.basename(fname))
        if specMatch == None: continue
        geom, teff, logg, mass, micro, feh, alpha, c, n, o, r, s = specMatch.groups()
        teff = float(teff)
        logg = float(logg)
        feh = float(feh)
        mass = float(mass)
        micro = float(micro)
        alpha = float(alpha)
        c = float(c)
        n = float(n)
        o = float(o)
        r = float(r)
        s = float(s)
        params.append((teff, logg, feh, alpha, c, n, o, r, s, mass, micro, geom))
        specFiles.append(fname)
    return specFiles, params

def insert2DB(dbName, params, specFiles, wave):
    conn= sqlite3.connect(dbName)
    conn.execute('insert into wave(wave) values (?)', (sqlite3.Binary(wave.tostring()), ))
    for param, specFile in zip(params, specFiles):
        print "@%s" % os.path.basename(specFile)
        spec = sqlite3.Binary(loadtxt(specFile).tostring())
        conn.execute('insert into marcs_tmp(teff, logg, feh, alpha, c, n, o, r, s, mass, micro, geometry, spec) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', tuple(list(param)+[spec]))
    conn.commit()
    conn.close()
    