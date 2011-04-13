import sqlite3
import gzip
import re

from glob import glob

def makeMarcsDB(dbName, tableName="marcs"):
    conn = sqlite3.connect(dbName)
    sqlStmt = "create table %s( id integer primary key,"
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
            "create table wave (wave npmap);" % tableName
    conn.executescript(sqlStmt)
    conn.commit()
    conn.close()
    
def getSpecList(path, separationChar = "_"):
    for fname in glob(path):