#==================================================================================================
#			HEADERS
#==================================================================================================
from __future__ import division
from struct import *
import numpy as np
import sys
import matplotlib as mpl
#mpl.use('Agg')		#No X11 mode
import os
import matplotlib.pyplot as plt
from pylab import *
import scipy.integrate as integ
import scipy.interpolate as interp
plt.close("all")


#==================================================================================================
#			VARIABLES
#==================================================================================================
#Data Fold
#datafold = "../../../box20VPH/"
datafold = "../data/"
#Codes Fold
codesfold = "./"


#==================================================================================================
#			FUNCTIONS
#==================================================================================================

#..................................................................................................
# Cutting the simulation box
#..................................................................................................
def slide( snapbase, snap, files, axis, coor, dx, sampling, type ):
    #Running code
    os.system( "%s/Cutter.out %s__%03d %d %d %f %f %d temp.tmp %d"%\
    (codesfold, snapbase, snap, files, axis, coor, dx, sampling, type) )
    data = np.loadtxt( "temp.tmp" )
    os.system( "rm temp.tmp" )
    
    return data
  
#..................................................................................................
# Saving data in ascii format
#..................................................................................................
def ascii( snapbase, snap, files, sampling, type ):
    #Running code
    os.system( "%s/Ascii.out %s__%03d %d %d temp.tmp %d"%\
    (codesfold, snapbase, snap, files, sampling, type) )
    data = np.loadtxt( "temp.tmp" )
    os.system( "rm temp.tmp" )
    
    return data
  
#..................................................................................................
# Finding distance to nearest halo in units of virial radius
#..................................................................................................
def gasdomains( X, Y, Z, sim, snapbase, snap, Lbox ):
    #Loading indexes
    index = np.loadtxt( "%s/%s/%s__%d.domain"%(datafold, sim, snapbase, snap) )
    #Loading catalogue of halos
    halos = np.loadtxt( "%s/%s/%s__%d.halo_catalog"%(datafold, sim, snapbase, snap) )
    
    #Calculating distances normalized with virial radius
    Ldomain = []
    for i in xrange( 0, len(X) ):
	#Position of domainant halo
	Xh = halos[index[i],0]; Yh = halos[index[i],1]; Zh = halos[index[i],2]
	#Virial radius
	Rvir = halos[index[i],4]
	#Distance
	Xrel = min( abs(X[i]-Xh), abs(X[i]-Xh+Lbox), abs(X[i]-Xh-Lbox) )
	Yrel = min( abs(Y[i]-Yh), abs(Y[i]-Yh+Lbox), abs(X[i]-Yh-Lbox) )
	Zrel = min( abs(Z[i]-Zh), abs(Z[i]-Zh+Lbox), abs(Z[i]-Zh-Lbox) )
	#Relative distance
	Ldomain.append( norm([Xrel,Yrel,Zrel])/Rvir )
	
    Ldomain = np.array(Ldomain)
    
    return Ldomain
  
