"""
==================================================================================================
NAME:		hdf52gadget.py
FUNCTION:	This script converts HDF5 outputs of gadget into the standar format gadget 1 or 2
ARGUMENTS:	snap_name, number of parts, output gadget format
USAGE:		python hdf52gadget.py snap_135 32 2
AUTHOR:		Sebastian Bustamante (macsebas33@gmail.com)
==================================================================================================
"""

#--------------------------------------------------------------------
#IMPORTS
#--------------------------------------------------------------------
import h5py
import numpy as np
import sys
import struct

#--------------------------------------------------------------------
#PARAMETERS
#--------------------------------------------------------------------
#Snapshot name
snapname = sys.argv[1]
#Number of hdf5 files
number_files = int(sys.argv[2])
#Format of output gadget file
gadgetfmt = int(sys.argv[3])

#--------------------------------------------------------------------
#CONVERTING
#--------------------------------------------------------------------
for i_file in xrange(number_files):
    #Opening hdf5 file
    fhdf = h5py.File('snap_135.%d.hdf5'%(i_file),'r')
    #Opening binary file for gadget file
    fbin = open('snap_135.%d'%(i_file),'wb')
    
    #................................................................
    #HEADER
    #................................................................
    #Loading header
    header = fhdf['/Header']
    #Writing header
    #Initial 
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('4B', *range(4)))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    
    #Particles in this file
    parts = header.attrs['NumPart_ThisFile']
    parts[2:] = 0	#Only DM and GAS particles!!
    N_part_file = sum(parts[2:])
    fbin.write(struct.pack('6i', *parts))
    #Masses of the particles
    mass = header.attrs['MassTable']
    fbin.write(struct.pack('6d', *mass))
    #Time
    time = header.attrs['Time']
    fbin.write(struct.pack('1d', time))
    #Redshift
    redshift = header.attrs['Redshift']
    fbin.write(struct.pack('1d', redshift))
    #Flag sfr
    flag_sfr = header.attrs['Flag_Sfr']
    fbin.write(struct.pack('1i', flag_sfr))
    #Flag feedback
    flag_feedback = header.attrs['Flag_Feedback']
    fbin.write(struct.pack('1i', flag_feedback))
    #Particles in the simulation
    parts_tot = header.attrs['NumPart_Total']
    parts_tot[2:] = 0	#Only DM and GAS particles!!
    fbin.write(struct.pack('6i', *parts_tot))
    #Flag cooling
    flag_cooling = header.attrs['Flag_Cooling']
    fbin.write(struct.pack('1i', flag_cooling))
    #Number of files
    num_files = header.attrs['NumFilesPerSnapshot']
    num_files = 1
    fbin.write(struct.pack('1i', num_files))
    #Box size
    boxL = header.attrs['BoxSize']
    fbin.write(struct.pack('1d', boxL))
    #Omega0
    Omega0 = header.attrs['Omega0']
    fbin.write(struct.pack('1d', Omega0))
    #OmegaLambda
    OmegaLambda = header.attrs['OmegaLambda']
    fbin.write(struct.pack('1d', OmegaLambda))
    #HubbleParam
    HubbleParam = header.attrs['HubbleParam']
    fbin.write(struct.pack('1d', HubbleParam))
    #Extra space
    fbin.write(struct.pack('96B', *xrange(96)))
    fbin.write(struct.pack('1i', 0))
    
    #................................................................
    #POSITIONS
    #................................................................
    fbin.write(struct.pack('1i', 0))
    #Saving name of this block
    block_name = [ord(l) for l in list("POS ")]
    fbin.write(struct.pack('4B', *block_name))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    
    for s_part in fhdf.keys()[1:3]:
	particles = fhdf['/%s'%(s_part)]
	positions = particles[u'Coordinates']
	for i_part in xrange(parts[ int(s_part[-1]) ]):
	    i_position = positions[i_part]
	    fbin.write(struct.pack('3f', *i_position))
	    
    fbin.write(struct.pack('1i', 0))
    
    #................................................................
    #VELOCITIES
    #................................................................
    fbin.write(struct.pack('1i', 0))
    #Saving name of this block
    block_name = [ord(l) for l in list("VEL ")]
    fbin.write(struct.pack('4B', *block_name))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    
    for s_part in fhdf.keys()[1:3]:
	particles = fhdf['/%s'%(s_part)]
	velocities = particles[u'Velocities']
	for i_part in xrange(parts[ int(s_part[-1]) ]):
	    i_velocities = velocities[i_part]
	    fbin.write(struct.pack('3f', *i_velocities))
	    
    fbin.write(struct.pack('1i', 0))
    
    #................................................................
    #IDs
    #................................................................
    fbin.write(struct.pack('1i', 0))
    #Saving name of this block
    block_name = [ord(l) for l in list("ID  ")]
    fbin.write(struct.pack('4B', *block_name))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    fbin.write(struct.pack('1i', 0))
    
    for s_part in fhdf.keys()[1:3]:
	particles = fhdf['/%s'%(s_part)]
	ids = particles[u'ParticleIDs']
	for i_part in xrange(parts[ int(s_part[-1]) ]):
	    i_ids = ids[i_part]
	    fbin.write(struct.pack('1i', 1))
	    
    fbin.write(struct.pack('1i', 0))
    
    fbin.close()