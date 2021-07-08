'''
DESCRIPTION:

usage:
python snapshotJoiner.py haloA haloB haloAB  0 0 0  0 0 0  0 0 0

where the zeros are meant to be relative positions in x y z,
relative velocities in vx vy vz, then angles of rotation
around the x, y and z axis.


Script that joins snapshots* and is able to write them, 
either from a set of text files, or by directly passing 
the necessary data to the function write_snapshot (for
details on the syntax, read the function).

A file named 'init.dat' will be created if no other name
is given.




*Modified to work with Python 3.8 and made from two
separate scripts (hence why there are two conventions
for writing the variable names and such).

Due to the changes made, it probably doesn't work with
Python 2.x anymore.


'''

#TODO make it more professional looking, adding "from argparse import ArgumentParser as parser"

import numpy as np
import pynbody as pnb

import sys
import os
import struct

# Creates a list containing the non-commented, non-empty lines
# of the input file for the header.
def process_input(file_):
    h_file = []
    input_ = open(file_, 'r')
    for line in input_:
        if line.find("#") != -1:
            continue
        elif line.find("\n") == 0:
            continue
        else:
            h_file.append(line.split('\t'))
    return h_file



def read_header_old(folder, n_part):
    h_file = process_input(folder + "/header.txt")
    h_data = []
    for j in n_part: # n_part
        h_data.append(int(j))
    for j in h_file[0][0:6]: # mass
        h_data.append(float(j))
    h_data.append(float(h_file[1][0])) # time
    h_data.append(float(h_file[2][0])) # redshift
    h_data.append(int(h_file[3][0])) # flag_sfr
    h_data.append(int(h_file[4][0])) # flag_feedback
    for j in n_part:
        h_data.append(int(j)) # n_part_total
    h_data.append(int(h_file[5][0])) # flag_coooling
    h_data.append(int(h_file[6][0])) # num_files
    h_data.append(float(h_file[7][0])) # box_size
    h_data.append(float(h_file[8][0])) # omega0
    h_data.append(float(h_file[9][0])) # omega_lambda
    h_data.append(float(h_file[10][0])) # hubble_param

    # blank, present in the header
    for i in np.arange(96):
        h_data.append(b'\x00')
    s = struct.Struct('iiiiii dddddd d d i i iiiiii i i dddd cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc')

    #h_data = [bytes(i, 'utf-8') for i in h_data]

    packed_data = s.pack(*h_data)
    return packed_data



#incorporated header
def read_header(n_part):
    h_data = []
    for j in n_part: # n_part
        h_data.append(int(j))
    for j in range(6): # mass table
        h_data.append(0.0)
    h_data.append(0.0) # time
    h_data.append(0.0) # redshift
    h_data.append(int(0)) # flag_sfr
    h_data.append(int(0)) # flag_feedback
    for j in n_part:
        h_data.append(int(j)) # n_part_total
    h_data.append(int(0)) # flag_coooling
    h_data.append(int(1)) # num_files
    h_data.append(0.0) # box_size
    h_data.append(0.0) # omega0
    h_data.append(0.0) # omega_lambda
    h_data.append(1.0) # hubble_param

    # blank, present in the header
    for i in np.arange(96):
        h_data.append(b'\x00')
    s = struct.Struct('iiiiii dddddd d d i i iiiiii i i dddd cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc')

    #h_data = [bytes(i, 'utf-8') for i in h_data]

    packed_data = s.pack(*h_data)
    return packed_data



def write_dummy(f, values_list):
    for i in values_list:
        dummy = [i]
        s = struct.Struct('i')
        d = s.pack(*dummy)
        f.write(d)



def write_block(f, block_data, data_type, block_name):
    write_dummy(f, [8])

    block_name_binary = [bytes(i, 'utf-8') for i in block_name]

    f.write(struct.pack('c' * 4, *block_name_binary))
    if(block_name == 'HEAD'):
        nbytes = 256
    else:
        fmt = data_type * len(block_data)
        nbytes = len(block_data) * 4
    write_dummy(f, [nbytes + 8, 8, nbytes]) 
    if(block_name == 'HEAD'):
        f.write(block_data) 
    else:
        f.write(struct.pack(fmt, *block_data))
    write_dummy(f, [nbytes])



def write_snapshot(n_part, folder=None, from_text=True, data_list=None, outfile='init.dat', old_header=False):
    N_gas = n_part[0]
    if(from_text and not folder):
        print ("error: can't call write_snapshot with from_text=True\n"
               "and without an input files folder.")
    if(not from_text):
        folder = os.getcwd()

    if old_header:
        header_data = read_header_old(folder, n_part)
    else:
        header_data = read_header(n_part)

    # Erasing the input file before opening it.
    with open(outfile, 'wb') as f:
        if(from_text):
            pos_data = np.fromfile(folder + "position.txt", sep='\t')
            vel_data = np.fromfile(folder + "velocity.txt", sep='\t')
            ID_data = np.fromfile(folder + "id.txt", dtype=int, sep='\t')
            mass_data = np.fromfile(folder + "masses.txt", sep='\t')
            if(N_gas > 0):
                U_data = np.fromfile(folder + "energy.txt", sep='\t')
                rho_data = np.fromfile(folder + "density.txt", sep='\t')
                smoothing_data = np.fromfile(folder + "smoothing.txt", sep='\t')
        else:
            pos_data = data_list[0]
            vel_data = data_list[1]
            ID_data = data_list[2]
            mass_data = data_list[3]
            if(N_gas > 0):
                U_data = data_list[4]
                rho_data = data_list[5]
                smoothing_data = data_list[6]

        write_block(f, header_data, None, 'HEAD')
        write_block(f, pos_data, 'f', 'POS ')
        write_block(f, vel_data, 'f', 'VEL ')
        write_block(f, ID_data, 'i', 'ID  ')
        write_block(f, mass_data, 'f', 'MASS')
        if(N_gas > 0):
            write_block(f, U_data, 'f', 'U   ')
            write_block(f, rho_data, 'f', 'RHO ')
            write_block(f, smoothing_data, 'f', 'HSML')



def rotation (vector, alpha=0, beta=0, gamma=0):
    """
    alpha: angle of rotation around the x axis
    beta: angle of rotation around the y axis
    gamma: angle of rotation around the z axis
    """
    
    vector = np.array(vector)

    #rotation matrix in x
    rAlpha = np.array([[1, 0, 0],
                       [0, np.cos(alpha), -np.sin(alpha)],
                       [0, np.sin(alpha), np.cos(alpha)]])
    
    #rotation matrix in y
    rBeta = np.array([[np.cos(beta), 0, -np.sin(beta)],
                      [0, 1, 0],
                      [-np.sin(beta), 0, np.cos(beta)]])

    #rotation matrix in z
    rGamma = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                       [np.sin(gamma), np.cos(gamma), 0],
                       [0, 0, 1]])

    rGeneral = np.matmul(np.matmul(rAlpha, rBeta), rGamma)

    return (np.matmul(rGeneral, np.array(vector)))



def join (snapshotZero, snapshotOne, output='init.dat', relativePos=[0.0, 0.0, 0.0],
          relativeVel=[0.0, 0.0, 0.0], rotationAngles=[0.0, 0.0, 0.0], shiftToCOM=True, writeNewSnapshot=True):
    """ 
    Joins snapshots and writes the result if writeNewSnapshot is True.

    The rotation will be applied to the second snapshot. Angles should
    be given in degrees.
    
    Rotation is applied using a for loop that transforms each vector
    with the function rotation(). It may be better to find a way to
    apply a rotation without using a for loop.
    """


    relativePos = [float(i) for i in relativePos]
    relativeVel = [float(i) for i in relativeVel]

    rotationAngles = np.radians([float(i) for i in rotationAngles])

    #standart families in gadget 2
    particleFamilies = ['gas', 'dm', 'disk', 'bulge', 'stars', 'bndry']


    #Loading snapshots
    snapshotZero = pnb.load(snapshotZero)
    snapshotOne = pnb.load(snapshotOne)


    #TODO solve the problem of joining two snapshots with non-standart families
    #Getting information and joining each type of particle
    nPart = []
    positions = []
    velocities = []
    masses = []
    energy = []
    rho = []
    smoothing = []


    for i in particleFamilies:
        existsInZero = i in [*snapshotZero.families()]
        existsInOne = i in [*snapshotOne.families()]

        if existsInZero and existsInOne:
            for j in snapshotZero.families():
                if j == i: #loop needed in order to iterate over the families
                    break

            nPart.append(len(snapshotZero[j]['iord']) + len(snapshotOne[j]['iord']))
            if np.any(rotationAngles): #applies rotation if any given angle is different of 0
                positions.append(np.concatenate((np.array(snapshotZero[j]['pos']),
                                                 np.array([rotation(h, *rotationAngles) for h in snapshotOne[j]['pos']]) + relativePos)))
                velocities.append(np.concatenate((np.array(snapshotZero[j]['vel']),
                                                 np.array([rotation(h, *rotationAngles) for h in snapshotOne[j]['vel']]) + relativeVel)))
            else:
                positions.append(np.concatenate((np.array(snapshotZero[j]['pos']), np.array(snapshotOne[j]['pos']) + relativePos)))
                velocities.append(np.concatenate((np.array(snapshotZero[j]['vel']), np.array(snapshotOne[j]['vel']) + relativeVel)))               
            masses.append(np.concatenate((np.array(snapshotZero[j]['mass']), np.array(snapshotOne[j]['mass']))))
            if i == 'gas':
                energy.append(np.concatenate((np.array(snapshotZero[j]['u']), np.array(snapshotOne[j]['u']))))
                rho.append(np.concatenate((np.array(snapshotZero[j]['rho']), np.array(snapshotOne[j]['rho']))))
                smoothing.append(np.concatenate((np.array(snapshotZero[j]['smooth']), np.array(snapshotOne[j]['smooth']))))
        
        elif existsInZero:
            for j in snapshotZero.families():
                if j == i:
                    break          
            nPart.append(len(snapshotZero[j]['iord']))
            positions.append(np.array(snapshotZero[j]['pos']))
            velocities.append(np.array(snapshotZero[j]['vel']))
            masses.append(np.array(snapshotZero[j]['mass']))
            if i == 'gas':
                energy.append(np.array(snapshotZero[j]['u']))
                rho.append(np.array(snapshotZero[j]['rho']))
                smoothing.append(np.array(snapshotZero[j]['smooth']))
        
        elif existsInOne:
            for j in snapshotOne.families():
                if j == i:
                    break
            nPart.append(len(snapshotOne[j]['iord']))
            if np.any(rotationAngles):
                positions.append(np.array([rotation(h, *rotationAngles) for h in snapshotOne[j]['pos']]) + relativePos)
                velocities.append(np.array([rotation(h, *rotationAngles) for h in snapshotOne[j]['vel']]) + relativeVel)
            else:
                positions.append(np.array(snapshotOne[j]['pos']) + relativePos)
                velocities.append(np.array(snapshotOne[j]['vel']) + relativeVel)                
            masses.append(np.array(snapshotOne[j]['mass']))
            if i == 'gas':
                energy.append(np.array(snapshotOne[j]['u']))
                rho.append(np.array(snapshotOne[j]['rho']))
                smoothing.append(np.array(snapshotOne[j]['smooth']))
        
        else:
            nPart.append(0)


    while len(nPart) < 6:
        nPart.append(0) #nPart needs to be a list with 6 objects


    dataList = []

    if ('gas' in [*snapshotZero.families()]) or ('gas' in [*snapshotZero.families()]):
        dataList = [np.concatenate(positions),
                    np.concatenate(velocities),
                    np.array(range(sum(nPart))) + 1,
                    np.concatenate(masses),
                    np.concatenate(energy),
                    np.concatenate(rho),
                    np.concatenate(smoothing)]
    
    else:
        dataList = [np.concatenate(positions),
                    np.concatenate(velocities),
                    np.array(range(sum(nPart))) + 1,
                    np.concatenate(masses)]




    #Shifting to center of mass
    if shiftToCOM:
        xCOM = sum(dataList[0][:, 0] * dataList[3]) / sum(dataList[0][:, 0])
        yCOM = sum(dataList[0][:, 1] * dataList[3]) / sum(dataList[0][:, 1])
        zCOM = sum(dataList[0][:, 2] * dataList[3]) / sum(dataList[0][:, 2])

        vxCOM = sum(dataList[1][:, 0] * dataList[3]) / sum(dataList[1][:, 0])
        vyCOM = sum(dataList[1][:, 1] * dataList[3]) / sum(dataList[1][:, 1])
        vzCOM = sum(dataList[1][:, 2] * dataList[3]) / sum(dataList[1][:, 2])

        dataList[0] = dataList[0] - [xCOM, yCOM, zCOM]
        dataList[1] = dataList[1] - [vxCOM, vyCOM, vzCOM]


    #Changing the shape for writting
    #the gadget format is composed of binary blocks of data which are written in a format like xyzxyz..xyz for positions and velocities
    dataList[0].shape = (-1, 1)
    dataList[1].shape = (-1, 1)

    if writeNewSnapshot:
        write_snapshot(n_part=nPart, from_text=False, outfile=output, data_list=dataList, old_header=False)
    
    else:
        return (nPart, dataList)



if __name__ == '__main__':

    join(*sys.argv[1:4], relativePos=sys.argv[4:7], relativeVel=sys.argv[7:10], rotationAngles=sys.argv[10:13])




