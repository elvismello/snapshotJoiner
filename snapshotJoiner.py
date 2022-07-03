'''
DESCRIPTION:

usage: python snapshotJoiner.py haloA haloB haloAB  0 0 0  0 0 0  0 0 0

where the zeros are meant to be relative positions in x y z, relative
velocities in vx vy vz, then angles of rotation around the x, y and z axis.


Script that joins snapshots and writes them by directly passing the necessary
data to the function write_snapshot.

'''


import numpy as np
import h5py
import snapwrite
import argparse


def rotation (vector, alpha=0.0, beta=0.0, gamma=0.0, returnMatrix=False, dtype='float64'):
    """
    alpha: angle of rotation around the x axis
    beta: angle of rotation around the y axis
    gamma: angle of rotation around the z axis
    """
    # It may be better to find a way to apply a rotation without using an external for loop.
    vector = np.array(vector)

    #rotation matrix in x
    rAlpha = np.array([[1,             0,              0],
                       [0, np.cos(alpha), -np.sin(alpha)],
                       [0, np.sin(alpha), np.cos(alpha)]], dtype=dtype)
    
    #rotation matrix in y
    rBeta = np.array([[np.cos(beta),  0, -np.sin(beta)],
                      [0,             1,             0],
                      [-np.sin(beta), 0, np.cos(beta)]], dtype=dtype)

    #rotation matrix in z
    rGamma = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                       [np.sin(gamma),  np.cos(gamma), 0],
                       [0,                          0, 1]], dtype=dtype)

    rGeneral = np.matmul(np.matmul(rGamma, rBeta), rAlpha, dtype=dtype)
    

    if not returnMatrix:
        return np.matmul(rGeneral, np.array(vector), dtype=dtype)
    else:
        return rGeneral



def join (snapshotZero, snapshotOne, output='init.dat',
          relativePos=[0.0, 0.0, 0.0], relativeVel=[0.0, 0.0, 0.0],
          rotationAngles=[0.0, 0.0, 0.0], shiftToCOM=True,
          writeNewSnapshot=True, outputFormat='gadget2',
          includeHaloZero=True, metallicity_in_everything=False):
    
    """ 
    Joins snapshots and writes the result as a new snapshot if
    writeNewSnapshot is True.

    The rotation will be applied to the second snapshot. Angles should
    be given in degrees.
    
    Rotation is applied using a for loop that transforms each vector
    with the function rotation().
    """


    relativePos = [float(i) for i in relativePos]
    relativeVel = [float(i) for i in relativeVel]

    rotationAngles = np.radians([float(i) for i in rotationAngles])

    #standard families in gadget 2
    particleTypes = ['PartType0', 'PartType1', 'PartType2', 'PartType3', 'PartType4', 'PartType5']


    #Loading snapshots
    snapshotZero = h5py.File(snapshotZero, 'r')
    snapshotOne = h5py.File(snapshotOne, 'r')
    print('Snapshots loaded!')

    #Getting information and joining each type of particle
    nPart = []
    positions = []
    velocities = []
    masses = []
    energy = []
    rho = []
    smoothing = []
    metallicity = []


    for i in particleTypes:
        print(f'Joining {i}...')

        if i != 'PartType1' or includeHaloZero: # skips writing the halo from the first snapshot if includeHaloZero is False
            existsInZero = i in [*snapshotZero.keys()]
        else:
            existsInZero = False
                
        existsInOne = i in [*snapshotOne.keys()]

        if existsInZero and existsInOne:
            # Possibly not needed
            #for j in snapshotZero.keys():
            #    if j == i: #loop needed in order to iterate over the families
            #        break

            nPart.append(len(snapshotZero[i]['ParticleIDs']) + len(snapshotOne[i]['ParticleIDs']))
            # Position and velocities
            if np.any(rotationAngles): #applies rotation if any given angle is different of 0
                positions.append(np.concatenate((np.array(snapshotZero[i]['Coordinates']),
                                                 np.array([rotation(h, *rotationAngles) for h in snapshotOne[i]['Coordinates']]) + relativePos)))
                velocities.append(np.concatenate((np.array(snapshotZero[i]['Velocities']),
                                                 np.array([rotation(h, *rotationAngles) for h in snapshotOne[i]['Velocities']]) + relativeVel)))
            else:
                positions.append(np.concatenate((np.array(snapshotZero[i]['Coordinates']), np.array(snapshotOne[i]['Coordinates']) + relativePos)))
                velocities.append(np.concatenate((np.array(snapshotZero[i]['Velocities']), np.array(snapshotOne[i]['Velocities']) + relativeVel)))               
            # Masses
            masses.append(np.concatenate((np.array(snapshotZero[i]['Masses']), np.array(snapshotOne[i]['Masses']))))
            # Gas Properties
            if i == 'PartType0':
                energy.append(np.concatenate((np.array(snapshotZero[i]['InternalEnergy']), np.array(snapshotOne[i]['InternalEnergy']))))
                rho.append(np.concatenate((np.array(snapshotZero[i]['Density']), np.array(snapshotOne[i]['Density']))))
                smoothing.append(np.concatenate((np.array(snapshotZero[i]['SmoothingLength']), np.array(snapshotOne[i]['SmoothingLength']))))
            # Metallicity
            if i != 'PartType1' and (metallicity_in_everything or i == 'PartType0'):
                metallicity.append(np.concatenate((np.array(snapshotZero[i]['Metallicity']), np.array(snapshotOne[i]['Metallicity']))))

        
        elif existsInZero:
            nPart.append(len(snapshotZero[i]['ParticleIDs']))
            positions.append(np.array(snapshotZero[i]['Coordinates']))
            velocities.append(np.array(snapshotZero[i]['Velocities']))
            masses.append(np.array(snapshotZero[i]['Masses']))
            if i == 'PartType0':
                energy.append(np.array(snapshotZero[i]['InternalEnergy']))
                rho.append(np.array(snapshotZero[i]['Density']))
                smoothing.append(np.array(snapshotZero[i]['SmoothingLength']))
            if i != 'PartType1' and (metallicity_in_everything or i == 'PartType0'):
                metallicity.append(np.array(snapshotZero[i]['Metallicity']))
        
        elif existsInOne:
            nPart.append(len(snapshotOne[i]['ParticleIDs']))
            if np.any(rotationAngles):
                positions.append(np.array([rotation(h, *rotationAngles) for h in snapshotOne[i]['Coordinates']]) + relativePos)
                velocities.append(np.array([rotation(h, *rotationAngles) for h in snapshotOne[i]['Velocities']]) + relativeVel)
            else:
                positions.append(np.array(snapshotOne[i]['Coordinates']) + relativePos)
                velocities.append(np.array(snapshotOne[i]['Velocities']) + relativeVel)                
            masses.append(np.array(snapshotOne[i]['Masses']))
            if i == 'PartType0':
                energy.append(np.array(snapshotOne[i]['InternalEnergy']))
                rho.append(np.array(snapshotOne[i]['Density']))
                smoothing.append(np.array(snapshotOne[i]['SmoothingLength']))
            if i != 'PartType1' and (metallicity_in_everything or i == 'PartType0'):
                metallicity.append(np.array(snapshotOne[i]['Metallicity']))
        
        else:
            nPart.append(0)


    while len(nPart) < 6:
        nPart.append(0) #nPart needs to be a list with 6 objects


    dataList = []

    if ('PartType0' in [*snapshotZero.keys()]) or ('PartType0' in [*snapshotZero.keys()]):
        dataList = [np.concatenate(positions),
                    np.concatenate(velocities),
                    np.array(range(sum(nPart))) + 1,
                    np.concatenate(masses),
                    np.concatenate(energy),
                    np.concatenate(rho),
                    np.concatenate(smoothing),
                    np.concatenate(metallicity)]

    else:
        dataList = [np.concatenate(positions),
                    np.concatenate(velocities),
                    np.array(range(sum(nPart))) + 1,
                    np.concatenate(masses)]




    #Shifting to center of mass
    if shiftToCOM:
        print('Shifting coordinates to center of mass...')
        xCOM  = sum(dataList[0][:, 0] * dataList[3]) / sum(dataList[3])
        yCOM  = sum(dataList[0][:, 1] * dataList[3]) / sum(dataList[3])
        zCOM  = sum(dataList[0][:, 2] * dataList[3]) / sum(dataList[3])

        vxCOM = sum(dataList[1][:, 0] * dataList[3]) / sum(dataList[3])
        vyCOM = sum(dataList[1][:, 1] * dataList[3]) / sum(dataList[3])
        vzCOM = sum(dataList[1][:, 2] * dataList[3]) / sum(dataList[3])

        dataList[0] = dataList[0] - [xCOM, yCOM, zCOM]
        dataList[1] = dataList[1] - [vxCOM, vyCOM, vzCOM]




    #Changing the shape for writting
    #the gadget format is composed of binary blocks of data which are written like xyzxyz..xyz for positions and velocities
    # snapwrite is able to write this format of lists into either hdf5 or gadget2 format
    dataList[0].shape = (-1, 1)
    dataList[1].shape = (-1, 1)

    if writeNewSnapshot:
        print('Writing snapshot...')
        snapwrite.write_snapshot(n_part=nPart, outfile=output, data_list=dataList, file_format=outputFormat)

    else:
        return nPart, dataList

    print('Done!')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='placeholder')

    parser.add_argument('snapshot0', help='The name of the first input file.')
    parser.add_argument('snapshot1', help='The name of the second input file.')
    parser.add_argument('-o', '--output', default='init.ic', help='The name of\
                        the output file')
    parser.add_argument('-rP', '--relative-position', nargs=3,
                        metavar=('X', 'Y', 'Z'), default=[0.0, 0.0, 0.0],
                        help='Relative position of the second snapshot with\
                         relation to the first. Must be given in terms of it\'s\
                         components in x, y and z')
    parser.add_argument('-rV', '--relative-velocity', nargs=3,
                        metavar=('vX', 'vY', 'vZ'), default=[0.0, 0.0, 0.0],
                        help='Relative velocity of the second snapshot with\
                         relation to the first. Must be given in terms of it\'s\
                         components in x, y and z')
    parser.add_argument('-r', '--rotation', nargs=3,
                        metavar=('angleX', 'angleY', 'angleZ'),
                        default=[0.0, 0.0, 0.0], help='Angles (in degrees) of\
                         the rotation to be applied to the second snapshot.\
                         Must be given in terms of rotations around the x\',\
                         y\' and z\' axis that pass by the origin of the second\
                         snapshot')
    parser.add_argument('--hdf5', action='store_true', help='Output initial\
                        conditions in HDF5')
    parser.add_argument('--noMainHalo', action='store_false', help='This will\
                         make the program skip the halo of dark matter in the\
                         first snapshot given.')
    parser.add_argument('--noCOMshift', action='store_false', help='This will\
                         skip the final shift into the center of mass.')

    args = parser.parse_args()

    if args.hdf5:
        outputFormat = 'hdf5'
    else:
        outputFormat = 'gadget2'

    join(args.snapshot0, args.snapshot1, args.output,
         relativePos=args.relative_position,
         relativeVel=args.relative_velocity,
         rotationAngles=args.rotation,
         outputFormat=outputFormat,
         includeHaloZero=args.noMainHalo,
         shiftToCOM=args.noCOMshift)




