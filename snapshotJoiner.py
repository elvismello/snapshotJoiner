'''
DESCRIPTION:

usage:
python snapshotJoiner.py haloA haloB haloAB  0 0 0  0 0 0  0 0 0

where the zeros are meant to be relative positions in x y z,
relative velocities in vx vy vz, then angles of rotation
around the x, y and z axis.


Script that joins snapshots and writes them by directly
passing the necessary data to the function write_snapshot.

'''


import numpy as np
import pynbody as pnb
import snapwrite


def rotation (vector, alpha=0, beta=0, gamma=0, returnMatrix=False):
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
                       [0, np.sin(alpha), np.cos(alpha)]])
    
    #rotation matrix in y
    rBeta = np.array([[np.cos(beta),  0, -np.sin(beta)],
                      [0,             1,             0],
                      [-np.sin(beta), 0, np.cos(beta)]])

    #rotation matrix in z
    rGamma = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                       [np.sin(gamma),  np.cos(gamma), 0],
                       [0,                          0, 1]])

    rGeneral = np.matmul(np.matmul(rAlpha, rBeta), rGamma)

    if not returnMatrix:
        return np.matmul(rGeneral, np.array(vector))
    else:
        return rGeneral



def join (snapshotZero, snapshotOne, output='init.dat',
          relativePos=[0.0, 0.0, 0.0], relativeVel=[0.0, 0.0, 0.0],
          rotationAngles=[0.0, 0.0, 0.0], shiftToCOM=True,
          writeNewSnapshot=True, outputFormat='gadget2',
          includeHaloZero=True):
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

        if i != 'dm' or includeHaloZero: # skips writing the halo from the first snapshot if includeHaloZero is False
            existsInZero = i in [*snapshotZero.families()]
        else:
            existsInZero = False
                
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
    # TODO make this part more readable
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
        snapwrite.write_snapshot(n_part=nPart, outfile=output, data_list=dataList, file_format=outputFormat)
    
    else:
        return nPart, dataList



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='placeholder')

    parser.add_argument('snapshot0', help='The name of the first input file.')
    parser.add_argument('snapshot1', help='The name of the second input file.')
    parser.add_argument('-o', '--output', default='init.ic', help='The name of\
                        the output file')
    parser.add_argument('-rP', '--relative-position', nargs='+',
                        default=[0.0, 0.0, 0.0], help='Relative position of\
                         the second snapshot with relation to the first.\
                         Must be given in terms of it\'s components in x, y \
                         and z')
    parser.add_argument('-rV', '--relative-velocity', nargs='+',
                        default=[0.0, 0.0, 0.0], help='Relative velocity of the\
                         second snapshot with relation to the first. Must be \
                         given in terms of it\'s components in x, y and z')
    parser.add_argument('-r', '--rotation', nargs='+',
                        default=[0.0, 0.0, 0.0], help='Angle of the rotation\
                         to be applied to the second snapshot. Must be given\
                         in terms of rotations around the x, y and z axis that\
                         pass by the origin of the second snapshot')
    parser.add_argument('--hdf5', action='store_true', help='Output initial\
                        conditions in HDF5')
    parser.add_argument('--noMainHalo', action='store_false', help='This will\
                         make the program skip the halo of dark matter in the\
                         first snapshot parsed.')

    args = parser.parse_args()

    if args.hdf5:
        outputFormat = 'hdf5'
    else:
        outputFormat = 'gadget2'

    join(args.snapshot0, args.snapshot1, args.output,
         relativePos=args.relative_position,
         relativeVel=args.relative_velocity,
         rotationAngles=args.rotation, outputFormat=outputFormat,
         includeHaloZero=args.noMainHalo)




