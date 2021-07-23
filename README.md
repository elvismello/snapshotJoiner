## snapshotJoiner
Self contained program written in python 3 that can join snapshots using Gadget2 binary format.

Based in the module "snapwrite" contained in [Rafael Ruggiero](https://github.com/ruggiero)'s IC generation programs, but modified to work in python 3.


## Usage
```
usage: snapshotJoiner.py [-h] [-o OUTPUT] [-rP RELATIVE_POSITION [RELATIVE_POSITION ...]]
                         [-rV RELATIVE_VELOCITY [RELATIVE_VELOCITY ...]] [-r ROTATION [ROTATION ...]] [--hdf5]
                         snapshot0 snapshot1

placeholder

positional arguments:


  snapshot0             The name of the first input file.
  
  snapshot1             The name of the second input file.

optional arguments:


  -h, --help            show this help message and exit
  
  
  -o OUTPUT, --output OUTPUT
                        The name of the output file
  
  
  -rP RELATIVE_POSITION [RELATIVE_POSITION ...], --relative-position RELATIVE_POSITION [RELATIVE_POSITION ...]
                        Relative position of the second snapshot with relation to the first. Must be given in terms
                        of it's components in x, y and z
                        
                        
  -rV RELATIVE_VELOCITY [RELATIVE_VELOCITY ...], --relative-velocity RELATIVE_VELOCITY [RELATIVE_VELOCITY ...]
                        Relative velocity of the second snapshot with relation to the first. Must be given in terms
                        of it's components in x, y and z
                        
                        
  -r ROTATION [ROTATION ...], --rotation ROTATION [ROTATION ...]
                        Angle of the rotation to be applied to the second snapshot. Must be given in terms of
                        rotations around the x, y and z axis that pass by the origin of the second snapshot
                        
                        
  --hdf5                Output initial conditions in HDF5
  ```
