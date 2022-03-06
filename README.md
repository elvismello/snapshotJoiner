## snapshotJoiner
Self contained program written in python 3 that can join snapshots using hdf5 format.

Based in the module "snapwrite" contained in [Rafael Ruggiero](https://ruggiero.github.io/)'s IC generation programs, but modified to work in python 3.


## Usage Examples

```
python3 snapshotJoiner.py cluster.hdf5 galaxy.hdf5 -o rampressure.hdf5 --relative-position 1000 100 0 --relative-velocity -100 0 0  --rotation 0 90 0
```

```
python3 snapshotJoiner.py clusterA.hdf5 clusterB.hdf5 -o rampressure.hdf5 --relative-position 1000 100 0
```

```
python3 snapshotJoiner.py galaxyA.hdf5 galaxyB.hdf5 -o rampressure.hdf5 --relative-position 1000 0 0 --rotation 0 90 0
```
