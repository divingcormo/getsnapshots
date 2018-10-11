get-snapshots
=============

   **Usage:** get-snapshots  [-h] [--csv ``CSV``] [--simulate] [--verbose]

   Options:
     --help, -h         Show this help message and exit
     --csv, -c          ``CSV`` can be one of [a|s0|old|predict|model], for S1/S0/pre-v36 snapshots or for prediction/modeling
     --simulate, -s     Dry run, don't write csv file.
     --verbose, -v      Verbose output for debugging.

This script produces "snapshots-a.csv"-like files. These files store
the following data:

* ``snapshots-a.csv`` contains qm14-fix-v36/v38 S1 snapshots.
* ``snapshots-s0.csv`` contains qm12 and qm14-v38 S0 snapshots.
* ``snapshots-old.csv`` contains snapshots from older potentials.
* ``snapshots-for-modeling.csv`` contains geometry data and
    absorption energies from from om2-opt snapshots, mostly from
    v38 trajectories.
* ``snapshots-for-prediction.csv`` contains geometries and absorption
    energies from om2-opt snapshots (v76 trajectories) and cro-opt
    (S1 OM2/MRCIS optimized) emission energies for selected
    snapshots.

These files are then used for data analysis with R or pandas.
