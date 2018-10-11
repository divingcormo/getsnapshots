==================================================================
Gather data from chemsh outputs and write csv file for R-modelling
==================================================================

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
