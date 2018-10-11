#!/home/wanko/bin/python3
# -*- coding: utf-8 -*-

"""Reads data from chemsh calculations and writes a csv file.

Structures are stored in atoms objects, which are documented in the pyrmsd
module."""

#  This file is part of the getsnapshots package.
#
#  Copyright 2018 Marius Wanko <marius.wanko@gmail.com>
#
#  getsnapshots is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
import glob
import argparse
import asyncio
import math
from collections import defaultdict, namedtuple
import shutil
import subprocess

import numpy as np
import tqdm

from getsnapshots.pyrmsd import read_c

__version__ = '1.2'

CONCUR_TASKS = 1

VARIABLES = tuple('Min Cycle Frame Omega M66cAcB M66cBcG M66cGsD \
            alpha alphap beta gamma phi tau \
            E16cAcB E16rot E16orie L199rot K70cDcE K70cGcD K70cBcG \
            I65cAcB I65rot Q109rot S111rot S62rot L119rot \
            hbE16acy hbI65acy bE16oh I197rot Q213rot \
            hbE215imi hbW239o64 hbW239o66 hbW239E16 hbE16W239 hbW240imi \
            hbW240E215 hbK70E148 hbK70imi hbW304imi hbW239acy hbQ213acy \
            hbK70W301 dW304Trp93 dW240Q42 M18cBcG M18cGsD hbW230Q109 \
            dW301Q109 hbW304K70 o66orie bb66dihe hbW260phe hbS146phe hbR95imi \
            hbW240n2 hbW230S111 hbW230o66 hbW230E16 hbW230o64 hbW230W239 \
            hbW301W230 hbW301W304 hbW301S111 hbW301o66 hbW304E148 \
            W301rmsd W304rmsd'.split())    # ordered, immutable

MIN_S1 = tuple("v36b v36c v38a v38b v38c v38d v38f v38g v38h v39d \
            v38e v38e1 v38e3 v38e5 v38e7\
            v38i v38j v38k v38l v38m v38n v38o v38p v38q v38r \
            v38s v38t v38u v38v v38w v38x v38y \
            v76a v76b v76c v76d v76e v76h v76i v76j v76k v76m v76n v76o".split())

MIN_S0 = tuple("a_qm12 v38s0 v76s0 a_s0-charges-fix-met105-lys-ile-glu-ser \
            a_s0-charges-fix-met105-lys-ile-glu-ser-hb".split())

MIN_OLD = tuple("a a_fix-lys-ile a_fix-lys-ile-glu a_fix-met105 \
            a_fix-met105-glu a_fix-met105-lys a_fix-met105-lys-glu \
            a_fix-met105-lys-ile a_fix-met105-lys-ile-glu \
            a_fix-met105-lys-ile-glu-ser a_fix-lys-ile-glu-ser \
            a_fix-met105-lys-ile-glu-ser-hb a_fix-lys-ile-glu-ser-hb".split())

"""
For GS, calculate both GS abs (Min=="*abs") and cro-opt emission (Min=="*s0").
The only difference is that Omega is obtained from
    ../6b-S0_om2-mrcis/snapshots-v*s0/cro-opt-*-mndo.out
in the first case and from
    ./snapshots-v*s0/cro-opt-*-mndo.out
in the latter case. The "cro-opt emission" is not trivial to interpret as only
the cro is actually relaxed in the excited state.
"""
MIN_MODELING = tuple("v38e v38e1 v38e3 v38e5 v38e7 \
            a_fix-met105-lys-ile-glu-ser-hb v36b v36c v38s0 v38abs v38a \
            v38b v38c v38d \
            v38f v38g v38h v39d v38i v38j v38k v38l v38m v38n \
            v38o v38p v38q v38r v38s v38t v38u v38v v38w v38x v38y".split())

MIN_PREDICT = tuple("v76abs v76s0 v76a v76b v76c v76d v76e v76h v76i v76j \
             v76k v76m v76n v76o".split())

FIRST_CYCLE, LAST_CYCLE = 2, 92
FIRST_FRAME, LAST_FRAME, SKIP_FRAME = 10, 200, 10
CHECK_CONVERGED = True
GET_SNAPSHOTS_SCRIPT = 'get-snapshots.csv.bash'


Variable = namedtuple('Variable', 'var_type ind dihe_shift ref_coords'.split())


def define_variables():
    """Reads geometric variable definitions from file GET_SNAPSHOTS_SCRIPT."""
    DEBUG = False
    variables = {}
    with open(GET_SNAPSHOTS_SCRIPT, 'r') as fp:
        lines = fp.readlines()
    var_type = None
    for NR, line in enumerate(lines):
        if line.startswith('d1='):
            var_type = 'hb'
            f = line.replace(';', ' ').split()
            indices = [int(s) for s in (f[1], f[2], f[5], f[6])]
            continue
        if line.startswith('var[') and not line.startswith('var[$'):
            line = line.replace('var[', '').replace(']=', ' ').replace('`', '')
            f = line.split()
            name = f[0]
            if name in 'Min Cycle Frame Omega'.split():
                continue
            dihe_shift, ref_coords = 0, [0, 0, 0]
            if var_type == 'hb':
                pass
            elif 'getdihedral' in f[1]:
                var_type = 'dihe'
                indices = [int(s) for s in (f[2], f[3], f[4], f[5])]
                if '360+$1' in line:
                    dihe_shift = 360
            elif 'getdistance' in f[1]:
                var_type = 'dist'
                indices = [int(s) for s in (f[2], f[3])]
            elif 'rmsd' in name:
                var_type = 'rmsd'
                f = line.replace('{', ' ').replace('=', ' ').split()
                indices = int(f[3])
                ref_coords = [float(s) for s in (f[7], f[9], f[11])]
            else:
                print("ERROR processing", GET_SNAPSHOTS_SCRIPT, ":")
                print("  no variable type reckognized in line", NR+1)
                print(line)
                exit(1)
            # Instantiate and add variable to dict
            if DEBUG:
                print("variable", name, "(", var_type, ")", indices,
                      dihe_shift, ref_coords)
            variables[name] = Variable(var_type, indices, dihe_shift, ref_coords)
            if var_type == 'hb':
                var_type = None     # reset for next line
    msg = 'Read {} variable definitions from "{}".'
    print(msg.format(len(variables), GET_SNAPSHOTS_SCRIPT))
    return variables


def get_value(v, atoms):
    """Reads a Variable's definition and obtains its value.

    Args:
        atoms: object that contains a list of atomic coordinates
        v: instance of Variable

    A variable is defined by the folling attributes:
        var_type:   One of 'hb', 'dist', 'rmsd', or 'dihedral'.
        ind:        A list of indices
        dihe_shift: 0 or 360, added to negative values. This is useful if
                    dihedrals around 180/-180 are populated.
        ref_coords: The reference point for calculating an RMSD.

    Returns:
        float: The variable's value.
    """
    if not isinstance(v, Variable):
        raise TypeError("{} is not an instance of Variable".format(v))
    if v.var_type == 'hb':
        dist_1 = get_distance(v.ind[0], v.ind[1], atoms)
        dist_2 = get_distance(v.ind[2], v.ind[3], atoms)
        return "{:.3f}".format(min(dist_1, dist_2))
    if v.var_type == 'dist':
        return "{:.3f}".format(get_distance(v.ind[0], v.ind[1], atoms))
    if v.var_type == 'dihe':
        dihe = get_dihedral(v.ind[0], v.ind[1], v.ind[2], v.ind[3], atoms)
        if dihe < 0:
            dihe += v.dihe_shift
        return "{:.2f}".format(dihe)
    if v.var_type == 'rmsd':
        # v.ind has type int not list here
        r_1 = [atoms[v.ind - 1][i] for i in range(3)]
        r_2 = v.ref_coords
        rmsd = math.sqrt(sum((a - b)**2 for a, b in zip(r_1, r_2)))
        return "{:.5f}".format(rmsd)
    raise ValueError("var_tpye {} not supported".format(v.var_type))


def get_dihedral(index1, index2, index3, index4, atoms):
    """Calculates the dihedral angle between four atoms.

    Args:
        index1...4 (int): atom indices (first atom: 1)
        atoms: list of atom coordinates (first atom: atoms[0])

    Returns:
        float: The dihedral angle in degrees.
    """
    xyz = slice(0, 3)
    x1 = np.array(atoms[index1-1][xyz])
    x2 = np.array(atoms[index2-1][xyz])
    x3 = np.array(atoms[index3-1][xyz])
    x4 = np.array(atoms[index4-1][xyz])
    u = x1 - x2
    v = x3 - x2
    w = x3 - x4
    uxv = np.cross(u, v)
    vxw = np.cross(v, w)
    uv_vw = float(uxv.dot(vxw))
    uvxvw = np.cross(uxv, vxw)
    uvxvw_v = float(uvxvw.dot(v))
    if uvxvw_v >= 0:
        sign = 1
    else:
        sign = -1
    if abs(uv_vw) < 1e-6:
        return 90.0 * sign
    norm = math.sqrt(uvxvw.dot(uvxvw))
    return sign * 180 / math.pi * math.atan2(norm, uv_vw)


def get_distance(index1, index2, atoms):
    """Returns an interatomic distance.

    Args:
        index1...2 atom indices (first atom: 1).
        atoms: list of atom coordinates (first atom: atoms[0]).
    """
    xyz = slice(0, 3)
    x1 = np.array(atoms[index1 - 1][xyz])
    x2 = np.array(atoms[index2 - 1][xyz])
    d = x2 - x1
    return math.sqrt(d.dot(d))


def converged(log):
    """Returns True if log is a chemsh log file with converged optimization."""
    with open(log, 'r') as fp:
        return 'Converged!' in fp.read()


def get_omega(mndo_out):
    """Extracts S0->S1 excitation energy from mndo (gugaci) output file.

    Args:
        mndo_out: The Filename of an mndo output.

    Returns:
        A float, if mndo_out is found, with the excitation energy in eV of the
        lowest electronic transition of a GUGACI calculation.
        Otherwise, if the input file is found, and mndo_out matches "om2-opt"
        or "cro-opt", tries to calculate it.
        Makes a backup of the input file (if non-existent), modifies it by
        adding GUGACI keywords for a single-point 2-root calculation.
        'NA' if mndo_out does not exist or the excitation energy could not be
        obtained.

    Raises:
        ValueError: If "State  2" appears more than once.
        NotImplementedError: If the reason for the missing excitation energy
            could not be determined.
    """
    DEBUG = False
    try:
        with open(mndo_out, 'r') as fp:
            lines = fp.read().splitlines()
    except FileNotFoundError:
        return 'NA'

    n_match = 0
    for line in lines:
        if 'State  2' in line:
            Omega = line.split()[5]
            n_match += 1
            if DEBUG:
                print("excitation energy from", mndo_out, ":", Omega)
            Omega = float(Omega)
    if n_match == 1:
        return Omega
    if n_match > 1:
        msg = "Warning: 'State  2' found several times in {}"
        raise ValueError(msg.format(mndo_out))
    assert n_match == 0
    if 'om2-opt' not in mndo_out and 'cro-opt' not in mndo_out:
        print("Warning: No excitation energy found in", mndo_out)
        return 'NA'

    # In om2-opt-*-mndo.out and cro-opt-*-mndo.out the
    # S1 excitation energy may be missing in the following cases:
    # - in mode predict if snapshot is not part of cro-opt sample,
    #   then calculate it for abs if om2-opt.
    # - in mode predict if snapshot is part of cro-opt sample,
    #   then calculate it as well.
    mndo_in = mndo_out.replace('mndo.out', 'mndo.in')
    assert mndo_in != mndo_out
    with open(mndo_in, 'r') as fp:
        fstr = fp.read()
    if 'iuvcd=2 ipop=2' in fstr:
        print("Warning: found iuvcd=2 and iroot=2 in", mndo_in,
              "but no excitation energy in", mndo_out)
        return 'NA'
    msg = "{} does not contain excited-state emission energy"\
        + " - trying to calculate it..."
    print(msg.format(mndo_out), end='')
    sys.stdout.flush
    if not os.path.exists(mndo_in + '.bak'):
        shutil.copy2(mndo_in, mndo_in + '.bak')

    if 'iroot=1 lroot=1' in fstr:
        # This is a GS cro-opt, modify mndo.in for excited-state calc.
        fstr = fstr.replace('iroot=1 lroot=1', 'iroot=2 lroot=2')
    elif 'iuvcd=' not in fstr:
        # This is an om2-opt, modify mndo.in for excited-state calc.
        fstr = fstr.replace('ktrial=11', 'iuvcd=2 ipop=2 ktrial=11')
        fstr = fstr.replace('numatm=4035', 'numatm=4035 +\nkci=5 '
                        + 'ici1=42 ici2=43 ioutci=1 movo=0 multci=1 '
                        + 'nciref=3 +\nmciref=0 levexc=1 '
                        + 'iroot=2 lroot=2')
    else:
        shutil.copy2(mndo_in + '.bak', mndo_in)
        msg = "Analysis of {} failed, no plan to handle this case."
        raise NotImplementedError(msg.format(mndo_in))

    fstr = fstr.replace('jop=-2', 'jop=-1')
    with open(mndo_in, 'w') as fp:
        fp.write(fstr)
    mndo_call = 'mndo99 < ' + mndo_in + ' > ' + mndo_out
    assert CONCUR_TASKS == 1
    subprocess.call(mndo_call, shell=True)  # must block
    subprocess.call('rm fort.*', shell=True)
    with open(mndo_out, 'r') as fp:
        lines = fp.read().splitlines()
    n_match = 0
    for line in lines:
        if 'State  2' in line:
            Omega = line.split()[5]
            n_match += 1
            Omega = float(Omega)
    if n_match != 1:
        print(' FAILURE')
        shutil.copy2(mndo_in + '.bak', mndo_in)
        return 'NA'
        exit()  # Or stop and investigate the problem?
    print(' SUCCESS')
    return Omega


async def read_one(Min, Cycle, Frame, variables, mode, semaphore, verbose):
    """Returns values for all variables of one MD snapshot.

    Args:
        Min, Cycle, Frame: These strings are used to determine the filenames
            from a chemsh calculation, but are also part of the return.
        variables (dict): Maps variable names (key) to Variable instances.
        mode (str): One of ('analyse', 'model', 'predict').
        semaphore: A Semaphore instance for the asyncio scheduler.
        verbose (boole): If True, prints text instead of a progress bar.

    Example:
        Min = "v76n" (charmm trajectory),
        Cycle = "13" (MD cycle),
        Frame = "200" (frame number from the charmm dcd file).
        mode = "model".

        The following chemsh files are expected in ./snapshots-v76n:
        om2-opt-13-200.sge, om2-opt-13-200.log, om2-opt-13-200-opt.c,
        om2-opt-13-200-mndo.in, om2-opt-13-200-mndo.out

    Raises:
        ValueError: If the excitation energy (Omega) could not be obtained in
        mode "analyse" or "model".
        ValueError: If the emission energy (Abs) could not be obtained in mode
        'model' or 'predict'.

    Returns:
        A dict(str) of variable names and their values (str). It contains all
        variables from the global VARIABLES.
        NO_DATA (False) if essential data is missing.

    """
    DEBUG = False
    NO_DATA = False         # return value if files are missing
    loop = asyncio.get_event_loop()

    if mode == 'analyse':
        opt = 'cro-opt'
    else:
        opt = 'om2-opt'     # evaluate variables from om2-opt
    i = str(Cycle) + '-' + str(Frame)
    pth = 'snapshots-' + str(Min)
    if 'abs' in str(Min):
        pth = 'snapshots-' + str(Min).replace('abs', 's0')
    log = pth + '/' + opt + '-' + i + '.log'
    sge = pth + '/' + opt + '-' + i + '.sge'
    mndo_out = pth + '/cro-opt-' + i + '-mndo.out'  # for Omega
    if 'abs' in str(Min):
        mndo_out = '../6b-S0_om2-mrcis/' + pth + '/cro-opt-' + i + '-mndo.out'
    abs_out = pth + '/om2-opt-' + i + '-mndo.out'   # for Abs
    optc = pth + '/' + opt + '-' + i + '-opt.c'

    # Check file stati
    async with semaphore:
        if DEBUG:
            print("   ", Min, Cycle, Frame, "checking log...")
        log_exists = await loop.run_in_executor(None, os.path.exists, log)
    if not log_exists:
        print("    Warning: found", sge, "but not", log, ".")
        return NO_DATA

    async with semaphore:
        if DEBUG:
            print("   ", Min, Cycle, Frame, "checking converged...")
        log_converged = await loop.run_in_executor(None, converged, log)
    if not log_converged:
        print("    Warning:", log, "not converged.")
        return NO_DATA

    async with semaphore:
        if DEBUG:
            print("   ", Min, Cycle, Frame, "checking optc...")
        optc_exists = await loop.run_in_executor(None, os.path.exists, optc)
    if not optc_exists:
        print("Warning:", optc, "missing although", log, "signals convergence!")
        return NO_DATA

    if mode in ('model', 'predict'):
        async with semaphore:
            Abs = await loop.run_in_executor(None, get_omega, abs_out)

    async with semaphore:
        if DEBUG:
            print("   ", Min, Cycle, Frame, "get_omega...")
        Omega = await loop.run_in_executor(None, get_omega, mndo_out)
    if Omega == 'NA' and not mode == 'predict':
        msg = 'Could not get Omega for {} {}-{} (required in mode {}).'
        raise ValueError(msg.format(Min, Cycle, Frame, mode))

    async with semaphore:
        if DEBUG or verbose:
            print("   ", Min, Cycle, Frame, "reading coordinates...")
        na, atoms, _ = await loop.run_in_executor(None, read_c, optc)

    # Now we have all raw data in memory. Time to evaluate the variable.
    # All values in snapshot are expected to be of type str.
    snapshot = defaultdict(lambda: 'NA')
    snapshot['Min'] = str(Min)
    snapshot['Cycle'] = str(Cycle)
    snapshot['Frame'] = str(Frame)
    if Omega == 'NA':
        snapshot['Omega'] = 'NA'
    else:
        snapshot['Omega'] = "{:.4f}".format(Omega)
    if mode in ('model', 'predict'):
        if Abs == 'NA':
            msg = "Missing Abs in {} {}-{} (needed in mode {})."
            raise ValueError(msg.format(Min, Cycle, Frame, mode))
        snapshot['Abs'] = "{:.4f}".format(Abs)
        if any(x in Min for x in 's0 qm12 qm15 abs'.split()):
            State = 'S0'
        else:
            State = 'S1'
        snapshot['State'] = State

    for var in VARIABLES:
        if var in ('Min', 'Cycle', 'Frame', 'Omega', 'Abs', 'State'):
            continue
        snapshot[var] = get_value(variables[var], atoms)

    return snapshot


def read_all(Mins, mode, verbose):
    """Asyncronuously retrieve all data required for the csv file.

    Args:
        Mins: An iterable of strings with MD trajectories to be processed.
        mode (str): One of ('analyse', 'model', 'predict').
        verbose (boole): If True, prints text instead of a progress bar.

    Returns:
        A list(defaultdict) of snapshots. The dict maps all variable names
        (str) in global VARIABLES to their values (str). The default value is
        'NA'.

    If global flag CONCUR_TASKS > 1, each snapshot forms a task in an asyncio
    event loop. This, however, prevents the posterior calculation of missing
    excitation energies. It is possible that no speed is gained if the OS
    limits concurrent file IO.
    """
    DEBUG = False
    variables = define_variables()
    tasks = {}
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    if mode in ('analyse', 'model'):
        opt = 'cro-opt'     # For "analyse" and "model", Omega is required,
    else:                   # consider only snapshots where cro-opt exists.
        opt = 'om2-opt'     # For "predict", all om2-opt are used (Omega
                            # exists only for a small sample).
    print("Preparing tasks...")
    for Min in Mins:
        pth = 'snapshots-' + str(Min)
        if 'abs' in str(Min):
            pth = 'snapshots-' + str(Min).replace('abs', 's0')
        # Search for existing log files. Cycle and Frame are str.
        logfiles = glob.glob(pth + '/' + opt + '-*.log')
        cycles = set(os.path.basename(f).split('-')[2] for f in logfiles)
        n_cycles = len(cycles)
        n_frames = 0
        for Cycle in cycles:
            # Search for existing log files.
            logfiles = glob.glob(pth + '/' + opt + '-' + Cycle + '-*.log')
            frames = set(os.path.basename(f).replace('.', '-').split('-')[3]
                         for f in logfiles)
            n_frames += len(frames)
            for Frame in frames:
                coro_obj = read_one(Min, Cycle, Frame, variables, mode,
                                    semaphore, verbose)
                task_name = "Min {} Cycle {} Frame {}".format(Min, Cycle, Frame)
                tasks[asyncio.ensure_future(coro_obj)] = task_name
                if DEBUG:
                    print("    scheduled task:", task_name)
                if DEBUG:
                    print("   ", task_name)
        if verbose:
            msg = "    {}: {} snapshots from {} cycles."
            print(msg.format(Min, n_frames, n_cycles))

    # Setup the monitor task
    monitor_obj = monitor(tasks, verbose)
    monitor_task = asyncio.ensure_future(monitor_obj)
    if DEBUG:
        print("scheduled monitor.")

    # Make common future and wait for completion
    common_future = asyncio.gather(*list(tasks.keys()), monitor_task)
    print("Searching for data...")
    data = loop.run_until_complete(common_future)
    print("data complete.")
    loop.close()

    # Clean and analyse data
    for i, snapshot in enumerate(data):
        if not isinstance(snapshot, dict):
            if snapshot == 'monitor closed':
                pass
            elif snapshot is False:
                pass
            else:
                msg = 'unexpected return from read_one: {!r}'.format(snapshot)
                raise TypeError(msg)
            del data[i]
    if DEBUG:
        for i, snapshot in enumerate(data):
            # print(i, snapshot, type(snapshot))
            line = ','.join([str(snapshot[v]) for v in VARIABLES])
            print(i, line)
    return data


async def monitor(tasks, verbose):
    "Produces a tqdm graphical progress bar for snapshot evaluation."
    loop = asyncio.get_event_loop()
    done_iter = asyncio.as_completed(tasks, loop=loop)
    if not verbose:
        done_iter = tqdm.tqdm(done_iter, total=len(tasks.values()))
    for task in done_iter:
        await task
    return 'monitor closed'


def main():
    """Parses CL arguments, calls read_all() and writes the csv file."""
    print("\n    --=== get-snapshots, version", __version__, "===--")

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            'Produce "snapshots-a.csv"-like file.\n'
            'snapshots-a.csv contains qm14-fix-v36/v38 S1 snapshots\n'
            'snapshots-s0.csv contains qm12 and qm14-v38 S0 snapshots\n'
            'snapshots-old.csv contains snapshots from older potentials'
        )
    )
    arg = parser.add_argument
    arg("--csv", "-c", default='a', help=(
                        '[a|s0|old|predict|model] for '
                        'S1/S0/pre-v36 snapshots or for prediction/modeling'))
    arg("--simulate", "-s", action="store_true", help=(
                        "Dry run, don't write csv file."))
    arg("--verbose", "-v", action="store_true", help=(
                        "verbose output for debugging."))

    args = parser.parse_args()

    mode = 'analyse'
    global VARIABLES    # modified in modes "predict" and "model"
    if args.csv == 'a':
        csv = 'snapshots-a.csv'
        Mins = MIN_S1
    elif args.csv == 's0':
        csv = 'snapshots-s0.csv'
        Mins = MIN_S0
    elif args.csv == 'old':
        csv = 'snapshots-old.csv'
        Mins = MIN_OLD
    elif args.csv == 'model':
        csv = 'snapshots-for-modeling.csv'
        Mins = MIN_MODELING
        mode = 'model'
        VARIABLES += ('State', 'Abs')
    elif args.csv == 'predict':
        csv = 'snapshots-for-prediction.csv'
        Mins = MIN_PREDICT
        mode = 'predict'
        VARIABLES += ('State', 'Abs')
    else:
        msg = "No rule for making csv file {}".format(args.csv)
        raise ValueError(msg)

    if args.simulate:
        print("just a simulation...")
        csv = '/dev/null'

    verbose = args.verbose

    assert os.path.exists(GET_SNAPSHOTS_SCRIPT)

    if verbose:
        nvar = len(VARIABLES)
        print(nvar, "variables:", VARIABLES)

    # Get snapshot data from files
    data = read_all(Mins, mode, verbose)

    # Write csv file
    csv_header = ','.join(VARIABLES) + '\n'
    with open(csv, 'w') as fp:
        fp.write(csv_header)
        for snapshot in data:
            line = ','.join(str(snapshot[x]) for x in VARIABLES)
            fp.write(line + '\n')

    # Final Report
    msg = 'Wrote {} with {} entries, {} variables.'
    print(msg.format(csv, len(data), len(VARIABLES)))


if __name__ == '__main__':
    main()
