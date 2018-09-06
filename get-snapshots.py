#!/home/wanko/bin/python3
# -*- coding: utf-8 -*-
"""Gather data from chemsh outputs and write csv file for R-modelling."""

import os
import glob
import argparse
import asyncio
import math
from collections import defaultdict, namedtuple
import shutil
import subprocess

import numpy as np
import tqdm

from pyrmsd import read_c

version = 1.1

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

MIN_MODELING = tuple("v38e v38e1 v38e3 v38e5 v38e7 \
            a_fix-met105-lys-ile-glu-ser-hb v36b v36c v38s0 v38a \
            v38b v38c v38d \
            v38f v38g v38h v39d v38i v38j v38k v38l v38m v38n \
            v38o v38p v38q v38r v38s v38t v38u v38v v38w v38x v38y".split())

MIN_PREDICT = tuple("v76s0 v76a v76b v76c v76d v76e v76h v76i v76j v76k v76m \
            v76n v76o".split())

FIRST_CYCLE, LAST_CYCLE = 2, 92
FIRST_FRAME, LAST_FRAME, SKIP_FRAME = 10, 200, 10
CHECK_CONVERGED = True
GET_SNAPSHOTS_SCRIPT = 'get-snapshots.csv.bash'


Variable = namedtuple('Variable', 'var_type ind dihe_shift ref_coords'.split())


def define_variables():
    DEBUG = False
    variables = {}
    with open(GET_SNAPSHOTS_SCRIPT, 'r') as fo:
        lines = fo.readlines()
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
    if not isinstance(v, Variable):
        raise TypeError("{} is not an instance of Variable".format(v))
    if v.var_type == 'hb':
        d1 = get_distance(v.ind[0], v.ind[1], atoms)
        d2 = get_distance(v.ind[2], v.ind[3], atoms)
        return "{:.3f}".format(min(d1, d2))
    if v.var_type == 'dist':
        return "{:.3f}".format(get_distance(v.ind[0], v.ind[1], atoms))
    if v.var_type == 'dihe':
        dihe = get_dihedral(v.ind[0], v.ind[1], v.ind[2], v.ind[3], atoms)
        if dihe < 0:
            dihe += v.dihe_shift
        return "{:.2f}".format(dihe)
    if v.var_type == 'rmsd':
        # v.ind has type int not list here
        r1 = [atoms[v.ind - 1][i] for i in range(3)]
        r2 = v.ref_coords
        rmsd = math.sqrt(sum((a - b)**2 for a, b in zip(r1, r2)))
        return "{:.5f}".format(rmsd)
    raise ValueError("var_tpye {} not supported".format(v.var_type))


def get_dihedral(i1, i2, i3, i4, atoms):
    # i1..i4 start from 1, atoms indices from 0.
    xyz = slice(0, 3)
    x1 = np.array(atoms[i1-1][xyz])
    x2 = np.array(atoms[i2-1][xyz])
    x3 = np.array(atoms[i3-1][xyz])
    x4 = np.array(atoms[i4-1][xyz])
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
    else:
        norm = math.sqrt(uvxvw.dot(uvxvw))
        return sign * 180 / math.pi * math.atan2(norm, uv_vw)


def get_distance(i1, i2, atoms):
    xyz = slice(0, 3)
    x1 = np.array(atoms[i1-1][xyz])
    x2 = np.array(atoms[i2-1][xyz])
    d = x2 - x1
    return math.sqrt(d.dot(d))


def converged(log):
    with open(log, 'r') as fo:
        return ('Converged' in fo.read())


def get_omega(mndo_out, chemsh_log):
    DEBUG = False
    try:
        with open(mndo_out, 'r') as fo:
            lines = fo.read().splitlines()
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
    if n_match == 0:
        if 'om2-opt' in mndo_out:
            mndo_in = mndo_out.replace('mndo.out', 'mndo.in')
            assert mndo_in != mndo_out
            with open(mndo_in, 'r') as fo:
                fs = fo.read()
            if not 'iuvcd=2 ipop=2' in fs:
                msg = "{} does not contain excited-state emission energy"\
                    + " - trying to calculate it..."
                print(msg.format(mndo_in), end='')
                if not os.path.exists(mndo_in + '.bak'):
                    shutil.copy2(mndo_in, mndo_in + '.bak')
                fs = fs.replace('ktrial=11', 'iuvcd=2 ipop=2 ktrial=11')
                fs = fs.replace('numatm=4035', 'numatm=4035 +\nkci=5 '\
                              + 'ici1=42 ici2=43 ioutci=1 movo=0 multci=1 '\
                              + 'nciref=3 +\nmciref=0 levexc=1 iroot=2 lroot=2')
                fs = fs.replace('jop=-2', 'jop=-1')
                with open(mndo_in, 'w') as fo:
                    fo.write(fs)
                mndo_call = 'mndo99 < ' + mndo_in + ' > ' + mndo_out
                assert CONCUR_TASKS == 1
                subprocess.call(mndo_call, shell=True)  # must block
                subprocess.call('rm fort.*', shell=True)
                with open(mndo_out, 'r') as fo:
                    lines = fo.read().splitlines()
                n_match = 0
                for line in lines:
                    if 'State  2' in line:
                        Omega = line.split()[5]
                        n_match += 1
                        Omega = float(Omega)
                if n_match == 1:
                    print(' SUCCESS')
                else:
                    print(' FAILURE')
                    shutil.copy2(mndo_in + '.bak', mndo_in)
                    return 'NA'
                    exit()
            else:
                print("Warning: found iuvcd=2 in", mndo_in, "but no excitation"\
                    + " energy in", mndo_out)
                return 'NA'
        else:
            print("Warning: No excitation energy found in", mndo_out)
            return 'NA'
    elif n_match > 1:
        print("Warning: 'State  2' found several times in", mndo_out)
    return Omega


async def read_one(Min, Cycle, Frame, variables, mode, semaphore, verbose):
    DEBUG = False
    NO_DATA = False         # return value if files are missing
    loop = asyncio.get_event_loop()

    if mode == 'analyse':
        opt = 'cro-opt'
    else:
        opt = 'om2-opt'     # evaluate variables from om2-opt
    i = str(Cycle) + '-' + str(Frame)
    pth = 'snapshots-' + str(Min)
    log = pth + '/' + opt + '-' + i + '.log'
    sge = pth + '/' + opt + '-' + i + '.sge'
    mndo_out = pth + '/cro-opt-' + i + '-mndo.out'  # for Omega
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
            Abs = await loop.run_in_executor(None, get_omega, abs_out, log)

    async with semaphore:
        if DEBUG:
            print("   ", Min, Cycle, Frame, "get_omega...")
        Omega = await loop.run_in_executor(None, get_omega, mndo_out, log)
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
    snapshot['Min'] = Min
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
        if 's0' in Min or 'qm12' in Min or 'qm15' in Min:
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
    "Return a list of defaultdict's with keys from VARIABLES, default: 'NA'."
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
        # Search for existing log files. Cycle and Frame are str.
        logfiles = glob.glob(pth + '/' + opt + '-*.log')
        cycles = set(os.path.basename(f).split('-')[2] for f in logfiles)
        n_cycles = len(cycles)
        n_frames = 0
        for Cycle in cycles:
            # Search for existing log files.
            logfiles = glob.glob(pth + '/' + opt + '-' + Cycle + '-*.log')
            frames = set(os.path.basename(f).replace('.', '-').split('-')[3] for f in logfiles)
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
    data = loop.run_until_complete(common_future)    # returns list of return values
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
    loop = asyncio.get_event_loop()
    done_iter = asyncio.as_completed(tasks, loop=loop)
    if not verbose:
        done_iter = tqdm.tqdm(done_iter, total=len(tasks.values()))
    for task in done_iter:
        await task
    return 'monitor closed'


def test():
    assert len(VARIABLES) == len(set(VARIABLES))
    Min = 'v38s'
    Cycle = 13
    Frame = 50
    assert tuple('v38s'.split())[0] == Min
    opt = 'cro-opt'
    i = str(Cycle) + '-' + str(Frame)
    pth = 'snapshots-' + str(Min)
    log = pth + '/' + opt + '-' + i + '.log'
    assert converged(log) is True
    optc = pth + '/' + opt + '-' + i + '-opt.c'
    na, atoms, _ = read_c(optc)
    variables = define_variables()
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    verbose = True
    mode = 'analyse'
    snapshot = loop.run_until_complete(read_one(
            Min, Cycle, Frame, variables, mode, semaphore, verbose))
    phi = get_value(variables['phi'], atoms)
    assert phi == "12.07"
    phi = snapshot['phi']
    assert phi == "12.07"
    hbW240E215 = get_value(variables['hbW240E215'], atoms)
    assert hbW240E215 == "1.760"
    hbW240E215 = snapshot['hbW240E215']
    assert hbW240E215 == "1.760"
    W304rmsd = get_value(variables['W304rmsd'], atoms)
    assert W304rmsd == "1.39187"
    W304rmsd = snapshot['W304rmsd']
    assert W304rmsd == "1.39187"
    L199rot = get_value(variables['L199rot'], atoms)
    assert L199rot == "-66.47"
    L199rot = snapshot['L199rot']
    assert L199rot == "-66.47"
    bb66dihe = get_value(variables['bb66dihe'], atoms)
    assert bb66dihe == "100.74"
    bb66dihe = snapshot['bb66dihe']
    assert bb66dihe == "100.74"
    E16cAcB = get_value(variables['E16cAcB'], atoms)
    assert E16cAcB == "190.78"
    E16cAcB = snapshot['E16cAcB']
    assert E16cAcB == "190.78"
    Omega = snapshot['Omega']
    assert Omega == "1.8418"
    # print(snapshot)
    return "Tests pass.\n"


def main():
    print("\n    --=== get-snapshots, version", version, "===--")

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Produce "snapshots-a.csv"-like file.\n'
                    'snapshots-a.csv contains qm14-fix-v36/v38 S1 snapshots\n'
                    'snapshots-s0.csv contains qm12 and qm14-v38 S0 snapshots\n'
                    'snapshots-old.csv contains snapshots from older potentials')
    parser.add_argument("--csv", "-c",
                                help='[a|s0|old|predict|model] for \
                                S1/S0/pre-v36 snapshots or for prediction/modeling',
                                default='a')
    parser.add_argument("--simulate", "-s",
                                help="Dry run, don't write csv file.",
                                action="store_true")
    parser.add_argument("--test", "-t",
                                help="Perform self test.",
                                action="store_true")
    parser.add_argument("--verbose", "-v",
                                help="verbose output for debugging.",
                                action="store_true")

    args = parser.parse_args()

    mode = 'analyse'
    global VARIABLES    # need to add some in modes "predict" and "model"
    if args.test:
        csv = 'snapshots-test.csv'
        Mins = tuple('v38s'.split())
    elif args.csv == 'a':
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

    if args.test:
        print("Selftest...")
        print(test())

    if args.simulate:
        print("just a simulation...")
        csv = '/dev/null'

    verbose = args.verbose

    if not os.path.exists(GET_SNAPSHOTS_SCRIPT):
        print("get-snapshots.csv.bash script not found")
        exit(1)

    if verbose:
        nvar = len(VARIABLES)
        print(nvar, "variables:", VARIABLES)

    # Get snapshot data from files
    data = read_all(Mins, mode, verbose)

    # Write csv file
    csv_header = ','.join(VARIABLES) + '\n'
    with open(csv, 'w') as fo:
        fo.write(csv_header)
        for snapshot in data:
            line = ','.join(str(snapshot[x]) for x in VARIABLES)
            fo.write(line + '\n')

    # Final Report
    msg = 'Wrote {} with {} entries, {} variables.'
    print(msg.format(csv, len(data), len(VARIABLES)))


if __name__ == '__main__':
    # Selftest
    if False:
        print("Selftest...")
        print(test())
        loop = asyncio.get_event_loop()
        loop.close()
    else:
        main()
