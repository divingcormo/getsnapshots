import pytest
import os
import sys

# If tests contains __init__.py (=> (sub)package) this should work but does not:
#from tests import test_pyrmsd

# Somehow Python needs to know where the package modules are.
# Use $PYTHONPATH (there must be a better way!!!)

# Need to have get-snapshots.csv.bash in cwd (hardcoded).
# All other data files are also in tests/ and subdirs -> cd tests
# No matter where os.chdir() appears in the code, it is executed at module time,
# i.e., before any test functions are called.
if not os.getcwd().endswith('/tests'):
    print("test module called from", os.getcwd(), "- changing to tests/")
    os.chdir('tests')
    #os.chdir('getsnapshots/tests')

# If this is not an installed package, append it to sys.path. Otherwise, Python
# won't find modules in ../
root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_path = root_path + '/src'
package_path = src_path + '/getsnapshots'
print("sys.path =",sys.path)                                                       # for debugging
assert package_path in sys.path or src_path in sys.path or root_path in sys.path   # for debugging
if not (package_path in sys.path or src_path in sys.path):
    sys.path.append(package_path)
    if os.path.exists('../setup.py'):
        print("Found ../setup.py but package was not in sys.path")

# Check if we are in a virtualenv.
if hasattr(sys, 'real_prefix'):
    print("Python is running in a virtualenv.")


# ----- Tests for pyrmsd.py ---------------------------------------------------

from getsnapshots.pyrmsd import read_xyz, read_c, read_cor, read_pdb, read_pair, get_rmsd

@pytest.fixture
def testfile():
    return "test01", "test02", "test03"

@pytest.fixture
def atoms_pair(testfile):
    files = (testfile[i] + ".pdb" for i in (0, 1))
    atoms1, atoms2 = read_pair(*files)
    return atoms1, atoms2

@pytest.fixture
def atoms_pair2(testfile):
    files = (testfile[i] + ".cor" for i in (0, 2))
    atoms1, atoms2 = read_pair(*files)
    return atoms1, atoms2

@pytest.fixture
def atoms_pair_fail(testfile):
    files = (testfile[0] + ".pdb", testfile[1] + ".xyz")
    atoms1, atoms2 = read_pair(*files)
    return atoms1, atoms2

def test_read_xyz(testfile):
    natoms, atoms, comment = read_xyz(testfile[0] + ".xyz")
    assert natoms == 27
    assert len(atoms) == natoms
    assert atoms[26] == [5.2146779, -0.1780056, -0.8875147, 'H']
    assert comment == "opt-FC_b3lyp_aug-cc-pVTZ, Energy = -723.7785131930"

def test_read_xyz_e01():
    with pytest.raises(FileNotFoundError):
        _, _, _ = read_xyz("non_existing.xyz")

def test_read_xyz_e02():
    with pytest.raises(ValueError):
        _, _, _ = read_xyz("test01.c")

def test_read_xyz_e03():
    natoms, atoms, comment = read_xyz("empty")
    assert natoms == 0
    assert atoms == []

def test_read_c(testfile):
    natoms, atoms, comment = read_c(testfile[0] + ".c")
    assert natoms == 4104
    assert atoms[119] == [0.4117886476578333, 14.754402106849028, -12.2876364981553, 'S']
    assert atoms[4103] == [9.655255693247717, 16.873287551257818, 11.485072379493625, 'H']
    assert comment == "chemsh coordinates from" + testfile[0] + ".c"

def test_read_c_e01():
    with pytest.raises(FileNotFoundError):
        _, _, _ = read_c("non_existing.xyz")

def test_read_c_e02(testfile):
        natoms, _, _ = read_c(testfile[0] + ".xyz")
        assert natoms == 0

def test_read_pdb(testfile):
    natoms, atoms, comment = read_pdb(testfile[0] + ".pdb", rm_link_atoms=True)
    assert natoms == 4107
    assert atoms[4103] == [10.669, 17.266, 11.050, 'H', 'SOLV', 'Wat', '461', 0.42]
    assert atoms[76] == [3.367, 10.945, -17.015, 'C', 'PROT', 'Glu', '10', 0.07]
    assert atoms[77] == [2.181, 10.890, -16.500, 'H', 'PROT', 'Glu', '10', 6.38]
    assert atoms[4106] == [18.503, 23.267, -21.028, 'O', 'ETOH', 'Eoh', '246', 14.73]
    assert atoms[945] == [1.515, 5.544, 1.527, 'N', 'PROT', 'CRO', '66', 3.00]

def test_read_pdb_num_resid(testfile):
    natoms, atoms, comment = read_pdb(testfile[0] + ".pdb", numeric_resid=True)
    assert atoms[76] == [3.367, 10.945, -17.015, 'C', 'PROT', 'Glu', 10, 0.07]

def test_read_pdb_e01():
    with pytest.raises(FileNotFoundError):
        natoms, atoms, comment = read_pdb("non_existing.cor")

def test_read_pdb_e02():
    natoms, _, _ = read_pdb("empty")
    assert natoms == 0

def test_read_cor(testfile):
    natoms, atoms, comment = read_cor(testfile[0] + ".cor", rm_link_atoms=True)
    assert natoms == 4104
    assert atoms[4103] == [9.65465, 16.87463, 11.48756, 'H', 'SOLV', 'Wat', '461', 0.417]
    assert atoms[76] == [3.28666, 11.06625, -16.96437, 'C', 'PROT', 'Glu', '10', 0.07]

def test_read_cor_e01():
    with pytest.raises(FileNotFoundError):
        natoms, atoms, comment = read_cor("non_existing.cor")

def test_read_cor_e02(testfile):
    with pytest.raises(NotImplementedError):
        natoms, _, _ = read_cor(testfile[0] + ".xyz")

def test_read_cor_e03(testfile):
    with pytest.raises(NotImplementedError):
        natoms, _, _ = read_cor(testfile[1] + ".cor")

def test_read_cor_e04():
    natoms, _, _ = read_cor("empty")
    assert natoms == 0

def test_read_pair(atoms_pair):
    assert len(atoms_pair[0]) == len(atoms_pair[1])
    assert atoms_pair[0][20] != atoms_pair[1][20]

def test_read_pair2(testfile):
    atoms1, atoms2 = read_pair(testfile[0] + ".c", testfile[1] + ".pdb")
    assert len(atoms1) == 4104
    assert len(atoms1[20]) == 4
    assert -4.48896 < atoms1[0][0] < -4.48894

def test_read_pair_e01(testfile):
    with pytest.raises(NotImplementedError):
        atoms1, atoms2 = read_pair(testfile[0] + ".c", "non-existent.foo")

def test_get_rmsd(atoms_pair):
    rmsd = get_rmsd(*atoms_pair)
    assert 0.64588 < rmsd[0] < 0.64590
    assert 0.426476 < rmsd[1]['PROT Gly 24'] < 0.426478

def test_get_rmsd_cor(atoms_pair2):
    atoms1 = atoms_pair2[0]
    gly = [atm[5] == 'Gly' for atm in atoms1]
    rmsd = get_rmsd(*atoms_pair2, selection=gly)
    assert any(gly)
    assert 0.367737 < rmsd[0] < 0.367739
    assert 0.424693 < rmsd[1]['PROT Gly 24'] < 0.424695

def test_get_rmsd_e01(atoms_pair_fail):
    rmsd = get_rmsd(*atoms_pair_fail)
    assert rmsd == (None, None)

def test_get_rmsd_e02():
    with pytest.raises(ValueError):
        rmsd = get_rmsd([], [])

def test_get_rmsd_e03(atoms_pair):
    with pytest.raises(ValueError):
        rmsd = get_rmsd(*atoms_pair, selection=[True for i in range(5)])

def test_get_rmsd_sel1(atoms_pair):
    rmsd = get_rmsd(*atoms_pair, selection=[True for i in atoms_pair[0]])
    assert 0.64588 < rmsd[0] < 0.64590

def test_get_rmsd_sel1(atoms_pair):
    rmsd = get_rmsd(*atoms_pair, selection=[False for i in atoms_pair[0]])
    assert rmsd[0] == 0
    assert len(rmsd[1]) == 0


# ---  Tests for get-snapshots.py ---------------------------------------------

# There is a conflict in the "inline" layout when module name and package name coincide.
# Here, import from getsnapshots fails, but copying it to snap2csv.py and importing from that
# works. Other potential problem: package name with underscore.

from getsnapshots.snap2csv import VARIABLES, CONCUR_TASKS, converged, get_value, define_variables, read_one, Variable, get_dihedral, read_all, get_omega, main

def test_getsnapshots_module():
    assert len(VARIABLES) == len(set(VARIABLES))

@pytest.fixture
def snap():
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    Min = 'v38s'
    Cycle = 13
    Frame = 50
    opt = 'cro-opt'
    i = str(Cycle) + '-' + str(Frame)
    pth = 'snapshots-' + str(Min)
    log = pth + '/' + opt + '-' + i + '.log'
    optc = pth + '/' + opt + '-' + i + '-opt.c'
    na, atoms, _ = read_c(optc)
    variables = define_variables()
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    verbose = True
    mode = 'analyse'
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    Snap = namedtuple('snap', 'log variables atoms shot Min Cycle Frame')
    return Snap(log, variables, atoms, snapshot, Min, Cycle, Frame)

@pytest.fixture
def snap_model(snap):
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, snap.Frame
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'model'
    verbose = False
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    Snap = namedtuple('snap', 'log variables atoms shot')
    return Snap(log, variables, atoms, shot=snapshot)

@pytest.fixture
def snap_predict(snap):
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, snap.Frame
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'predict'
    verbose = True
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    Snap = namedtuple('snap', 'log variables atoms shot')
    return Snap(log, variables, atoms, shot=snapshot)

def test_converged(snap):
    assert converged(snap.log) is True

def test_get_value(snap):
    variables = snap.variables
    atoms = snap.atoms
    phi = get_value(variables['phi'], atoms)
    assert phi == "12.07"
    hbW240E215 = get_value(variables['hbW240E215'], atoms)
    assert hbW240E215 == "1.760"
    W304rmsd = get_value(variables['W304rmsd'], atoms)
    assert W304rmsd == "1.39187"
    L199rot = get_value(variables['L199rot'], atoms)
    assert L199rot == "-66.47"
    bb66dihe = get_value(variables['bb66dihe'], atoms)
    assert bb66dihe == "100.74"
    E16cAcB = get_value(variables['E16cAcB'], atoms)
    assert E16cAcB == "190.78"
    # fake Variable type
    fields = 'var_type ind dihe_shift ref_coords'.split()
    values = ['hb', [2, 3], 0, [0, 0, 0]]
    fake_variable = dict(zip(fields, values))
    with pytest.raises(TypeError):
        get_value(fake_variable, atoms)
    # real Variable of unsupported type
    bad_variable = Variable('foo', 1, 0, [0, 0, 0])
    with pytest.raises(ValueError):
        get_value(bad_variable, atoms)

def test_get_dihe():
    # dihe of 90 deg and <90 deg
    atoms = [[-1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 1, 1]]
    dihe = get_dihedral(1, 2, 3, 4, atoms)
    assert dihe == 90
    atoms = [[-1, 0, 0], [0, 0, 0], [0, 1, 0], [-0.1, 1, 1]]
    dihe = get_dihedral(1, 2, 3, 4, atoms)
    assert 80 < dihe < 90

def test_read_one_analyse(snap):
    snapshot = snap.shot
    phi = snapshot['phi']
    assert phi == "12.07"
    hbW240E215 = snapshot['hbW240E215']
    assert hbW240E215 == "1.760"
    W304rmsd = snapshot['W304rmsd']
    assert W304rmsd == "1.39187"
    L199rot = snapshot['L199rot']
    assert L199rot == "-66.47"
    bb66dihe = snapshot['bb66dihe']
    assert bb66dihe == "100.74"
    E16cAcB = snapshot['E16cAcB']
    assert E16cAcB == "190.78"
    Omega = snapshot['Omega']
    assert Omega == "1.8418"

def test_read_one_model(snap_model):
    snapshot = snap_model.shot
    assert 292.16 < float(snapshot['K70cBcG']) < 292.18
    assert 189.27 < float(snapshot['K70cGcD']) < 189.29
    assert 40.98 < float(snapshot['I65cAcB']) < 41.0
    assert 163.32 < float(snapshot['I65rot']) < 163.34
    assert 28.34 < float(snapshot['alpha']) < 28.36
    assert 179.27 < float(snapshot['beta']) < 179.29
    assert 33.44 < float(snapshot['gamma']) < 33.46
    assert 190.80 < float(snapshot['E16cAcB']) < 190.82
    assert 144.86 < float(snapshot['E16rot']) < 144.88
    assert 28.22 < float(snapshot['E16orie']) < 28.24
    assert 167.04 < float(snapshot['M18cBcG']) < 167.06
    assert 62.14 < float(snapshot['M18cGsD']) < 62.16
    assert 2.251 < float(snapshot['Abs']) < 2.253

def test_read_one_predict(snap_predict):
    snapshot = snap_predict.shot
    assert 5.28 < float(snapshot['phi']) < 5.30
    assert 14.61 < float(snapshot['tau']) < 14.63
    assert 77.47 < float(snapshot['M66cAcB']) < 77.49
    assert 1.06 < float(snapshot['bE16oh']) < 1.08
    assert 2.89 < float(snapshot['hbE16acy']) < 2.91
    assert 1.37 < float(snapshot['hbE16W239']) < 1.39
    assert 1.585 < float(snapshot['hbW239acy']) < 1.587
    assert 3.209 < float(snapshot['hbQ213acy']) < 3.211
    assert 1.980 < float(snapshot['hbW239E16']) < 1.982
    assert 1.648 < float(snapshot['hbW304imi']) < 1.650
    assert 1.670 < float(snapshot['hbW260phe']) < 1.672
    assert 1.644 < float(snapshot['hbS146phe']) < 1.646
    assert 1.640 < float(snapshot['hbR95imi']) < 1.642
    assert 291.12 < float(snapshot['S111rot']) < 291.14

def test_read_one_e01(snap):
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = 'v38abs', snap.Cycle, snap.Frame
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'predict'
    verbose = True
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    assert snapshot == False

def test_read_one_e02(snap):
    # om2-opt-13-100.log not converged
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, 100
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'model'
    verbose = True
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    assert snapshot == False

def test_read_one_e03(snap):
    # om2-opt-13-200-opt.c does not exist
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, 200
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'model'
    verbose = True
    snapshot = loop.run_until_complete(read_one(
        Min, Cycle, Frame, variables, mode, semaphore, verbose))
    assert snapshot == False

def test_read_one_e04(snap):
    # om2-opt-13-150-mndo.out does not exist
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, 150
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'model'
    verbose = True
    with pytest.raises(ValueError):
        snapshot = loop.run_until_complete(read_one(
            Min, Cycle, Frame, variables, mode, semaphore, verbose))

def test_read_one_e05(snap):
    # cro-opt-13-150-mndo.out does not exist
    import asyncio
    from collections import namedtuple      # Use dataclass in Python 3.7
    loop = asyncio.get_event_loop()
    semaphore = asyncio.Semaphore(CONCUR_TASKS)
    Min, Cycle, Frame = snap.Min, snap.Cycle, 150
    log, variables, atoms = snap.log, snap.variables, snap.atoms
    mode = 'predict'
    verbose = True
    with pytest.raises(ValueError):
        snapshot = loop.run_until_complete(read_one(
            Min, Cycle, Frame, variables, mode, semaphore, verbose))

def test_read_all():
    Mins = ['v38s0']
    verbose = True
    mode = 'model'
    data = read_all(Mins, mode, verbose)
    assert len(data) == 1
    line = ','.join(str(data[0][x]) for x in VARIABLES)
    assert line == ('v38s0,13,20,2.2471,115.93,181.66,201.93,18.13,190.69,'
    '173.09,92.79,8.77,-1.57,186.05,157.15,-13.20,31.79,71.28,197.15,'
    '292.88,-69.98,289.59,271.02,292.38,70.67,-37.87,1.518,3.171,1.030,'
    '167.63,-59.35,2.766,1.956,1.834,3.103,4.222,3.777,2.748,2.591,5.719,'
    '1.566,3.722,5.389,1.826,3.289,3.122,-176.43,-68.79,3.880,5.849,4.003,'
    '195.73,82.46,1.642,1.638,1.655,3.163,1.785,1.878,6.776,5.982,3.574,'
    '1.873,1.778,4.137,3.110,6.606,0.22655,0.21945')

def test_read_all_e01():
    # fails for missing 13-150-mndo.out
    import asyncio
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    Mins = ['v38s']
    verbose = True
    for mode in ('analyse', 'model', 'predict'):
        with pytest.raises(ValueError):
            read_all(Mins, mode, verbose)

def test_get_omega():
    import subprocess
    import os
    assert get_omega('foo.out') == 'NA'
    # If not om2-opt or cro-opt, get_omega should just return 'NA'
    assert get_omega("snapshots-v38s/pre-opt-13-200.sge") == 'NA'
    # Same if mndo.out lacks 'State  2' but 'iuvcd=2 ipop=2' is found in mndo.in
    assert get_omega("failed_om2-opt-mndo.out") == 'NA'
    # "om2/cro-opt-*-mndo.in" is not allowed as mndo_out
    with pytest.raises(AssertionError):
        get_omega("snapshots-v38s/om2-opt-13-200-mndo.in")
    # Otherwise, get_omega tries to calculate missing omega. Fails if no mndo_exe
    assert get_omega("original_om2-opt-mndo.out") == 'NA'
    # Assure that mndo.in was restored and delete bak
    assert os.path.exists('original_om2-opt-mndo.in.bak')
    diff_call = 'diff original_om2-opt-mndo.in original_om2-opt-mndo.in.bak'
    diff = subprocess.call(diff_call, shell=True)
    assert diff == 0
    subprocess.call('rm original_om2-opt-mndo.in.bak', shell=True)

#def test_main():
