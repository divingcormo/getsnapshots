import pytest
import os

# Somehow Python needs to know where the package modules are.
# Use $PYTHONPATH (there must be a better way!!!)

# Need to have get-snapshots.csv.bash in cwd (hardcoded).
# All other data files are also in tests/ and subdirs -> cd tests
# No matter where os.chdir() appears in the code, it is executed at module time,
# i.e., before any test functions are called.
if not os.getcwd().endswith('tests'):
    print("test module called from", os.getcwd(), "- changing to tests/")
    os.chdir('tests')


from getsnapshots.pyrmsd import read_xyz

# @pytest.fixture
# def atoms1

def test_read_xyz():
    natoms, atoms, comment = read_xyz("test01.xyz")
    assert natoms == 27
    assert atoms[26] == [5.2146779, -0.1780056, -0.8875147, 'H']
    assert comment == "opt-FC_b3lyp_aug-cc-pVTZ, Energy = -723.7785131930"

def test_read_xyz_e01():
    with pytest.raises(FileNotFoundError):
        _, _, _ = read_xyz("non_existing.xyz")

def test_read_xyz_e02():
    with pytest.raises(ValueError):
        _, _, _ = read_xyz("test01.c")
