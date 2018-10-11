#!/home/wanko/bin/python3
# -*- coding: utf-8 -*-
"""Functions for atomic coordinate files and calculating the RMS displacement.

Supported import file formats: xyz, pdb, and chemsh coordinate files.

main() implements a shell command to calculate the RMSD between two coordinate
files.

This module stores atomic structures in atoms objects. An atoms object is a
list of lists:

[[x, y, z, element, optional fields],...]

The optional fields are (CHARMM nomenclature): segid, resname, resid.
"""
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

from math import sqrt
import sys
import os
import argparse
from collections import defaultdict

DEBUG = False


def read_xyz(filename):
    """Reads xyz atom coordinate files.

    Expected file format:
    ::
        n_atoms
        comment
        element_symbol x y z [optional fields (ignored)]

    Args:
        filename (str): The xyz file to be read.

    Returns:
        tuple(int, atoms, str): Number of atoms, atoms object, comment

        comment: 2nd line of the xyz file.
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError("{} not found".format(filename))
    with open(filename, 'r') as file_obj:
        lines = file_obj.readlines()
    atoms = []
    natoms, header_natoms = 0, 0
    header_comment = ""
    for NR, line in enumerate(lines, start=1):
        fields = line.split()
        NF = len(fields)
        if NR == 1:
            assert NF >= 1
            header_natoms = int(fields[0])
        elif NR == 2:
            header_comment = line.strip()
        elif NF >= 4:
            natoms += 1
            atom = [float(fields[i]) for i in range(1, 4)]
            atom.append(str(fields[0]).title())
            atoms.append(atom)
        elif 0 < NF < 4:
            print("WARNING: ignoring line", NR, ":", line, file=sys.stderr)
    # Some checks.
    if natoms != header_natoms:
        print("WARNING: found", natoms, "coordinates in", filename,
              "whereas header indicates", header_natoms, "atoms.", file=sys.stderr)
    if natoms == 0:
        print("WARNING: no atom coordinates found in file", filename, file=sys.stderr)
    return natoms, atoms, header_comment


def read_c(filename):
    """Reads chemsh coordinate files.

    Atomic coordinates are assumed in bohr and are converted to angstrom!
    Expected file format:
    ::
        ...
        block = coordinates records = (header_na)
        element_symbol x y z [optional fields (ignored)]
        block|EOF
        ...

    Args:
        filename (str): The file to be read.

    Returns:
        tuple(int, atoms, str): Number of atoms, atoms object, comment

        comment: Refers to the imported file.
    """
    ANGSTROM = 0.52917720859
    natoms = 0
    header_natoms = 0
    atoms = []
    if not os.path.isfile(filename):
        raise FileNotFoundError("{} not found".format(filename))
    with open(filename, 'r') as file_obj:
        lines = file_obj.readlines()
    header = True
    for line in lines:
        if line.startswith('block = coordinates records'):
            header = False
            fields = line.split()
            header_natoms = int(fields[-1])
            continue
        if header:
            continue    # skip this line
        if line.startswith('block'):
            break       # job done, ignore rest of file
        # read coordinates
        fields = line.split()
        NF = len(fields)
        if NF >= 4:
            natoms += 1
            atom = [float(fields[i])*ANGSTROM for i in range(1, 4)]
            atom.append(str(fields[0]).title())
            atoms.append(atom)
        elif 0 < NF < 4:
            print("WARNING: ignoring line:", line, file=sys.stderr)
    if natoms != header_natoms:
        print("WARNING: found", natoms, "coordinates in", filename,
              "whereas header indicates", header_natoms, "atoms.",
              file=sys.stderr)
    if natoms == 0:
        print("WARNING: no atom coordinates found in file", filename,
              file=sys.stderr)
    return natoms, atoms, 'chemsh coordinates from' + str(filename)


def read_cor(filename, **kwargs):
    """Reads CHARMM coordinate files.

    Args:
        filename (str): The file to be read.

    Keyword Args:
        rm_link_atoms (boole): If True, skips link atoms (name Q*).

    Returns:
        tuple(int, atoms, str): Number of atoms, atoms object, comment

        comment: Extracted from file header.
    """
    natoms, header_natoms = 0, 0
    nlink = 0
    NR = 0
    header_comment = ""
    atoms = []
    rm_link_atoms = kwargs.get('rm_link_atoms', False)
    if not os.path.isfile(filename):
        raise FileNotFoundError("{} not found".format(filename))
    with open(filename, 'r') as file_obj:
        lines = file_obj.readlines()
    for NR, line in enumerate(lines, start=1):
        fields = line.split()
        NF = len(fields)
        if line.startswith('*'):
            if NR == 2:
                header_comment = line[1:].strip()
            continue
        if NF == 1:
            header_natoms = int(line.strip())
            continue
        # Must be an atom, but assert correct number of fields.
        if NF != 10:
            msg = "read_cor: col separation not implemented\n"
            msg += "read_cor: cannot handle line {NR}:\n{line}"
            raise NotImplementedError(msg.format_map(vars()))
        atom = [float(fields[4]), float(fields[5]), float(fields[6]),
                str(fields[3])[0], fields[7].strip(), fields[2].title(),
                fields[8].strip(), float(fields[9].strip())]
        if atom[3] == 'Q' and rm_link_atoms:
            nlink += 1
            continue
        if atom[5] in ["Hoh", "Tip3", "Tip"]:
            atom[5] = "Wat"     # normalize resid for water
        natoms += 1
        atoms.append(atom)
    if natoms != header_natoms - nlink:
        print("WARNING: found", natoms, "coordinates in", filename,
              "whereas header indicates", header_natoms, "atoms.", file=sys.stderr)
    if natoms == 0:
        print("WARNING: no atom coordinates found in file", filename,
              file=sys.stderr)
    return natoms, atoms, header_comment


def read_pdb(filename, chain="A", conformation="A", rm_link_atoms=False,
             numeric_resid=False):
    """Reads Protein Data Bank and CHARMM PDB files.

    Args:
        filename (str): PDB file to be read.
        chain: Only the specified chain (and chain=" " atoms) is read.
        conformation: In superpositions, only the specified conformation is
            considered.
        rm_link_atoms: If True, skip link atoms (name Q*), default: False.
        numeric_resid: If True, resids are converted to int (default: str).

    Returns:
        tuple(int, atoms, str): Number of atoms, atoms object, comment

        comment: Refers to the imported file.
    """

    # Assign variables by column ranges (consider only ATOM/HETATM lines).
    RECORD = slice(0, 6)        # REMARK, ATOM, HETATOM, etc.
    INDEX = slice(6, 11)        # Atom serial number (not used)
    NAME = slice(12, 16)        # Atom name (type in charmm)
    ALTLOC = 16                 # Alternative location indicator
    RESNAME = slice(17, 21)     # Residue name (CHARMM uses the 4th letter)
    CHAINID = 21                # Chain identifier
    RESID = slice(22, 26)       # Residue sequence number
    ICODE = 26                  # Code for insertion of residues (not used)
    CRD = slice(30, 54)         # Coordinates (use X, Y, Z slices instead)
    X = slice(30, 38)
    Y = slice(38, 46)
    Z = slice(46, 54)
    OCCUPANCY = slice(54, 60)   # Currently not checked
    BFACTOR = slice(60, 66)     # Temperature (beta) factor
    SEGID = slice(72, 76)       # charmm segment ID (empty in databank files)
    ELEMENT = slice(76, 78)     # Element symbol
    CHARGE = slice(78, 80)      # Charge on the atom (currently ignored)

    na = 0
    nlink = 0
    n_lines_66 = 0
    n_lines_78 = 0
    n_lines_atm = 0
    n_lines_chain = 0
    n_lines_altloc = 0
    atoms = []

    if not os.path.isfile(filename):
        raise FileNotFoundError("{} not found".format(filename))
    with open(filename, 'r') as file_obj:
        lines = file_obj.readlines()

    for line in lines:
        # Require all fields up to B-factor.
        if len(line) < 66:
            continue
        n_lines_66 += 1

        # Check if the line contains an atom and whether to include it.
        record = line[RECORD].strip().lower()
        if record not in ("atom", "hetatm"):
            continue
        n_lines_atm += 1
        if not line[CHAINID].lower() in (chain.lower(), ' '):
            continue
        n_lines_chain += 1
        altloc = line[ALTLOC].strip().lower()
        if altloc not in ['', conformation.lower()]:
            continue
        n_lines_altloc += 1
        if len(line) < 78:
            # Extract element from atom type (only correct for single-letter
            # symbols).
            element = line[NAME].strip()[0]
        else:
            element = line[ELEMENT].strip()
            n_lines_78 += 1
        if rm_link_atoms and element == "Q":
            nlink += 1
            continue
        resname = line[RESNAME].strip().title()
        if line[SEGID].strip() != "":
            segid = line[SEGID].strip()
        elif resname in ["Hoh", "Tip3", "Tip"]:
            resname, segid = "Wat", "SOLV"
        elif resname in ["Eoh"]:
            segid = "ETOH"
        elif resname in ["Moh"]:
            segid = "MEOH"
        elif resname in ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln',
                         'Glu', 'Gly', 'His', 'Hsd', 'Hse', 'Hsp',
                         'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro',
                         'Ser', 'Thr', 'Trp', 'Tyr', 'Val']:
            segid = "PROT"
        elif resname in ['Cro', 'Cr8', 'Crq', 'Nrq', 'Cfy',
                         'Csy', 'Gyc', 'Gys', 'Iey', 'Imd']:
            resname, segid = "CRO", "PROT"
        else:
            segid = "NA"
        resid = line[RESID].strip().lstrip('0')
        if numeric_resid:
            resid = int(resid)

        # Add atom to list.
        na += 1
        atom = [float(line[X].strip()), float(line[Y].strip()),
                float(line[Z].strip()), element, segid, resname,
                resid, float(line[BFACTOR].strip())]
        atoms.append(atom)

    # Debug output.
    if DEBUG:
        if nlink > 0:
            print("read_cor:", nlink, "link atom(s) found.", file=sys.stderr)
        print("# lines with 66+ chars:", n_lines_66)
        print("# lines with 78+ chars:", n_lines_78)
        print("# lines with atoms:", n_lines_atm)
        print("# lines with atoms in chain " + chain + ":", n_lines_chain)
        print("# lines with atoms in conformation " + conformation + ":",
              n_lines_altloc)

    return na, atoms, "Data exctracted from {filename}".format_map(vars())


def get_rmsd(atoms1, atoms2, selection=None):
    """Returns the total and residue-resolved RMSD's between two structures.

    Args:
        atoms1/atoms2 (atoms object): See module doc.
        selection (list(boole)): List of same length. Unselected atoms are
            ignored.

    Returns:
        tuple(float, dict): The total RMSD, residue RMSD's.

        The dict uses as keys a string with the (space separated) segid,
        resname, resid from atoms, for example, "PROT Glu 16".
    """
    na1 = len(atoms1)
    na2 = len(atoms2)
    if na1 != na2:
        print("get_rmsd: different number of atoms:", na1, na2,
              file=sys.stderr)
        return None, None
    if na1 == 0:
        print("ERROR: get_rmsd called with empty atom objects")
        raise ValueError
    if selection:
        len_sel = len(selection)
        if len_sel != na1:
            msg = "get_rmsd: selection list has wrong length ({}, not {})"
            raise ValueError(msg.format(len_sel, na1))
    else:
        selection = [True for i in range(na1)]

    # Store residue RMSD.
    res_rmsd = defaultdict(lambda: 0, {})
    res_nsel = defaultdict(lambda: 0, {})
    rmsd, nsel = 0, 0
    # lambda version
    for i in range(na1):
        if selection[i]:
            nsel += 1
            crd1 = list(atoms1[i][j] for j in range(3))
            crd2 = list(atoms2[i][j] for j in range(3))
            key = " ".join(atoms1[i][4:7])
            squared_displacement = sum(map(lambda x, y: (x-y)**2, crd1, crd2))
            rmsd += squared_displacement
            res_rmsd[key] += squared_displacement
            res_nsel[key] += 1
    for key, val in res_rmsd.items():
        if res_nsel[key] != 0:
            res_rmsd[key] = sqrt(val/res_nsel[key])
    if nsel == 0:
        print("WARNING: atom selection is empty!", file=sys.stderr)
        return 0, res_rmsd
    assert sum(res_nsel.values()) == nsel
    return sqrt(rmsd/nsel), res_rmsd


def read_pair(filename1, filename2):
    "Reads coordinates from two files, returns two atoms objects."
    pair = []
    for filename in (filename1, filename2):
        if filename.endswith('.cor'):
            _, atoms, _ = read_cor(filename, rm_link_atoms=True)
        elif filename.lower().endswith('.pdb'):
            _, atoms, _ = read_pdb(filename, rm_link_atoms=True)
        elif filename.endswith('.c'):
            _, atoms, _ = read_c(filename)
        elif filename.endswith('.xyz'):
            _, atoms, _ = read_xyz(filename)
        else:
            msg = "Filetype of {} not recognized"
            raise NotImplementedError(msg.format(filename))
        pair.append(atoms)
    return pair[0], pair[1]


def main():
    "Parses command line arguments and prints RMSD(s)."
    parser = argparse.ArgumentParser(description='Calculate RMSD between two structures.')
    parser.add_argument("file1", help='Name of first file. Supported: xyz, chemsh, charmm, PDB.')
    parser.add_argument("file2", help='Name of second file.')
    parser.add_argument("-r",    help="Calculate also RMSD for each resid (file1)", action="store_true")
    args = parser.parse_args()

    atoms1, atoms2 = read_pair(args.file1, args.file2)
    rmsd, res_rmsd = get_rmsd(atoms1, atoms2)
    if not rmsd:
        print("No RMSD calculated.", file=sys.stderr)
        exit()
    if args.r:
        print("segid  residue    RMSD\n=============================")
        for res, val in res_rmsd.items():
            print("{:15} {:9.6f}".format(res, val))
        print("\nTotal RMSD:")
    print(rmsd)


if __name__ == "__main__":
    main()
