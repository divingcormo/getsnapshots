#!/home/wanko/bin/python3
from math import sqrt
import sys
#import re
import argparse
from collections import defaultdict

def read_xyz(filename):
    """Read xyz files.
    
    Expected format: 
        n_atoms
        comment
        element_symbol x y z [optional fields (ignored)]
    
    Returned:
    number of atoms (int), atoms object, comment line (str)
    atoms object (list of lists): [[x, y, z, element, optional fields],...]
        optional fields (charmm): segid, resname, resid, B-factor
    comment line: 2nd line of the xyz file
    """
    na = 0
    NR = 0
    atoms = []
    try:
        with open(filename, 'r') as file_obj:
            for line in file_obj:
                NR += 1
                fields = line.split()
                NF = len(fields)
                if NR == 1:
                    assert NF >= 1
                    header_na = int(fields[0])
                elif NR == 2:
                    header_comment = line.strip()
                elif NF >= 4:
                    na += 1
                    # I could also use a dict here
                    atom_obj = [
                        float(fields[1]), 
                        float(fields[2]), 
                        float(fields[3]), 
                        str(fields[0]).title()]
                    atoms.append(atom_obj)
                elif NF < 4 and NF > 0:
                    print("WARNING: ignoring line", NR, ":", line, file=sys.stderr)
    except FileNotFoundError:
        print("read_xyz: File " + filename + " not found.", file=sys.stderr)
        return 0, [], ''
    # Some checks.
    if na != header_na:
        print("WARNING: found", na, "coordinates in", filename,\
            "whereas header indicates", header_na, "atoms.", file=sys.stderr)
    if na == 0:
        print("WARNING: no atom coordinates found in file", filename, file=sys.stderr)
    #else:
    #    print "Unknown error with:", filename
    return na, atoms, header_comment

def read_c(filename):
    """Read chemsh coordinate files.
    
    Expected format: 
        ...
        block = coordinates records = (header_na)
        element_symbol x y z [optional fields (ignored)]
        block|EOF
        ...
    Atomic coordinates are assumed in bohr and are converted to angstrom!
    
    Returned:
    number of atoms (int), atoms object, comment line (str)
    atoms object (list of lists): [[x, y, z, element, optional fields],...]
        optional fields (charmm): segid, resname, resid, B-factor
    comment line: 2nd line of the xyz file
    """
    angstrom = 0.52917720859
    na = 0
    header_na = 0
    atoms = []
    try:
        with open(filename, 'r') as file_obj:
            ignore = True
            for line in file_obj:
                if line.startswith('block = coordinates records'):
                    ignore = False
                    fields = line.split()
                    header_na = int(fields[-1])
                    continue
                if not ignore:
                    if not line.startswith('block'):
                        # read coordinates
                        fields = line.split()
                        NF = len(fields)
                        if NF >= 4:
                            na += 1
                            # I could also use a dict here
                            atom_obj = [
                                float(fields[1])*angstrom, 
                                float(fields[2])*angstrom, 
                                float(fields[3])*angstrom, 
                                str(fields[0]).title()]
                            atoms.append(atom_obj)
                        elif NF < 4 and NF > 0:
                            print("WARNING: ignoring line:", line, file=sys.stderr)
                    else:
                        # job finished
                        break
    except FileNotFoundError:
        print("read_c: File " + filename + " not found.", file=sys.stderr)
        return 0, [], ''
    # Some checks.
    if na != header_na:
        print("WARNING: found", na, "coordinates in", filename,\
            "whereas header indicates", header_na, "atoms.", file=sys.stderr)
    if na == 0:
        print("WARNING: no atom coordinates found in file", filename, file=sys.stderr)
    #else:
    #    print "Unknown error with:", filename
    return na, atoms, 'chemsh coordinates from' + str(filename)

def read_cor(filename, **kwargs):
    """Read charmm coordinate files.
    
    Options:
        rm_link_atoms: if True, skip link atoms (type QQ*), default: False
    
    Returned:
    number of atoms (int), atoms object, comment line (str)
    atoms object (list of lists): [[x, y, z, element, optional fields],...]
        optional fields (charmm): segid, resname, resid, B-factor
    """
    na = 0
    nlink = 0
    NR = 0
    atoms = []
    rm_link_atoms = kwargs.get('rm_link_atoms', False)
    try:
        with open(filename, 'r') as file_obj:
            for line in file_obj:
                NR += 1
                fields = line.split()
                NF = len(fields)
                if line.startswith('*'):
                    if NR == 2:
                        header_comment = line[1:].strip()
                elif NF == 1:
                    header_na = int(line.strip())
                else:
                    # Must be an atom.
                    na += 1
                    # I could also use a dict for atom_obj.
                    # If there are 10 fields use field separator,
                    # otherwise figure out columns for each field.
                    if NF == 10:
                        atom_obj = [
                            float(fields[4]),
                            float(fields[5]),
                            float(fields[6]),
                            str(fields[3])[0],
                            fields[7].strip(),
                            fields[2].title(),
                            fields[8].strip()]
                        if atom_obj[3] == 'Q':
                            if rm_link_atoms:
                                nlink += 1
                                na -= 1
                            else:
                                atoms.append(atom_obj)
                        else:
                            atoms.append(atom_obj)
                    else:
                        print("read_cor: col separation not implemented", file=sys.stderr)
                        print("read_cor: cannot handle line {na}:\n{line}".format_map(vars()), file=sys.stderr)
                        return 0, [], ''
    except FileNotFoundError:
        print("read_cor: File " + filename + " not found.", file=sys.stderr)
        return 0, [], ''
    # Some checks.
    #if nlink > 0:
    #    print("read_cor:",nlink,"link atom(s) have been removed.", file=sys.stderr)
    if na != header_na - nlink:
        print("WARNING: found", na, "coordinates in", filename,\
            "whereas header indicates", header_na, "atoms.", file=sys.stderr)
    if na == 0:
        print("WARNING: no atom coordinates found in file", filename, file=sys.stderr)
    #else:
    #    print("Unknown error with:", filename)
    return na, atoms, header_comment

def read_pdb(filename, chain="A", conformation="A", rm_link_atoms=False):
    """Read Protein Data Bank and CHARMM PDB files.
    
    Options:
        chain: Consider only the specified chain (and chain=" " atoms)
        conformation: In superpositions, consider only the specified conformation
        rm_link_atoms: If True, skip link atoms (type Q*), default: False
    
    Returned:
    number of atoms (int), atoms object, comment line (str)
    atoms object (list of lists): [[x, y, z, element, optional fields],...]
        optional fields (charmm): segid, resname, resid, B-factor
    """
    na = 0
    nlink = 0
    n_lines_66 = 0
    n_lines_78 = 0
    n_lines_atm = 0
    n_lines_chain = 0
    n_lines_altloc = 0
    atoms = []

    # Debug
    #print("read_pdb: filename is", filename, ", chain:", chain, "conf:", conformation,
    #    ", rm_link_atoms:", rm_link_atoms)
    
    # If more kwargs are needed, replace them all by **kwargs
    #rm_link_atoms = kwargs.get('rm_link_atoms', False)

    # Assign variables by column ranges (consider only ATOM/HETATM)
    RECORD = slice(0,6)
    INDEX = slice(6,11)     # Atom serial number
    NAME = slice(12,16)     # Atom name (type in charmm)
    ALTLOC = 16             # Alternative location indicator
    RESNAME = slice(17,20)  # Residue name
    CHAINID = 21            # Chain identifier
    RESID = slice(22,26)    # Residue sequence number
    ICODE = 26              # Code for insertion of residues
    CRD = slice(30,54)
    X = slice(30,38)
    Y = slice(38,46)
    Z = slice(46,54)
    OCCUPANCY = slice(54,60)
    BFACTOR = slice(60,66)  # Temperature (beta) factor
    SEGID = slice(72,76)    # charmm segment ID (empty in databank files)
    ELEMENT = slice(76,78)  # Element symbol
    CHARGE = slice(78,80)   # Charge on the atom (currently ignored)

    try:
        with open(filename, 'r') as file_obj:
            for line in file_obj:
                
                # Require all fields up to B-factor
                if len(line) < 66:
                    continue
                n_lines_66 += 1
                    
                # Check if the line contains an atom and whether to include it
                record = line[RECORD].strip().lower()
                if record != "atom" and record != "hetatm":
                    continue
                n_lines_atm += 1
                if line[CHAINID].lower() in (chain.lower(), ' ') == False:
                    continue
                n_lines_chain += 1
                altloc = line[ALTLOC].strip().lower()
                if altloc not in ['', conformation.lower()]:
                    continue
                n_lines_altloc += 1
                if len(line) < 78:
                    # Extract element from atom type
                    # (only correct for single-letter symbols)
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
                elif resname in ["Hoh", "Tip3"]:
                    resname = "W"
                    segid = "SOLV"
                elif resname in ["eoh"]:
                    segid = "ETOH"
                elif resname in ["moh"]:
                    segid = "MEOH"
                elif resname in ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln',
                        'Glu', 'Gly', 'His', 'Hsd', 'Hse', 'Hsp', 'Ile', 'Leu',
                        'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr',
                        'Val']:
                    segid = "PROT"
                elif resname in ['Cro', 'Cr8', 'Crq', 'Nrq', 'Cfy',
                        'Csy', 'Gyc', 'Gys', 'Iey', 'Imd']:
                    resname = "CRO"
                    segid = "PROT"
                else:
                    segid = "NA"
                    
                # Add atom to list
                na += 1
                atom_obj = [
                    float(line[X].strip()),
                    float(line[Y].strip()),
                    float(line[Z].strip()),
                    element,
                    segid,
                    resname,
                    line[RESID].strip(),
                    float(line[BFACTOR].strip())]
                atoms.append(atom_obj)
                
    except FileNotFoundError:
        print("read_cor: File " + filename + " not found.", file=sys.stderr)
        return 0, [], ''
    
    # Some checks for debugging.
    if False:
        if nlink > 0:
            print("read_cor:",nlink,"link atom(s) found.", file=sys.stderr)
        print("# lines with 66+ chars:", n_lines_66)
        print("# lines with 78+ chars:", n_lines_78)
        print("# lines with atoms:", n_lines_atm)
        print("# lines with atoms in chain " + chain + ":", n_lines_chain)
        print("# lines with atoms in conformation " + conformation + ":", n_lines_altloc)
        
    return na, atoms, "Data exctracted from {filename}".format_map(vars())

def get_rmsd (atoms1, atoms2, selection=None):
    """Return the RMSD between two structures.
    
    atoms1, atoms2, and the optional list selection must have the same length.
    atoms1, atoms2 (list of lists): [[x, y, z, element, optional fields],...]
        optional fields (charmm): segid, resname, resid
    selection: list of logicals
    """
    na1 = len(atoms1)
    na2 = len(atoms2)
    if na1 != na2:
        print("get_rmsd: different number of atoms:", na1, na2, file=sys.stderr)
        return None, None
    if na1 == 0:
        print("ERROR: get_rmsd called with empty atom objects")
        raise ValueError
    if selection:
        nasel = len(selection)
        if na1 != nasel:
            print("get_rmsd: selection list has wrong length:", nasel, file=sys.stderr)
            return None, None
    else:
        selection = []
        # SyntaxError (why?):
        #selection[i] = True for i in range(0,na1)
        for i in range(0,na1):
            selection.append(True)
        #print("Including all atoms.")
    rmsd = 0
    # Store residue RMSD
    res_rmsd = defaultdict(lambda: 0, {})
    nsel = 0
    # lambda version
    for i in range(0, na1):
        if selection[i]:
            nsel += 1
            crd1 = list(atoms1[i][j] for j in range(0,3))
            crd2 = list(atoms2[i][j] for j in range(0,3))
            x = sum(map(lambda x, y: (x-y)**2, crd1, crd2))
            rmsd += x
            res_rmsd[" ".join(atoms1[i][4:7])] += x
    for key, val in res_rmsd.items():
        if nsel != 0:
            res_rmsd[key] = sqrt(val/nsel)
    if nsel == 0:
        print("ERROR: atom selection is empty!")
        raise ValueError
    return sqrt(rmsd/nsel), res_rmsd

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate RMSD between two structures.')
    parser.add_argument("file1", help='Name of first file. Supported: xyz, chemsh, charmm, PDB.')
    parser.add_argument("file2", help='Name of second file.')
    parser.add_argument("-r", help="Calculate also RMSD for each resid (file1)", action="store_true")

    args = parser.parse_args()

    #print("\nReading .cor and opt.c files...", file=sys.stderr)
    filename = args.file1
    if filename.endswith('.cor'):
        na1, atoms1, comment = read_cor(filename, rm_link_atoms=True)
    elif filename.lower().endswith('.pdb'):
        na1, atoms1, comment = read_pdb(filename, rm_link_atoms=True)
    elif filename.endswith('.c'):
        na1, atoms1, comment = read_c(filename)
    elif filename.endswith('.xyz'):
        na1, atoms1, comment = read_xyz(filename)
    else:
        print("Filetype of ",filename,"not recognized", file=sys.stderr)
        exit()
    filename = args.file2
    if filename.endswith('.cor'):
        na2, atoms2, comment = read_cor(filename, rm_link_atoms=True)
    elif filename.lower().endswith('.pdb'):
        na2, atoms2, comment = read_pdb(filename, rm_link_atoms=True)
    elif filename.endswith('.c'):
        na2, atoms2, comment = read_c(filename)
    elif filename.endswith('.xyz'):
        na2, atoms2, comment = read_xyz(filename)
    else:
        print("Filetype of ",filename,"not recognized", file=sys.stderr)
        exit()
    
    #print "na1: {} na2: {}".format(na1, na2)
    #print "\nCalculate RMSD..." 
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
    
    
