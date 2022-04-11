import numpy as np
import biobox as bb
import sys, os
import argparse
import re
import pandas as pd

def get_reorder_resid(fasta, order, dimer):
    """
    Get the idx order to restructure an input PDB based on a fasta file (i.e. split dimers into two chains, and restructure from there)
    This also creates a new verbose fasta file based on the updated topology
    :param fasta: input verbose fasta file to reorder based on dimer domain sites
    :param order: order to place domains based on all_verbose.fasta (e.g. D4 D3 D2 D1 D5 etc.)
    :param dimer: dimerisation indices from the fasta file to base chain designations off
    """
    dimer = dimer - 1 # start counting domains from 1 (so need to count from 0 here)
    order = order - 1
    fasta_chainA = []
    fasta_chainB = []
    reorder_residA = []
    reorder_residB = []
    fasta_record = 0 # current fasta length (i.e. current PDB structure prior to reordering)
    shift= 0 # linker shift (i.e. based on linker presence in fasta)
    chainA = True # are we transversing chain A right now?
    # for now (for CSF1R) assume order is always from top to bottom - but this will need to be addressed later
    for i in range(len(fasta)): # can't move through via order and shift based on linker (screws up ordering) - #TODO fix this

        if fasta[i][0][:7] == ">linker" :
                if chainA:
                    fasta_chainA.append([fasta[i][0] + "_A", fasta[i][1]])
                    # shift by current length of ordering
                    reorder_residA.append(np.arange(len(fasta[i][1])) + fasta_record)
                    fasta_record += len(fasta[i][1])
                    fasta_chainB.append([fasta[i][0] + "_B", ""])
                else:
                    fasta_chainB.append([fasta[i][0] + "_B", fasta[i][1]])
                    reorder_residB.append(np.arange(len(fasta[i][1])) + fasta_record)
                    fasta_record += len(fasta[i][1])
                    fasta_chainA.append([fasta[i][0] + "_A", ""])
                shift += 1

        elif np.any(i-shift == dimer):

            # need a way of managing fasta linker counts
            #TODO future proof this against fasta of weird length (maybe doesn't matter during intermediate stages) - basically want nice way of restructuring fasta

            fasta_len = len(fasta[i][1]) # assumes homodimer of same length each side
            f_half0 = fasta[i][1][:int(fasta_len/2)]; f_half1 = fasta[i][1][int(fasta_len/2):]

            if chainA:
                # add verbose fasta
                fasta_chainA.append([fasta[i][0] + "_A", f_half0])
                fasta_chainB.append([fasta[i][0] + "_B", f_half1])
                # add idx order for restructuring PDB
                reorder_residA.append(np.arange(fasta_len/2) + fasta_record)
                fasta_record += fasta_len/2
                reorder_residB.append(np.arange(fasta_len/2) + fasta_record)
                fasta_record += fasta_len/2
                chainA = False # crossing point
            else:
                fasta_chainA.append([fasta[i][0] + "_A", f_half1])
                fasta_chainB.append([fasta[i][0] + "_B", f_half0])
                reorder_residB.append(np.arange(fasta_len/2) + fasta_record)
                fasta_record += fasta_len/2
                reorder_residA.append(np.arange(fasta_len/2) + fasta_record)
                fasta_record += fasta_len/2
                chainA = True # crossing point

        else: # we're not on a dimer domain
            if chainA:
                 fasta_chainA.append([fasta[i][0] + "_A", fasta[i][1]])
                 reorder_residA.append(np.arange(len(fasta[i][1])) + fasta_record)
                 fasta_record += len(fasta[i][1])
                 fasta_chainB.append([fasta[i][0] + "_B", ""])
            else:
                 fasta_chainB.append([fasta[i][0] + "_B", fasta[i][1]])
                 reorder_residB.append(np.arange(len(fasta[i][1])) + fasta_record)
                 fasta_record += len(fasta[i][1])
                 fasta_chainA.append([fasta[i][0] + "_A", ""])


    # Now use the new resid positions to rejuggle the pdb metadata
    new_idx = np.concatenate((np.concatenate(reorder_residA).ravel(), np.concatenate(reorder_residB).ravel()))

    return new_idx, [reorder_residA, reorder_residB], [fasta_chainA, fasta_chainB]

def insert_resid(M, fasta):
    """
    Insert resid at missing positions to match 
    Note: the intended use of this is for linker insertion, not domain (assembly should be used for that), so please ensure you're not inserting missing domains here.
    :params M: Biobox molecule to insert resid into
    :params fasta: list of fastas for chain A and B. There must be at least one missing linker for one of the chains - they should have the same lengths though, just some empty list postions
    """

    def _insert_seq_(M, pos, seq):
        """
        Insert sequence into positions given by pos in molecule M
        There is a heavy assumption here that everything has chain A as a designation
        :params M: Molecule object to insert to
        :params pos: positions to insert into (first resid to be shifted by seq essentially)
        :param seq: sequence to insert
        """

        ### METADATA ON AMINO ACIDS TO INSERT ###
        cols = M.data.columns
        AA_knowledge = {'G': 'GLY', 'A': 'ALA', 'L': 'LEU', 'M': 'MET', 'F': 'PHE',
                        'W': 'TRP', 'K': 'LYS', 'Q': 'GLN', 'E': 'GLU', 'S': 'SER',
                        'P': 'PRO', 'V': 'VAL', 'I': 'ILE', 'C': 'CYS', 'Y': 'TYR',
                        'H': 'HIS', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'T': 'THR'}

        atomnames = { "ARG" : ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
                      "VAL" : ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
                      "GLU" : ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
                      "THR" : ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
                      "ALA" : ["N", "CA", "C", "O", "CB"],
                      "TYR" : ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
                      "TRP" : ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                      "ILE" : ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"], 
                      "SER" : ["N", "CA", "C", "O", "CB", "OG"],
                      "GLY" : ["N", "CA", "C", "O"],
                      "HIS" : ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
                      "PRO" : ["N", "CA", "C", "O", "CB", "CG", "CD"],
                      "LEU" : ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
                      "PHE" : ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                      "CYS" : ["N", "CA", "C", "O", "CB", "SG"],
                      "LYS" : ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
                      "ASP" : ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
                      "ASN" : ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
                      "MET" : ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"], 
                      "GLN" : ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"]}
 
        # first, shift resid positions after pos by value of pos
        for i, p in enumerate(pos):
            # current values
            # subsequent p should consider the shifted values
            shift_val = np.asarray(M.data.loc[M.data["resid"] >= p, "resid"]) + len(seq[i])
            # now shift
            M.data.loc[M.data["resid"] >= p, "resid"] = shift_val 

        # now insert
        for i, p in enumerate(pos):
            # first build dummy Molecule to insert
            D = np.asarray([['', '', '', '', '', '', '', '', '', '']])
            for cnt, s in enumerate(seq[i]):
                atoms = atomnames[AA_knowledge[s]]
                arr = np.asarray([["ATOM"] * len(atoms), np.arange(len(atoms)), atoms, 
                      [AA_knowledge[s]] * len(atoms), ["A"]*len(atoms), [p + cnt] * len(atoms),
                      [1.0] * len(atoms), [-1.0] * len(atoms), [""] * len(atoms), [1.0] * len(atoms)])
                D = np.concatenate((D, arr.T), axis=0)

            D = D[1:] # remove placeholder
            # create dataframe from this
            df = pd.DataFrame(D, columns = cols)
            df = df.astype({cols[0]: str, cols[1]: int, cols[2]: str, cols[3]: str, cols[4]: str,
                            cols[5]: int, cols[6]: float, cols[7]: float, cols[8]: str, cols[9]: float})
            # create coordinates
            cwd = np.zeros((1, df.shape[0], 3))
        
            # get the index just before the insert position, and use as index basis for insertion
            idx = M.data[M.data["resid"] == p - 1].iloc[-1].name + 1
            M.data = pd.concat([
                     M.data.iloc[:idx],
                     df,
                     M.data.iloc[idx:]
                     ])
            # insert coordinates
            M.coordinates = np.concatenate((np.concatenate((M.coordinates[:, :idx],
                                             cwd), axis=1),
                                             M.coordinates[:, idx:]), axis=1)

            M.data.reset_index(drop=True, inplace=True)
            M.data["index"] = np.arange(len(M.data))

        return M

    chainA_fasta = fasta[0]
    chainB_fasta = fasta[1]

    # get the length of chainA so we know how much to shift B insert positions by
    chainA_len = np.sum([len(x[1]) for x in chainA_fasta])

    # move through len of fasta file to know resid positions to insert between
    insert_positions = []
    insert_seq = []
    current_pos = 1 # resid numbering
    for i, f in enumerate(chainA_fasta):
        if f[0][:7] == ">linker":
            
            # check whether insert position is in A or B to know exact position
            if f[1] == "": # insert position in chain A
                insert_positions.append(current_pos)
                insert_seq.append(chainB_fasta[i][1]) # sequence to insert at this position
                current_pos += len(chainB_fasta[i][1])
            else:
                insert_positions.append(current_pos + chainA_len)
                insert_seq.append(f[1])
                current_pos += len(f[1])

        else:
            current_pos += len(f[1])

    return _insert_seq_(M, insert_positions, insert_seq)

def create_loopfile(fasta, outname="loopfile"):
    """
    Create loopfile, e.g.:
    LOOP   287  291    0  0.0  
    LOOP  1205 1212    0  0.0  
    Based on input fasta (tells rosetta where to add loops)
    See https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/loops-file
    """
    chainA_fasta = fasta[0]
    chainB_fasta = fasta[1]

    # get the length of chainA so we know how much to shift B loop positions by
    chainA_len = np.sum([len(x[1]) for x in chainA_fasta])

    # move through len of fasta file to know resid positions to insert between
    current_pos = 1 # resid numbering
    insert_positions = []
    for i, f in enumerate(chainA_fasta):
        if f[0][:7] == ">linker":

            # check whether insert position is in A or B to know exact position
            if f[1] == "": # insert position in chain A
                insert_positions.append([current_pos, current_pos + len(chainB_fasta[i][1])])
                current_pos += len(chainB_fasta[i][1])
            else:
                insert_positions.append([current_pos + chainA_len, current_pos + chainA_len + len(f[1])])
                current_pos += len(f[1])

        else:
            current_pos += len(f[1])

    with open(outname, "w") as f:
        for p in insert_positions:
            f.write("LOOP  %4i %4i    0  0.0"%(p[0], p[1]))

# get new resid idx based on needed reordering for linker fill in
# need to know inital structure of input data (i.e. dimerisation domains etc.)
parser = argparse.ArgumentParser(description='Input parameters for restructuring resid needed for loopmodel etc.')
parser.add_argument('-s','--pdb', help='<Required> Input PDB structure to be reordered', required=True)
parser.add_argument('-d', '--dimer', nargs='+', required=True, help='<Required> List of PDB domains that either dimerise or are participate in LBD (in order from EC to CT), starting from index 1 from domain 1. This needs to match the dimer input from the initial scaffold building, which determined how the PDB would be built at the time (essentially we are undoing that phase here).')
parser.add_argument('-o', '--order', nargs='+', required=False, default=0, help='Order of PDB domains you need relative to verbose fasta (e.g. 2 1 4 3 means reorder as domain 2, then 1 etc.). Default of 0 means D1 D2 D3 etc. order of domains (just splitting by chains instead)')
parser.add_argument('-f','--fasta', help='<Required> Verbose fasta file created in the previous assembly stage', required=True)

args = parser.parse_args()

# testing
pdb = str(args.pdb)
dimer_domains = np.asarray(args.dimer).astype(int)
fasta_file= str(args.fasta)

#TODO how to deal with heterodimers?
fasta = []
with open(fasta_file, "r") as f:
    lines = f.readlines()
    for i in range(int(len(lines)/2)):
        fasta.append([lines[2*i][:-1], lines[2*i+1][:-1]]) # remove \n at end of line
    # dimer domains
    if args.order == 0:
        cnt = 1
        order = []
        for line in lines:
            if line[:7] != ">linker" and line[0] == ">":
                order.append(cnt)
                cnt +=1
            else:
                continue
        order = np.asarray(order)
    else:
        order = np.asarray(args.dimer).astype(int)

new_idx, reorder_resid, fasta_newchain = get_reorder_resid(fasta, order, dimer_domains)
# now load in PDB and reorder metadata, and add in missing linkers in the correct place
A = bb.Molecule(pdb)
A.reorder_resid(new_idx)

# insert resid at position in pandas dataframe
S = insert_resid(A, fasta_newchain)
S.write_pdb(pdb)

# create loopfile for remodel
create_loopfile(fasta_newchain)
