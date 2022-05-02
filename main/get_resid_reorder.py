import numpy as np
import biobox as bb
import sys, os
import argparse
import re
import pandas as pd

def get_reorder_resid(fasta, order, dimer, ligand=0):
    """
    Get the idx order to restructure an input PDB based on a fasta file (i.e. split dimers into two chains, and restructure from there)
    This also creates a new verbose fasta file based on the updated topology
    :param fasta: input verbose fasta file to reorder based on dimer domain sites
    :param order: order to place domains based on all_verbose.fasta (e.g. D4 D3 D2 D1 D5 etc.)
    :param dimer: dimerisation indices from the fasta file to base chain designations off
    :param ligand: Domain that a ligand is attached to (will treat the fasta differently if so) - should be start of fasta by definition from previous stages.
    """
    dimer = dimer - 1 # start counting domains from 1 (so need to count from 0 here)
    order = order - 1
    ligand -= 1
    fasta_chainA = []
    fasta_chainB = []
    reorder_residA = []
    reorder_residB = []
    reorder_resid = []
    fasta_record = 0 # current fasta length (i.e. current PDB structure prior to reordering)
    shift= 0 # linker shift (i.e. based on linker presence in fasta)
    chainA = True # are we transversing chain A right now?

    # first thing, split the dimerisation domains into two - this doesn't really work through if we have a ligand attached
    fasta_new = []
    cnt = 0
    for f in fasta:
        if f[0][:7] == ">linker" :
            fasta_new[-1].append(f) # join the linker onto the end of the last domain
        else:
            if np.any(dimer == cnt):
                # split this domain
                if cnt == ligand:
                    # ligand present (at start of fasta) - remove this first. Assume the ligand is a dimer (if it isn't, we need to reassess)
                    # also what if ligand on two non dimerising domains?
                    start = f[1][:12] # matching sequence for ligand using first 12 AA
                    ligand_len = [match.end() for match in re.finditer(start, f[1])][1] - 12
                    f_nonligand = f[1][2*ligand_len:] # skip over both monomeric units in ligand
                    f_split = int(len(f_nonligand) / 2)
                    fasta_new.append([f[0]+"_A", f[1][:f_split + ligand_len*2]]) # include the ligand here
                    fasta_new.append([f[0]+"_B", f_nonligand[f_split:]])
                else:
                    f_split = int(len(f[1]) / 2)
                    fasta_new.append([f[0]+"_A", f[1][:f_split]])
                    fasta_new.append([f[0]+"_B", f[1][f_split:]])
                cnt += 1
            else:
                fasta_new.append(f)
                cnt += 1

    # create an intermediate fasta file without all the periphral content to assist with knowing the resid reordering
    fasta_condensed = []
    for f in fasta_new:
        fasta_condensed.append(f[1])
        if len(f) == 3:
            fasta_condensed[-1] += f[2][1]

    # as part of the reordering, linkers between domains need to be moved with them
    # maybe throw up an error if there is a linker between two domains that are being moved into different positions?
    # since we're splitting domains, we move to chain B after the halfway point of order
    for cnt, o in enumerate(order):
        f = fasta_new[o]

        # the reordering is based on everything in the input PDB up till now - information stored in the condensed fasta
        length_so_far = len("".join(fasta_condensed[:o]))
        reorder_resid.append(np.arange(length_so_far, length_so_far + len(fasta_condensed[o]))) # this covers domains and linkers
        if chainA:
            fasta_chainA.append([f[0] + "_A", f[1]])

            #reorder_residA.append(np.arange(len(fasta[1])) + fasta_record)
            #fasta_record += len(fasta[1])
            if len(f) == 3: # a linker is present
                fasta_chainA.append([f[2][0] + "_A", f[2][1]])

            # check if we need to switch to chain B
            if cnt+1 == int(len(order)/2):
                chainA = False
        else:
            fasta_chainB.append([f[0] + "_B", f[1]])

            if len(f) == 3: # a linker is present
                fasta_chainB.append([f[2][0] + "_B", f[2][1]])

    # for now (for CSF1R) assume order is always from top to bottom - but this will need to be addressed later
    #for i in range(len(fasta)): # can't move through via order and shift based on linker (screws up ordering) - #TODO fix this

    #    if fasta[i][0][:7] == ">linker" :
    #            if chainA:
    #                fasta_chainA.append([fasta[i][0] + "_A", fasta[i][1]])
    #                # shift by current length of ordering
    #                reorder_residA.append(np.arange(len(fasta[i][1])) + fasta_record)
    #                fasta_record += len(fasta[i][1])
    #                fasta_chainB.append([fasta[i][0] + "_B", ""])
    #            else:
    #                fasta_chainB.append([fasta[i][0] + "_B", fasta[i][1]])
    #                reorder_residB.append(np.arange(len(fasta[i][1])) + fasta_record)
    #                fasta_record += len(fasta[i][1])
    #                fasta_chainA.append([fasta[i][0] + "_A", ""])
    #            shift += 1

    #    elif np.any(i-shift == dimer):

    #        # need a way of managing fasta linker counts
    #        #TODO future proof this against fasta of weird length (maybe doesn't matter during intermediate stages) - basically want nice way of restructuring fasta

    #        fasta_len = len(fasta[i][1]) # assumes homodimer of same length each side
    #        f_half0 = fasta[i][1][:int(fasta_len/2)]; f_half1 = fasta[i][1][int(fasta_len/2):]

    #        if chainA:
    #            # add verbose fasta
    #            fasta_chainA.append([fasta[i][0] + "_A", f_half0])
    #            fasta_chainB.append([fasta[i][0] + "_B", f_half1])
    #            # add idx order for restructuring PDB
    #            reorder_residA.append(np.arange(fasta_len/2) + fasta_record)
    #            fasta_record += fasta_len/2
    #            reorder_residB.append(np.arange(fasta_len/2) + fasta_record)
    #            fasta_record += fasta_len/2
    #            chainA = False # crossing point
    #        else:
    #            fasta_chainA.append([fasta[i][0] + "_A", f_half1])
    #            fasta_chainB.append([fasta[i][0] + "_B", f_half0])
    #            reorder_residB.append(np.arange(fasta_len/2) + fasta_record)
    #            fasta_record += fasta_len/2
    #            reorder_residA.append(np.arange(fasta_len/2) + fasta_record)
    #            fasta_record += fasta_len/2
    #            chainA = True # crossing point

    #    else: # we're not on a dimer domain
    #        if chainA:
    #             fasta_chainA.append([fasta[i][0] + "_A", fasta[i][1]])
    #             reorder_residA.append(np.arange(len(fasta[i][1])) + fasta_record)
    #             fasta_record += len(fasta[i][1])
    #             fasta_chainB.append([fasta[i][0] + "_B", ""])
    #        else:
    #             fasta_chainB.append([fasta[i][0] + "_B", fasta[i][1]])
    #             reorder_residB.append(np.arange(len(fasta[i][1])) + fasta_record)
    #             fasta_record += len(fasta[i][1])
    #             fasta_chainA.append([fasta[i][0] + "_A", ""])


    # Now use the new resid positions to rejuggle the pdb metadata
    #new_idx = np.concatenate((np.concatenate(reorder_residA).ravel(), np.concatenate(reorder_residB).ravel()))
    new_idx = np.concatenate(reorder_resid)

    return new_idx, [fasta_chainA, fasta_chainB]

def insert_resid(M, fasta_newchain):
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

    chainA_fasta = fasta_newchain[0]
    chainB_fasta = fasta_newchain[1]

    # get the length of chainA so we know how much to shift B insert positions by
    chainA_len = np.sum([len(x[1]) for x in chainA_fasta])

    # move through len of fasta file to know resid positions to insert between
    insert_positions = []
    insert_seq = []
    current_posA = 1 # resid numbering
    current_posB = 1 # need different numbering in case the sequences are different
    shift=0
    for cnt, f in enumerate(chainA_fasta): # loop through whatever chain is longer (or the same length)
        if f[0][:7] == ">linker":

            # check whether same is also true of chain B, if not we need to insert here!
            # check whether insert position is in A or B to know exact position
            if chainB_fasta[cnt-shift][0][:7] != ">linker": # insert here!
                insert_positions.append(chainA_len + current_posB)
                insert_seq.append(f[1])
                #current_pos += len(f[1])
                shift += 1
            #else:
                #current_pos += len(chainB_fasta[cnt][1])
            current_posA += len(f[1])
            current_posB += len(f[1])

        #    if f[1] == "": # insert position in chain A
        #        insert_positions.append(current_pos)
        #        insert_seq.append(chainB_fasta[i][1]) # sequence to insert at this position
        #        current_pos += len(chainB_fasta[i][1])
        #    else:
        #        insert_positions.append(current_pos + chainA_len)
        #        insert_seq.append(f[1])
        #        current_pos += len(f[1])

        elif chainB_fasta[cnt-shift][0][:7] == ">linker": # no linker on chain A but linker on chain B
            insert_positions.append(current_posA)
            insert_seq.append(chainB_fasta[cnt-shift][1])
            current_posA += len(chainB_fasta[cnt-shift][1])
            current_posB += len(chainB_fasta[cnt-shift][1])
            chainA_len += len(chainB_fasta[cnt-shift][1]) # add to this total length also since this extends the length of A
            shift -= 1
            # this effectively triggers a skip of chainA since we not accounting for any indexing drift there, so we need to add in the length of the missing domain we're skipping over
            current_posA += len(f[1])
            current_posB += len(chainB_fasta[cnt-shift][1])
        else:
            current_posA += len(f[1])
            current_posB += len(chainB_fasta[cnt-shift][1])

    return _insert_seq_(M, insert_positions, insert_seq), insert_positions, insert_seq

def create_loopfile(insert_positions, insert_seq, outname="loopfile"):
    """
    Create loopfile, e.g.:
    LOOP   287  291    0  0.0  
    LOOP  1205 1212    0  0.0  
    Based on input fasta (tells rosetta where to add loops)
    See https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/loops-file
    """
    #chainA_fasta = fasta[0]
    #chainB_fasta = fasta[1]

    # get the length of chainA so we know how much to shift B loop positions by
    #chainA_len = np.sum([len(x[1]) for x in chainA_fasta])

    # move through len of fasta file to know resid positions to insert between
    #current_pos = 1 # resid numbering
    #insert_positions = []
    #for i, f in enumerate(chainA_fasta):
    #    if f[0][:7] == ">linker":

    #        # check whether insert position is in A or B to know exact position
    #        if f[1] == "": # insert position in chain A
    #            insert_positions.append([current_pos, current_pos + len(chainB_fasta[i][1])])
    #            current_pos += len(chainB_fasta[i][1])
    #        else:
    #            insert_positions.append([current_pos + chainA_len, current_pos + chainA_len + len(f[1])])
    #            current_pos += len(f[1])

    #    else:
    #        current_pos += len(f[1])

    inserts = []
    for cnt, val in enumerate(insert_positions):
        inserts.append([val, val + len(insert_seq[cnt])])

    with open(outname, "w") as f:
        for p in inserts:
            f.write("LOOP  %4i %4i    0  0.0\n"%(p[0], p[1]))

# get new resid idx based on needed reordering for linker fill in
# need to know inital structure of input data (i.e. dimerisation domains etc.)
parser = argparse.ArgumentParser(description='Input parameters for restructuring resid needed for loopmodel etc.')
parser.add_argument('-s','--pdb', help='<Required> Input PDB structure to be reordered', required=True)
parser.add_argument('-d', '--dimer', nargs='+', required=True, help='<Required> List of PDB domains that either dimerise or are participate in LBD (in order from EC to CT), starting from index 1 from domain 1. This needs to match the dimer input from the initial scaffold building, which determined how the PDB would be built at the time (essentially we are undoing that phase here).')
parser.add_argument('-o', '--order', nargs='+', required=False, default=0, help='Order of PDB domains you need relative to verbose fasta (e.g. 2 1 4 3 means reorder as domain 2, then 1 etc.). Default of 0 means D1 D2 D3 etc. order of domains in verbose_fasta')
parser.add_argument('-f','--fasta', help='<Required> Verbose fasta file created in the previous assembly stage', required=True)
parser.add_argument('-l', '--ligand', default=0, required=False, help='Is there a ligand present? If so, which domain is it attached to?')

args = parser.parse_args()

# To run this correctly then, one needs to take a look at all_verbose.fasta to base their decisions on order and dimerisation
pdb = str(args.pdb)
# dimer domain means domains where we still need to split, with respect to the input all_verbose.fasta
dimer_domains = np.asarray(args.dimer).astype(int)
fasta_file= str(args.fasta)
ligand=int(args.ligand)

# order needs to correspond with the number of domains (NOT linkers) present in the all_verbose.fasta so we know how to reorder
# it also needs to know that the dimer array will split some of the domains in two (so you need to shift the new ordering after that point)
# You need to take the second half of the dimer domain as well to pair with what follows (e.g. 4 instead of 3 if the third domain is a dimer)
# e.g. order = [4, 5, 6, 3, 1, 2] if verbose_fasta looks like:
#>D2
#DEFLFTPRVETATETA
#>linker_D3
#W
#>D3
#ISLVTALHLVLGLSAVLGLLLL
#>D1
#SAYLNLSSEQNLIQEVTVGEGLNLKVMVEAYPGLQGFNWTYLGPFSDHQPEPKLANATTKDTYRHTFTLSLPRLKPSEAGRYSFLARNPGGWRALTFELTLRYPPEVSVIWTFINGSGTLLCAASGYPQPNVTWLQCSGHTDRCDEAQVLQVWDDPYPEVLSQEPFHKVTVQSLLTVETLEHNQTYECRAHNSVGSGSWAFIPSAYLNLSSEQNLIQEVTVGEGLNLKVMVEAYPGLQGFNWTYLGPFSDHQPEPKLANATTKDTYRHTFTLSLPRLKPSEAGRYSFLARNPGGWRALTFELTLRYPPEVSVIWTFINGSGTLLCAASGYPQPNVTWLQCSGHTDRCDEAQVLQVWDDPYPEVLSQEPFHKVTVQSLLTVETLEHNQTYECRAHNSVGSGSWAFIP
#>linker_D2
#ISAGAHTHPP
#>D2
#DEFLFTPRVETATETA
#>linker_D3
#W
#>D3
#ISLVTALHLVLGLSAVLGLLLL
# D1 is the dimerisation domain remaining

# alternatively, say all_verbose.fasta looks like:
#>D13
#EGICRNRVTNNVKDVTKLVANLPKDYMITLKYVPGMDVLPSHCWISEMVVQLSDSLTDLLDKFSNISEGLSNYSIIDKLVNIVDDLVECVKENSSKDLKKSFKSPEPRLFTPEEFFRIFNRSIDAFKDFVVASETSDCVVEGICRNRVTNNVKDVTKLVANLPKDYMITLKYVPGMDVLPSHCWISEMVVQLSDSLTDLLDKFSNISEGLSNYSIIDKLVNIVDDLVECVKENSSKDLKKSFKSPEPRLFTPEEFFRIFNRSIDAFKDFVVASETSDCVVEPSPPSIHPGKSDLIVRVGDEIRLLCTDPGFVKWTFEILDETNENKQNEWITEKAEATNTGKYTCTNKHGLSNSIYVFVRDPAKLFLVDRSLYGKEDSDTLVRCPLTDPEVTSYSLKGCQGKPLPKDLRFIPDPKAGIMIKSVKRAYHRLCLHCSVDQEGKSVLSEKFILKVRPAFKAVPVVSVSKASYLLREGEEFTVTCTIKDVSSSVYSTWKRENSQTKLQEKYNSWHHGDFNYERQATLTISSARVNDSGVFMCYANNTFGSANVTTTLEVVEPSPPSIHPGKSDLIVRVGDEIRLLCTDPGFVKWTFEILDETNENKQNEWITEKAEATNTGKYTCTNKHGLSNSIYVFVRDPAKLFLVDRSLYGKEDSDTLVRCPLTDPEVTSYSLKGCQGKPLPKDLRFIPDPKAGIMIKSVKRAYHRLCLHCSVDQEGKSVLSEKFILKVRPAFKAVPVVSVSKASYLLREGEEFTVTCTIKDVSSSVYSTWKRENSQTKLQEKYNSWHHGDFNYERQATLTISSARVNDSGVFMCYANNTFGSANVTTTLEVV
#>linker_D13
#DKG
#>D45
#FINIFPMINTTVFVNDGENVDLIVEYEAFPKPEHQQWIYMNRTFTDKWEDYPKSENESNIRYVSELHLTRLKGTEGGTYTFLVSNSDVNAAIAFNVYVNTKPEILTYDRLVNGMLQCVAAGFPEPTIDWYFCPGCSASVLPVDVQTLNSSGPPFGKLVVQSSIDSSAFKHNGTVECKAYNDVGKTSAYFNFAFINIFPMINTTVFVNDGENVDLIVEYEAFPKPEHQQWIYMNRTFTDKWEDYPKSENESNIRYVSELHLTRLKGTEGGTYTFLVSNSDVNAAIAFNVYVNTKPEILTYDRLVNGMLQCVAAGFPEPTIDWYFCPGCSASVLPVDVQTLNSSGPPFGKLVVQSSIDSSAFKHNGTVECKAYNDVGKTSAYFNFA
#>linker_D7
#RVETATETAW
#>D7
#ISLVTALHLVLGLSAVLGLLLLISLVTALHLVLGLSAVLGLLLL
#
# 1 O - O 2
#    \ /
#  4 OO 3
#    ||
#  5 OO 6 
# Numbers indicate domain trajectory in terms of counting
#
# They are all dimerisation domains, and you want to fill in just the linkers. Order would look like [1, 4, 5, 2, 3, 6] (note 1 will have a ligand attached if present)

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
        order = np.asarray(args.order).astype(int)

new_idx, fasta_newchain = get_reorder_resid(fasta, order, dimer_domains, ligand)
# now load in PDB and reorder metadata, and add in missing linkers in the correct place
A = bb.Molecule(pdb)
A.reorder_resid(new_idx)

# insert resid at position in pandas dataframe
S, insert_positions, insert_seq = insert_resid(A, fasta_newchain)
S.write_pdb(pdb)

# create loopfile for remodel
create_loopfile(insert_positions, insert_seq) #fasta_newchain)

# write a new all_verbose fasta file with the new fasta
with open("new_verbose.fasta", "w") as f:
    for n in fasta_newchain[0]:
        f.write(n[0] + "\n")
        f.write(n[1] + "\n")
    for n in fasta_newchain[1]:
        f.write(n[0] + "\n")
        f.write(n[1] + "\n")        
