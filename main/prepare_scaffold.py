"""
Notes on usage
The linker file should only contain resid you want to remove
You should not include linkers at the start and end of the construct (e.g. D1 in EC and CT domains)
If there are any additional linker resid you desire (i.e. no structural data), you can add them in additional file
Cut PDB outputs contain no linker data (to be added in remodelling later)
"""

import sys, os
import numpy as np
import argparse
import biobox as bb

parser = argparse.ArgumentParser(description='Input parameters for scaffold prep')
parser.add_argument('-s','--pdb', nargs='+', help='<Required> List of PDB domains (in order from EC to CT) used for construction', required=True)
parser.add_argument('-d', '--dimer', nargs='+', required=False, default=1, help='List of PDB domains that either dimerise or are participate in LBD (in order from EC to CT), starting from index 1 from domain 1. While not required, this is crucial for the generated constraint files')
parser.add_argument('-a', '--linker_avoid', nargs='+', default=0, required=False, help='Domains to avoid saving linkers of in fasta file prior to loop reconstruction starting from index 1. Value equal 0 (default) assumes no domains will avoid linker cutting. Note you can still remove linkers with the linker position file, and potentially readd them with the extra linker option.')
parser.add_argument('-l', '--linker_positions', required=True, help='<Required> Text file containing the position of the linkers to remove from the input domains. Each line in the file must correspond to an input domain, and four (eight for dimer) columns for start/end points of cutting region (example: 1 3 X X means resid 1-3 will be identified as starting linkers and removed at the start of the protein, while nothing will be removed at the end.) This file must have the same number of lines as domains. Resid positions must be with respect to input data.')
parser.add_argument('-x', '--extra_linker', required=False, default=None, help='Text file containing the domain index and desired extra linker residues for START and END of domain, X = no extra residues (e.g. 2 KS X on a line corresponds to domain 2 getting KS extra linker at the start BEFORE the current linkers, and X (nothing in this case) AFTER the current linkers at the end of the domain). This can be applied to any domain.')

# notification on which are dimerisation/LBD domains
args = parser.parse_args()

pdb = args.pdb
dimer_domains = np.asarray(args.dimer).astype(int)
linker_avoid = np.asarray(args.linker_avoid).astype(int)
linkers = str(args.linker_positions)
extra_linkers = args.extra_linker

def pdb2fasta(A, name):
    """
    Convert a PDB to fasta 
    :param A: BioBox molecule object
    :param name: Description of fasta sequence (to be written as >name)
    """

    CA_idx = A.atomselect("*", "*", "CA", get_index=True)[1]
    resname = A.data["resname"][CA_idx].reset_index()["resname"]

    AA_knowledge = {"GLY": "G", "ALA": "A", "LEU": "L", "MET": "M", "PHE": "F",
                    "TRP": "W", "LYS": "K", "GLN": "Q", "GLU": "E", "SER": "S",
                    "PRO": "P", "VAL": "V", "ILE": "I", "CYS": "C", "TYR": "Y",
                    "HIS": "H", "ARG": "R", "ASN": "N", "ASP": "D", "THR": "T"}

    AA_convert = ""
    for i in resname:
        AA_convert += AA_knowledge[i]

    return [">%s"%(name), AA_convert]

def cut_linker(name, line):
    """
    Remove the linker atoms based on the input text file
    Dimer or not is decided by line length
    linker file should only contain atoms you actually want to remove...
    you should try and have the removed linker on only one side of the protein each time (to simplify things)
    e.g. D1 has no removed linkers, D2 has top and bottom, D3 has only bottom, etc
    :param name: name of domain to import
    :param line: line of linker line to use as basis to remove linkers
    """

    A = bb.Molecule(name)

    pdb_name = name.split('.')[0] # get name without .pdb suffix

    rem = np.array([])

    line = np.asarray(line.split(" "))
    # convert values of X (indicating nothing is to be cut) into zeros to avoid int issues
    # THis does mean only X is accepted for no slicing
    line[line == "X"] = 0
    line = line.astype(int)
    fasta_l_tmp = []
    fasta_l = []

    # check if dimer or not
    if len(line) == 4:
        dimer = False
    elif len(line) == 8:
        dimer = True
    else:
        raise Exception("ERROR: Number of atom index locations for %s in linker file are not monomeric (4 locations) or dimeric (8 locations)"%
                         (name))

    # first check if all zeros (i.e. no changes are being made)
    if np.sum(line) == 0:
        B = A
        fasta_d = pdb2fasta(B, pdb_name) # fasta of the domain
        fasta = [[], fasta_d, []]
    else: # otherwise, run as normal
        if dimer:
            for i in range(4):
                if line[i*2] == line[i*2 + 1]:
                    fasta_l_tmp.append([[], []]) # blank array to make comparison easier later
                else:
                    links = np.arange(line[i*2], line[i*2 + 1]+1)
                    rem = np.concatenate((rem, links))
                    link_idx = A.atomselect("*", links, "*", get_index=True)[1]
                    fasta_l_tmp.append( pdb2fasta( A.get_subset(link_idx), "linker_%s"%(pdb_name)) )     
    
            # don't need repeat linkers for the dimers (as it'll need to be same both sides)
            # for homodimer anyway - eventually for heterodimers we'll need to split this
            # that might just need to be a flag to introduce different behaviour, because otherwise it's too complicated to split hetero/homo
            if len(fasta_l_tmp[0][1]) > len(fasta_l_tmp[2][1]):
                fasta_l.append(fasta_l_tmp[0])
            else:
                fasta_l.append(fasta_l_tmp[2])
            if len(fasta_l_tmp[1][1]) > len(fasta_l_tmp[3][1]):
                fasta_l.append(fasta_l_tmp[1])
            else:
                fasta_l.append(fasta_l_tmp[3])

        else:
            for i in range(2):
                if line[i*2] == line[i*2 + 1]:
                    fasta_l.append([])
                else:
                    links = np.arange(line[i*2], line[i*2 + 1]+1)
                    rem = np.concatenate((rem, links))
                    link_idx = A.atomselect("*", links, "*", get_index=True)[1]
                    fasta_l.append( pdb2fasta( A.get_subset(link_idx), "linker_%s"%(pdb_name)) )    
    
        ######

        keep_idx = A.atomignore("*", rem, "*", get_index=True)[1]
        B = A.get_subset(keep_idx)

        # need the fastas to be correctly placed around domains
        fasta_d = pdb2fasta(B, pdb_name) # fasta of the domain
    
        # now create final fasta sequence for verboseness
        fasta = [fasta_l[0], fasta_d, fasta_l[1]]

    return B, fasta # return domain without links, and cooresonding fasta file 

def linker_constraint(pdb0, pdb1, fasta_parts, fasta_total):
    """
    The constraints will be obtained and then applied where no linkers are being added 
    and where we skip a domain. x0 is determined by the distance between C termini and N termini
    CA of the domains being constrained, while tol is determined by the corresponding linker length
    in the other side of the pre-built scaffold for this region (i.e. mirror image)
    """
    # e.g. AtomPair CA 195 CA  810 FLAT_HARMONIC 48.7 1  8 # D: 40.7, R: 4 
    # AtomPair CA 809 CA 1307 FLAT_HARMONIC 40.0 1 18 # D: 22, R: 9
    # FLAT_HARMONIC x0 sd tol 
    # Zero in the range of x0 - tol to x0 + tol. Harmonic with width parameter sd outside that range. Basically, a HARMONIC potential (see above) split at x0 with a 2*tol length region of zero inserted.
    # so sd=1 always, x0= expectec distanced (ang), tol the allowed stretching
    
    # WARNING: consider the geometry of the domains with respect to each other
    # Do they need to be pre-orientated prior to processing? 
    # Otherwise, will one dimer domain be flipped to accomdate the bioinformatic input?

    #
    #  O
    #   \
    #    O O
    #    I |
    #    O O
    #   / 
    #  O 
    #
    # I = constraint being defined (last of chain A in top dimer, first of chain B in bottom dimer)


    A = bb.Molecule(pdb0); B = bb.Molecule(pdb1)

    A.guess_chain_split() # split into A and B chains so we can access the monomeric unit we want
    last_resid = A.data["resid"][A.data["chain"] == "A"][-1]
    start_cwdA = A.atomselect("A", last_resid, "CA") # x y z of last CA coordinate

    # This time we need the second chain as we're moving back along the dimer
    B.guess_chain_split()
    first_resid = B.data["resid"][B.data["chain"] == "B"][0] # and its the firest resid of this chain
    end_cwdB = B.atomselect("B", first_resid, "CA")

    # can't use the coordinates - need to just use linker length and average size
    # cwd will be used in the domain constraint



def domain_constraint(pdbs, fasta_parts, fasta_total):
    """
    This will need to be tested (maybe with the VEGFR example?)
    """
    return 0

f = open(linkers, 'r')
lines = f.readlines()
f.close()

if extra_linkers is not None:
    f = open(extra_linkers, 'r')
    extras = f.readlines()
    f.close()

    extras = np.array(' '.join(np.asarray(extras)).split(' '))
    extra_domains = extras[::3].astype(int) # split information about domains and linker residiues
    mask = np.ones(extras.size, dtype=bool); mask[::3] = 0
    extra_resid = extras[mask].reshape(len(extra_domains), 2)

domain_count = 1
fasta_verbose = []
all_fasta = ""
for cnt, p in enumerate(pdb):

    pdb_name = p.split('.')[0]

    # cut linkers
    M, fasta_M = cut_linker(p, lines[cnt][:-1]) # remove \n from end of lines
    M.write_pdb(pdb_name + "_cut.pdb")

    if np.any(linker_avoid == domain_count):
 
        fasta_tmp = fasta_M[1]

        fasta_verbose.append([])
        fasta_verbose.append(fasta_tmp)
        fasta_verbose.append([])    

        # write all fasta
        fasta_tmp = fasta_tmp[1]    
    else:
        fasta_verbose = fasta_verbose + fasta_M
        # for all fasta
        fasta_tmp = fasta_M[1][1] # skip linkers

    if extra_linkers is not None and np.any(extra_domains == domain_count):
        d_idx = extra_domains == domain_count #domain idx of where to add extra linker

        extraSTART = extra_resid[d_idx][0][0]; extraEND = extra_resid[d_idx][0][1][:-1]

        if extraSTART != 'X':
            if len(fasta_verbose[-3]) == 0:
                fasta_verbose[-3] = ['>linker_%s'%(pdb_name), extraSTART]
                all_fasta += extraSTART
            else:
                fasta_verbose[-3][1] = extraSTART + fasta_verbose[-3][1]
                all_fasta += extraSTART + fasta_verbose[-3][1]
        else: 
            all_fasta += fasta_verbose[-3][1]

        all_fasta += fasta_tmp

        if extraEND != 'X':
            if len(fasta_verbose[-1]) == 0:
                fasta_verbose[-1] = ['>linker_%s'%(pdb_name), extraEND]
                all_fasta += extraEND
            else:
                fasta_verbose[-1] = fasta_verbose[-1][1] + extraEND 
                all_fasta += fasta_verbose[-1][1] + extraEND 
        else:
            all_fasta += fasta_verbose[-1][1]

    else: # write all fasta without the extras
        if len(fasta_verbose[-3]) != 0:
            all_fasta += fasta_verbose[-3][1]

        all_fasta += fasta_tmp

        if len(fasta_verbose[-1]) != 0:
            all_fasta += fasta_verbose[-1][1]

    domain_count += 1

# now, we want to move through and add linker regions if they are adjacent
fasta_verbose = [x for x in fasta_verbose if x != []]
fasta_verbose_collapsed = []
cnt = 0
for i in range(len(fasta_verbose)-1):
    txt1 = fasta_verbose[cnt][0]; txt2 = fasta_verbose[cnt+1][0]
    if "linker" in txt1 and "linker" in txt2:
        name1 = txt1.split("_")[1]; name2 = txt2.split("_")[1]
        new_fasta = [">linker_%s%s"%(name1, name2), fasta_verbose[cnt][1] + fasta_verbose[cnt+1][1]]
        fasta_verbose_collapsed.append(new_fasta)
        cnt += 2
    else:
        fasta_verbose_collapsed.append(fasta_verbose[cnt])
        cnt += 1

    if cnt >= len(fasta_verbose) - 1:
        fasta_verbose_collapsed.append(fasta_verbose[cnt])
        break

# write out the fasta file (both verbose and all)
# first verbose
with open("all_verbose.fasta", "w") as f:
    for fasta in fasta_verbose_collapsed:
        f.write(fasta[0] + "\n")
        f.write(fasta[1] + "\n")

# then all
with open("all.fasta", "w") as f:
    f.write(">all\n")
    f.write(all_fasta)
    
# now we need to identify and add constraints
# at each dimerisation point, we cross the mirror plane of the receptor
# the linker is built on the new side we're occupying, a constraint is applied on the other
# We need to know between which domains (and specfiically residues) constraints must be applied
# We need to know how long the linker is for tolerence
# the size of the in between domain not being constructed determines the actual distance
# For that value - you will need to calculate the difference between the start and end atoms of the domain not being built
# the count for the residue will be in the new structure (i.e. renumbered based on the fasta...)
# if there are multiple "hallucinated" domains, you need multiple domains in the distance calculation

# three cases:
# 1) two consecutive dimer domains (distance determined by linker)
# 2) two non-consecutive (distance determined by domain)
# 3) only 1 dimerisation site (no constraint necessary... at least in this round)

# condense down to strings so we can match for constraints
fasta_strings = np.asarray(fasta_verbose_collapsed)[:,0]

if len(dimer_domains) == 1:
    pass

elif len(dimer_domains) > 1:

    for d in range(len(dimer_domains) - 1): #  don't need constraint for last dimer domain

        domain0 = dimer_domains[d]; domain1 = dimer_domains[d+1]
        pdb_name0 = pdb[dimer_domains[d]-1].split('.')[0]
        pdb_name1 = pdb[dimer_domains[d+1]-1].split('.')[0]
        fasta_idx0 = [s for s, cnt in enumerate(fasta_strings) if pdb_name0 in cnt]
        fasta_idx1 = [s for s, cnt in enumerate(fasta_strings) if pdb_name1 in cnt]

        # remove start and end linkers if they exist
        if "linker" in fasta_verbose_collapsed[fasta_idx0[0]][0]:
            fasta_idx0 = fasta_idx0[1]
        else:
            fasta_idx0 = fasta_idx0[0]
        if "linker" in fasta_verbose_collapsed[fasta_idx1[-1]][0]:
            fasta_idx1 = fasta_idx1[-2]
        else:
            fasta_idx1 = fasta_idx1[-1]
        # the range between these two fasta_idx are what we're interested in for the constraint definition     

        if domain0 + 1 == domain1:
            # constraint determined by linkers
            # for pdb reference, indexing starts at zero, not 1 like the domain count
            constraint = linker_constraint(pdb[dimer_domains[d]-1], pdb[dimer_domains[d+1]-1], fasta_verbose_collapsed[fasta_idx0:fasta_idx1+1], all_fasta)
        
        elif domain0 + 1 < domain1:
            # constraint determined by hallucinated domains
            # no -1 because of python indexing
            constraint = domain_constraint(pdb[dimer_domains[d]-1 : dimer_domains[d+1]], fasta_verbose_collapsed[fasta_idx0:fasta_idx1+1], all_fasta)
        
        else:
            raise Exception("ERROR: There is something strange about your domain index input, did you put the -d indexes in ascending order?")

else:
    raise Exception("ERROR: You have indicated no dimerisation domains are present (empty -d field), you do not need this software if this is the case.") 
    




# for domains - how do we reintegrate linkers later if fasta is jumbled?

# the chunks of the pdb domains in the overall PDB may need to have be rejumbled each time
# dependingon what mp domain assembly actually does....

# don't want to be adding linkers at the ends (i.e. without them being between domains)
# linker file criticil information - needs to be chain A START, END, chain B START END