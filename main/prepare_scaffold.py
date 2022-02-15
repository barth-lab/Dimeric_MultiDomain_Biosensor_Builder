#!/usr/bin/env python3
"""
Notes on usage
The linker file should only contain resid you want to remove
You should not include linkers at the start and end of the construct (e.g. D1 in EC and CT domains)
If there are any additional linker resid you desire (i.e. no structural data), you can add them in additional file
Cut PDB outputs contain no linker data (to be added in remodelling later)
Note there can be NO gaps in the input files (e.e. space at end of a line)
"""

import sys, os
import numpy as np
import argparse
import biobox as bb
import re

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

def linker_max_span(fasta):
    """
    Get the maximum distance of a linker span based on input linker sequence
    see interesting article on loop stretch and span (to determine constraint length here) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3628373/
    here, the maximum span (lmax) is given by (n = number of AA, gamma = 6.046 Ang., delta = 3.46 Ang.)
    This more or less coresponds to the tolerance of the harmonic potential
    :param fasta: Input fasta (str) of linker
    """

    gamma = 6.046 # Ang.
    delta = 3.46 # Ang.
    n = len(fasta) # length of linker

    if n%2 == 0: # even linker length
        lmax = gamma * (n/2 -1) + delta 
    else: # odd linker length
        lmax = gamma * (n - 1) / 2    
    
    return lmax

def linker_span(fasta):
    """
    Get the most likely span distance (i.e. not fully stretched) of a linker based on statistical data and its length
    See here for distrubutions: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3628373/
    Note that it has been found span length is independent of the nature of local anchor points on domains or global topology
    Here, it was also found that loop span is independent of the number of residues, and is likely entirely a local folding problem
    It is captured by a simple Maxwell-Boltzmann distribution
    Note that the authors did note that very long loops between domains may have some local perturbation caused by the structure
    :param fasta: Input fasta (str) of linker
    """

    lmin = 3.8 # Ang.
    lmode=12.7 # Ang.
    lmax = linker_max_span(fasta)

    #def _maxw_(size = None, lmode=12.7): # Ang. (estimated from Gaussian kernel density estimation of all 20 000 loops checked in paper)
    #    """Generates size samples of maxwell"""
    #    vx = np.random.normal(size=size, scale=lmode / np.sqrt(2))
    #    vy = np.random.normal(size=size, scale=lmode / np.sqrt(2))
    #    vz = np.random.normal(size=size, scale=lmode / np.sqrt(2))
    #    return np.sqrt(vx*vx + vy*vy + vz*vz)

    #mdata = _maxw_(100000)
    #h, bins = np.histogram(mdata, bins = 101, range=(lmin, lmax))

    # rather than sample a boltzman dist., instead just use the mode (since it is seemingly independent), but constrain by lmax and lmin
    if lmax < lmin:
        l = lmin
    elif lmax < lmode:
        l = lmax
    else:
        l = lmode

    return l

def get_linker_stats(linker_fasta):
    """
    Get the span length (x0 for constraint, l) and maximum span length (tol for constraint, lmax) of linker
    :param linker_fasta: string of linker sequence
    """
    lmin = 3.8 # Ang.
 
    lmax = linker_max_span(linker_fasta)
    l = linker_span(linker_fasta) # x0

    if l == lmin:
        tol = 0 # too small to constrict (i.e. 1 or 2 amino acid linkers)
    else:
        tol = lmax - l + lmin # shrink to 1 AA length minimum

    return l, tol

def linker_constraint(pdb0, pdb1, fasta_parts, fasta_total):
    """
    The constraints will be obtained and then applied where no linkers are being added 
    and where we skip a domain. x0 is determined by the distance between C termini and N termini
    CA of the domains being constrained, while tol is determined by the corresponding linker length
    in the other side of the pre-built scaffold for this region (i.e. mirror image)

    Linker constraint is based on middle "chunk" in fasta_parts (i.e. this is only ever the linker between two domains)
    """
    # e.g. AtomPair CA 195 CA  810 FLAT_HARMONIC 48.7 1  8 # D: 40.7, R: 4 
    # FLAT_HARMONIC x0 sd tol 
    # Zero in the range of x0 - tol to x0 + tol. Harmonic with width parameter sd outside that range. Basically, a HARMONIC potential (see above) split at x0 with a 2*tol length region of zero inserted.
    # so sd=1 always, x0= expectec distanced (ang), tol the allowed stretching

    #
    #  O
    #   \
    #    O O
    #    I |
    #    O O
    #   / 
    #  O 
    #
    # I = constraint being defined between (last of chain A in top dimer, first of chain B in bottom dimer)
 
    # get linker constraint data
    x0, tol = get_linker_stats(fasta_parts[1][1])

    # Now we need to identify between which atoms we will place these constraints - remember we need the cut structure (i.e. no linkers in domains)
    A = bb.Molecule(pdb0.split('.')[0] + '_cut.pdb'); B = bb.Molecule(pdb1.split('.')[0] + '_cut.pdb')

    A.guess_chain_split() # split into A and B chains so we can access the monomeric unit we want
    A_chainA = A.atomselect("A", "*", "*", get_index=True)[1]
    A_sub = A.get_subset(A_chainA)
    fasta_d0 = pdb2fasta(A_sub, "tmp")[1] 

    # This time we need the second chain as we're moving back along the dimer
    B.guess_chain_split()
    B_chainB = B.atomselect("B", "*", "*", get_index=True)[1]
    B_sub = B.get_subset(B_chainB)
    fasta_d1 = pdb2fasta(B_sub, "tmp")[1] 

    # find the location of these in the larger fasta - we need the end of the first chain in domain0, and the start of the second chain in domains
    # in case of slight mismatch in AA, check for size of lise
    A_end = [match.end() for match in re.finditer(fasta_d0, fasta_total)]
    B_start = [match.start() for match in re.finditer(fasta_d1, fasta_total)]
    if len(A_end) == 2 or len(A_end) == 1:
        loc0 = A_end[0] - 1 # -1 because this returns first element AFTER the match substring ends
    else:
        raise Exception("ERROR: There are %i repeat units for one of your domains we're calculating constraints between"%(len(A_end)))
    
    if len(B_start) == 2:
        loc1 = B_start[1]
    elif len(B_start) == 1:
        loc1 = B_start[0]
    else:
        raise Exception("ERROR: There are %i repeat units for one of your domains we're calculating constraints between"%(len(A_end)))

    return ["AtomPair CA %i CA %i FLAT_HARMONIC %f 1 %f\n"%(loc0, loc1, x0, tol)]

def domain_constraint(pdbs, fasta_parts, fasta_total):
    """
    Apply constraints across domain gaps.
    """

    #
    #  O
    #   \
    #    O O
    #    I  \
    #    I   O
    #    I  /
    #    O O
    #   /
    #  O
    #
    # I = constraint being defined between (last of chain A in top dimer, first of chain B in bottom dimer)

    # In principle, we could have any number of domains and linkers between the two dimer sites
    cnt = 1 # count for fasta_parts array
    x0 = 0; tol = 0 # prep x0 and tol for final returned distance information
    for p in pdbs[1:-1]:
        P = bb.Molecule(p.split('.')[0] + '_cut.pdb')

        no_chains = P.guess_chain_split()[0]
        if no_chains != 1:
            raise Exception("ERROR: Intermediate domain %s detected to be a %i-mer. Please make sure there are no structural breaks in input data, or dimers in non -d flagged domains."%(p, no_chains))

        # first get distance between start and end of domain 
        first_resid = P.data["resid"].iloc[0]
        last_resid = P.data["resid"].iloc[-1]
        start_cwd = P.atomselect("*", [first_resid], "CA") # x y z of last CA coordinate
        last_cwd = P.atomselect("*", [last_resid], "CA")
        dist = np.linalg.norm(start_cwd - last_cwd)

        # now grab linker stats
        l, lmax = get_linker_stats(fasta_parts[cnt][1])

        # distance, x0, between domains is sum of linker span and domain sizes
        x0 += dist + l
        # tol is sum of lmax stats (slightly manipulated to accound for short linkers)
        tol += lmax

        cnt += 1

    # Finally, determine anchor points for constraint - see linker_constraint function for more commented information
    A = bb.Molecule(pdbs[0].split('.')[0] + '_cut.pdb'); B = bb.Molecule(pdbs[-1].split('.')[0] + '_cut.pdb')

    A.guess_chain_split()
    A_chainA = A.atomselect("A", "*", "*", get_index=True)[1]
    A_sub = A.get_subset(A_chainA)
    fasta_d0 = pdb2fasta(A_sub, "tmp")[1] 

    B.guess_chain_split()
    B_chainB = B.atomselect("B", "*", "*", get_index=True)[1]
    B_sub = B.get_subset(B_chainB)
    fasta_d1 = pdb2fasta(B_sub, "tmp")[1] 

    A_end = [match.end() for match in re.finditer(fasta_d0, fasta_total)]
    B_start = [match.start() for match in re.finditer(fasta_d1, fasta_total)]
    if len(A_end) == 2 or len(A_end) == 1:
        loc0 = A_end[0] 
    else:
        raise Exception("ERROR: There are %i repeat units for one of your domains we're calculating constraints between"%(len(A_end)))
    
    if len(B_start) == 2:
        loc1 = B_start[1] + 1 # + 1 to avoid pythonic indexing
    elif len(B_start) == 1:
        loc1 = B_start[0] + 1
    else:
        raise Exception("ERROR: There are %i repeat units for one of your domains we're calculating constraints between"%(len(A_end)))

    return ["AtomPair CA %i CA %i FLAT_HARMONIC %f 1 %f\n"%(loc0, loc1, x0, tol)]

##### INPUT ARGUMENTS #####
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

##### MAIN RUNCODE ######

##### DEFINE FASTA FILES #####

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

##### DEFINE CONSTRAINTS #####

# condense down to strings so we can match for constraints
fasta_strings = np.asarray(fasta_verbose_collapsed)[:,0]

constraint = []
if len(dimer_domains) == 1:
    pass
elif len(dimer_domains) > 1:

    for d in range(len(dimer_domains) - 1): #  don't need constraint for last dimer domain

        domain0 = dimer_domains[d]; domain1 = dimer_domains[d+1]
        pdb_name0 = pdb[dimer_domains[d]-1].split('.')[0]
        pdb_name1 = pdb[dimer_domains[d+1]-1].split('.')[0]
        fasta_idx0 = [s for s, cnt in enumerate(fasta_strings) if pdb_name0 in cnt]
        fasta_idx1 = [s for s, cnt in enumerate(fasta_strings) if pdb_name1 in cnt]

        #TODO: Put an exception here to stop people running linkers on the N and C termini of the protein
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
            constraint_linker = linker_constraint(pdb[dimer_domains[d]-1], pdb[dimer_domains[d+1]-1], fasta_verbose_collapsed[fasta_idx0:fasta_idx1+1], all_fasta)
            constraint += constraint_linker
        elif domain0 + 1 < domain1:
            # constraint determined by hallucinated domains and linkers
            # no -1 because of python indexing
            constraint_domain = domain_constraint(pdb[dimer_domains[d]-1 : dimer_domains[d+1]], fasta_verbose_collapsed[fasta_idx0:fasta_idx1+1], all_fasta)
            constraint += constraint_domain
        else:
            raise Exception("ERROR: There is something strange about your domain index input, did you put the -d indexes in ascending order?")

else:
    raise Exception("ERROR: You have indicated no dimerisation domains are present (empty -d field), you do not need this software if this is the case.") 
    
# now write the constraints to file
with open("cst", "w") as f:
    for c in constraint:
        f.write(c)

# Second stage domain progression (i.e. after the first initial models have been built), will require their own domain approach
# TM_pose_number is which domain (starting from 1) that contains the TM region for rosetta

##### NOTES ####
    
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

# for domains - how do we reintegrate linkers later if fasta is jumbled?

# the chunks of the pdb domains in the overall PDB may need to have be rejumbled each time
# dependingon what mp domain assembly actually does....

# don't want to be adding linkers at the ends (i.e. without them being between domains)
# linker file criticil information - needs to be chain A START, END, chain B START END