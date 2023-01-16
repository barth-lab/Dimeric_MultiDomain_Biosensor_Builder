#!/usr/bin/env python3
"""
Create a fasta file with the reordered structure so we can add new domains trivially for a second assembly stage
"""
import sys, os
import numpy as np
import argparse
import biobox as bb
import re

# don't need to alter the PDB structures actually, we just need to restructure the fasta based on the input domains and linker information (like stage 1)
# THe main kind of "innovation" at this stage, is the reordering of the fasta file for all.fasta. The PDB domains themselves can be left immutable.
# So what input data do we want? Proably the previous verbose_fasta can be used to determine how we will restructure the fasta.

def pdb2fasta(A):
    """
    Convert a PDB to fasta 
    :param A: BioBox molecule object
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

    return AA_convert

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

    lmin = 3.0 # Ang.
    lmode= 12.7 # Ang.
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
    lmin = 3.0 # Ang.

    lmax = linker_max_span(linker_fasta)
    l = linker_span(linker_fasta) # x0

    if l == lmin:
        # to bring in line with previous work, increase by lmin
        tol = lmin # too small to constrict (i.e. 1 or 2 amino acid linkers)
    else:
        tol = lmax - l # shrink to 1 AA length minimum

    return l, tol

def linker_constraint_D0(fasta_verbose, linker):
    """
    Constraints are determined by the verbose fasta built already, with the length of the linker determined by get_linker_stats
    on the predetermined linker that's going to be added
    This linker constraint method has a different procotol for the domain added before what is prebuilt, and the one afterwards (linker_constraint_D1)
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
    x0, tol = get_linker_stats(linker)

    # remove \n from fasta verbose to make matching easier
    fasta_verbose = [f[:-1] for f in fasta_verbose]

    # the very last seq is repeated earlier on - the seq between that last domain and the end of the first "chain" 
    last_seq_chainA = np.where(np.asarray(fasta_verbose) == fasta_verbose[-1])[0][0]
    # now we want to move from there to the start of the end of the first chain in the dimer domain BEFORE our added domain of interest
    # first, what is the domain before added domain of interest?
    domain_added = np.where(np.asarray(fasta_verbose) == fasta_verbose[1])[0][1] 
    if fasta_verbose[domain_added - 3][:7] == ">linker":
        domain_added_before = domain_added - 4
    else: # no linker, so shift is different
        domain_added_before = domain_added - 2
    # now find actual shift between the start of the added domain (1)
    resid_distance = 0
    for f in fasta_verbose:
        if f.startswith(">"):
            continue
        elif f == fasta_verbose[domain_added_before]: # by definition, this should be a dimerisation domain - it assumes symmetry (homodimer...)
            resid_distance += int(len(fasta_verbose[domain_added_before]) / 2)
            break
        else:
            resid_distance += len(f)

    return "AtomPair CA 1 CA %i FLAT_HARMONIC %f 1 %f\n"%(resid_distance, x0, tol)

##### INPUT ARGUMENTS #####
parser = argparse.ArgumentParser(description='Input parameters for inserting domains into scaffold')
parser.add_argument('-s','--pdb', nargs='+', help='<Required> List of PDB domains used for construction. Need either 2 or 3 inputs, either DX D or DX D DZ, where DX is the domain closest to the TM to be added, D is the assembled structure from the previous round, and DZ the "highest" domain yet to be added.', required=True)
parser.add_argument('-d', '--dimer', nargs='+', required=False, default=1, help='List of PDB domains that either dimerise or are participate in LBD (in order from EC to CT), starting from index 1 from domain 1. This should correspond to the previous input (i.e. all_verbose.fasta)')
parser.add_argument('-p', '--add_domains', nargs='+', help='Where, with respect to the original all_verbose.fasta, are the new domains being added (starting from domain 1), i.e. if you are inserting the second domain from the first assembly stage in, here you should have -p 2')

# how does this fit in if you need a third stage of assembly?

# notification on which are dimerisation/LBD domains
args = parser.parse_args()

pdb = args.pdb
dimer_domains = np.asarray(args.dimer).astype(int)
domain_positions = np.asarray(args.add_domains).astype(int)

domain_count = 1
fasta_verbose = []
all_fasta = ""

f = open("../input_scaffold/all_verbose.fasta", 'r')
fasta_verbose_old = f.readlines()
f.close()

fasta_verbose_old_condensed = []
# just the domains version
linker = False
for f in fasta_verbose_old:
    if f.startswith(">linker"):
        linker = True
        continue
    elif f.startswith(">D"):
        linker=False
        continue
    elif linker:
        continue
    else:
        fasta_verbose_old_condensed.append(f)

# we also need to reorder the resid to match the fasta

#### Different behaviour if 1 or 2 domains are being inserted. Cover the easier case first (1 domain insertion)
if len(pdb) == 2:
    # load in the infoformation from the previous verbose_fasta
    # its the dimerisation domain underneath the added domain where we need to split
    # if you need to insert domains higher up... what do? Need to always move from N to C termini

    # find starting point for building fasta
    loc = np.where(np.asarray(fasta_verbose_old) == fasta_verbose_old_condensed[domain_positions[-1] - 1])[0][0]

    # first add in the new domain (which will be the start of your new fasta)
    for line in fasta_verbose_old[loc-1:]: # from the fasta title
        
        if np.any(line == np.array(fasta_verbose_old_condensed)[dimer_domains[-1] - 1]):
            # split line, half needs to go here, half at the other end of the pre-existing chain
            #line = line[:-1] # no \n on the final line
            line = line[:int(len(line)/2)] + "\n"
            fasta_verbose.append(line)
        else:
            fasta_verbose.append(line)

        if not line.startswith('>'):
            all_fasta += line[:-1] # remove \n
        
    # now complete the other, pre-existing chain
    for line in fasta_verbose_old:
        # treat differently for that dimer domain after the added domain
        if np.any(line == np.array(fasta_verbose_old_condensed)[dimer_domains[-1] - 1]):
            #line = line[:-1]
            line = line[:int(len(line)/2)] + "\n"
            fasta_verbose.append(line)
        else:
            fasta_verbose.append(line)

        if not line.startswith('>'):
            all_fasta += line[:-1] # remove \n
        
# write out the fasta file (both verbose and all)
# first verbose
with open("all_verbose.fasta", "w") as f:
    for fasta in fasta_verbose:
        f.write(fasta)

# then all
with open("all.fasta", "w") as f:
    f.write(">all\n")
    f.write(all_fasta)

##### DEFINE CONSTRAINTS #####
# now we need the domain constraint file 
# this is slightly different than before - now we need the linker constraint in the space not being filled
# for cheap and easy constraints - find what comes above the domain from all_verbose, plug into linker file
# also need to know location within all.fasta too ofc

constraint = []
# the constraint here is defined slightly differently - above for the first domain, and below for the second
if len(pdb) == 2:
    # first find linker between domains
    # there should be two instances of this domain remember, use that to your advantage
    domain_matches = np.where(np.asarray(fasta_verbose) == fasta_verbose[0])
    # the linker is given by the second domain_matches position minus 1
    linker = fasta_verbose[domain_matches[0][1] - 1][:-1] # remove \n
    link_constraint = linker_constraint_D0(fasta_verbose, linker) # works if domain being added is BEFORE prebuilt structure
    constraint.append(link_constraint)

with open("cst", "w") as f:
    for c in constraint:
        f.write(c)

# Finally, reorder the input PDBs to match the fasta (specifically the cluster centers)
if len(pdb) == 2:
    # ignore the first domain because that is what's being added
    cluster_fasta = fasta_verbose[2:]
    # ignore the next linker if it's present also
    if cluster_fasta[0].startswith(">linker"):
        cluster_fasta = fasta_verbose[4:]
    else:
        pass

    # remove \n and >
    cluster_fasta = [f[:-1] for f in cluster_fasta]; fasta_bool = [not f.startswith(">") for f in cluster_fasta]
    cluster_fasta = np.asarray(cluster_fasta)[fasta_bool]

    # load in PDB to reorder
    A = bb.Molecule(pdb[1])
    A_fasta = pdb2fasta(A)

    # now, based on cluster fasta and the known topology of the dimer, get the necessary reorder resid
    # find where we have the other chain in the dimer (that is already built and doesn't need rejuggling)
    loc = np.where(np.asarray(cluster_fasta) == cluster_fasta[-1])[0][0]  
    # everything up to and including loc needs to be placed first - for hopping over the dimer, every resid after can be grouped up
    # the resid up to the restart of the second chain are actually second 
    resid_cnt = []
    # join up the sequence after loc, which should be immutable from the input (unless second domain is being added, then other half of chain needs mainpulating) 
    remaining_fasta = ''.join(cluster_fasta[loc+1:])
    for cnt, f in enumerate(cluster_fasta):
        if cnt <= loc:
            resid_loc = [match.span() for match in re.finditer(f, A_fasta)][-1]
            resid_cnt.append(np.arange(resid_loc[0], resid_loc[1]))
            # remove these residues from A_fasta - actually A_fasta needs to be immutable w.r.t. PDB
            #fasta_list = [char for char in A_fasta]
            #del fasta_list[resid_loc[0]:resid_loc[1]]
            #A_fasta = ''.join(fasta_list)
    # join on the rest of the sequence
    resid_loc = [match.span() for match in re.finditer(remaining_fasta, A_fasta)][0]
    resid_cnt.append(np.arange(resid_loc[0], resid_loc[1]))

    # Now use the new resid positions to rejuggle the pdb metadata
    new_idx = np.concatenate(resid_cnt)
    A.data["chain"] = "A" # set all to A so the reordering works as expected
    A.reorder_resid(new_idx)
    A.write_pdb(pdb[1])

# remember, you need to also delete any TER lines that aren't the last one

#TODO: Figure out also how to deal with domains being added elsewhere (filling in the gaps, i.e. the ones we need to end on)
#TODO: Deal with heterodimers
