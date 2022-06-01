"""
Simple script to split the chains of a molecule such that idealise/relax can run correctly
"""
import biobox as bb
import sys, os

M = bb.Molecule(sys.argv[1])
M.guess_chain_split()
M.write_pdb(sys.argv[1])
