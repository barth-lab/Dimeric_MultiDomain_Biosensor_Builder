"""
Get the average interfacial energy (energy of binding) for a complex
"""

import numpy as np
import sys, os

def get_interface_score(filename):
    """
    Get the interfacial score from specific filename
    """
    data = np.loadtxt(filename, skiprows=2, usecols)


