"""
Functions for dealing with molecule names.
"""

import re
from .helpers import issubset

# a list of valid elements
valid_elements = set([
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
    'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
    'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
    'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
    'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
    'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    'Uut', 'Uuq', 'Uup', 'Uuh', 'Uuo'
    ])

def is_valid_element(tocheck):
    return tocheck in valid_elements

def is_valid_molecule(tocheck):
    return issubset(get_elements(tocheck), valid_elements)

def get_elements(molecule):
    els = re.compile('[A-Z][a-z]{0,2}')
    
    return set(els.findall(molecule))

def valence_long2short(ion):
    """
    Converts, e.g., Ca++ to Ca+2
    """
    charge = re.compile('[+-]+$')
    
    chg = charge.findall(ion)
    if chg:
        out = chg[0][0]
        n = len(chg[0])
        if n > 1:
            out = f"{out}{n}"
    
        return charge.sub(out, ion)
    else:
        return ion
    
def valence_short2long(ion):
    """
    Converts, e.g., Ca+2 to Ca++
    """
    charge = re.compile('[+-][0-9]+$')
    
    chg = charge.findall(ion)
    if chg:
        out = chg[0][0]
        
        if len(chg[0]) > 1:
            n = int(chg[0][1:])
            out *= n    
        return charge.sub(out, ion)
    else:
        return ion

def decompose_molecule(molecule, n=1):
    """
    Returns the chemical constituents of the molecule, and their number.

    Parameters
    ----------
    molecule : str
        A molecule in standard chemical notation with optional valence 
        in +N or -N format, e.g. 'CO2', 'HCO3' or 'B(OH)4', Ca+2.
    n : int
        number of times the molecule is repeated.
    
    Returns
    -------
    All elements in molecule with their associated counts and valence : dict
    """
    if isinstance(n, str):
        if n == '':
            n = 1
        n = int(n)
    
    # define regexs
    parens = re.compile('\(([A-z0-9()]+)\)([0-9]+)?')
    stoich = re.compile('([A-Z][a-z]?)([0-9]+)?')
    valence = re.compile('.*([+-][0-9]{1})')

    ps = parens.findall(molecule)  # find subgroups in parentheses
    rem = parens.sub('', molecule)  # get remainder
    
    if len(ps) > 0:
        for s, ns in ps:
            comp = decompose_molecule(s, ns)
        for k, v in comp.items():
            comp[k] = v * n
    else:
        comp = {}
        
    for e, ns in stoich.findall(rem):
        if e not in comp:
            comp[e] = 0
        if ns == '':
            ns = 1 * n
        else:
            ns = int(ns) * n
        comp[e] += ns
    
    # valence information, if present:
    if valence.match(rem):
        comp['valence'] = valence.findall(rem)[0]
    
    return comp

def SW():
    """
    Returns dict containing composition of standard seawater for use with make_PHREEQC_input.

    From Dickson, Sabine and Christian (2007) 'Best Practices in Ocean Acidification Research'
    Chapter 5, Table 3.
    """
    sw = {
        'Cl': 0.54586,
        'S(6)': 0.02824,
        'Br': 0.00084,
        'F': 0.00007,
        'Na': 0.46906,
        'Mg': 0.05282,
        'Ca': 0.01028,
        'K': 0.01021,
        'Sr': 0.00009,
        'B': 0.00042,
        'C(6)': 0.00204,
        'unit': 'mol/kgs'
    }

    return sw