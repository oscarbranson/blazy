"""
Functions for handling iPHREEQC inputs and outputs.
"""

import os
import re
import platform
import pandas as pd
import pkg_resources as pkgrs
from .parser import database

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

default_output = """SELECTED_OUTPUT
    -pH
    -temperature
    -alkalinity
    -ionic_strength
    -totals Cl Na Mg K B Ca C S(6)
    -m OH- H+
    # molality outputs
    -m B(OH)4- B(OH)3 CaB(OH)4+ MgB(OH)4+ NaB(OH)4 B3O3(OH)4- B4O5(OH)4-2  # boron
    -m HCO3- CO3-2 CO2 CaCO3 MgCO3 BaCO3 SrCO3 NaCO3- KCO3- # carbon
    -m SO4-2 HSO4-  # S
    # activities
    -a OH- H+
    -a B(OH)4- B(OH)3 CaB(OH)4+ MgB(OH)4+ NaB(OH)4 B3O3(OH)4- B4O5(OH)4-2  # boron
    -si Calcite Aragonite
    """ 

def get_database_path(database_name='pitzer'):
    database_dir = os.path.join(pkgrs.resource_filename('blazy', 'resources'), 'database')
    return os.path.join(database_dir, database_name.replace('.dat', '') + '.dat')

def make_solution(inputs, n=1):
    inp = [f"SOLUTION {int(n):d}"]
    for k, v in inputs.items():
        if isinstance(v, str):
            inp.append(f'    {k:20s}{v:s}')
        else:
            inp.append(f'    {k:20s}{float(v):.8e}')
    return '\n'.join(inp) + '\n'

# def make_default_output(database, inputs):

def input_str(inputs, outputs=None):
    """
    Generate an input for calculating PHREEQC solutions.

    Parameters
    ----------
    inputs : dict, or list of dicts
        Where each key is a valid PHREEQC input key, and each 
        value is its value.
        
        If a dict of dicts, multiple solutions are specified for 
        calculation with the names of the input dicts.
    
    outputs : array-like or string
        A full output string or a list of output lines.
    """
    # inputs
    if isinstance(inputs, dict):
        inputs = [inputs]
    solutions = []
    for n, v in enumerate(inputs):
        solutions.append(make_solution(v, n))

    # outputs
    output = []
    if outputs is None:
        # if not specified, use default (defined at top ^)
        output.append(default_output)
    elif isinstance(outputs, str):
        # if it's a string
        output.append(outputs)
    else:
        # if it's a list
        output += outputs
    
    return '\n'.join(solutions) + '\n' + '\n'.join(output) + '\nEND'

def phreeqfind():
    """
    Return the path to the default libiphreeqc install location on your system.
    """
    system = platform.system()
    if system in ['Linux']:
        phreeq_path = '/usr/local/lib/libiphreeqc.so'
    elif system in ['Darwin']:
        phreeq_path = '/usr/local/lib/libiphreeqc.dylib'
    elif system in ['Windows']:
        return None
    
    if os.path.exists(phreeq_path):
        return phreeq_path
    else:
        raise ValueError(f'Cannot find libiphreeqc at the default loation ({phreeq_path}). Have you installed it, or do you have a custom install? If so, please specify `phreeq_path` manually to point at the correct location.')


def run_phreeqc(input_string, database=None, phreeq_path=None, output_file=False, parse_output=True):
    """
    Run input string in phreeqc with specified database.

    Parameters
    ----------
    input_string : str
        Valid phreeqc input string with SELECTED_OUTPUT.
    database : str
        Name of an included database to use (e.g. 'pitzer'), or 
        a complete path to a different phreeqc database (e.g. './path/to/pitzer.dat')
    phreeq_path : str
        Path to iphreeqc shared library. Defaults to '/usr/local/lib/libiphreeqc.so',
        which should work for standard installs on Linux machines

    Returns
    -------
    pandas.Series of calculated species
    """
    if phreeq_path is None:
        phreeq_path = phreeqfind()

    if database is None:
        print('No database specified  :  using pitzer')
        database = get_database_path()
    elif not os.path.exists(database):
        database = os.path.join(get_database_path(database))

    if not os.path.exists(database):
        raise ValueError(f"Can't phreeqc database: {database}\n   Please check that it exists.")

    phreeqc = phreeqc_mod.IPhreeqc(phreeq_path)
    phreeqc.load_database(database)
    if output_file:
        phreeqc.set_output_file_on()
    phreeqc.run_string(input_string)
    out = phreeqc.get_selected_output_array()
    if parse_output:
        return output_parser(out)
    else:
        return out

def output_parser(phreeqc_out):
    out = pd.DataFrame(phreeqc_out[1:], columns=phreeqc_out[0])

    interpreter = {
        'molality (mol/kgw)': re.compile('m_(.*)\(mol/kgw\)'),
        'total (mol/kgw)': re.compile('^([^m_]+)\(mol/kgw\)'),
        'log10(activity)': re.compile('la_(.*)'),
        'log10(saturation)': re.compile('si_(.*)')
    }

    index = []
    for c in out.columns:
        done = False
        for s, pat in interpreter.items():
            if pat.match(c):
                cclean = pat.findall(c)[0]
                index.append((s, cclean))
                done = True
        if not done:
            index.append(('general', c))
        
    out.columns = pd.MultiIndex.from_tuples(index)
    
    return out