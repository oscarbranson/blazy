import os
from glob import glob
import pandas as pd
import pkg_resources as pkgrs

def issubset(target, reference):
    """
    Returns True is all elements of target are present in reference.

    Both inputs must be iterable.
    """
    return set(target).issubset(set(reference))

def load_reference_data(file=None):
    if file is None:
        file = pkgrs.resource_filename('blazy', os.path.join('resources', 'test_data', 'SolutionKs.csv'))
    
    return pd.read_csv(file, comment='#')

def list_databases(show_headers=True, silent=False):
    dbases = glob(pkgrs.resource_filename('blazy', os.path.join('resources', 'database', '*.dat')))

    if not silent:
        for d in sorted(dbases):
            if show_headers:
                print(get_database_header(d))
            else:
                print(os.path.basename(d))
            print()
    
    return dbases

def get_database_header(d):
    out = [os.path.basename(d)]
    with open(d, 'r', errors='replace') as f:
        header = True
        while header:
            line = f.readline()
            if line.startswith('#'):
                out.append(f'   {line[1:].lstrip().rstrip()}')
            else:
                header = False
    return '\n'.join(out)