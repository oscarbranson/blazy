"""
Monte-Carlo uncertainties in phreeqc solution chemistry.
"""

# from .phreeq import input_str, run_phreeqc
# import uncertainties as un
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing as mp
from tqdm.autonotebook import tqdm
from .io import run_phreeqc
from . import istarmap

# Monte Carlo functions

def mc_input_dfs(df, N=1000, uncertainty_id='_std', distribution=None, outputs=None, db=None):
    if distribution is None:
        distribution = stats.norm

    cols = [c for c in df.columns if uncertainty_id not in c]  # all column names
    uncs = [c for c in df.columns if uncertainty_id in c]  # all uncertainty columns
    col_no_unc = [c for c in cols if c + uncertainty_id not in uncs]  # all columns without uncertainties
    col_with_unc = [c for c in cols if c not in col_no_unc]

    if not col_with_unc:
        raise ValueError(f"None of your columns contain uncertainties (have '{uncertainty_id}' in the name).")

    for i, r in df.iterrows():
        out = pd.DataFrame(columns=cols, index=range(N))

        for c in col_no_unc:
            out.loc[:, c] = r.loc[c]
        for c in col_with_unc:
            if all(~np.isnan([r.loc[c], r.loc[c + uncertainty_id]])):
                out.loc[:, c] = distribution(r.loc[c], r.loc[c + uncertainty_id]).rvs(N)
            else:
                out.loc[:, c] = r.loc[c]

        yield out, outputs, db

def concat_mc_results(mc_dfs):
    for i,o in enumerate(mc_dfs):
        o.index = pd.MultiIndex.from_product([[i], o.index], names=['sample', 'iteration'])
    return pd.concat(mc_dfs)

def make_and_run_input(inputs, outputs, db):
    input_str = db.generate_SOLUTIONS(inputs) + '\n' + outputs + '\nEND'
    out = run_phreeqc(input_str, parse_output=True, database=db.name)
    return out

def run_mc(inputs, N, database, targets=None, output_totals=True, output_molalities=True, output_activities=True, output_phases=True, phase_targets=None, allow_HCO_phases=True, drop_OH_species=True, uncertainty_id='_std', distribution=None):
    inputs = database.check_inputs(inputs, uncertainty_id=uncertainty_id)
    if targets is None:
        targets = database.get_target_elements(inputs, drop_OH=drop_OH_species, uncertainty_id=uncertainty_id)
    outputs = database.generate_SELECTED_OUTPUT(targets, totals=output_totals,
                                                molalities=output_molalities, activities=output_activities,
                                                phases=output_phases, phase_targets=phase_targets, 
                                                allow_HCO=allow_HCO_phases)
    with mp.Pool(mp.cpu_count()) as pool:
        out = list(tqdm(pool.istarmap(make_and_run_input, mc_input_dfs(df=inputs, N=N, uncertainty_id=uncertainty_id, 
                                                                       distribution=distribution, outputs=outputs, db=database)), total=len(inputs), desc='Running MC'))

    return concat_mc_results(out)

def calc_mc_quantiles(mc_output, CI=0.95, quantiles=None):
    if quantiles is None:
        quantiles = [0.5 - CI / 2, 0.5, 0.5 + CI / 2]
    out = []
    for i, g in mc_output.groupby(level=0):
        out.append(g.quantile(quantiles))
    N = len(out)
    out = pd.concat(out)
    out.index = pd.MultiIndex.from_product([range(N), quantiles], names=['sample', 'quantile'])
    
    return out

# def mc_input_str(input_dict, N, outputs=None):
#     """
#     Generates phreeqc input string for Monte-Carlo uncertainties
    
#     Parameters
#     ----------
#     input_dict : dict 
#         Containing `{entry: value}` pairs. The `value` may be either:
#             - string, float or int : the value will be used for all iterations.
#             - tuple : the first value will be taken as the mean, the second as 
#               the standard deviation
#             - uncertainties.core.Variable object : the nominal_value and std_dev
#               will be used.
#             - A scipy `rv_frozen` distribution object
#             - a custom object that contains a `.rvs()` method to generate a draw
#               for each iteration.
#         In practice, all of the former are converted into objects with a 'rvs()'
#         method before performing the iteration.
#     N : int
#         The number of monte-carlo iterations to generate.
#     outputs : str or list
#         a complete phreeqc output string, or a list of lines of an output string.

#     Returns
#     -------
#     generator : where each iteration yields a new random dict drawn from the inputs.
#     """

#     return input_str(mc_input_dicts(input_dict=input_dict, N=N, outputs=outputs))

# class dummy_str:
#     def __init__(self, string):
#         self.string = string
    
#     def rvs(self):
#         return self.string

# class dummy_numeric:
#     def __init__(self, num):
#         self.num = num
    
#     def rvs(self):
#         return self.num

# def mc_input_dicts(input_dict, N):
#     """
#     Generates phreeqc input dicts for Monte-Carlo uncertainties
    
#     Parameters
#     ----------
#     input_dict : dict 
#         Containing `{entry: value}` pairs. The `value` may be either:
#             - string, float or int : the value will be used for all iterations.
#             - tuple : the first value will be taken as the mean, the second as 
#               the standard deviation
#             - uncertainties.core.Variable object : the nominal_value and std_dev
#               will be used.
#             - A scipy `rv_frozen` distribution object
#             - a custom object that contains a `.rvs()` method to generate a draw
#               for each iteration.
#         In practice, all of the former are converted into objects with a 'rvs()'
#         method before performing the iteration.
#     N : int
#         The number of monte-carlo iterations to generate.

#     Returns
#     -------
#     generator : where each iteration yields a new random dict drawn from the inputs.
#     """
    
#     # construct input dict of objects with .rvs() methods
#     dists = {}
#     for k, v in input_dict.items():
#         if isinstance(v, str):
#             dists[k] = dummy_str(v)
#         elif isinstance(v, (float, int)):
#             dists[k] = dummy_numeric(v)
#         elif isinstance(v, un.core.Variable):
#             dists[k] = stats.norm(v.nominal_value, v.std_dev)
#         elif isinstance(v, tuple):
#             dists[k] = stats.norm(v[0], v[1])
#         elif hasattr(v, 'rvs'):
#             dists[k] = v
#         else:
#             raise ValueError('Entry for {k} is invalid. See function doc for valid entry types.')
    
#     for i in range(N):
#         yield {k: v.rvs() for k, v in dists.items()}

