"""
Functions for parsing iPHREEQC databases.
"""

import re
import os
import itertools
import warnings
import pandas as pd
from glob import glob
import pkg_resources as pkgrs

from ..helpers import issubset, get_database_header
from ..chemistry import get_elements, valid_elements
from .io import make_solution

# make sure warnings are always shown
warnings.filterwarnings('always', category=UserWarning)

def reacsplit(reac):
    """
    Splits a reaction into its components.
    
    For example:
    reacsplit('Ca+2 + B(OH)3 + H2O = CaB(OH)4+ + H+')
    > [['Ca+2', 'B(OH)3', 'H2O'], ['CaB(OH)4+', 'H+']]
    
    list of reactants in the first item, list of products in the second
    
    Parameters
    ----------
    reac : str
        The reaction to be separated, in phreeqc notation.

    Returns
    -------
    list : containing [[left, hand, side], [right, hand, side]]
    """
    sides = reac.split('=')
    return [[c for c in side.split(' ') if c not in ['+', '']] for side in sides]

def remove_stoich(species):
    """
    Removes stoichiometric multipler from the start of a species name.

    For example:
    2B(OH)3  ->  B(OH)3

    Parameters
    ----------
    species : str
        The species to be treated
    
    Returns
    -------
    str : the input species with stoichiometric multiplier removed.
    """
    srm = re.compile('^[0-9]+')
    return srm.sub('', species)


class datParser:
    """
    Class for loading and parsing PHREEQC databases.

    Particularly useful for identifying outputs that are
    relevant to your solution constituents.

    Parameters
    ----------
    database : path
        Name of database or path to phreeqc database.
    """
    def __init__(self, database, silent=False):
        """
        Class for loading and parsing PHREEQC databases.

        Particularly useful for identifying outputs that are
        relevant to your solution constituents.

        Parameters
        ----------
        database : path
            Name of database or path to phreeqc database.
        """
        self.path = self._database_path_handler(database)
        if not silent:
            print('Using ' + get_database_header(self.path))
        self.name = os.path.basename(database).split('.')[0]
        
        self.db = self.load(self.path)
        self.sections = self.find_sections()

        self.get_SOLUTION_MASTER_SPECIES()

        self._exempt_inputs = [
            'temperature', 'temp', 
            'pressure', 'press',
            'pH', 'pe', 
            'redox', 'O(0)', '-water'
            'density', 'unit', 'units'
        ]

    def _database_path_handler(self, database):
        """
        Convenience function for checking/getting the path to the database.

        Parameters
        ----------
        database : str
            Either path to a database, or the name of a database

        Returns
        -------
        str : Checked path to the database.
        """
        if os.path.exists(database):
            return database
        else:
            database = database.replace('.dat', '')
        
        dbase_path = pkgrs.resource_filename('blazy', os.path.join('resources','database'))
        valid_dbases = glob(os.path.join(dbase_path, '*.dat'))
        dbase_names = sorted([os.path.basename(d).split('.')[0] for d in valid_dbases])

        if database in dbase_names:
            return pkgrs.resource_filename('blazy', os.path.join('resources','database', database.replace('.dat','') + '.dat'))
        
        raise ValueError(f"The database '{database}' does not exist. Please provide a complete path to a PHREEQC database, or use one of: [{', '.join(dbase_names)}]")
    
    def _targets_handler(self, targets):
        """
        Convenience function for checking the format of 'targets'
        """
        if targets is None:
            return targets
        if isinstance(targets, str):
            targets = [targets]
        return set(targets)

    def load(self, database, keep_comments=False):
        """
        Load selected database into memory.

        Parameters
        ----------
        database : str
            Path to valid PHREEQC database (.dat).
        keep_comments : bool
            If False, all comments (lines starting with #) are removed.

        Returns
        -------
        list : All lines in the database.
        """
        with open(database, 'r', errors='replace') as f:
            db = [line.rstrip() for line in f]
        db= [line for line in db if line != '']  # remove empty lines
        if not keep_comments:
            db = [line for line in db if line[0] != '#']
        return db

    def find_sections(self):
        """
        Finds all section headings in database and records them with line numbers.

        Result

        Parameters
        ----------
        None

        Returns
        -------
        dict : containing {SECTION: (startline, endline)} entries
        """
        sectionhead = re.compile('[A-Z_]+$')
        headind = {}
        ilast = None
        headlast= None
        for i, line in enumerate(self.db):
            if sectionhead.match(line):
                if ilast is not None:
                    headind[headlast] = (ilast, i)
                headlast = line.strip()
                ilast = i        
        return headind

    def get_section(self, section, remove_comments=True):
        """

        Parameters
        ----------
        section : str
            The name of the section in the database, for example 'SOLUTION_SPECIES'.
        remove_comments : bool
            If True, in-line comments (anything preceeding '#') are removed.

        Returns
        -------
        generator : yields each line of the section
        """
        if section not in self.sections:
            secnames = ', '.join(self.sections.keys())
            raise KeyError(f"'{section}' is not a valid section. Please choose one of: {secnames}")
        lines = self.sections[section]

        for i in range(lines[0]+1, lines[1]):
            if remove_comments:
                yield self.db[i].split('#')[0]
            else:
                yield self.db[i]
    
    def parse_SOLUTION_SPECIES(self, remove_comments=True):
        """
        Turns tabbed database entries into a dict of {entry: [lines]}.
        
        Each entry is a title, followed by a number of tab-inset lines that
        are associated with that entry.
        
        Parameters
        ----------
        section : str
            The name of the section in the database, for example 'SOLUTION_SPECIES'.
        remove_comments : bool
            If True, in-line comments (anything preceeding '#') are removed.

        Returns
        -------
        dict : containing each entry in the section and the lines it contains.
        
        """
        section = 'SOLUTION_SPECIES'
        out = {}
        entry = None
        for p in self.get_section(section=section, remove_comments=remove_comments):
            if '=' in p:
                entry = p.lstrip()
                out[entry] = []
            else:
                out[entry].append(p.lstrip())
        
        return out

    # def parse_entries(self, section='SOLUTION_SPECIES', remove_comments=True):
    #     """
    #     Turns tabbed database entries into a dict of {entry: [lines]}.
        
    #     Each entry is a title, followed by a number of tab-inset lines that
    #     are associated with that entry.
        
    #     Parameters
    #     ----------
    #     section : str
    #         The name of the section in the database, for example 'SOLUTION_SPECIES'.
    #     remove_comments : bool
    #         If True, in-line comments (anything preceeding '#') are removed.

    #     Returns
    #     -------
    #     dict : containing each entry in the section and the lines it contains.
        
    #     """
    #     out = {}
    #     entry = None
    #     for p in self.get_section(section=section, remove_comments=remove_comments):
    #         if not p.startswith(('\t', ' ')):
    #             entry = p
    #             out[entry] = []
    #         else:
    #             out[entry].append(p.lstrip())
        
    #     return out

    def get_SOLUTION_MASTER_SPECIES(self):
        """
        Generate dictionaries of valid 'element' and 'master species' names.

        Creates four dictionaries:
        element_2_master : converts from 'element name' to 'master species'
        master_2_element : converts from 'master species' to 'element name'
        master_nocharge_2_element : converts from 'master species' with the valence information removed to 'element name'
        element_2_master_nocharge : converts from 'element name' to 'master species' with the valence information removed

        Returns
        -------
        None
        """
        solution_species = self.get_section('SOLUTION_MASTER_SPECIES', remove_comments=True)
        self.master_species_table = {}
        for s in solution_species:
            sp = s.split()
            self.master_species_table[sp[0]] = sp[1:]

        novalence = re.compile('[+-][0-9]?')
        self.element_2_master = {k: v[0] for k, v in self.master_species_table.items()}
        self.element_2_master_nocharge = {k: novalence.sub('', v) for k, v in self.element_2_master.items()}
        self.master_2_element = {v: k for k, v in self.element_2_master.items()}
        self.master_nocharge_2_element = {novalence.sub('', v): k for k, v in self.element_2_master.items()}
    
    def list_valid_species(self):
        """
        Prints a list of valid master species in the database.
        """
        headers = ['Species Name', 'Species', 'Alk', 'Gram Formula', 'Weight']
        pad = 3
        L0 = len(headers[0])
        Ls = [len(h) for h in headers[1:]]

        for k, v in self.master_species_table.items():
            if len(k) > L0:
                L0 = len(k)
            for i, vi in enumerate(v):
                if len(vi) > Ls[i]:
                    Ls[i] = len(vi)
        
        print(f'{headers[0]:>{L0}}' + ''.join([f'{h:>{L + pad}}' for h, L in zip(headers[1:], Ls)]))
        for k, v in self.master_species_table.items():
            print(f'{k:>{L0}}' + ''.join([f'{h:>{L + pad}}' for h, L in zip(v, Ls)]))
        
        print('\nFor more info, see https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-57.html')

    def get_SOLUTION_SPECIES(self, targets=None, include_foreign=True):
        """
        Returns a list of solution species present in the database that are relevant to the target elements.

        Parameters
        ----------
        targets : str or list
            An element or list of elements of interest. If None, all species are returned.
        include_foreign : bool
            If True, include species containing any of targets *and* other elements
            not in targets.

        Returns
        -------
        set : calculated species in the database containing the target elements.
        """
        targets = self._targets_handler(targets=targets)
        solution_species = self.parse_SOLUTION_SPECIES(remove_comments=True)
        srm = re.compile('^[0-9]+')  # pattern for removing stoichiometric multiplier from speciues
        out_species = set()

        if targets is None:
            for s in solution_species.keys():
                _, prod = reacsplit(s)  # only look at reaction products
                for p in prod:
                    out_species.add(srm.sub('', p))
        else:
            for s in solution_species.keys():
                _, prod = reacsplit(s)  # only look at reaction products
                for p in prod:
                    pels = get_elements(p)
                    if include_foreign:
                        if pels.intersection(targets):
                            out_species.add(srm.sub('', p))
                    else:
                        if issubset(pels, targets) and pels.intersection(targets):
                            out_species.add(srm.sub('', p))

        # this isn't quite right - should get all possible species given solution compositions,
        # *then* filter by targets. At the moment returns species containing irrelevant elements.
        
        return out_species.difference(['', None])

    def parse_PHASES(self, remove_comments=True):
        """
        Turns tabbed database entries into a dict of {entry: [lines]}.
        
        Each entry is a title, followed by a number of tab-inset lines that
        are associated with that entry.
        
        Parameters
        ----------
        section : str
            The name of the section in the database, for example 'SOLUTION_SPECIES'.
        remove_comments : bool
            If True, in-line comments (anything preceeding '#') are removed.

        Returns
        -------
        dict : containing each entry in the section and the lines it contains.
        
        """
        section = 'PHASES'
        out = {}
        entry = None
        for p in self.get_section(section=section, remove_comments=remove_comments):
            if not p.startswith(('\t', ' ')):
                entry = p.split()[0]  # split because some have an additional number in the database
                out[entry] = []
            else:
                out[entry].append(p.lstrip())
        
        return out

    def get_PHASES(self, targets=None, allow_HCO=True):
        """
        Returns all phases that contain target elements.

        Parameters
        ----------
        targets : str or list
            An element or list of elements of interest. If None, all phases are returned.
        allow_HCO : bool
            If True, also return species that contain the target elements
            and one or more of H, C and O.

        Returns
        -------
        set : calculated phases in the database containing the target elements.
        """
        
        self.phases = self.parse_PHASES()

        targets = self._targets_handler(targets)

        if targets is None:
            return self.phases

        possible_phases = set()
        for p, i in self.phases.items():
            for t in targets:
                phase_formula = reacsplit(i[0])[0][0]
                if t in phase_formula:
                    possible_phases.add((p, phase_formula))
        
        if allow_HCO:
            ftargets = targets.union(['H', 'C', 'O'])
        out_phases = set()
        
        for p, f in possible_phases:
            if issubset(get_elements(f), ftargets):
                out_phases.add(p)

        
        return out_phases
    
    def list_valid_phases(self):
        """
        Prints a list of valid phases in the database.
        """
        headers = ['Phase', 'Formula']
        pad = 3
        L0 = len(headers[0])
        Ls = [len(h) for h in headers[1:]]
        
        phases = self.get_PHASES()
        
        L = max([L0] + [len(k) for k in phases]) + 3
        
        print(f'{headers[0]:>{L}}' + ''.join([f'{h:>{L + pad}}' for h, L in zip(headers[1:], Ls)]))
        for k, v in phases.items():
            formula = v[0].split(' =') [0]
            print(f'{k:>{L}}' +  f'   {formula:>{pad}}')

    def check_inputs(self, inputs, remove_failures=True, uncertainty_id='_std'):
        """
        Checks the validity of an input dictionary, replacing names where necessary.

        Parameters
        ----------
        input_dict : dict
            A dictionary used to create an input string, consisting of
            {element_name: value}.

        Returns
        -------
        dict : a checked input dict with any required modification.
        """
        if isinstance(inputs, dict):
            inputs = pd.DataFrame(inputs, index=[0])
        elif isinstance(inputs, list):
            inputs = pd.DataFrame(index=[0], columns=inputs)


        retain = []
        final_names = []
        msg_subs = []
        msg_rems = []
        changed = False
        for k in inputs.columns:
            k0 = k
            if uncertainty_id in k:
                k = k.replace(uncertainty_id, '')
            if k in self._exempt_inputs or (k[0] == '-') or k in self.element_2_master:
                # these are generic options that won't be in the database
                # there are a lot of abbreviated options that start with a dash
                retain.append(k0)
                final_names.append(k0)
            else:
                if k in self.master_2_element:
                    retain.append(k0)
                    kn = k0.replace(k, self.master_2_element[k])
                    final_names.append(kn)
                    msg_subs.append(f"   - {k0} --> {kn}")
                    changed = True
                elif k in self.master_nocharge_2_element:
                    retain.append(k0)
                    kn = k0.replace(k, self.master_nocharge_2_element[k])
                    final_names.append(kn)
                    msg_subs.append(f"   - {k0} --> {kn}")
                    changed = True
                elif remove_failures:
                    msg_rems.append(f"   - {k0}")
                    changed = True
                else:
                    raise ValueError(f"{k} is not a valid element or species name for the {self.name} database.")
        if changed:
            if len(msg_rems + msg_subs) > 0:
                msg = (f"\n\nThere were columns in your inputs which aren't valid in the {self.name} database." + 
                       "\nWe have either tried to substitute them with valid inputs or removed them.")
                if len(msg_subs) > 0:
                    msg += f"\nInvalid columns which we were able to substitute:\n" + '\n'.join(msg_subs)
                if len(msg_rems) > 0:
                    msg += f"\nInvalid columns which we have removed:\n" + '\n'.join(msg_rems)
                msg += (
                    "\nPlease make sure these substitions/removals make sense!" + 
                    "\nIf they don't please manually specify links between column names\nand valid species in your chosen database." +
                    "\n\nCreate a lookup dictionary linking column names to valid species names {'column_name': 'species_name'}," + 
                    "\nthen prepare your data for input using the `.select_inputs()` function." +
                    f"\n\nTo see a list of valid species names for {self.name}.dat, use the `.list_valid_species()` function."
                )

                warnings.warn(msg)

            inputs = inputs.loc[:, retain]
            inputs.columns = final_names

        return inputs

    def select_inputs(self, inputs, column_lookup, uncertainty_id='_std'):
        select_columns = []
        colnames = []
        for c in column_lookup.keys():
            select_columns.append(c)
            colnames.append(column_lookup[c])
            if c + uncertainty_id in inputs.columns:
                select_columns.append(c + uncertainty_id)
                colnames.append(column_lookup[c] + uncertainty_id)

        selected = inputs.loc[:, select_columns]
        selected.columns = colnames

        selected = self.check_inputs(selected, uncertainty_id=uncertainty_id)
        
        return selected

    def get_target_elements(self, inputs, drop_OH=True, uncertainty_id='_std'):
        """
        Gets a list of all elements in the input file.

        Parameters
        ----------
        input_dict : dict
            A dictionary used to create an input string, consisting of
            {element_name: value}.
        
        Returns
        -------
        set : containing names of all elements in the input
        """
        inputs = self.check_inputs(inputs)
        
        targets = set()
        for k in inputs.columns:
            if k in self._exempt_inputs or (k[0] == '-') or uncertainty_id in k:
                continue
            targets.update(get_elements(self.element_2_master_nocharge[k]))
        if drop_OH:
            targets = targets.difference({'H', 'O'})  # get rid of OH ions

        return targets.intersection(valid_elements)


    def generate_SELECTED_OUTPUT(self, targets, totals=True, molalities=True, activities=True, phases=True, phase_targets=None, allow_HCO=True):
        """
        Generates a SELECTED_OUTPUT string containing all relevant species and phases.

        Parameters
        ----------
        targets : str or list
            An element or list of elements of interest.
        totals : bool
            If True, the totals of the input elements are included in the output.
        molalities : bool
            If True, the molalities of all viable species are included in the output. 
        activities : bool
            If True, the log10(activities) of all viable species are included in the output.
        phases : bool
            If True, the log10(saturatio) of all viable phases are included in the output.
        phase_targets : str or list
            An element or list of elements of interest for PHASES, if different from `targets`.
        allow_HCO : bool
            If True, also return species that contain the target elements
            and one or more of H, C and O.

        Returns
        -------
        str : SELECTED_OUTPUT string for including in PHREEQC input.
        """
        
        if phase_targets is None:
            phase_targets = targets

        outstr = [
            'SELECTED_OUTPUT', 
            '    -pH',
            '    -temperature',
            '    -alkalinity',
            '    -ionic_strength'
        ]
        
        if totals:
            outstr.append('    -totals ' + ' '.join(targets))
        
        if molalities or activities:
            species = self.get_SOLUTION_SPECIES(targets).union(['H+', 'OH-'])
        
        if molalities:
            outstr.append('    -m ' + ' '.join(species))
        
        if activities:
            outstr.append('    -a ' + ' '.join(species))
        
        if phases:
            phases = self.get_PHASES(phase_targets, allow_HCO=allow_HCO)
            outstr.append('    -si ' + ' '.join(phases))
        
        return '\n'.join(outstr)
    
    def generate_SOLUTIONS(self, inputs):
        inputs = self.check_inputs(inputs)

        solutions = []
        for n, v in inputs.iterrows():
            solutions.append(make_solution(v, n))
        return solutions
    
    def add_EQUILIBRIUM_PHASES(self, phases):
        """
        Add EQUILIBRIUM_PHASES to the PHREEQC input.

        Parameters
        ----------
        phases : list of tuples
            Of [(name, target_log10SI, *options)] of the phases to be added to the
            input, where *options are valid entries for EQUILIBRIUM_PHASES lines.
        """
        if not isinstance(phases[0][0], str):
            return [self.add_EQUILIBRIUM_PHASES(p) for p in phases]

        valid_phases = []
        for p in phases:
            if p[0] not in self.get_PHASES():
                print(f'{p[0]} is not in the databse; removed.')
            else:
                valid_phases.append(map(str,p))
        
        lines = ['  '.join(p) for p in valid_phases]
        
        return 'EQUILIBRIUM_PHASES\n   ' + '\n   '.join(lines)



    def make_PHREEQC_input(self, inputs, targets=None, output_totals=True, output_molalities=True, output_activities=True, output_phases=True, phase_targets=None, allow_HCO_phases=True, drop_OH_species=True, uncertainty_id='_std', equilibrium_phases=None):
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
        database : str
            Which database the input_string is for. Required for input
            checking and automatic output generation.
        equilibrium_phases : list of tuples
            Of [(name, target_log10SI, *options)] of the phases to be added to the
            input, where *options are valid entries for EQUILIBRIUM_PHASES lines.
            Can also be a list containing lists of tuples, where each entry in that
            list is applied to each solution separately.
        """
        inputs = self.check_inputs(inputs)
        
        if targets is None:
            targets = self.get_target_elements(inputs, drop_OH=drop_OH_species, uncertainty_id=uncertainty_id)

        # inputs
        solutions = self.generate_SOLUTIONS(inputs)

        # outputs
        output = self.generate_SELECTED_OUTPUT(targets, totals=output_totals, molalities=output_molalities, activities=output_activities, phases=output_phases, phase_targets=phase_targets, allow_HCO=allow_HCO_phases)
        
        if equilibrium_phases is None:
            self.input_str = '\n'.join(solutions) + '\n' + output + '\nEND'
        else:
            eqps = self.add_EQUILIBRIUM_PHASES(equilibrium_phases)

            if isinstance(eqps, str):
                self.input_str = '\n'.join([s + '\n' + eqps + '\n' + output + '\nEND' for s in solutions])
            else:
                inp = []
                for sol, eq in itertools.product(solutions, eqps):
                    inp.append(sol + '\n' + eq + '\n' + output  + '\nEND')
                self.input_str = '\n'.join(inp)
    
        return self.input_str