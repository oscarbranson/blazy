"""
Functions for parsing iPHREEQC databases.
"""

import re
import os
import warnings
from glob import glob
import pkg_resources as pkgrs

from ..helpers import issubset
from ..chemistry import get_elements

class database:
    def __init__(self, database):
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
        self.name = os.path.basename(database).split('.')[0]
        
        self.db = self.load(self.path)
        self.sections = self.find_sections()

        self.get_SOLUTION_MASTER_SPECIES()

        self._exempt_inputs = [
            'temperature', 'temp', 
            'pressure', 'press',
            'pH', 'pe', 
            'redox'
            'density', 'unit'
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
    
    @staticmethod
    def _reacsplit(reac):
        """
        Splits a reaction into its components.
        
        For example:
        self._reacsplit('Ca+2 + B(OH)3 + H2O = CaB(OH)4+ + H+')
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

    def get_species(self, targets, section='SOLUTION_SPECIES', remove_comments=True):
        """
        Returns all species that containing the elements in 'target'.

        Parameters
        ----------
        targets : str or array-like
        section : str
        remove_comments : bool

        Returns
        -------
        set : all of the lines that contain the target elements.

        """
        
        targets = self._targets_handler(targets=targets)
        
        active_lines = set()
        for t in targets:
            for line in self.get_section(section=section, remove_comments=remove_comments):
                if t in line:
                    active_lines.add(line)
                
        return active_lines

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
        solution_species = self.get_section('SOLUTION_MASTER_SPECIES')
        out_species = {}
        for s in solution_species:
            sp = s.split()
            out_species[sp[0]] = sp[1]

        novalence = re.compile('[+-][0-9]?')
        self.element_2_master = out_species
        self.element_2_master_nocharge = {k: novalence.sub('', v) for k, v in self.element_2_master.items()}
        self.master_2_element = {v: k for k, v in self.element_2_master.items()}
        self.master_nocharge_2_element = {novalence.sub('', v): k for k, v in self.element_2_master.items()}

    def get_SOLUTION_SPECIES(self, targets, allow_HCO=True):
        """
        Returns a list of solution species present in the database that are relevant to the target elements.

        Parameters
        ----------
        targets : str or list
            An element or list of elements of interest.
        allow_HCO : bool
            If True, also return species that contain the target elements
            and one or more of H, C and O.

        Returns
        -------
        list : calculated species in the database containing the target elements.
        """
        solution_species = self.get_species(targets, 'SOLUTION_SPECIES')
        targets = self._targets_handler(targets=targets)
        if allow_HCO:
            ftargets = targets.union(['H', 'C', 'O'])
        else:
            ftargets = targets

        out_species = set()
        for s in solution_species:
            _, prod = self._reacsplit(s)  # only look at reaction products
            for p in prod:
                pels = get_elements(p)
                if issubset(pels, ftargets) and pels.intersection(targets):
                    out_species.add(p)
        
        return out_species
    
    def parse_PHASES(self):
        """
        Parse the PHASES part of the database into a dictionary

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        phases = {}
        phase = ''
        for p in self.get_section('PHASES'):
            if '\t' not in p:  # phase names don't have tabs
                phase = p
                phases[p] = []
            else:
                phases[phase].append(p.lstrip())
        
        self.phases = phases

    def get_PHASES(self, targets, allow_HCO=True):
        """
        Returns all phases that contain target elements.

        Parameters
        ----------
        targets : str or list
            An element or list of elements of interest.
        allow_HCO : bool
            If True, also return species that contain the target elements
            and one or more of H, C and O.

        Returns
        -------
        set : calculated phases in the database containing the target elements.
        """
        
        if not hasattr(self, 'phases'):
            self.parse_PHASES()

        targets = self._targets_handler(targets)

        possible_phases = set()
        for p, i in self.phases.items():
            for t in targets:
                phase_formula = self._reacsplit(i[0])[0][0]
                if t in phase_formula:
                    possible_phases.add((p, phase_formula))
        
        if allow_HCO:
            ftargets = targets.union(['H', 'C', 'O'])
        out_phases = set()
        
        for p, f in possible_phases:
            if issubset(get_elements(f), ftargets):
                out_phases.add(p)
        
        return out_phases

    def check_input_dict(self, input_dict):
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
        for k, v in input_dict.items():
            if k in self._exempt_inputs:
                # these are generic options that won't be in the database
                continue
            if k[0] == '-':
                # there are a lot of abbreviated options that start with a dash
                continue
            if k not in self.element_2_master:
                kn = None
                if k in self.master_2_element:
                    kn = self.master_2_element[k]
                elif k in self.master_nocharge_2_element:
                    kn = self.master_nocharge_2_element[k]
                else:
                    raise ValueError(f"{k} is not a valid element or species name for the {self.name} database.")
                msg = f"{k} is not a valid input for the {self.name} database, but {kn} is. We've swapped {k} for {kn}."
                warnings.warn(msg)
                if kn is not None:
                    input_dict.pop(k)
                    input_dict[kn] = v
        
        return input_dict

    def get_target_elements(self, input_dict):
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
        input_dict = self.check_input_dict(input_dict)
        targets = set()
        for k, v in input_dict.items():
            if k in self._exempt_inputs:
                continue
            targets.update(get_elements(self.element_2_master_nocharge[k]))
        
        return targets


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
            species = self.get_SOLUTION_SPECIES(targets)
        
        if molalities:
            outstr.append('    -m ' + ' '.join(species))
        
        if activities:
            outstr.append('    -a ' + ' '.join(species))
        
        if phases:
            phases = self.get_PHASES(phase_targets, allow_HCO=allow_HCO)
            outstr.append('    -si ' + ' '.join(phases))
        
        return '\n'.join(outstr)