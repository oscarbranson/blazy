"""
Functions for parsing iPHREEQC databases.
"""

import re
import os
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
        
        self.load(self.path)
        self.find_sections()
    
    def _database_path_handler(self, database):
        if os.path.exists(database):
            return database
        
        dbase_path = pkgrs.resource_filename('blazy', os.path.join('resources','database'))
        valid_dbases = glob(os.path.join(dbase_path, '*.dat'))
        dbase_names = sorted([os.path.basename(d).split('.')[0] for d in valid_dbases])

        if database in dbase_names:
            return pkgrs.resource_filename('blazy', os.path.join('resources','database', database + '.dat'))
        
        raise ValueError(f"The database '{database}' does not exist. Please provide a complete path to a PHREEQC database, or use one of: [{', '.join(dbase_names)}]")
    
    def _targets_handler(self, targets):
        if isinstance(targets, str):
            targets = [targets]
        return set(targets)

    def load(self, database, keep_comments=False):
        with open(database, 'r', errors='replace') as f:
            db = [line.rstrip() for line in f]
        db= [line for line in db if line != '']  # remove empty lines
        if not keep_comments:
            db = [line for line in db if line[0] != '#']
        self.db = db

    def find_sections(self):
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
        self.sections = headind
    
    @staticmethod
    def _reacsplit(reac):
        """
        Splits a reaction into its components.
        
        For example:
        self._reacsplit('Ca+2 + B(OH)3 + H2O = CaB(OH)4+ + H+')
        > [['Ca+2', 'B(OH)3', 'H2O'], ['CaB(OH)4+', 'H+']]
        
        list of reactants in the first item, list of products in the second
        
        """
        sides = reac.split('=')
        return [[c for c in side.split(' ') if c not in ['+', '']] for side in sides]

    def get_section(self, section, remove_comments=True):
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
        
        targets = self._targets_handler(targets=targets)
        
        active_lines = set()
        for t in targets:
            for line in self.get_section(section=section, remove_comments=remove_comments):
                if t in line:
                    active_lines.add(line)
                
        return active_lines
    
    def get_SOLUTION_SPECIES(self, targets, allow_HCO=True):
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
    
    def get_SOLUTION_MASTER_SPECIES(self, targets):
        solution_species = self.get_species(targets, 'SOLUTION_MASTER_SPECIES')
        targets = self._targets_handler(targets)

        out_species = set()
        for s in solution_species:
            pels = get_elements(s.split()[0])
            if issubset(pels, targets):
                out_species.add(s.split()[0])
        
        return out_species

    def parse_PHASES(self):
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

    def generate_SELECTED_OUTPUT(self, targets, totals=True, molalities=True, activities=True, phases=True, phase_targets=None, allow_HCO=True):
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