"""
The main user interaction class
"""
import numpy as np
import pandas as pd
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

from .parser import datParser
from ..chemistry import get_elements
from .io import phreeqfind, output_parser
from .montecarlo import run_mc

class iphreeqc:
    def __init__(self, database='pitzer', iphreeqc_path=None):
        self.db = datParser(database)
        
        if iphreeqc_path is None:
            iphreeqc_path = phreeqfind()
        self.iphreeqc_path = iphreeqc_path

        self.make_PHREEQC_input = self.db.make_PHREEQC_input

        # if isinstance(solutions, dict):
        #     self.solutions = pd.Series(solutions)
        # else:
        #     self.solutions = solutions
    
    def _spawn(self, output_file=False):
        self.phreeqc = phreeqc_mod.IPhreeqc(self.iphreeqc_path)
    
    def _load_database(self):
        self.phreeqc.load_database(self.db.path)
    
    def _run(self, input_string):
        self.phreeqc.run_string(input_string)

    def _kill(self):
        self.phreeqc.destroy_iphreeqc()
        del self.phreeqc

    def _getoutput(self):
        return self.phreeqc.get_selected_output_array()

    def run_phreeqc(self, input_string, keepalive=False):
        if not hasattr(self, 'phreeqc'):
            self._spawn()
            self._load_database()

        self._run(input_string)

        out = output_parser(self._getoutput()).replace(-999.999, np.nan)

        if not keepalive:
            self._kill()
        
        return out

    def change_database(self, database):
        self.db = datParser(database)
        if hasattr(self, 'phreeqc'):
            self._load_database()

    # def make_input_string(self, inputs, targets=None, output_totals=True, output_molalities=True, output_activities=True, output_phases=True, phase_targets=None, allow_HCO_phases=True, drop_OH_species=True, uncertainty_id='_std'):
        # return self.db.make_PHREEQC_input(inputs=inputs, targets=targets, output_totals=output_totals, output_molalities=output_molalities, output_activities=output_activities, output_phases=output_phases, phase_targets=phase_targets, allow_HCO_phases=allow_HCO_phases, drop_OH_species=drop_OH_species, uncertainty_id=uncertainty_id)
    
    def run(self, inputs, targets=None, output_totals=True, output_molalities=True, output_activities=True, output_phases=True, phase_targets=None, allow_HCO_phases=True, drop_OH_species=True, uncertainty_id='_std'):
        inputs = self.db.check_inputs(inputs, uncertainty_id=uncertainty_id)
        
        inputs = inputs.loc[:, [c for c in inputs.columns if uncertainty_id not in c]]

        input_string = self.db.make_PHREEQC_input(inputs=inputs, targets=targets, output_totals=output_totals, output_molalities=output_molalities, output_activities=output_activities, output_phases=output_phases, phase_targets=phase_targets, allow_HCO_phases=allow_HCO_phases, drop_OH_species=drop_OH_species, uncertainty_id=uncertainty_id)

        return self.run_phreeqc(input_string)

    def run_mc(self, inputs, N, targets=None, output_totals=True, output_molalities=True, output_activities=True, output_phases=True, phase_targets=None, allow_HCO_phases=True, drop_OH_species=True, uncertainty_id='_std', distribution=None):
        return run_mc(inputs=inputs, N=N, database=self.db, targets=targets, output_totals=output_totals, output_molalities=output_molalities, output_activities=output_activities, output_phases=output_phases, phase_targets=phase_targets, allow_HCO_phases=allow_HCO_phases, drop_OH_species=drop_OH_species, uncertainty_id=uncertainty_id)

    def list_valid_species(self):
        print(f'Valid Species for {self.db.name}.dat:\n')
        self.db.list_valid_species()
        print(
            "\nCreate a lookup dictionary linking column names to valid species names {'column_name': 'species_name'}," + 
            "\nthen prepare your data for input using the `.select_inputs()` function."
        )

    def select_inputs(self, inputs, column_lookup, uncertainty_id='_std'):
        return self.db.select_inputs(inputs=inputs, column_lookup=column_lookup, uncertainty_id=uncertainty_id)