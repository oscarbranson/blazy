import os
import unittest
from glob import glob
from importlib import resources

from blazy.phreeqc.parser import datParser

class TestDatabaseInteraction(unittest.TestCase):

    def test_file_paths(self):
        # where the datatbase should be
        with resources.as_file(resources.files('blazy').joinpath('resources', 'database', 'pitzer.dat')) as dbase_path:
            dbase_path_str = str(dbase_path)

            # where the database is
            db = datParser(database='pitzer')
            self.assertEqual(db.path, dbase_path_str)

            db = datParser(database='pitzer.dat')
            self.assertEqual(db.path, dbase_path_str)

            db = datParser(database=dbase_path_str)
            self.assertEqual(db.path, dbase_path_str)

    def test_all_databases(self):
        with resources.as_file(resources.files('blazy').joinpath('resources', 'database')) as dbase_path:
            databases = glob(os.path.join(str(dbase_path), '*.dat'))
            print()

            for d in databases:
                print(f'\nChecking {os.path.basename(d)}')
            db = datParser(d)
            print(f'  -> {db.path}')

            print('  -> get_SOLUTION_MASTER_SPECIES')
            db.get_SOLUTION_MASTER_SPECIES()
            print('       ' + ', '.join(db.element_2_master.keys()))

            print('  -> get_target_elements')
            targets = db.get_target_elements(db.master_nocharge_2_element)
            print('       ' + ', '.join(targets))

            print('  -> get_SOLUTION_SPECIES')
            species = db.get_SOLUTION_SPECIES()
            print('       ' + ', '.join(species))

            if 'PHASES' in db.sections:
                print('  -> get_PHASES')
                phases = db.get_PHASES()
                print('       ' + ', '.join(phases))
            else:
                print('  -X skipped get_PHASES (no PHASES in database)')

if __name__ == '__main__':
    unittest.main()