import os
import unittest
from glob import glob

import blazy
from blazy.phreeqc.parser import datParser

database_dir = os.path.join(os.path.dirname(blazy.__file__), 'resources', 'database')

class TestDatabaseInteraction(unittest.TestCase):

    def test_file_paths(self):
        # where the datatbase should be
        dbase_path = os.path.join(database_dir, 'pitzer.dat')

        # where the database is
        db = datParser(database='pitzer')
        self.assertEqual(db.path, dbase_path)

        db = datParser(database='pitzer.dat')
        self.assertEqual(db.path, dbase_path)

        db = datParser(database=dbase_path)
        self.assertEqual(db.path, dbase_path)

    def test_all_databases(self):
        databases = glob(os.path.join(database_dir, '*.dat'))
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