import tempfile
import shutil
import os

from tests import QcipToolsTestCase
from qcip_tools.chemistry_files import gaussian


class GaussianTestCase(QcipToolsTestCase):
    """Gaussian stuffs"""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.input_file = os.path.join(self.temp_dir, 'input.com')

        with open(self.input_file, 'w') as f:
            with open('tests/tests_files/gaussian_input.com') as fx:
                f.write(fx.read())

        self.fchk_file = os.path.join(self.temp_dir, 'file.fchk')

        with open(self.fchk_file, 'w') as f:
            with open('tests/tests_files/gaussian_fchk.fchk') as fx:
                f.write(fx.read())

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_input_files(self):
        """Test the behavior of gaussian input file class"""

        text = 'some random text'

        # test with existing input
        gi1 = gaussian.Input()

        with open(self.input_file) as f:
            gi1.read(f)

        self.assertIsNotNone(gi1.molecule)
        self.assertEqual(len(gi1.molecule), 3)
        self.assertEqual(gi1.molecule.multiplicity, 2)
        self.assertEqual(gi1.molecule.charge, -1)
        self.assertEqual(gi1.molecule.number_of_electrons(), 11)  # = 2 (H and H) + 8 (O) + 1 (charge)

        self.assertEqual(len(gi1.options_dict), 3)
        self.assertIn('mem', gi1.options_dict)
        self.assertIn('nprocshared', gi1.options_dict)
        self.assertIn('chk', gi1.options_dict)

        self.assertEqual(len(gi1.other_blocks), 1)  # one block for the electric field

        # write it
        other_file = os.path.join(self.temp_dir, 'other.com')
        gi1.title = text

        with open(other_file, 'w') as f:
            gi1.write(f)

        self.assertTrue(os.path.exists(other_file))

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk=other\n' in content)  # chk has changed based on the file name
            self.assertTrue('%mem=3gb' in content)
            self.assertTrue('%nprocshared=4' in content)
            self.assertTrue(text in content)
            self.assertTrue('{} {}'.format(gi1.molecule.charge, gi1.molecule.multiplicity) in content)
            self.assertTrue(gi1.molecule.output_atoms() in content)
            self.assertTrue(gi1.other_blocks[0][0] in content)

        # write it with a custom chk
        custom_chk = 'another_one'
        with open(other_file, 'w') as f:
            gi1.write(f, chk=custom_chk)

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk={}\n'.format(custom_chk) in content)  # chk is custom

        # write with initial chk
        with open(other_file, 'w') as f:
            gi1.write(f, original_chk=True)

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk=calculation_0002\n' in content)  # chk remains the same as the original

    def test_fchk_file(self):

        fi = gaussian.FCHK()

        with open(self.fchk_file) as f:
            fi.read(f)

        self.assertEqual(fi.calculation_title, 'Ground state of LiH_mini with gen')
        self.assertEqual(fi.calculation_method, 'RHF')
        self.assertEqual(fi.calculation_basis, 'Gen')
        self.assertEqual(fi.calculation_type, 'SP')

        already_parsed = 6  # number of atoms, positions, weights, number of electron, charge and multiplicity
        self.assertEqual(len(fi.chunks_parsed), already_parsed)

        # get an array:
        self.assertTrue('Nuclear spins' in fi)
        self.assertTrue('Nuclear spins' not in fi.chunks_parsed)
        self.assertEqual(len(fi.get('Nuclear spins')), 2)
        self.assertEqual(len(fi.chunks_parsed), already_parsed + 1)
        self.assertTrue('Nuclear spins' in fi.chunks_parsed)  # OK, stored.
        self.assertTrue(fi.get('Nuclear spins')[0], 3)

        # get a single number
        self.assertTrue('Number of basis functions' in fi)
        self.assertTrue('Number of basis functions' not in fi.chunks_parsed)
        self.assertEqual(fi.get('Number of basis functions'), 6)
        self.assertTrue('Number of basis functions' in fi.chunks_parsed)

        # fool the storage system to make sure that it use it:
        fi.chunks_parsed['Number of basis functions'] = 8
        self.assertEqual(fi.get('Number of basis functions'), 8)
        self.assertEqual(fi.property('Number of basis functions'), 8)
        self.assertEqual(fi['Number of basis functions'], 8)

        # test molecule (conversion from a.u. to angstrom):
        self.assertAlmostEqual(fi.molecule[0].position[2], 0.04791742, places=3)
        self.assertAlmostEqual(fi.molecule[1].position[2], -1.45865742, places=3)
