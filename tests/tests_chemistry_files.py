import tempfile
import shutil
import os
import random

from tests import QcipToolsTestCase
from qcip_tools import math as qcip_math
from qcip_tools.chemistry_files import gaussian, dalton, helpers


class GaussianTestCase(QcipToolsTestCase):
    """Gaussian stuffs"""

    def setUp(self):

        self.temp_dir = tempfile.mkdtemp()

        self.input_file = os.path.join(self.temp_dir, 'input.com')

        with open(self.input_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_input.com')) as fx:
                f.write(fx.read())

        self.fchk_file = os.path.join(self.temp_dir, 'file.fchk')

        with open(self.fchk_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_fchk.fchk')) as fx:
                f.write(fx.read())

        self.log_file = os.path.join(self.temp_dir, 'file.log')

        with open(self.log_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_output.log')) as fx:
                f.write(fx.read())

        self.cube_file = os.path.join(self.temp_dir, 'file.cub')

        with open(self.cube_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_cube.cub')) as fx:
                f.write(fx.read())

        self.cube_file2 = os.path.join(self.temp_dir, 'file2.cub')

        with open(self.cube_file2, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_cube_2.cub')) as fx:
                f.write(fx.read())

        self.cube_file3 = os.path.join(self.temp_dir, 'file3.cub')

        with open(self.cube_file3, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/gaussian_cube_3.cub')) as fx:
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

        self.assertEqual(fi.title, 'Ground state of LiH_mini with gen')
        self.assertEqual(fi.calculation_method, 'RHF')
        self.assertEqual(fi.basis_set, 'Gen')
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

    def test_log_file(self):

        fo = gaussian.Output()

        with open(self.log_file) as f:
            fo.read(f)

        self.assertTrue(fo.link_called(1))  # the real beginning of the program
        self.assertTrue(fo.link_called(202))  # Link 202 deals with geometry
        self.assertFalse(fo.link_called(9998))  # Does not exists (apart from 9999 and 1, all link have the form xxx)

        # test molecule
        self.assertEqual(len(fo.molecule), 3)

        symbols = ['O', 'H', 'H']
        for index, a in enumerate(fo.molecule):
            self.assertEqual(a.symbol, symbols[index])

        self.assertEqual(fo.molecule.charge, 0)
        self.assertEqual(fo.molecule.multiplicity, 1)

        # test find string:o
        to_find = 'RESTRICTED RIGHTS LEGEND'
        line_found = 30
        self.assertEqual(fo.search(to_find), line_found)
        self.assertEqual(fo.search(to_find, line_start=line_found), line_found)
        self.assertEqual(fo.search(to_find, line_end=line_found), -1)
        self.assertEqual(fo.search(to_find, line_start=line_found + 1), -1)
        self.assertEqual(fo.search(to_find, in_link=1), line_found)
        self.assertEqual(fo.search(to_find, in_link=101), -1)

    def test_cube_file(self):
        """Test the cube files"""

        fc = gaussian.Cube()

        with open(self.cube_file) as f:
            fc.read(f)

        self.assertEqual(fc.records_per_direction, [29, 29, 29])
        self.assertEqual(fc.data_per_record, 2)
        self.assertEqual(fc.MOs, [2, 3])
        self.assertEqual(fc.number_of_data(), 29 * 29 * 29 * 2)
        self.assertEqual(fc.cube_type, 'MO')

        # test molecule
        self.assertEqual(len(fc.molecule), 2)

        symbols = ['Li', 'H']
        for index, a in enumerate(fc.molecule):
            self.assertEqual(a.symbol, symbols[index])

        # test slice operation
        MO_1 = fc.slice(0)
        self.assertEqual(MO_1.data_per_record, 1)
        self.assertEqual(MO_1.MOs[0], 2)
        self.assertEqual(MO_1.records.shape, (29, 29, 29, 1))

        MO_2 = fc.slice(1)
        self.assertEqual(MO_2.data_per_record, 1)
        self.assertEqual(MO_2.MOs[0], 3)
        self.assertEqual(MO_2.records.shape, (29, 29, 29, 1))

        random_place = random.randrange(0, 28), random.randrange(0, 28), random.randrange(0, 28)
        self.assertEqual(MO_1.records[random_place][0], fc.records[random_place][0])
        self.assertEqual(MO_2.records[random_place][0], fc.records[random_place][1])

        # test mathematical operations
        square_subtract = (MO_1 - MO_2) ** 2

        fc2 = gaussian.Cube()

        with open(self.cube_file2) as f:
            fc2.read(f)

        self.assertArrayAlmostEqual(square_subtract.records, fc2.records, places=5)

        square_subtract.title = fc2.title
        square_subtract.subtitle = fc2.subtitle
        square_subtract.cube_type = 'density'
        out = square_subtract.to_string()

        with open(self.cube_file2) as f:
            ref = f.read()
            self.assertEqual(out, ref)

    def test_cube_charge_transfer(self):
        """Test the computation of the charge transfer"""

        fc = gaussian.Cube()

        with open(self.cube_file) as f:
            fc.read(f)

        MO_1 = fc.slice(0)
        MO_2 = fc.slice(1)

        diff_of_square = MO_2 ** 2 - MO_1 ** 2
        ct = diff_of_square.compute_charge_transfer()

        # those results are checked against the original implementation of D. Jacquemin.
        self.assertArrayAlmostEqual(ct.charge, 0.8836, places=4)
        self.assertAlmostEqual(ct.distance, 2.4217, places=4)
        self.assertArrayAlmostEqual(ct.vector, [.0, .0, -2.4217], places=4)

    def test_sum_density(self):
        """Test the sum of density.

        Note that due to the low resolution of the cube file, results are approximates"""

        sets = {
            'almost_everything': qcip_math.BoundingSet(),
            'almost_nothing': qcip_math.BoundingSet(),
            'lithium_atom': qcip_math.BoundingSet(),
            'hydrogen_atom': qcip_math.BoundingSet()
        }

        fc = gaussian.Cube()

        with open(self.cube_file3) as f:
            fc.read(f)

        sets['almost_everything'].include(qcip_math.AABoundingBox(
            fc.origin * 0.52917165, maximum=(fc.origin + fc.records_per_direction * fc.increments) * 0.52917165))

        sets['almost_nothing'].include(qcip_math.AABoundingBox(fc.origin * 0.52917165, size=[.1, .1, .1]))
        sets['lithium_atom'].include(qcip_math.BoundingSphere(fc.molecule[0].position, radius=1.28))
        sets['hydrogen_atom'].include(qcip_math.BoundingSphere(fc.molecule[1].position, radius=0.31))

        results = fc.sum_density_of_sets(sets)
        self.assertAlmostEqual(results['almost_everything'], fc.molecule.number_of_electrons(), delta=.5)
        self.assertAlmostEqual(results['almost_nothing'], .0)
        self.assertAlmostEqual(results['lithium_atom'], 2.25, delta=.5)
        self.assertAlmostEqual(results['hydrogen_atom'], .25, delta=.5)

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.input_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gaussian.Input)

        with open(self.fchk_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gaussian.FCHK)

        with open(self.log_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gaussian.Output)

        with open(self.cube_file2) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gaussian.Cube)


class DaltonTestCase(QcipToolsTestCase):
    """Dalton stuffs"""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.input_mol_file = os.path.join(self.temp_dir, 'input.mol')

        with open(self.input_mol_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/dalton_molecule.mol')) as fx:
                f.write(fx.read())

        self.input_dal_file = os.path.join(self.temp_dir, 'input.dal')

        with open(self.input_dal_file, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/dalton_input.dal')) as fx:
                f.write(fx.read())

        self.output_archive = os.path.join(self.temp_dir, 'output.tar.gz')

        with open(self.output_archive, 'wb') as f:
            with open(os.path.join(self.test_directory, 'tests_files/dalton_archive.tar.gz'), 'rb') as fx:
                f.write(fx.read())

        self.output_log = os.path.join(self.temp_dir, 'output.out')

        with open(self.output_log, 'w') as f:
            with open(os.path.join(self.test_directory, 'tests_files/dalton_output.out'), 'r') as fx:
                f.write(fx.read())

    def test_input_mol(self):

        fm = dalton.MoleculeInput()

        with open(self.input_mol_file) as f:
            fm.read(f)

        self.assertEqual(fm.basis_set, 'd-aug-cc-pVDZ')
        self.assertEqual(fm.title, 'no-field')

        # test molecule
        self.assertEqual(len(fm.molecule), 3)

        symbols = ['O', 'H', 'H']
        for index, a in enumerate(fm.molecule):
            self.assertEqual(a.symbol, symbols[index])

        self.assertEqual(fm.molecule.charge, -2)
        self.assertEqual(fm.molecule.multiplicity, 1)
        self.assertEqual(fm.molecule.number_of_electrons(), 12)

        # test generation with no symmetry
        new_input = os.path.join(self.temp_dir, 'new_mol.mol')
        with open(new_input, 'w') as f:
            fm.write(f, nosym=True)

        with open(new_input) as f:
            content = f.readlines()
            self.assertTrue('Angstrom' in content[4])
            self.assertTrue('Charge=-2' in content[4])
            self.assertTrue('Nosymmetry' in content[4])
            self.assertTrue('0.115904592' in content[6])
            self.assertTrue('Charge=1.0 Atoms=1' in content[7])

        # test generation in atomic units
        new_input = os.path.join(self.temp_dir, 'new_mol.mol')
        with open(new_input, 'w') as f:
            fm.write(f, in_angstrom=False)

        with open(new_input) as f:
            content = f.readlines()
            self.assertTrue('Charge=-2' in content[4])
            self.assertTrue('Nosymmetry' not in content[4])
            self.assertTrue('Angstrom' not in content[4])
            self.assertAlmostEqual(float(content[6].split()[-1]), 0.115904592 * 0.52917165, places=5)
            self.assertTrue('Charge=1.0 Atoms=1' in content[7])

        # test generation with grouping
        new_input = os.path.join(self.temp_dir, 'new_mol.mol')
        with open(new_input, 'w') as f:
            fm.write(f, group_atoms=True)

        with open(new_input) as f:
            content = f.readlines()
            self.assertTrue('Nosymmetry' not in content[4])
            self.assertTrue('Charge=-2' in content[4])
            self.assertTrue('Angstrom' in content[4])
            self.assertTrue('0.115904592' in content[6])
            self.assertTrue('Charge=1.0 Atoms=2' in content[7])  # hydrogens are grouped

    def test_output_archive(self):
        """Test archive output"""

        fa = dalton.ArchiveOutput()

        with open(self.output_archive, 'rb') as f:
            fa.read(f)

            # test molecule
            self.assertEqual(len(fa.molecule), 3)

            symbols = ['O', 'H', 'H']
            for index, a in enumerate(fa.molecule):
                self.assertEqual(a.symbol, symbols[index])

            with self.assertRaises(FileNotFoundError):
                fa.get_file('whatever')

    def test_output(self):
        """Test the bare dalton output"""

        fl = dalton.Output()

        with open(self.output_log) as f:
            fl.read(f)

        # test molecule
        self.assertEqual(len(fl.molecule), 3)

        symbols = ['O', 'H', 'H']
        for index, a in enumerate(fl.molecule):
            self.assertEqual(a.symbol, symbols[index])

        # test find string:
        to_find = '@    Occupied SCF orbitals'
        line_found = 515
        self.assertEqual(fl.search(to_find), line_found)
        self.assertEqual(fl.search(to_find, line_start=line_found), line_found)
        self.assertEqual(fl.search(to_find, line_end=line_found), -1)
        self.assertEqual(fl.search(to_find, line_start=line_found + 1), -1)
        self.assertEqual(fl.search(to_find, in_section='SIRIUS'), line_found)
        self.assertEqual(fl.search(to_find, in_section='CC'), -1)

        # compare with the archive output:
        fa = dalton.ArchiveOutput()

        with open(self.output_archive, 'rb') as f:
            fa.read(f)

            for index, a in enumerate(fa.molecule):
                self.assertArrayAlmostEqual(a.position, fl.molecule[index].position)

        # test get inputs
        dal, mol = fl.get_inputs()

        for index, a in enumerate(mol.molecule):
            self.assertEqual(a.symbol, symbols[index])
            self.assertArrayAlmostEqual(a.position, fl.molecule[index].position)  # it is the same molecule

        self.assertTrue('.STATIC' in dal['WAVE FUNCTION']['CCQR'])  # it is probably the good dal as well

    def test_input_dal(self):
        """Test .dal files"""

        fi = dalton.Input()
        with open(self.input_dal_file) as fx:
            fi.read(fx)

        # check if reading is done properly (and if "__contains__()" works)
        available_modules = ['DALTON INPUT', 'INTEGRALS', 'WAVE FUNCTIONS']
        for m in available_modules:
            self.assertTrue(m[:5] in fi.modules)
            self.assertTrue(m in fi)

        self.assertFalse('RESPONSE' in fi)

        self.assertTrue('DIRECT' in fi['DALTON INPUT'].input_cards)
        self.assertTrue('.DIRECT' in fi['DALTON INPUT'])
        self.assertTrue('.CC' in fi['WAVE FUNCTIONS'])
        self.assertTrue('.STATIC' in fi['WAVE FUNCTIONS']['CCLR'])

        self.assertFalse('.DIRECT' in fi['INTEGRALS'])
        self.assertFalse('.DIRECT' in fi['WAVE FUNCTIONS']['CCLR'])

        self.assertTrue('SCF INPUT' in fi['WAVE FUNCTIONS'])

        self.assertTrue(len(fi['DALTON INPUT']['.DIRECT'].parameters) == 0)
        self.assertTrue(len(fi['WAVE FUNCTIONS']['CCLR']['.FREQUE'].parameters) == 2)
        self.assertEqual(fi['WAVE FUNCTIONS']['CCLR']['.FREQUE'].parameters[0], '9')

        # check if file is identical to input
        with open(self.input_dal_file) as fx:
            self.assertEqual(str(fi), fx.read())

        # fire exceptions
        with self.assertRaises(Exception):
            fi['PROPERTIES'] = dalton.InputModule('WAVE FUNCTIONS', 0)  # key and name divergence

        with self.assertRaises(Exception):
            fi['DALTON INPUT']['DIRECT'] = dalton.InputCard('DIRECT')  # input card should start with "."

        with self.assertRaises(Exception):
            fi['XX'] = dalton.InputModule('XX', 0)  # XX is not an allowed module name

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR'] = dalton.InputModule('CCLR', 1)  # module already exists

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR']['XX'] = dalton.InputModule('XX', 1)  # no subsubmodules

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR']['.XX'] = dalton.InputCard('YY')  # key and name divergence

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.output_log) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), dalton.Output)

        with open(self.input_mol_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), dalton.MoleculeInput)

        with open(self.input_dal_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), dalton.Input)

        with open(self.output_archive) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), dalton.ArchiveOutput)
