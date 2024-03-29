import os
import random
import io
import argparse
import numpy
# import unittest
from unittest.mock import MagicMock, patch

from tests import QcipToolsTestCase, factories
from qcip_tools import molecule as qcip_molecule, atom as qcip_atom, derivatives, bounding, \
    chemistry_files
from qcip_tools.chemistry_files import ChemistryFile, gaussian, dalton, helpers, xyz, gamess, chemistry_datafile, pdb, \
    crystal


class HelpersTestCase(QcipToolsTestCase):
    """test the helpers"""

    def setUp(self):
        self.input_file = self.copy_to_temporary_directory('gaussian_input.com')
        self.input_file_2 = self.copy_to_temporary_directory('dalton_molecule.mol')

    def test_argparse_action(self):
        """Test the helper action"""

        parser = argparse.ArgumentParser()
        parser.add_argument('-a', action=helpers.create_open_chemistry_file_action())
        parser.add_argument('-b', action=helpers.create_open_chemistry_file_action(must_be=[gaussian.Input]))
        parser.add_argument('-c', action=helpers.create_open_chemistry_file_action(), nargs='*')  # list

        # 1. Open
        args = parser.parse_args(['-a', self.input_file])
        self.assertEqual(type(args.a), gaussian.Input)

        args = parser.parse_args(['-a', self.input_file_2])
        self.assertEqual(type(args.a), dalton.MoleculeInput)

        args = parser.parse_args(['-c', self.input_file, self.input_file_2])
        self.assertEqual(type(args.c[0]), gaussian.Input)
        self.assertEqual(type(args.c[1]), dalton.MoleculeInput)

        # 2. File does not exists:
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):  # ... Please shut up.
            with self.assertRaises(SystemExit) as cm:
                parser.parse_args(['-a', 'xxx'])
            self.assertNotEqual(cm.exception.code, 0)

        # 3. Open the right type:
        args = parser.parse_args(['-b', self.input_file])
        self.assertEqual(type(args.b), gaussian.Input)  # ok

        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit) as cm:
                parser.parse_args(['-b', self.input_file_2])
            self.assertNotEqual(cm.exception.code, 0)

        # open with identifier
        args = parser.parse_args(['-a', 'GAUSSIAN_INP:' + self.input_file])
        self.assertEqual(type(args.a), gaussian.Input)  # ok

        args = parser.parse_args(['-b', 'GAUSSIAN_INP:' + self.input_file])
        self.assertEqual(type(args.b), gaussian.Input)  # ok

        args = parser.parse_args(['-c', 'GAUSSIAN_INP:' + self.input_file, 'DALTON_MOL:' + self.input_file_2])
        self.assertEqual(type(args.c[0]), gaussian.Input)
        self.assertEqual(type(args.c[1]), dalton.MoleculeInput)  # ... ok :)

        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit) as cm:
                parser.parse_args(['-a', 'DALTON_MOL:' + self.input_file])  # not a dalton mol
            self.assertNotEqual(cm.exception.code, 0)

        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit) as cm:
                parser.parse_args(['-b', 'DALTON_MOL:' + self.input_file_2])  # not a correct possibility
            self.assertNotEqual(cm.exception.code, 0)

    def test_extra_chemistry_files(self):
        """Test extra chemistry file"""

        name = 'x.tmp'
        path = self.temporary_directory + '/' + name
        into = 'xxx'

        class DummyChemFile(chemistry_files.ChemistryFile, chemistry_files.WithIdentificationMixin):

            @classmethod
            def possible_file_extensions(cls):
                return [name.split('.')[-1]]

            @classmethod
            def attempt_identification(cls, f):
                return name in f.name

            def __init__(self):
                self.content = ''

            def read(self, f):
                self.content = f.read()

        helpers.EXTRA_CHEMISTRY_FILES.append(DummyChemFile)

        with open(path, 'w') as f:
            f.write(into)

        with open(path) as f:
            x = helpers.open_chemistry_file(f)
            self.assertIsInstance(x, DummyChemFile)
            self.assertEqual(x.content, into)


class ChemistryFileTestCase(QcipToolsTestCase):
    """base stuffs"""

    def test_add_property(self):

        f = ChemistryFile()

        # use existing callback
        self.assertTrue(f.has_property('file_type'))
        self.assertEqual(f.property('file_type'), f.file_type)

        # define custom property
        property_name = 'custom_property'
        self.assertFalse(f.has_property(property_name))

        with self.assertRaises(Exception):
            f.property(property_name)

        @ChemistryFile.define_property(property_name)
        def get_(obj, **kwargs):
            return obj.file_type

        self.assertTrue(f.has_property(property_name))
        self.assertEqual(f.property(property_name), f.file_type)

        # test heritage
        class T(ChemistryFile):
            file_type = 'XXX'

        ft = T()

        self.assertTrue(ft.has_property(property_name))
        self.assertEqual(ft.property(property_name), T.file_type)


class GaussianTestCase(QcipToolsTestCase):
    """Gaussian stuffs"""

    def setUp(self):
        self.input_file = self.copy_to_temporary_directory('gaussian_input.com')
        self.fchk_file = self.copy_to_temporary_directory('gaussian_fchk.fchk')
        self.fchk_file_v3 = self.copy_to_temporary_directory('gaussian_fchk_v3.fchk')
        self.fchk_one_atom = self.copy_to_temporary_directory('oneatom.fchk')
        self.fchk_ECP = self.copy_to_temporary_directory('gaussian_ECP.fchk')
        self.log_file = self.copy_to_temporary_directory('gaussian_output.log')

        self.cube_file = self.copy_to_temporary_directory('gaussian_cube.cub')
        self.cube_file2 = self.copy_to_temporary_directory('gaussian_cube_2.cub')
        self.cube_file3 = self.copy_to_temporary_directory('gaussian_cube_3.cub')

        self.basis_set1 = self.copy_to_temporary_directory('STO-3G.gbs')
        self.basis_set2 = self.copy_to_temporary_directory('cc-pVDZ.gbs')

    def test_input_files(self):
        """Test the behavior of gaussian input file class"""

        text = 'some random text'

        # test with existing input
        gi1 = gaussian.Input()
        self.assertFalse(gi1.from_read)

        with open(self.input_file) as f:
            gi1.read(f)

        self.assertTrue(gi1.from_read)

        self.assertIsNotNone(gi1.molecule)
        self.assertEqual(len(gi1.molecule), 3)
        self.assertEqual(gi1.molecule.multiplicity, 2)
        self.assertEqual(gi1.molecule.charge, -1)
        self.assertEqual(gi1.molecule.number_of_electrons(), 11)  # = 2 (H and H) + 8 (O) + 1 (charge)

        symbols = ['O', 'H', 'H']
        self.assertEqual(len(symbols), len(gi1.molecule))

        for index, a in enumerate(gi1.molecule):
            self.assertEqual(a.symbol, symbols[index])

        self.assertEqual(len(gi1.options), 3)
        self.assertIn('mem', gi1.options)
        self.assertIn('nprocshared', gi1.options)
        self.assertIn('chk', gi1.options)

        self.assertEqual(len(gi1.other_blocks), 1)  # one block for the electric field

        # write it
        other_file = os.path.join(self.temporary_directory, 'other.com')
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

        #  test creation:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)

        inp_ = gaussian.Input.from_molecule(m, title='test', input_card='#P B3LYP', options={'mem': '1GB'})
        self.assertEqual(inp_.title, 'test')
        self.assertEqual(inp_.molecule, m)
        self.assertEqual(inp_.input_card[0], '#P B3LYP')
        self.assertIn('mem', inp_.options)

    def test_fchk_file(self):

        fi = gaussian.FCHK()
        self.assertFalse(fi.from_read)

        with open(self.fchk_file) as f:
            fi.read(f)

        self.assertTrue(fi.from_read)

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
        self.assertEqual(fi.get('Number of basis functions'), 8)
        self.assertEqual(fi['Number of basis functions'], 8)

        # test molecule (conversion from a.u. to angstrom):
        self.assertAlmostEqual(fi.molecule[0].position[2], 0.04791742, places=3)
        self.assertAlmostEqual(fi.molecule[1].position[2], -1.45865742, places=3)

        # test energy
        energies = fi.property('computed_energies')
        self.assertIn('HF', energies)
        self.assertIn('SCF/DFT', energies)
        self.assertIn('total', energies)

        # one atom fchk (fix #35)
        fi = gaussian.FCHK()
        with open(self.fchk_one_atom) as f:
            fi.read(f)

        # ECP fchk (fix #36)
        fi = gaussian.FCHK()
        with open(self.fchk_ECP) as f:
            fi.read(f)

        self.assertNotEqual(fi.frozen_electrons, 0)
        self.assertEqual(fi.molecule.number_of_electrons(), fi.get('Number of electrons') + fi.frozen_electrons)

    def test_fchk_file_v3(self):
        """FCHK file with v3 are a bit different, because it includes characters"""

        fi = gaussian.FCHK()
        self.assertFalse(fi.from_read)

        with open(self.fchk_file_v3) as f:
            fi.read(f)

        self.assertTrue(fi.from_read)
        self.assertEqual(fi.get('Route'), '#p hf/cc-pVDZ polar=(DCSHG,cubic) cphf=rdfreq nosym')
        self.assertEqual(fi.get('Full Title'), 'Water gamma (cc-pVDZ)')

    def test_log_file(self):

        fo = gaussian.Output()
        self.assertFalse(fo.from_read)

        with open(self.log_file) as f:
            fo.read(f)

        self.assertTrue(fo.from_read)

        self.assertTrue(fo.chunk_exists(1))  # the real beginning of the program
        self.assertTrue(fo.chunk_exists(202))  # Link 202 deals with geometry
        self.assertFalse(fo.chunk_exists(9998))  # Does not exist (apart from 9999 and 1, all links have the form xxx)

        # test molecule
        symbols = ['O', 'H', 'H']
        self.assertEqual(len(symbols), len(fo.molecule))

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
        self.assertEqual(fo.search(to_find, into=1), line_found)
        self.assertEqual(fo.search(to_find, into=101), -1)

        # test energy
        energies = fo.property('computed_energies')
        self.assertIn('HF', energies)
        self.assertIn('MP2', energies)
        self.assertIn('total', energies)

    def test_cube_file(self):
        """Test the cube files"""

        fc = gaussian.Cube()
        self.assertFalse(fc.from_read)

        with open(self.cube_file) as f:
            fc.read(f)

        self.assertTrue(fc.from_read)

        self.assertEqual(fc.records_per_direction, [29, 29, 29])
        self.assertEqual(fc.data_per_record, 2)
        self.assertEqual(fc.MOs, [2, 3])
        self.assertEqual(fc.number_of_data(), 29 * 29 * 29 * 2)
        self.assertEqual(fc.cube_type, 'MO')

        # test molecule
        self.assertEqual(len(fc.molecule), 2)

        symbols = ['Li', 'H']
        self.assertEqual(len(symbols), len(fc.molecule))

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

        self.assertArraysAlmostEqual(square_subtract.records, fc2.records, places=5)

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

        # those results are checked against the original implementation of D. Jacquemin (like really!)
        self.assertArraysAlmostEqual(ct.charge, 0.882013, places=4)
        self.assertAlmostEqual(ct.distance, 2.42609, places=4)
        self.assertArraysAlmostEqual(ct.vector, [.0, .0, -2.42609], places=4)

    def test_sum_density(self):
        """Test the sum of density.

        Note that due to the low resolution of the cube file, results are approximates"""

        sets = {
            'almost_everything': bounding.BoundingSet(),
            'almost_nothing': bounding.BoundingSet(),
            'lithium_atom': bounding.BoundingSet(),
            'hydrogen_atom': bounding.BoundingSet()
        }

        fc = gaussian.Cube()

        with open(self.cube_file3) as f:
            fc.read(f)

        sets['almost_everything'].include(bounding.AABoundingBox(
            fc.origin * 0.52917165, maximum=(fc.origin + fc.records_per_direction * fc.increments) * 0.52917165))

        sets['almost_nothing'].include(bounding.AABoundingBox(fc.origin * 0.52917165, size=[.1, .1, .1]))
        sets['lithium_atom'].include(bounding.BoundingSphere(fc.molecule[0].position, radius=1.28))
        sets['hydrogen_atom'].include(bounding.BoundingSphere(fc.molecule[1].position, radius=0.31))

        results = fc.sum_density_of_sets(sets)
        self.assertAlmostEqual(results['almost_everything'], fc.molecule.number_of_electrons(), delta=.5)
        self.assertAlmostEqual(results['almost_nothing'], .0)
        self.assertAlmostEqual(results['lithium_atom'], 2.25, delta=.5)
        self.assertAlmostEqual(results['hydrogen_atom'], .25, delta=.5)

    def test_basis_set(self):
        """Test basis set"""

        # test reading
        b1 = gaussian.BasisSet()

        with open(self.basis_set1) as f:
            b1.read(f)

        self.assertTrue(b1.from_read)

        b2 = gaussian.BasisSet()

        with open(self.basis_set2) as f:
            b2.read(f)

        self.assertTrue(b2.from_read)

        atoms = ['H', 'C', 'O']

        for a in atoms:
            self.assertTrue(a in b1, msg=a)
            self.assertTrue(a in b2, msg=a)

        self.assertEqual(str(b1['H']), 'H [3s|1s]')
        self.assertEqual(str(b1['C']), 'C [6s3p|2s1p]')
        self.assertEqual(str(b1['O']), 'O [6s3p|2s1p]')

        self.assertEqual(str(b2['H']), 'H [4s1p|2s1p]')
        self.assertEqual(str(b2['C']), 'C [17s4p1d|3s2p1d]')
        self.assertEqual(str(b2['O']), 'O [17s4p1d|3s2p1d]')

        # try to fire exceptions:
        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\n****')  # Empty basis set
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nX\n****')  # Too short
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nQ 0\n S 1 1.0\n 10.0 1.0\n****')  # atom does not exists
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0 0\n S 1 1.0\n 10.0 1.0\n****')  # too much arguments for atom
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n S 1 1.0 1.0\n 10.0 1.0\n****')  # Too much argument for basis function
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n Q 1 1.0\n 10.0 1.0\n****')  # Not a shell
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n S 2 1.0\n 10.0 1.0\n****')  # Not the good number of primitives
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n S 1 1.0\n 10.0 1.0 1.0 1.0 1.0\n****')  # Wrong definition of primitive
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n S 1 1.0\n 10.0 1.0 1.0\n****')  # p coefficient for s-type shell
            b_test.read(o)

        with self.assertRaises(gaussian.BasisSetFormatError):
            b_test = gaussian.BasisSet()
            o = io.StringIO('****\nH 0\n S 1 1.0\n 10.0 1.0\nSP 1 1.0\n10.0 1.0\n****')  # Not p coefficient for SP
            b_test.read(o)

        # test writing
        other_file = os.path.join(self.temporary_directory, 'other.gbs')

        with open(other_file, 'w') as f:
            b1.write(f)

        b1p = gaussian.BasisSet()

        with open(other_file, 'r') as f:
            b1p.read(f)

        self.assertEqual(str(b1p['H']), 'H [3s|1s]')
        self.assertEqual(str(b1p['C']), 'C [6s3p|2s1p]')
        self.assertEqual(str(b1p['O']), 'O [6s3p|2s1p]')

        for a in atoms:
            a_1 = b1.basis_set[a]
            a_2 = b1p.basis_set[a]

            for bl_1, bl_2 in zip(a_1.shells(), a_2.shells()):
                self.assertEqual(len(a_1[bl_1]), len(a_2[bl_2]))
                for b_1, b_2 in zip(a_1[bl_1], a_2[bl_2]):
                    self.assertEqual(len(b_1), len(b_2))
                    for p_1, p_2 in zip(b_1.primitives, b_2.primitives):
                        self.assertAlmostEqual(p_1.exponent, p_2.exponent, places=5)
                        self.assertAlmostEqual(p_1.contraction_coefficient, p_2.contraction_coefficient, places=5)
                        self.assertAlmostEqual(p_1.p_coefficient, p_2.p_coefficient, places=5)

        #  test writing, only for given atoms
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)

        with open(other_file, 'w') as f:
            b1.write(f, for_molecule=m)

        b1p = gaussian.BasisSet()

        with open(other_file, 'r') as f:
            b1p.read(f)

        self.assertTrue('H' in b1p)
        self.assertTrue('O' in b1p)
        self.assertFalse('C' in b1p)

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
        self.input_mol_file = self.copy_to_temporary_directory('dalton_molecule.mol')
        self.input_mol_file2 = self.copy_to_temporary_directory('dalton_atombasis.mol')
        self.input_dal_file = self.copy_to_temporary_directory('dalton_input.dal')
        self.output_archive = self.copy_to_temporary_directory('dalton_archive.tar.gz')
        self.output_log = self.copy_to_temporary_directory('dalton_output.out')

    def test_input_mol(self):

        fm = dalton.MoleculeInput()
        self.assertFalse(fm.from_read)

        with open(self.input_mol_file) as f:
            fm.read(f)

        self.assertTrue(fm.from_read)

        self.assertEqual(fm.basis_set, 'd-aug-cc-pVDZ')
        self.assertEqual(fm.title, 'no-field')

        # test molecule
        self.assertEqual(len(fm.molecule), 3)

        symbols = ['O', 'H', 'H']
        self.assertEqual(len(symbols), len(fm.molecule))

        for index, a in enumerate(fm.molecule):
            self.assertEqual(a.symbol, symbols[index])

        self.assertEqual(fm.molecule.charge, -2)
        self.assertEqual(fm.molecule.multiplicity, 1)
        self.assertEqual(fm.molecule.number_of_electrons(), 12)

        # test atombasis input
        fm2 = dalton.MoleculeInput()
        self.assertFalse(fm2.from_read)

        with open(self.input_mol_file2) as f:
            fm2.read(f)

        self.assertTrue(fm2.from_read)
        self.assertEqual(fm2.molecule[0].extra['basis_set'], 'lanl2tz')
        self.assertEqual(fm2.molecule[1].extra['basis_set'], '6-311G*')
        self.assertEqual(fm2.molecule[0].extra['ecp'], 'lanl2tz')

        # test generation with no symmetry
        new_input = os.path.join(self.temporary_directory, 'new_mol.mol')
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
        new_input = os.path.join(self.temporary_directory, 'new_mol.mol')
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
        new_input = os.path.join(self.temporary_directory, 'new_mol.mol')
        with open(new_input, 'w') as f:
            fm.write(f, group_atoms=True)

        with open(new_input) as f:
            content = f.readlines()
            self.assertTrue('Nosymmetry' not in content[4])
            self.assertTrue('Charge=-2' in content[4])
            self.assertTrue('Angstrom' in content[4])
            self.assertTrue('0.115904592' in content[6])
            self.assertTrue('Charge=1.0 Atoms=2' in content[7])  # hydrogens are grouped

        # test creation:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)

        inp_ = dalton.MoleculeInput.from_molecule(m, title='test', basis_set='aug-cc-pVDZ')
        self.assertEqual(inp_.title, 'test')
        self.assertEqual(inp_.basis_set, 'aug-cc-pVDZ')
        self.assertEqual(inp_.molecule, m)

    def test_output_archive(self):
        """Test archive output"""

        fa = dalton.ArchiveOutput()
        self.assertFalse(fa.from_read)

        with open(self.output_archive, 'rb') as f:
            fa.read(f)

            self.assertTrue(fa.from_read)

            # test molecule
            self.assertEqual(len(fa.molecule), 3)

            symbols = ['O', 'H', 'H']
            self.assertEqual(len(symbols), len(fa.molecule))

            for index, a in enumerate(fa.molecule):
                self.assertEqual(a.symbol, symbols[index])

            with self.assertRaises(FileNotFoundError):
                fa.get_file('whatever')

            # test energy
            energies = fa.property('computed_energies')
            self.assertIn('CCS/SCF', energies)
            self.assertIn('total', energies)

    def test_output(self):
        """Test the bare dalton output"""

        fl = dalton.Output()
        self.assertFalse(fl.from_read)

        with open(self.output_log) as f:
            fl.read(f)

        self.assertTrue(fl.from_read)

        # test molecule
        self.assertEqual(len(fl.molecule), 3)

        symbols = ['O', 'H', 'H']
        self.assertEqual(len(symbols), len(fl.molecule))

        for index, a in enumerate(fl.molecule):
            self.assertEqual(a.symbol, symbols[index])

        # test find string:
        self.assertTrue(fl.chunk_exists('SIRIUS'))
        self.assertFalse(fl.chunk_exists('XXXX'))

        to_find = '@    Occupied SCF orbitals'
        line_found = 515
        self.assertEqual(fl.search(to_find), line_found)
        self.assertEqual(fl.search(to_find, line_start=line_found), line_found)
        self.assertEqual(fl.search(to_find, line_end=line_found), -1)
        self.assertEqual(fl.search(to_find, line_start=line_found + 1), -1)
        self.assertEqual(fl.search(to_find, into='SIRIUS'), line_found)
        self.assertEqual(fl.search(to_find, into='CC'), -1)

        # test energy
        energies = fl.property('computed_energies')
        self.assertIn('SCF/DFT', energies)
        self.assertIn('total', energies)

        # compare with the archive output:
        fa = dalton.ArchiveOutput()

        with open(self.output_archive, 'rb') as f:
            fa.read(f)

            for index, a in enumerate(fa.molecule):
                self.assertArraysAlmostEqual(a.position, fl.molecule[index].position)

            energies_archive = fa.property('computed_energies')
            self.assertAlmostEqual(energies_archive['CCS/SCF'], energies['SCF/DFT'])

        # test get inputs
        dal, mol = fl.get_inputs()

        for index, a in enumerate(mol.molecule):
            self.assertEqual(a.symbol, symbols[index])
            self.assertArraysAlmostEqual(a.position, fl.molecule[index].position)  # it is the same molecule

        self.assertTrue('.STATIC' in dal['WAVE FUNCTION']['CCQR'])  # it is probably the good dal as well

    def test_input_dal(self):
        """Test .dal files"""

        fi = dalton.Input()
        self.assertFalse(fi.from_read)

        with open(self.input_dal_file) as fx:
            fi.read(fx)

        self.assertTrue(fi.from_read)

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
            fi['DALTON INPUT']['DIRECT'] = dalton.InputCard()  # input card should start with "."

        with self.assertRaises(Exception):
            fi['DALTON INPUT']['.XX'] = dalton.InputModule()  # module should not start with "."

        with self.assertRaises(Exception):
            fi['XX'] = dalton.InputModule(0, 'XX')  # XX is not an allowed (level 0) module name

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR'] = dalton.InputModule(1, 'CCLR')  # module already exists

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR']['XX'] = dalton.InputModule(level=2)  # no subsubmodules

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR']['XX'] = dalton.InputModule(level=1)  # no subsubmodules, I said!

        with self.assertRaises(Exception):
            fi['PROPERTIES'] = dalton.InputModule(0, 'WAVE FUNCTIONS')  # key and name divergence

        with self.assertRaises(Exception):
            fi['WAVE FUNCTIONS']['CCLR']['.XX'] = dalton.InputCard('YY')  # key and name divergence

        # test file creation
        fix = dalton.Input()

        fix['DALTON'] = dalton.InputModule()  # "DALTON" is 5 characters, so it's good
        fix['DALTON']['.DIRECT'] = dalton.InputCard()
        fix['DALTON']['.RUN WAVE FUNCTIONS'] = dalton.InputCard()
        fix['WAVE F'] = dalton.InputModule()  # "WAVE F" is 6 characters, so it is larger than the required 5 ones
        fix['WAVE F']['.DFT'] = dalton.InputCard(parameters=['B3LYP'])

        with open(os.path.join(self.tests_files_directory, 'dalton_input_2.dal'), 'r') as fx:
            self.assertEqual(str(fix), fx.read())

        # test update
        file_before = fix.to_string()
        self.assertFalse('SCF INPUT' in fix['WAVE F'])

        fix.update('**WAVE FUNC\n*SCF INPUT\n.THRESH\n1.0D-6')
        file_after = fix.to_string()
        self.assertNotEqual(file_before, file_after)
        self.assertTrue('SCF INPUT' in fix['WAVE F'])
        self.assertTrue('.THRESH' in fix['WAVE F']['SCF INPUT'])
        self.assertEqual(fix['WAVE F']['SCF INPUT']['.THRESH'].parameters[0], '1.0D-6')

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


class XYZTestCase(QcipToolsTestCase):
    """XYZ stuffs"""

    def setUp(self):
        self.xyz_file = self.copy_to_temporary_directory('xyz_molecule.xyz')

    def test_xyz(self):
        """Test the behavior of the xyz"""

        fx = xyz.File()
        self.assertFalse(fx.from_read)

        with open(self.xyz_file) as f:
            fx.read(f)

        self.assertTrue(fx.from_read)

        # test molecule
        symbols = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H']
        self.assertEqual(len(symbols), len(fx.molecule))

        for index, a in enumerate(fx.molecule):
            self.assertEqual(symbols[index], a.symbol)

        # test writing
        other_xyz = os.path.join(self.temporary_directory, 'n.xyz')

        with open(other_xyz, 'w') as f:
            fx.write(f)

        with open(other_xyz) as f:
            fx2 = xyz.File()
            fx2.read(f)

            self.assertEqual(len(fx2.molecule), len(fx.molecule))

            for index, a in enumerate(fx2.molecule):
                self.assertEqual(a.symbol, fx.molecule[index].symbol)
                self.assertArraysAlmostEqual(a.position, fx.molecule[index].position)

        # test creation:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)

        xyz_ = xyz.File.from_molecule(m, title='test')
        self.assertEqual(xyz_.title, 'test')
        self.assertEqual(xyz_.molecule, m)

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.xyz_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), xyz.File)


class GAMESSTestCase(QcipToolsTestCase):
    """GAMESS stuffs"""

    def setUp(self):
        self.input_file = self.copy_to_temporary_directory('gamess_input.inp')
        self.input_file_2 = self.copy_to_temporary_directory('gamess_input_2.inp')
        self.output_file = self.copy_to_temporary_directory('gamess_output.log')

    def test_input(self):
        """Test the behavior of the input"""

        fi = gamess.Input()
        self.assertFalse(fi.from_read)

        with open(self.input_file) as f:
            fi.read(f)

        self.assertTrue(fi.from_read)

        # test modules
        modules = ['basis', 'contrl', 'scf', 'tdhfx', 'system', 'data']
        for m in modules:
            self.assertTrue(m in fi, msg='{} not in!'.format(m))

        basis_module = fi['basis']
        self.assertEqual(len(basis_module.options), 2)
        self.assertTrue('gbasis' in basis_module)

        control_module = fi['contrl']
        self.assertEqual(len(control_module.options), 2)
        self.assertTrue('scftyp' in control_module)
        self.assertTrue('runtyp' in control_module)

        # test molecule
        symbols = ['C', 'H', 'H', 'H', 'H']
        self.assertEqual(len(symbols), len(fi.molecule))

        for index, a in enumerate(fi.molecule):
            self.assertEqual(symbols[index], a.symbol)

        # test writing
        other_input = os.path.join(self.temporary_directory, 'u.inp')

        with open(other_input, 'w') as f:
            fi.write(f)

        with open(other_input) as f:
            fi2 = gamess.Input()
            fi2.read(f)

            self.assertEqual(len(fi.modules), len(fi2.modules))
            self.assertEqual(len(fi2.molecule), len(fi.molecule))

            for index, a in enumerate(fi2.molecule):
                self.assertEqual(a.symbol, fi.molecule[index].symbol)
                self.assertArraysAlmostEqual(a.position, fi.molecule[index].position)

        # test with a more tricky input (comment and too long line)
        fi2 = gamess.Input()

        with open(self.input_file_2) as f:
            fi2.read(f)

        self.assertTrue('force' in fi2)
        force_module = fi2['force']
        self.assertEqual(len(force_module.options), 4)

        options = ['method', 'vibanl', 'temp(1)', 'sclfac']
        for o in options:
            self.assertTrue(o in force_module, msg=o)

        self.assertNotIn('!', str(force_module))
        self.assertIn('\n', str(force_module))

        # test creation:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)
        mod = '$CONTRL SCFTYP=RHF RUNTYP=ENERGY $END'

        inp_ = gamess.Input.from_molecule(m, title='test', modules=[gamess.InputModule.from_string(mod)])
        self.assertEqual(inp_.title, 'test')
        self.assertEqual(inp_.molecule, m)
        self.assertEqual(len(inp_.modules), 1)
        self.assertTrue('contrl' in inp_)

    def test_output(self):
        """Test the behavior of the output"""

        fo = gamess.Output()
        self.assertFalse(fo.from_read)

        with open(self.output_file) as f:
            fo.read(f)

        self.assertTrue(fo.from_read)

        # test molecule
        symbols = ['O', 'H', 'H']
        self.assertEqual(len(symbols), len(fo.molecule))

        for index, a in enumerate(fo.molecule):
            self.assertEqual(symbols[index], a.symbol)

        # test find string:
        self.assertTrue(fo.chunk_exists('RHF CALCULATION'))
        self.assertFalse(fo.chunk_exists('XXXX'))

        to_find = 'TOTAL NUMBER OF BASIS SET SHELLS'
        line_found = 231
        self.assertEqual(fo.search(to_find), line_found)
        self.assertEqual(fo.search(to_find, line_start=line_found), line_found)
        self.assertEqual(fo.search(to_find, line_end=line_found), -1)
        self.assertEqual(fo.search(to_find, line_start=line_found + 1), -1)
        self.assertEqual(fo.search(to_find, into='SETTING UP THE RUN'), line_found)
        self.assertEqual(fo.search(to_find, into='RHF CALCULATION'), -1)

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.input_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gamess.Input)

        with open(self.output_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), gamess.Output)


class ChemistryDatafileTestCase(QcipToolsTestCase):
    """Chemistry data file stuffs"""

    def setUp(self):
        self.fchk_file = self.copy_to_temporary_directory('gaussian_fchk.fchk')
        self.chemistry_data_file = self.copy_to_temporary_directory('chemistry_data_file.hdf5')

    def test_datafile(self):
        """Test the behavior of the datafile"""

        fx = gaussian.FCHK()
        self.assertFalse(fx.from_read)

        with open(self.fchk_file) as f:
            fx.read(f)

        fd = chemistry_datafile.ChemistryDataFile()
        fd.title = 'test'
        fd.molecule = fx.molecule
        fd.spacial_dof = 3 * len(fx.molecule)
        fd.trans_plus_rot_dof = 5
        fd.derivatives['F'] = {'static': factories.FakeElectricDipole()}
        fd.derivatives['dD'] = {
            '1064nm': factories.FakePolarizabilityTensor(isotropy_factor=2),
            0.04: factories.FakePolarizabilityTensor(isotropy_factor=1.5),
        }
        fd.derivatives['FFF'] = {
            'static': factories.FakeFirstHyperpolarizabilityTensor()}
        fd.derivatives['G'] = derivatives.Tensor(
            'G', components=numpy.zeros((3 * len(fx.molecule),)), spacial_dof=3 * len(fx.molecule))
        fd.derivatives[''] = derivatives.Tensor('', components=numpy.array((1.05,)))

        # test writing
        other_input = os.path.join(self.temporary_directory, 'u.hdf5')

        with open(other_input, 'wb') as f:
            fd.write(f)

        # test reading
        fde = chemistry_datafile.ChemistryDataFile()
        self.assertFalse(fde.from_read)
        with open(other_input, 'rb') as f:
            fde.read(f)

        self.assertTrue(fde.from_read)

        self.assertEqual(fde.title, fd.title)
        self.assertEqual(fde.version, fd.version)
        self.assertEqual(fde.molecule.charge, fd.molecule.charge)
        self.assertEqual(fde.molecule.multiplicity, fd.molecule.multiplicity)
        self.assertEqual(len(fde.molecule), len(fd.molecule))

        for i, a in enumerate(fde.molecule):
            self.assertEqual(a.symbol, fd.molecule[i].symbol)
            self.assertArraysAlmostEqual(a.position, fd.molecule[i].position)

        self.assertArraysAlmostEqual(
            fde.derivatives['F']['static'].components, fd.derivatives['F']['static'].components)
        self.assertArraysAlmostEqual(
            fde.derivatives['dD']['1064nm'].components, fd.derivatives['dD']['1064nm'].components)
        self.assertArraysAlmostEqual(
            fde.derivatives['dD'][0.04].components, fd.derivatives['dD'][0.04].components)
        self.assertArraysAlmostEqual(
            fde.derivatives['FFF']['static'].components, fd.derivatives['FFF']['static'].components)
        self.assertArraysAlmostEqual(
            fde.derivatives['G'].components, fd.derivatives['G'].components)
        self.assertArraysAlmostEqual(
            fde.derivatives[''].components, fd.derivatives[''].components)  # yeah, we can store energy as well

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.chemistry_data_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), chemistry_datafile.ChemistryDataFile)


class PDBTestCase(QcipToolsTestCase):
    """PDB stuffs"""

    def setUp(self):
        self.pdb_file = self.copy_to_temporary_directory('6g7h.pdb')

    def test_pdb(self):
        """Test the behavior of the pdb file"""

        # test read
        fi = pdb.File()
        self.assertFalse(fi.from_read)

        with open(self.pdb_file) as f:
            fi.read(f)

        self.assertTrue(fi.from_read)

        dx = {
            'pdb_resSeq': '   5',
            'pdb_name': ' N  ',
            'pdb_altLoc': ' ',
            'pdb_resName': 'THR',
            'pdb_chainId': 'A',
            'pdb_iCode': ' ',
            'pdb_occupancy': 1.,
            'pdb_tempFactor': 79.36,
            'pdb_charge': '  ',
        }

        self.assertDictEqual(dx, fi.molecule[0].extra)

        # test write
        other_input = os.path.join(self.temporary_directory, 'tmp.pdb')
        with open(other_input, 'w') as f:
            fi.write(f)

        fxe = pdb.File()
        with open(other_input) as f:
            fxe.read(f)

        self.assertEqual(len(fxe.molecule), len(fi.molecule))
        self.assertEqual(fxe.TERs, fi.TERs)
        self.assertEqual(fxe.molecule[0].position, fi.molecule[0].position)
        self.assertDictEqual(fxe.molecule[0].extra, fi.molecule[0].extra)

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.pdb_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), pdb.File)


class CrystalTestCase(QcipToolsTestCase):
    """Crystal stuffs"""

    def setUp(self):
        self.output_file = self.copy_to_temporary_directory(
            'properties/electrical_derivatives/crystal_output_static.out')

    def test_primitive_cell(self):
        fo = crystal.Output()
        self.assertFalse(fo.from_read)

        with open(self.output_file) as f:
            fo.read(f)

        self.assertArraysAlmostEqual(fo.lattice_abc, [14.4597, 14.4597, 14.4597])
        self.assertArraysAlmostEqual(fo.lattice_albega, [90, 90, 90])
        self.assertEqual(len(fo.molecule), 300)

    def test_file_recognition(self):
        """Test that the helper function recognise file as it is"""

        with open(self.output_file) as f:
            self.assertIsInstance(helpers.open_chemistry_file(f), crystal.Output)
