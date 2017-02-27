import unittest
import numpy

from qcip_tools import atom as qcip_atom, molecule as qcip_molecule, math as qcip_math


class MoleculeTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_molecule_creation(self):
        """Test the behavior of the Molecule class"""

        # empty atom
        m1 = qcip_molecule.Molecule()
        self.assertEqual(m1.mass(), 0.0)
        self.assertEqual(m1.charge, 0.0)
        self.assertEqual(m1.number_of_electrons(), 0)
        self.assertEqual(m1.formula(), '')

        # H2 molecule
        atom_list = [qcip_atom.Atom(symbol='H', position=[-.5, 0, 0]), qcip_atom.Atom(symbol='H', position=[.5, 0, 0])]
        m2 = qcip_molecule.Molecule(atom_list=atom_list)
        self.assertEqual(m2.mass(), sum(a.mass for a in atom_list))
        self.assertEqual(m2.charge, 0)
        self.assertEqual(m2.number_of_electrons(), 2)
        self.assertEqual(m2.formula(), 'H2')

        self.assertTrue(numpy.array_equal(m2.center_of_mass(), [0, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.center_of_charges(), [0, 0, 0]))

        # Test shifted index:
        self.assertTrue(numpy.array_equal(m2.atom(1).position, [-.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.atom(2).position, [.5, 0, 0]))

        with self.assertRaises(KeyError):
            m2.atom(0)
            m2.atom(-1)
            m2.atom(3)

        # test add atom:
        self.assertEqual(len(m2), 2)
        m2.insert(qcip_atom.Atom(symbol='H', position=[0, .25, 0]))
        self.assertEqual(len(m2), 3)
        self.assertTrue(numpy.array_equal(m2.atom(3).position, [0, .25, 0]))
        self.assertEqual(m2.number_of_electrons(), 3)

        # test MO
        m2.charge = 1
        self.assertEqual(m2.number_of_electrons(), 2)
        self.assertEqual(m2.MO_to_str(1), 'HOMO')
        self.assertEqual(m2.MO_to_str(2), 'LUMO')
        self.assertEqual(m2.MO_to_str(3), 'LUMO+1')

        m2.charge = -1
        self.assertEqual(m2.number_of_electrons(), 4)
        self.assertEqual(m2.MO_to_str(1), 'HOMO-1')
        self.assertEqual(m2.MO_to_str(2), 'HOMO')
        self.assertEqual(m2.MO_to_str(3), 'LUMO')
        self.assertEqual(m2.MO_to_str(4), 'LUMO+1')

        # test remove atom
        m2.remove_atom(3)
        self.assertEqual(len(m2), 2)
        self.assertTrue(numpy.array_equal(m2.atom(1).position, [-.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.atom(2).position, [.5, 0, 0]))

        # test translate
        m2.translate([1, 0, 0])
        self.assertTrue(numpy.array_equal(m2.atom(1).position, [.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.atom(2).position, [1.5, 0, 0]))

        m2.translate_to_center_of_mass()
        self.assertTrue(numpy.array_equal(m2.atom(1).position, [-.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.atom(2).position, [.5, 0, 0]))

        # H2O molecule:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m3 = qcip_molecule.Molecule(atom_list=atom_list)
        self.assertEqual(m3.mass(), sum(a.mass for a in atom_list))
        self.assertEqual(m3.charge, 0)
        self.assertEqual(m3.number_of_electrons(), 10)
        self.assertEqual(m3.formula(), 'H2O')

        # distances and bonds
        distances = m3.distances()
        self.assertEqual(distances.shape, (3, 3))
        self.assertEqual(distances[0, 0], .0)
        self.assertEqual(distances[0, 1], qcip_math.distance(atom_list[0].position, atom_list[1].position))
        self.assertEqual(distances[0, 2], qcip_math.distance(atom_list[0].position, atom_list[2].position))
        self.assertEqual(distances[1, 2], qcip_math.distance(atom_list[1].position, atom_list[2].position))
        self.assertEqual(distances[1, 0], distances[0, 1])
        self.assertEqual(distances[2, 0], distances[0, 2])
        self.assertEqual(distances[1, 2], distances[2, 1])

        bonds = m3.bonds()
        self.assertEqual(len(bonds), 2)
        self.assertTrue(bonds[0].length == bonds[1].length == distances[0, 1] == distances[0, 2])

    def test_atom_groups(self):
        """Test the behavior or atom groups"""

        # H2O molecule:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[0, 0, -.115]),
            qcip_atom.Atom(symbol='H', position=[0, .767, .460]),
            qcip_atom.Atom(symbol='H', position=[0, -.767, .460])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)

        groups = qcip_molecule.AtomsGroups(m)

        groups.add_group_from_string('A: 1,2')  # defined with shifted indexes
        self.assertEqual(len(groups.groups), 1)
        self.assertEqual(groups.groups[0].name, 'A')
        self.assertEqual(len(groups.groups[0].atom_list), 2)
        self.assertEqual(groups.groups[0].atom_list, [0, 1])  # shifted indexes are not used internally

        # wildcard
        groups.add_group_from_string('B: 2-*')
        self.assertEqual(len(groups.groups), 2)
        self.assertEqual(groups.groups[-1].name, 'B')
        self.assertEqual(len(groups.groups[-1].atom_list), 2)
        self.assertEqual(groups.groups[-1].atom_list, [1, 2])  # two atoms selected by *

        # ranging:
        groups.add_group_from_string('C: 3-1')  # ranging backwards
        self.assertEqual(len(groups.groups), 3)
        self.assertEqual(groups.groups[-1].name, 'C')
        self.assertEqual(len(groups.groups[-1].atom_list), 3)
        self.assertEqual(groups.groups[-1].atom_list, [0, 1, 2])

        groups.add_group_from_string('D: 1-3')
        self.assertEqual(len(groups.groups), 4)
        self.assertEqual(groups.groups[-1].name, 'D')
        self.assertEqual(len(groups.groups[-1].atom_list), 3)
        self.assertEqual(groups.groups[-1].atom_list, [0, 1, 2])

        groups.add_group_from_string('E: 2-3,3')
        self.assertEqual(len(groups.groups), 5)
        self.assertEqual(groups.groups[-1].name, 'E')
        self.assertEqual(len(groups.groups[-1].atom_list), 2)
        self.assertEqual(groups.groups[-1].atom_list, [1, 2])  # remove doubles

        # some invalid stuffs
        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('E: 2,3')  # already exists

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: 1 : 2,3')  # too many ":"

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: *')  # invalid

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: ')  # invalid

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: a')  # invalid

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: 4,3,2')  # invalid initial index

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: 3,2,4')  # invalid initial index

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: 1-4')  # invalid final index

        with self.assertRaises(qcip_molecule.GroupError):
            groups.add_group_from_string('x: 1-2-3')  # invalid range
