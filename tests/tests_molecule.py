import numpy
import copy
import math

from tests import QcipToolsTestCase

from qcip_tools import atom as qcip_atom, molecule as qcip_molecule, math as qcip_math


class MoleculeTestCase(QcipToolsTestCase):

    def setUp(self):
        pass

    def test_molecule_creation(self):
        """Test the behavior of the Molecule class"""

        # empty atom
        m1 = qcip_molecule.Molecule()
        self.assertEqual(m1.mass(), 0.0)
        self.assertEqual(m1.charge, 0.0)
        self.assertEqual(m1.number_of_electrons(), 0)
        self.assertEqual(m1.formula(), '')
        self.assertEqual(m1.multiplicity, 1)

        # H2 molecule
        atom_list = [
            qcip_atom.Atom(symbol='H', position=[-.5, 0, 0]),
            qcip_atom.Atom(symbol='H', position=[.5, 0, 0])
        ]

        m2 = qcip_molecule.Molecule(atom_list=atom_list)
        self.assertEqual(m2.mass(), sum(a.mass for a in atom_list))
        self.assertEqual(m2.charge, 0)
        self.assertEqual(m2.number_of_electrons(), 2)
        self.assertEqual(m2.formula(), 'H2')
        self.assertEqual(m2.multiplicity, 1)

        self.assertTrue(numpy.array_equal(m2.center_of_mass(), [0, 0, 0]))

        # Test shifted index:
        self.assertTrue(numpy.array_equal(m2.atom(1).position, [-.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m2.atom(2).position, [.5, 0, 0]))

        with self.assertRaises(qcip_molecule.ShiftedIndexError):
            m2.atom(0)

        with self.assertRaises(qcip_molecule.ShiftedIndexError):
            m2.atom(-1)

        with self.assertRaises(qcip_molecule.ShiftedIndexError):
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
        self.assertEqual(m3.multiplicity, 1)

        # atoms()
        self.assertEqual(m3.atoms(atomic_number_in=[1]), [1, 2])
        self.assertEqual(m3.atoms(atomic_number_not_in=[1]), [0])
        self.assertEqual(m3.atoms(atomic_number_not_in=[1, 8]), [])
        self.assertEqual(m3.atoms(symbol_in=['H']), [1, 2])
        self.assertEqual(m3.atoms(symbol_not_in=['H']), [0])
        self.assertEqual(m3.atoms(symbol_not_in=['H'], atomic_number_not_in=[8]), [])

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

        # bounding box
        b1 = m3.bounding_box()

        angle_HOH = qcip_math.angle(m3.atom(2).position, m3.atom(1).position, m3.atom(3).position)
        self.assertArraysAlmostEqual(b1.origin, [0, -.767, -.115])
        # size of the box should be the projection of O-H on Y and H to H distance for X.
        self.assertArraysAlmostEqual(b1.size, [
            0,
            qcip_math.distance(m3.atom(2).position, m3.atom(3).position),
            qcip_math.distance(m3.atom(1).position, m3.atom(2).position) * math.cos(numpy.radians(angle_HOH / 2))
        ])

        extra_space = .5
        b2 = m3.bounding_box(extra_space=extra_space)  # additional space
        self.assertArraysAlmostEqual(b2.origin, [i - extra_space for i in b1.origin])
        self.assertArraysAlmostEqual(b2.size, [i + 2 * extra_space for i in b1.size])
        self.assertArraysAlmostEqual(b2.maximum(), [i + extra_space for i in b1.maximum()])

        # spin and multiplicities
        self.assertEqual(m3.square_of_spin_angular_moment(), .0)  # singlet
        m3.multiplicity = 3
        self.assertEqual(m3.square_of_spin_angular_moment(), 2)  # triplet
        m3.charge = -1
        m3.multiplicity = 2
        self.assertEqual(m3.square_of_spin_angular_moment(), .75)  # doublet

    def test_molecule_modifications(self):
        """Test the manipulation of the atoms"""

        # H2 molecule:
        atom_list = [
            qcip_atom.Atom(symbol='H', position=[-.5, 0, 0]),
            qcip_atom.Atom(symbol='H', position=[.5, 0, 0])
        ]

        m1 = qcip_molecule.Molecule(atom_list=atom_list)
        self.assertTrue(m1.linear())

        # test translate
        m1.translate_self(*[1, 0, 0])
        self.assertTrue(numpy.array_equal(m1.atom(1).position, [.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m1.atom(2).position, [1.5, 0, 0]))

        m1.translate_self_to_center_of_mass()
        self.assertTrue(numpy.array_equal(m1.atom(1).position, [-.5, 0, 0]))
        self.assertTrue(numpy.array_equal(m1.atom(2).position, [.5, 0, 0]))

        # H2O molecule:
        atom_list = [
            qcip_atom.Atom(symbol='O', position=[1, .0, -.115]),
            qcip_atom.Atom(symbol='H', position=[1, -.767, .460]),
            qcip_atom.Atom(symbol='H', position=[1, .767, .460])
        ]

        m2 = qcip_molecule.Molecule(atom_list=atom_list)
        m2.translate_self_to_center_of_mass()  # X has the largest moment of inertia, then Z, then Y
        m2_copy = copy.deepcopy(m2)

        self.assertFalse(m2.linear())

        # set to inertia axes
        m2.set_to_inertia_axes()

        for u in zip(m2, m2_copy):
            good_position = u[1].position[1], u[1].position[2], u[1].position[0]
            self.assertArraysAlmostEqual(u[0].position, good_position)

        # reorient
        m2_copy = copy.deepcopy(m2)
        m2.rotate_around_axis_self(numpy.array([0, 0, 1.]), 90, in_degree=True)  # rotation of 90° around Z

        for i in range(3):  # X and -Y where exchanged
            good_position = -m2_copy[i].position[1], m2_copy[i].position[0], m2_copy[i].position[2]
            self.assertArraysAlmostEqual(m2[i].position, good_position)

    def test_list_of_atoms(self):
        """Test the molecule.list_of_atoms() function"""

        # create a molecule like 1--2--3
        #                           |
        #                           4
        #                           |
        #                           5

        atom_list = [
            qcip_atom.Atom(symbol='C', position=[0, 0, 0]),
            qcip_atom.Atom(symbol='C', position=[1.5, 0, 0]),
            qcip_atom.Atom(symbol='C', position=[3, 0, 0]),
            qcip_atom.Atom(symbol='C', position=[1.5, 1.5, 0]),
            qcip_atom.Atom(symbol='C', position=[1.5, 3, 0])
        ]

        m = qcip_molecule.Molecule(atom_list=atom_list)
        self.assertFalse(m.linear())

        # starting with 2, excluding 3 and 1, should get 2, 4 and 5
        self.assertEqual(m.list_of_atoms(2, exclude_atoms=[1, 3]), [2, 4, 5])

        # starting with 4, excluding 2, should get 4 and 5
        self.assertEqual(m.list_of_atoms(4, exclude_atoms=[2]), [4, 5])

        # starting with 3, excluding 2, should get 3
        self.assertEqual(m.list_of_atoms(3, exclude_atoms=[2]), [3])

        # starting with 2, excluding 1 and 5, should get 2, 3, 4
        self.assertEqual(m.list_of_atoms(2, exclude_atoms=[1, 5]), [2, 3, 4])

        # and the raises
        with self.assertRaises(ValueError):
            m.list_of_atoms(6, exclude_atoms=[1, 2])  # non-existing atom

        with self.assertRaises(ValueError):
            m.list_of_atoms(0, exclude_atoms=[1, 2])  # shifted indexes starts at 1 !

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
