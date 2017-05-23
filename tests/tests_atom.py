import numpy

from tests import QcipToolsTestCase
from qcip_tools import atom as qcip_atom


class AtomTestCase(QcipToolsTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_atom_creation(self):
        """Test the behavior of the atom class"""

        # creation via symbol:
        a1 = qcip_atom.Atom(symbol='He')

        self.assertEqual(a1.atomic_number, 2)
        self.assertEqual(a1.symbol, 'He')
        self.assertEqual(a1.num_of_electrons, 2)
        self.assertEqual(a1.mass_number, 4)
        self.assertTrue(numpy.array_equal(a1.position, numpy.zeros(3)))
        self.assertEqual(a1.number_of_core_electrons(), 2)
        self.assertEqual(a1.number_of_valence_electrons(), 0)

        # creation via atomic number:
        a2 = qcip_atom.Atom(atomic_number=3)

        self.assertEqual(a2.atomic_number, 3)
        self.assertEqual(a2.symbol, 'Li')
        self.assertEqual(a2.num_of_electrons, 3)
        self.assertEqual(a2.mass_number, 7)
        self.assertTrue(numpy.array_equal(a2.position, numpy.zeros(3)))
        self.assertEqual(a2.number_of_core_electrons(), 2)
        self.assertEqual(a2.number_of_valence_electrons(), 1)

        # test with nothing:
        with self.assertRaises(Exception):
            qcip_atom.Atom()

        # test with non-existing atoms:
        with self.assertRaises(ValueError):
            qcip_atom.Atom(symbol='xxx')
        with self.assertRaises(ValueError):
            qcip_atom.Atom(atomic_number=-2)
        with self.assertRaises(ValueError):
            qcip_atom.Atom(atomic_number=185)

        # test with position:
        coordinates = [0, 1, 2]
        a3 = qcip_atom.Atom(symbol='Be', position=coordinates)
        self.assertEqual(a3.number_of_core_electrons(), 2)
        self.assertEqual(a3.number_of_valence_electrons(), 2)
        self.assertTrue(numpy.array_equal(a3.position, coordinates))

        with self.assertRaises(ValueError):
            qcip_atom.Atom(symbol='Be', position=[2, 2])  # incorrect size for position

        # test representation:
        r = str(a3)
        self.assertTrue(a3.symbol in r)

        # test functions:
        self.assertEqual(a3.mass_number, 9)
        self.assertEqual(a3.charge(), 0)
        self.assertEqual(a3.number_of_protons(), 4)
        self.assertEqual(a3.number_of_electrons(), 4)
        self.assertEqual(a3.number_of_neutrons(), 5)

        # change charge:
        a3.num_of_electrons = 3
        self.assertEqual(a3.charge(), 1)
        self.assertEqual(a3.number_of_protons(), 4)
        self.assertEqual(a3.number_of_electrons(), 3)

        # special test for core/valence electrons
        atms = [('Cl', 10, 7), ('Ne', 10, 0), ('I', 36, 17)]

        for atm in atms:
            ax = qcip_atom.Atom(symbol=atm[0])
            self.assertEqual(ax.number_of_core_electrons(), atm[1], atm[0])
            self.assertEqual(ax.number_of_valence_electrons(), atm[2], atm[0])
