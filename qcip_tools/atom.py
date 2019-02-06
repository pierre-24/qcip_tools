import numpy
import mendeleev

from qcip_tools import transformations

AtomicNumberToSymbol = {
    0: 'Xx',  # dummy atom
    1: 'H',
    2: 'He',
    3: 'Li',
    4: 'Be',
    5: 'B',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    10: 'Ne',
    11: 'Na',
    12: 'Mg',
    13: 'Al',
    14: 'Si',
    15: 'P',
    16: 'S',
    17: 'Cl',
    18: 'Ar',
    19: 'K',
    20: 'Ca',
    21: 'Sc',
    22: 'Ti',
    23: 'V',
    24: 'Cr',
    25: 'Mn',
    26: 'Fe',
    27: 'Co',
    28: 'Ni',
    29: 'Cu',
    30: 'Zn',
    31: 'Ga',
    32: 'Ge',
    33: 'As',
    34: 'Se',
    35: 'Br',
    36: 'Kr',
    37: 'Rb',
    38: 'Sr',
    39: 'Y',
    40: 'Zr',
    41: 'Nb',
    42: 'Mo',
    43: 'Tc',
    44: 'Ru',
    45: 'Rh',
    46: 'Pd',
    47: 'Ag',
    48: 'Cd',
    49: 'In',
    50: 'Sn',
    51: 'Sb',
    52: 'Te',
    53: 'I',
    54: 'Xe',
    55: 'Cs',
    56: 'Ba',
    57: 'La',
    58: 'Ce',
    59: 'Pr',
    60: 'Nd',
    61: 'Pm',
    62: 'Sm',
    63: 'Eu',
    64: 'Gd',
    65: 'Tb',
    66: 'Dy',
    67: 'Ho',
    68: 'Er',
    69: 'Tm',
    70: 'Yb',
    71: 'Lu',
    72: 'Hf',
    73: 'Ta',
    74: 'W',
    75: 'Re',
    76: 'Os',
    77: 'Ir',
    78: 'Pt',
    79: 'Au',
    80: 'Hg',
    81: 'Tl',
    82: 'Pb',
    83: 'Bi',
    84: 'Po',
    85: 'At',
    86: 'Rn',
    87: 'Fr',
    88: 'Ra',
    89: 'Ac',
    90: 'Th',
    91: 'Pa',
    92: 'U'
}

SymbolToAtomicNumber = dict((b, a) for a, b in AtomicNumberToSymbol.items())


class Atom(transformations.MutableTranslatable):
    """
    Create an atom. You must gives either symbol or atomic_number.

    The object is mutable.

    :param atomic_number: the atomic number
    :type atomic_number: int
    :param symbol: the atom symbol (starting with an upercase letter)
    :type symbol: str
    :param position: atom position
    :type position: list|numpy.ndarray
    :param mass: define the mass (for an isotope, for example)
    :type mass: float
    """

    def __init__(self, atomic_number=None, symbol=None, position=None, mass=None):

        if symbol is not None:

            if symbol not in SymbolToAtomicNumber:
                raise ValueError('{} is not a chemical element'.format(symbol))

            self.symbol = symbol
            self.atomic_number = SymbolToAtomicNumber[symbol]

        elif atomic_number is not None:

            if atomic_number < 0 or atomic_number > 92:
                raise ValueError('{} is not an allowed Z'.format(atomic_number))

            self.atomic_number = atomic_number
            self.symbol = AtomicNumberToSymbol[self.atomic_number]

        else:
            raise Exception('either atomic_number or symbol must be defined')

        mendeleev_element = mendeleev.element(self.symbol)

        self.num_of_electrons = self.atomic_number
        self.mass_number = mendeleev_element.mass_number
        self.mass = mendeleev_element.atomic_weight if mass is None else mass

        if position is not None:
            if len(position) != 3:
                raise ValueError('len of position must be 3')

            self.position = numpy.array(position)
        else:
            self.position = numpy.zeros(3)

    def __str__(self):
        return '{} @ ({:.3f}, {:.3f}, {:.3f})'.format(self.symbol, *self.position)

    def is_dummy(self):
        """Return if an atom is dummy

        :rtype: bool
        """
        return self.symbol in ['X', 'Xx']

    def number_of_electrons(self):
        """
        :return: the number of electrons in the atoms
        :rtype: int
        """
        return self.num_of_electrons

    def number_of_core_electrons(self):
        """Return the number of core electrons

        :rtype: int
        """

        prev = 0
        possible_values = [0, 2, 10, 18, 36, 54, 86]

        for core_electrons in possible_values[1:]:

            if core_electrons > self.num_of_electrons:
                return prev

            prev = core_electrons

    def number_of_valence_electrons(self):
        """Get the number of valence electrons

        :rtype: int
        """

        return self.num_of_electrons - self.number_of_core_electrons()

    def number_of_protons(self):
        """
        :return: the number of protons
        :rtype: int
        """
        return self.atomic_number

    def number_of_neutrons(self):
        """
        :return: The number of neutrons
        :rtype: int
        """
        return self.mass_number - self.atomic_number

    def charge(self):
        """Get atomic charge.

        :return: Charge
        :rtype: int
        """
        return self.number_of_protons() - self.num_of_electrons

    def _apply_transformation_self(self, transformation):
        """Appply transformation to atom

        :param transformation: the transformation
        :type transformation: numpy.ndarray
        """

        self.position = transformation.dot([*self.position, 1.])[:3]
