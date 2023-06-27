import numpy
# import mendeleev

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
    92: 'U',
    93: 'Np',
    94: 'Pu',
    95: 'Am',
    96: 'Cm',
    97: 'Bk',
    98: 'Cf',
    99: 'Es'
}

ATOMIC_WEIGHTS = {  # from https://iupac.qmul.ac.uk/AtWt/
    'H': 1.008,
    'He': 4.003,
    'Li': 6.940,
    'Be': 9.012,
    'B': 10.810,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'F': 18.998,
    'Ne': 20.180,
    'Na': 22.990,
    'Mg': 24.305,
    'Al': 26.982,
    'Si': 28.085,
    'P': 30.974,
    'S': 32.060,
    'Cl': 35.450,
    'Ar': 39.950,
    'K': 39.098,
    'Ca': 40.078,
    'Sc': 44.956,
    'Ti': 47.867,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.693,
    'Cu': 63.546,
    'Zn': 65.380,
    'Ga': 69.723,
    'Ge': 72.630,
    'As': 74.922,
    'Se': 78.971,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.468,
    'Sr': 87.620,
    'Y': 88.906,
    'Zr': 91.224,
    'Nb': 92.906,
    'Mo': 95.950,
    'Tc': 97.000,
    'Ru': 101.070,
    'Rh': 102.905,
    'Pd': 106.420,
    'Ag': 107.868,
    'Cd': 112.414,
    'In': 114.818,
    'Sn': 118.710,
    'Sb': 121.760,
    'Te': 127.600,
    'I': 126.904,
    'Xe': 131.293,
    'Cs': 132.905,
    'Ba': 137.327,
    'La': 138.905,
    'Ce': 140.116,
    'Pr': 140.908,
    'Nd': 144.242,
    'Pm': 145.000,
    'Sm': 150.360,
    'Eu': 151.964,
    'Gd': 157.250,
    'Tb': 158.925,
    'Dy': 162.500,
    'Ho': 164.930,
    'Er': 167.259,
    'Tm': 168.934,
    'Yb': 173.045,
    'Lu': 174.967,
    'Hf': 178.486,
    'Ta': 180.948,
    'W': 183.840,
    'Re': 186.207,
    'Os': 190.230,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.967,
    'Hg': 200.592,
    'Tl': 204.380,
    'Pb': 207.200,
    'Bi': 208.980,
    'Po': 209.000,
    'At': 210.000,
    'Rn': 222.000,
    'Fr': 223.000,
    'Ra': 226.000,
    'Ac': 227.000,
    'Th': 232.038,
    'Pa': 231.036,
    'U': 238.029,
    'Np': 237.000,
    'Pu': 244.000,
    'Am': 243.000,
    'Cm': 247.000,
    'Bk': 247.000,
    'Cf': 251.000,
    'Es': 252.000,
    'Fm': 257.000,
    'Md': 258.000,
    'No': 259.000,
    'Lr': 262.000,
    'Rf': 267.000,
    'Db': 270.000,
    'Sg': 269.000,
    'Bh': 270.000,
    'Hs': 270.000,
    'Mt': 278.000,
    'Ds': 281.000,
    'Rg': 281.000,
    'Cn': 285.000,
    'Nh': 286.000,
    'Fl': 289.000,
    'Mc': 289.000,
    'Lv': 293.000,
    'Ts': 293.000,
    'Og': 294.000,
}

DUMMY_SYMBOLS = ['Xx', 'X', 'x']

SymbolToAtomicNumber = dict((b, a) for a, b in AtomicNumberToSymbol.items())

_mendeleev_cache = {}


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

            if symbol not in SymbolToAtomicNumber and symbol not in DUMMY_SYMBOLS:
                raise ValueError('{} is not a chemical element'.format(symbol))

            self.symbol = symbol if symbol not in DUMMY_SYMBOLS else 'Xx'
            self.atomic_number = SymbolToAtomicNumber[symbol]

        elif atomic_number is not None:

            if atomic_number > 92:
                raise ValueError('{} is not an allowed Z'.format(atomic_number))

            self.atomic_number = atomic_number if atomic_number > 0 else 0
            self.symbol = AtomicNumberToSymbol[self.atomic_number]

        else:
            raise Exception('either atomic_number or symbol must be defined')

        self.num_of_electrons = 0
        self.mass = 0

        if self.atomic_number > 0:
            self.num_of_electrons = self.atomic_number
            self.mass = ATOMIC_WEIGHTS[self.symbol] if mass is None else mass

        if position is not None:
            if len(position) != 3:
                raise ValueError('len of position must be 3')

            self.position = numpy.array(position)
        else:
            self.position = numpy.zeros(3)

        self.extra = {}

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
