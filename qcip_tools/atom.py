import numpy

# `Definition` is given as follow :
#   0. atomic number, Z (int) ;
#   1. mass number, A (int) ;
#   2. relative atomic mass (float)
# (by definition, the most stable isotope is defined, for radioactive element, it's more an arbitrary choice)
# Obtained via http://ciaaw.org/pubs/TSAW%202005.pdf (Pure Appl. Chem., Vol. 78, No. 11, pp. 2051–2066, 2006)

Definition = {
    'X': [0, 0, 0.0],
    'H': [1, 1, 1.00794],
    'He': [2, 4, 4.002602],
    'Li': [3, 7, 6.941],
    'Be': [4, 9, 9.012182],
    'B': [5, 11, 10.811],
    'C': [6, 12, 12.0107],
    'N': [7, 14, 14.0067],
    'O': [8, 16, 15.9994],
    'F': [9, 19, 18.9984032],
    'Ne': [10, 20, 20.1797],
    'Na': [11, 23, 22.98976928],
    'Mg': [12, 24, 24.305],
    'Al': [13, 27, 26.9815386],
    'Si': [14, 28, 28.0855],
    'P': [15, 31, 30.973762],
    'S': [16, 32, 32.065],
    'Cl': [17, 35, 35.453],
    'Ar': [18, 40, 39.948],
    'K': [19, 39, 39.0983],
    'Ca': [20, 40, 40.078],
    'Sc': [21, 45, 44.955912],
    'Ti': [22, 48, 47.867],
    'V': [23, 51, 50.9415],
    'Cr': [24, 52, 51.9961],
    'Mn': [25, 55, 54.938045],
    'Fe': [26, 56, 55.845],
    'Co': [27, 59, 58.933195],
    'Ni': [28, 59, 58.6934],
    'Cu': [29, 64, 63.546],
    'Zn': [30, 65, 65.409],
    'Ga': [31, 70, 69.723],
    'Ge': [32, 73, 72.64],
    'As': [33, 75, 74.9216],
    'Se': [34, 79, 78.96],
    'Br': [35, 80, 79.904],
    'Kr': [36, 84, 83.798],
    'Rb': [37, 85, 85.4678],
    'Sr': [38, 88, 87.62],
    'Y': [39, 89, 88.90585],
    'Zr': [40, 91, 91.224],
    'Nb': [41, 93, 92.90638],
    'Mo': [42, 96, 95.94],
    'Tc': [43, 97, 96.9064],
    'Ru': [44, 101, 101.07],
    'Rh': [45, 103, 102.9055],
    'Pd': [46, 106, 106.42],
    'Ag': [47, 108, 107.8682],
    'Cd': [48, 112, 112.411],
    'In': [49, 115, 114.818],
    'Sn': [50, 119, 118.71],
    'Sb': [51, 122, 121.76],
    'Te': [52, 128, 127.6],
    'I': [53, 127, 126.90447],
    'Xe': [54, 131, 131.293],
    'Cs': [55, 133, 132.9054519],
    'Ba': [56, 137, 137.327],
    'La': [57, 139, 138.90547],
    'Ce': [58, 140, 140.116],
    'Pr': [59, 141, 140.90765],
    'Nd': [60, 144, 144.242],
    'Pm': [61, 147, 146.9151],
    'Sm': [62, 150, 150.36],
    'Eu': [63, 152, 151.964],
    'Gd': [64, 157, 157.25],
    'Tb': [65, 159, 158.92535],
    'Dy': [66, 163, 162.5],
    'Ho': [67, 165, 164.9332],
    'Er': [68, 167, 167.259],
    'Tm': [69, 169, 168.93421],
    'Yb': [70, 173, 173.04],
    'Lu': [71, 175, 174.967],
    'Hf': [72, 178, 178.49],
    'Ta': [73, 181, 180.94788],
    'W': [74, 184, 183.84],
    'Re': [75, 181, 180.94788],
    'Os': [76, 190, 190.23],
    'Ir': [77, 192, 192.217],
    'Pt': [78, 195, 195.084],
    'Au': [79, 197, 196.966569],
    'Hg': [80, 201, 200.59],
    'Tl': [81, 204, 204.3833],
    'Pb': [82, 207, 207.2],
    'Bi': [83, 209, 208.984],
    'Po': [84, 209, 209],
    'At': [85, 210, 210],
    'Rn': [86, 222, 222],
    'Fr': [87, 223, 223],
    'Ra': [88, 226, 226.0254],
    'Ac': [89, 227, 227],
    'Th': [90, 232, 232.03806],
    'Pa': [91, 231, 231.03588],
    'U': [92, 238, 238.02891]
}

AtomicNumberToSymbol = {
    0: 'X',  # dummy atoms
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

# defined using `Cordeo et al. Covalent radii revisited, Dalton Trans., 2008, 2832–2838`
# [currently only for single bond]
CovalentRadii = {
    'H': [0.31],
    'He': [0.28],
    'Li': [1.28],
    'Be': [0.96],
    'B': [0.84],
    'C': [0.76],
    'N': [0.71],
    'O': [0.66],
    'F': [0.57],
    'Ne': [0.58],
    'Na': [1.66],
    'Mg': [1.41],
    'Al': [1.21],
    'Si': [1.11],
    'P': [1.07],
    'S': [1.05],
    'Cl': [1.02],
    'Ar': [1.06],
    'K': [2.03],
    'Ca': [1.76],
    'Sc': [1.7],
    'Ti': [1.6],
    'V': [1.53],
    'Cr': [1.39],
    'Mn': [1.61],
    'Fe': [1.52],
    'Co': [1.50],
    'Ni': [1.24],
    'Cu': [1.32],
    'Zn': [1.22],
    'Ga': [1.22],
    'Ge': [1.2],
    'As': [1.19],
    'Se': [1.2],
    'Br': [1.2],
    'Kr': [1.16],
    'Rb': [2.2],
    'Sr': [1.95],
    'Y': [1.9],
    'Zr': [1.75],
    'Nb': [1.64],
    'Mo': [1.54],
    'Tc': [1.47],
    'Ru': [1.46],
    'Rh': [1.42],
    'Pd': [1.39],
    'Ag': [1.45],
    'Cd': [1.44],
    'In': [1.42],
    'Sn': [1.39],
    'Sb': [1.39],
    'Te': [1.38],
    'I': [1.39],
    'Xe': [1.4],
    'Cs': [2.44],
    'Ba': [2.15],
    'La': [2.07],
    'Ce': [2.04],
    'Pr': [2.03],
    'Nd': [2.01],
    'Pm': [1.99],
    'Sm': [1.98],
    'Eu': [1.98],
    'Gd': [1.96],
    'Tb': [1.94],
    'Dy': [1.92],
    'Ho': [1.92],
    'Er': [1.89],
    'Tm': [1.9],
    'Yb': [1.87],
    'Lu': [1.87],
    'Hf': [1.75],
    'Ta': [1.7],
    'W': [1.62],
    'Re': [1.51],
    'Os': [1.44],
    'Ir': [1.41],
    'Pt': [1.36],
    'Au': [1.36],
    'Hg': [1.32],
    'Tl': [1.45],
    'Pb': [1.46],
    'Bi': [1.48],
    'Po': [1.4],
    'At': [1.5],
    'Rn': [1.5],
    'Fr': [2.6],
    'Ra': [2.21],
    'Ac': [2.15],
    'Th': [2.06],
    'Pa': [2],
    'U': [1.96]
}

# colors defined by JMol (http://jmol.sourceforge.net/jscolors/)
# [converted from hexadecimal to floating point]
Colors = {
    'X': [.75, .75, .75],
    'H': [1.0000, 1.0000, 1.0000],
    'He': [0.8510, 1.0000, 1.0000],
    'Li': [0.8000, 0.5020, 1.0000],
    'Be': [0.7608, 1.0000, 0.0000],
    'B': [1.0000, 0.7098, 0.7098],
    'C': [0.5647, 0.5647, 0.5647],
    'N': [0.1882, 0.3137, 0.9725],
    'O': [1.0000, 0.0510, 0.0510],
    'F': [0.5647, 0.8784, 0.3137],
    'Ne': [0.7020, 0.8902, 0.9608],
    'Na': [0.6706, 0.3608, 0.9490],
    'Mg': [0.5412, 1.0000, 0.0000],
    'Al': [0.7490, 0.6510, 0.6510],
    'Si': [0.9412, 0.7843, 0.6275],
    'P': [1.0000, 0.5020, 0.0000],
    'S': [1.0000, 1.0000, 0.1882],
    'Cl': [0.1216, 0.9412, 0.1216],
    'Ar': [0.5020, 0.8196, 0.8902],
    'K': [0.5608, 0.2510, 0.8314],
    'Ca': [0.2392, 1.0000, 0.0000],
    'Sc': [0.9020, 0.9020, 0.9020],
    'Ti': [0.7490, 0.7608, 0.7804],
    'V': [0.6510, 0.6510, 0.6706],
    'Cr': [0.5412, 0.6000, 0.7804],
    'Mn': [0.6118, 0.4784, 0.7804],
    'Fe': [0.8784, 0.4000, 0.2000],
    'Co': [0.9412, 0.5647, 0.6275],
    'Ni': [0.3137, 0.8157, 0.3137],
    'Cu': [0.7843, 0.5020, 0.2000],
    'Zn': [0.4902, 0.5020, 0.6902],
    'Ga': [0.7608, 0.5608, 0.5608],
    'Ge': [0.4000, 0.5608, 0.5608],
    'As': [0.7412, 0.5020, 0.8902],
    'Se': [1.0000, 0.6314, 0.0000],
    'Br': [0.6510, 0.1608, 0.1608],
    'Kr': [0.3608, 0.7216, 0.8196],
    'Rb': [0.4392, 0.1804, 0.6902],
    'Sr': [0.0000, 1.0000, 0.0000],
    'Y': [0.5804, 1.0000, 1.0000],
    'Zr': [0.5804, 0.8784, 0.8784],
    'Nb': [0.4510, 0.7608, 0.7882],
    'Mo': [0.3294, 0.7098, 0.7098],
    'Tc': [0.2314, 0.6196, 0.6196],
    'Ru': [0.1412, 0.5608, 0.5608],
    'Rh': [0.0392, 0.4902, 0.5490],
    'Pd': [0.0000, 0.4118, 0.5216],
    'Ag': [0.7529, 0.7529, 0.7529],
    'Cd': [1.0000, 0.8510, 0.5608],
    'In': [0.6510, 0.4588, 0.4510],
    'Sn': [0.4000, 0.5020, 0.5020],
    'Sb': [0.6196, 0.3882, 0.7098],
    'Te': [0.8314, 0.4784, 0.0000],
    'I': [0.5804, 0.0000, 0.5804],
    'Xe': [0.2588, 0.6196, 0.6902],
    'Cs': [0.3412, 0.0902, 0.5608],
    'Ba': [0.0000, 0.7882, 0.0000],
    'La': [0.4392, 0.8314, 1.0000],
    'Ce': [1.0000, 1.0000, 0.7804],
    'Pr': [0.8510, 1.0000, 0.7804],
    'Nd': [0.7804, 1.0000, 0.7804],
    'Pm': [0.6392, 1.0000, 0.7804],
    'Sm': [0.5608, 1.0000, 0.7804],
    'Eu': [0.3804, 1.0000, 0.7804],
    'Gd': [0.2706, 1.0000, 0.7804],
    'Tb': [0.1882, 1.0000, 0.7804],
    'Dy': [0.1216, 1.0000, 0.7804],
    'Ho': [0.0000, 1.0000, 0.6118],
    'Er': [0.0000, 0.9020, 0.4588],
    'Tm': [0.0000, 0.8314, 0.3216],
    'Yb': [0.0000, 0.7490, 0.2196],
    'Lu': [0.0000, 0.6706, 0.1412],
    'Hf': [0.3020, 0.7608, 1.0000],
    'Ta': [0.3020, 0.6510, 1.0000],
    'W': [0.1294, 0.5804, 0.8392],
    'Re': [0.1490, 0.4902, 0.6706],
    'Os': [0.1490, 0.4000, 0.5882],
    'Ir': [0.0902, 0.3294, 0.5294],
    'Pt': [0.8157, 0.8157, 0.8784],
    'Au': [1.0000, 0.8196, 0.1373],
    'Hg': [0.7216, 0.7216, 0.8157],
    'Tl': [0.6510, 0.3294, 0.3020],
    'Pb': [0.3412, 0.3490, 0.3804],
    'Bi': [0.6196, 0.3098, 0.7098],
    'Po': [0.6706, 0.3608, 0.0000],
    'At': [0.4588, 0.3098, 0.2706],
    'Rn': [0.2588, 0.5098, 0.5882],
    'Fr': [0.2588, 0.0000, 0.4000],
    'Ra': [0.0000, 0.4902, 0.0000],
    'Ac': [0.4392, 0.6706, 0.9804],
    'Th': [0.0000, 0.7294, 1.0000],
    'Pa': [0.0000, 0.6314, 1.0000],
    'U': [0.0000, 0.5608, 1.0000]
}


class Atom:
    """
    Create an atom. You must gives either symbol or atomic_number.

    :param atomic_number: the atomic number
    :type atomic_number: int
    :param symbol: the atom symbol (starting with an upercase letter)
    :type symbol: str
    :param position: atom position
    :type position: list|numpy.ndarray
    """

    def __init__(self, atomic_number=None, symbol=None, position=None):

        if symbol is not None:

            if symbol not in Definition.keys():
                raise ValueError('{} is not a chemical element'.format(symbol))

            self.symbol = symbol
            self.atomic_number = Definition[symbol][0]

        elif atomic_number is not None:

            if atomic_number < 0 or atomic_number > 92:
                raise ValueError('{} is not an allowed Z'.format(atomic_number))

            self.atomic_number = atomic_number
            self.symbol = AtomicNumberToSymbol[self.atomic_number]

        else:
            raise Exception('either atomic_number or symbol must be defined')

        self.num_of_electrons = self.atomic_number
        self.mass_number = Definition[self.symbol][1]
        self.mass = Definition[self.symbol][2]

        if position is not None:
            if len(position) != 3:
                raise ValueError('len of position must be 3')

            self.position = numpy.array(position)
        else:
            self.position = numpy.zeros(3)

    def __str__(self):
        return '{} @ ({:.3f}, {:.3f}, {:.3f})'.format(self.symbol, *self.position)

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
