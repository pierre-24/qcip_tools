import math

import numpy

from qcip_tools import atom as qcip_atom, math as qcip_math


class Bond:
    """Define a bond between two atoms"""

    atom_1 = None
    atom_2 = None
    type = 'None'
    length = 0.0
    index1 = 0
    index2 = 0
    vector = []

    def __init__(self, atom_1, atom_2, bond_type='single', indexes=None):
        """Create a bond.

        :param atom_1: atom
        :type atom_1:  qcip_tools.atom.Atom
        :param atom_2: atom
        :type atom_2:  qcip_tools.atom.Atom
        :param bond_type: single, double, triple, ...
        :type bond_type: str
        :param indexes: index of the corresponding atoms in the molecule, if any
        :type indexes: list
        """

        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.type = bond_type

        if indexes and type(indexes) is list and len(indexes) == 2:
            self.index1 = indexes[0]
            self.index2 = indexes[1]

        self.vector = numpy.array(self.atom_2.position) - numpy.array(self.atom_1.position)
        self.length = numpy.linalg.norm(self.vector)

    def __str__(self):
        return 'bond between {} and {}'.format(self.atom_1, self.atom_2)


class Molecule:
    """
    Class to define a molecule (basically a list of atoms).

    .. note ::
        Each time a function requires a ``shifted_index``, indexing starts at 1 instead of 0.
    """

    def __init__(self, atom_list=None, charge=0.0):
        """Create a molecule

        :param atom_list: list of atoms
        :type atom_list: list of qcip_tools.atom.Atom
        :param charge: global charge of the molecule
        :type charge: float
        """

        self.atom_list = atom_list if atom_list else []
        self.symbols_contained = []
        self.charge = charge

        for a in self.atom_list:
            if a.symbol not in self.symbols_contained:
                self.symbols_contained.append(a.symbol)

    def __str__(self):
        return self.formula()

    def __iter__(self):
        """
        Allow a molecule to act as an iterable
        :return: an iterable of the atom list
        """

        for a in self.atom_list:
            yield a

    def __getitem__(self, item):
        return self.atom_list[item]

    def __len__(self):
        return len(self.atom_list)

    def __contains__(self, item):
        if 0 <= item < len(self.atom_list):
            return True
        else:
            return False

    def insert(self, atom, position=None):
        """
        Add atom in the list.

        :param atom: an atom (raise Exception if it's not the case)
        :type atom: qcip_tools.atom.Atom
        :param position: position in the list.
        :type position: int
        """

        if type(atom) is not qcip_atom.Atom:
            raise TypeError(atom)

        if position is None:
            self.atom_list.append(atom)
        else:
            self.atom_list.insert(position, atom)

        if atom.symbol not in self.symbols_contained:
            self.symbols_contained.append(atom.symbol)

    def remove_atom(self, shifted_index):
        """
        Remove an atom in the list
        :param shifted_index: index of the atom
        """

        if len(self) >= shifted_index > 0:
            self.atom_list.pop(shifted_index - 1)
        else:
            raise Exception('The atom shifted_index (' + str(shifted_index) + ') is out of range')

    def atom(self, shifted_index):
        """
        Get the atom defined by index

        :param shifted_index: Index of the atom in the list
        :type shifted_index: int
        :rtype: qcip_atom.Atom
        """

        if len(self) >= shifted_index > 0:
            return self.atom_list[shifted_index - 1]
        else:
            raise KeyError(shifted_index)

    def number_of_electrons(self):
        """
        :return: The number of electrons in the molecule
        """

        num = 0
        for a in self:
            num += a.number_of_electrons()

        return num - self.charge

    def MO_to_str(self, mo_level):
        """Get the corresponding Molecular Orbital (MO) as HOMO-x/LUMO+x.

        .. note ::
            Assume closed-shell.

        :param mo_level: MO number (starting at 1)
        :type mo_level: int
        :rtype: str
        """

        homo_level = math.ceil(self.number_of_electrons() / 2)

        if mo_level <= homo_level:
            r = 'HOMO'
        else:
            r = 'LUMO'

        diff = mo_level - homo_level

        if diff not in [0, 1]:
            if diff > 1:
                diff -= 1

            r = '{}{:+d}'.format(r, diff)

        return r

    def mass(self):
        """Get the mass of the molecule

        :rtype: float
        """
        mass = 0
        for a in self:
            mass += a.mass
        return mass

    def center_of_mass(self):
        """Position of the center of mass

        :rtype: numpy.ndarray
        """

        com = numpy.array([0, 0, 0])
        for a in self:
            com = com + a.mass * a.position
        return 1 / self.mass() * com

    def center_of_charges(self):
        """Position of the center of charge

        :rtype: numpy.ndarray
        """

        com = numpy.array([0, 0, 0])
        all_atomic_numbers = 0
        for a in self:
            com = com + a.atomic_number * a.position
            all_atomic_numbers += a.atomic_number
        return 1 / all_atomic_numbers * com

    def translate(self, coordinates):
        """
        Translate the molecule (each of its atom).

        :param coordinates: new position of the molecule
        :type coordinates: list|numpy.ndarray
        """

        coordinates = numpy.array(coordinates, float)
        for a in self:
            a.position = a.position + coordinates

    def translate_to_center_of_mass(self):
        """
        Translate the molecule to the center of mass
        """
        self.translate(-self.center_of_mass())

    def translate_to_center_of_charges(self):
        """
        Translate the molecule to the center of charge
        """
        self.translate(-self.center_of_charges())

    def distances(self):
        """
        Retrieve the distance matrix between all atoms in the molecule.

        :rtype: numpy.ndarray
        """

        distance_matrix = numpy.zeros((len(self), len(self)))
        for x, a1 in enumerate(self):
            for y, a2 in enumerate(self):
                if y < x:
                    distance_matrix[x, y] = distance_matrix[y, x] = qcip_math.distance(a1.position, a2.position)
                else:
                    break

        return distance_matrix

    def connectivities(self, threshold=.1):
        """
        Get connectivity atom by atom, using the VdW radii.

        :param threshold: threshold for acceptance (because bond are sometimes a little larger than the sum of VdWs)
        :type threshold: float
        :return: the connectivity, as a dictionary per atom index (starting with one !)
        :rtype: dict
        """
        distances = self.distances()
        connectivities = {}
        for i, a1 in enumerate(self):
            tmp = []
            covalent_radius1 = qcip_atom.CovalentRadii[a1.symbol][0]  # assume single bond
            for j, a2 in enumerate(self):
                if i is not j:  # an atom cannot bond with himself
                    covalent_radius2 = qcip_atom.CovalentRadii[a2.symbol][0]  # assume single bond
                    if distances[i, j] < (covalent_radius1 + covalent_radius2 + threshold):
                        tmp.append(j + 1)  # there is a bond
            connectivities[i + 1] = tmp

        return connectivities

    def bonds(self, threshold=.1):
        """Get a list of the bond in the molecule.

        :param threshold:  threshold (see `self. connectivites()`)
        :type threshold: float
        :rtype: list
        """

        connectivity = self.connectivities(threshold=threshold)
        bonds_handler = []
        for i in range(0, len(self)):
            atm = self.atom(i + 1)
            if connectivity[i + 1]:
                for index in connectivity[i + 1]:
                    if index > i:
                        bonds_handler.append(Bond(atm, self.atom(index), 'single', [i + 1, index]))
        return bonds_handler

    def formula(self, enhanced=False):
        """Get the formula, based on self.formula_holder and self. symbols_contained.
        Implementation note : logically, it can write the atoms in whatever order. But in chemistry,
        the usual way is to start in this order C,H,N,O,X ... Where X is the other(s) atom(s). So
        does this. "enhanced" use "X_{i}" to indicate indice i instead of just "Xi".
        """

        all_atoms = ''.join(a.symbol for a in self)

        first_atoms = ['C', 'N', 'O', 'H']
        s = ''
        for symbol in first_atoms:
            if symbol in self.symbols_contained:
                s += symbol
                num_s = all_atoms.count(symbol)

                if num_s > 1:
                    if enhanced:
                        s += '_{'
                    s += str(num_s)
                    if enhanced:
                        s += '}'

        for symbol in self.symbols_contained:
            num_s = all_atoms.count(symbol)

            if symbol not in first_atoms and num_s:
                s += symbol
                if num_s > 1:
                    if enhanced:
                        s += '_{'
                    s += str(num_s)
                    if enhanced:
                        s += '}'

        if s == 'OH2':
            s = 'H2O'

        if s == 'OH_{2}':
            s = 'H_{2}O'

        return s

    def output_atoms(self, use_z_instead_of_symbol=False):
        """List atoms.

        :param use_z_instead_of_symbol: put atomic number instead of symbol
        :type use_z_instead_of_symbol: bool
        :return: formatted string `<symbol or number> <x> <y> <z>`
        :rtype: str
        """

        atoms_string = ''
        for a in self:
            atoms_string += '{:5} {:12.8f} {:12.8f} {:12.8f}\n'.format(
                a.symbol if not use_z_instead_of_symbol else a.atomic_number, *a.position)

        return atoms_string

    def output_as_xyz(self, title='', use_z_instead_of_symbol=False):
        """Return an XYZ content.

        :param title: optional name of the molecule
        :type title: str
        :param use_z_instead_of_symbol: put atomic number instead of symbol
        :type use_z_instead_of_symbol: bool
        :return: Molecule formatted for an xyz file
        :rtype: str
        """

        rstr = str(len(self)) + '\n'
        rstr += title + '\n'
        rstr += self.output_atoms(use_z_instead_of_symbol)
        return rstr

    def output_as_dalton(self, nosym=False):
        """Return atoms as in a Dalton molecule file.

        :param nosym: set nosym !
        :type nosym: bool
        :rtype: str
        """

        rstr = 'Atomtypes={}  Angstrom{}\n'.format(len(self), ' Nosymmetry' if nosym else '')

        for a in self:
            rstr += 'Charge={:.1f} Atoms={}\n'.format(a.atomic_number, 1)
            rstr += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(a.symbol, *a.position)

        return rstr

    def moments_of_inertia(self):
        """Get the moment of inertia tensor.

        :rtype: numpy.ndarray
        """

        tensor = numpy.zeros((3, 3))

        coordinates, masses, r = \
            numpy.zeros((len(self), 3)), \
            numpy.zeros(len(self)), \
            numpy.zeros(len(self))

        for index, atm in enumerate(self):
            coordinates[index] = atm.coordinates
            r[index] = numpy.linalg.norm(atm.coordinates)
            masses[index] = atm.mass

        for i in [0, 1, 2]:
            for j in [0, 1, 2]:
                tmp = 0
                for k in range(0, len(self)):
                    delta = i == j
                    tmp += masses[k] * (delta * r[k] ** 2 - coordinates[k][i] * coordinates[k][j])
                tensor[i, j] = tmp

        return tensor

    def principal_axes(self):
        """Return the moment of inertia
        (from https://github.com/moorepants/DynamicistToolKit/blob/master/dtk/inertia.py).

        :return: The principal moment of inertia (sorted from lowest to largest) and the rotation matrix
        """

        inertia_tensor = self.moments_of_inertia()
        Ip, C = numpy.linalg.eig(inertia_tensor)
        indices = numpy.argsort(Ip)
        Ip = Ip[indices]
        C = C.T[indices]
        return Ip, C

    def bounding_box(self, extra_space=0.0):
        """
        Note: may be a good idea to perform a `geomoperation.set_to_inertia_axis()` before, it may reduce the
        volume of the bounding bow

        :return: (origin, size) tuple
        :rtype: tuple
        """

        low = numpy.array([0, 0, 0])
        high = numpy.array([0, 0, 0])
        for atm in self:
            coo = atm.position
            for i in [0, 1, 2]:
                if coo[i] < low[i]:
                    low[i] = coo[i]
                if coo[i] > high[i]:
                    high[i] = coo[i]

        extra_space_array = numpy.zeros(3)
        if type(extra_space) is float:
            extra_space_array = numpy.array([extra_space, extra_space, extra_space])
        elif type(extra_space) is numpy.ndarray:
            extra_space_array = extra_space
        size = high - low + 2 * extra_space_array
        return low - extra_space_array, size


class GroupError(Exception):
    pass


class AtomsGroups:
    def __init__(self, molecule):
        """

        :param molecule: geometry
        :type molecule: Molecule
        """

        self.molecule = molecule
        self.groups = []

    def add_group(self, name, atom_list):
        """Create a group of atoms

        :param name: name of the group
        :type name: str
        :param atom_list: list of atoms in the group
        :type atom_list: list
        """

        self.groups.append(AtomsGroup(name, atom_list))

    def add_group_from_string(self, group_str):
        """Add a group from a string corresponding to

        `` name : [[number(-number|*)?], ...]``

        :param group_str: the string that defines the group
        :type group_str: str
        """

        if ':' not in group_str:
            raise GroupError('no name!')
        li = group_str.split(':')

        if len(li) != 2:
            raise GroupError('too much groups')

        group_name = li[0].strip()

        for g in self.groups:
            if g.name == group_name:
                raise GroupError('{} already defines a group'.format(group_name))

        subgroups_list = li[1].split(',')

        group_def = []

        for subggroup in subgroups_list:
            subggroup = subggroup.strip()

            if '-' in subggroup:  # it's a range
                s = subggroup.split('-')
                if len(s) != 2:
                    raise GroupError('{} is not a range definition'.format(subggroup))
                try:
                    first_atom = int(s[0]) - 1
                except ValueError:
                    raise GroupError('{} is not an atom definition'.format(s[0]))

                if first_atom < 0 or first_atom > len(self.molecule) - 1:
                    raise GroupError('atom {} is not a valid atom index'.format(s[0]))

                if s[1] == '*':
                    last_atom = len(self.molecule) - 1
                else:
                    try:
                        last_atom = int(s[1]) - 1
                    except ValueError:
                        raise GroupError('{} is not an atom definition'.format(s[1]))
                    if last_atom < 0 or last_atom > len(self.molecule) - 1:
                        raise GroupError('atom {} is not a valid atom index'.format(s[1]))

                if last_atom < first_atom:
                    last_atom, first_atom = first_atom, last_atom

                group_def.extend(range(first_atom, last_atom + 1))

            else:  # it's a single atom
                try:
                    atom = int(subggroup) - 1
                except ValueError:
                    raise GroupError('{} is not an atom definition'.format(subggroup))
                if atom < 0 or atom > len(self.molecule) - 1:
                    raise GroupError('atom {} is not a valid atom index'.format(subggroup))

                group_def.append(atom)

        self.add_group(group_name, list(set(group_def)))  # uniqueness


class AtomsGroup:
    def __init__(self, name, atom_list):
        """Create a group of atoms

        :param name: name of the group
        :type name: str
        :param atom_list: list of atoms in the group
        :type atom_list: list of int
        """

        self.name = name
        self.atom_list = atom_list

    def __contains__(self, item):
        return item in self.atom_list
