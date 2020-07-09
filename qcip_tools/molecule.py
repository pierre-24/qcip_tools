import math

import numpy
import mendeleev
import random

from qcip_tools import bounding, transformations, symmetry
from qcip_tools import atom as qcip_atom, ValueOutsideDomain


class Bond:
    """Define a bond between two atoms"""

    atom_1 = None
    atom_2 = None
    type = 'None'
    length = 0.0
    index1 = 0  # shifted_index
    index2 = 0  # shifted_index
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


class ShiftedIndexError(ValueOutsideDomain):
    def __init__(self, val, mol):
        super().__init__(val, 1, len(mol))


class Molecule(transformations.MutableTranslatable, transformations.MutableRotatable):
    """
    Class to define a molecule (basically a list of atoms).

    This object is mutable.

    .. note ::

        Each time a function requires a ``shifted_index``, indexing starts at 1 instead of 0.

    """

    def __init__(self, atom_list=None, charge=0):
        """Create a molecule

        :param atom_list: list of atoms
        :type atom_list: list of qcip_tools.atom.Atom
        :param charge: global charge of the molecule
        :type charge: float
        """

        self.atom_list = []
        self.symbols_contained = []
        self.charge = charge

        self.point_group = symmetry.PointGroup.C_n(1)

        if atom_list is not None:
            for a in atom_list:
                self.insert(a)

        self.multiplicity = 1 if self.number_of_electrons() % 2 == 0 else 2

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

    def _apply_transformation_self(self, transformation):
        """Appply transformation to the molecule (broadcast to each atom)

        :param transformation: the transformation
        :type transformation: numpy.ndarray
        """

        for a in self.atom_list:
            a._apply_transformation_self(transformation)

    def insert(self, atom, position=None):
        """
        Add atom in the list.

        :param atom: an atom (raise Exception if it's not the case)
        :type atom: qcip_tools.atom.Atom|int|str
        :param position: position in the list.
        :type position: int
        """

        if type(atom) is int:
            atom = qcip_atom.Atom(atomic_number=atom)
        elif type(atom) is str:
            atom = qcip_atom.Atom(symbol=atom)
        elif type(atom) is not qcip_atom.Atom:
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
            raise ShiftedIndexError(shifted_index, self)

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
            raise ShiftedIndexError(shifted_index, self)

    def atoms(self, **kwargs):
        """Get a list of atom (index) based on some criterion"""

        atoms = []

        symbols_in = kwargs.get('symbol_in', None)
        symbols_not_in = kwargs.get('symbol_not_in', None)
        atomic_number_in = kwargs.get('atomic_number_in', None)
        atomic_number_not_in = kwargs.get('atomic_number_not_in', None)

        for index, a in enumerate(self):
            if symbols_in and a.symbol not in symbols_in:
                continue
            if symbols_not_in and a.symbol in symbols_not_in:
                continue
            if atomic_number_in and a.atomic_number not in atomic_number_in:
                continue
            if atomic_number_not_in and a.atomic_number in atomic_number_not_in:
                continue

            atoms.append(index)

        return atoms

    def list_of_atoms(self, shifted_index, exclude_atoms=None):
        """Make a list of atoms following the connectivity and excluding some if any.

        .. note:;
            It may seems obvious, but the first atom IS included in the list.

        :param shifted_index: starting atom
        :type shifted_index: int
        :param exclude_atoms: stop following atom (listed as shifted_indexes) if on one of those
        :type exclude_atoms: list of int
        :return: list of atoms (as shifted_indexes)
        :rtype: list of int
        """

        def lookup(current_atom, exclude_atoms, added_atom_indexes):
            for atom_index in connectivity[current_atom]:
                if atom_index not in added_atom_indexes:
                    if atom_index in exclude_atoms:
                        continue
                    added_atom_indexes.append(atom_index)

                    lookup(atom_index, exclude_atoms, added_atom_indexes)

        if exclude_atoms is None:
            return [a for a in range(len(self))]

        if shifted_index < 1 or shifted_index > len(self):
            raise ShiftedIndexError(shifted_index, self)

        connectivity = self.connectivities()

        atom_indexes = [shifted_index]
        lookup(shifted_index, exclude_atoms, atom_indexes)  # recursive

        return atom_indexes

    def position_matrix(self, with_='Z', randomise_dummy=True):
        """Get position matrix

        :param with_: extra first column, either ``Z``, ``masss`` or ``charge``
        :type with_: str
        :param randomise_dummy: use a random number, outside of range, for dummy atoms
        :type randomise_dummy: bool
        :rtype: numpy.ndarray
        """

        m = []
        for a in self:
            i = .0
            if with_ == 'Z':
                i = a.atomic_number
            elif with_ == 'mass':
                i = a.mass
            elif with_ == 'charge':
                i = a.charge()

            if randomise_dummy and a.is_dummy():
                i = random.randrange(0, 200.)

            m.append([i, *a.position])

        return numpy.vstack(m)

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

        m = self.position_matrix(with_='mass')
        return numpy.einsum('i,ij->j', m[:, 0], m[:, 1:]) / m[:, 0].sum()

    def center_of_charges(self):
        """Position of the center of charge

        :rtype: numpy.ndarray
        """

        m = self.position_matrix(with_='charge')
        if m[:, 0].sum() == .0:
            raise RuntimeError('charge is null')

        return numpy.einsum('i,ij->j', m[:, 0], m[:, 1:]) / m[:, 0].sum()

    def translate_self_to_center_of_mass(self):
        """
        Translate the molecule to the center of mass
        """
        self.translate_self(*[-i for i in self.center_of_mass()])

    def translate_self_to_center_of_charges(self):
        """
        Translate the molecule to the center of charge
        """
        self.translate_self(*[-i for i in self.center_of_charges()])

    def moments_of_inertia(self):
        """Get the moment of inertia tensor.

        :rtype: numpy.ndarray
        """

        m = self.position_matrix(with_='mass')
        w = m[:, 0]
        m[:, 1:] -= self.center_of_mass()
        c = m[:, 1:]

        t = numpy.zeros((3, 3))
        for i in range(3):
            t[i, i] = numpy.sum(w * (c[:, (i + 1) % 3] ** 2 + c[:, (i + 2) % 3] ** 2))

        for i, j in [(0, 1), (0, 2), (1, 2)]:
            x = -numpy.sum(w * c[:, i] * c[:, j])
            t[i, j] = t[j, i] = x

        return t

    def principal_axes(self):
        """Return the moment of inertia
        (from https://github.com/moorepants/DynamicistToolKit/blob/master/dtk/inertia.py).

        :return: The principal moment of inertia (sorted from lowest to largest) and the rotation matrix
        :rtype: tuple
        """

        inertia_tensor = self.moments_of_inertia()
        e, t = numpy.linalg.eigh(inertia_tensor)
        indices = numpy.argsort(e)
        e = e[indices]
        t = t.T[indices]
        return e, t

    def set_to_inertia_axes(self):
        """Translate and rotate molecule to its principal inertia axes
        """

        self.translate_self_to_center_of_mass()
        _, rot = self.principal_axes()
        r = numpy.eye(4)
        r[:3, :3] = rot
        self._apply_transformation_self(r)

    def distances(self):
        """
        Retrieve the distance matrix between all atoms in the molecule.

        .. note::
            The array indexes are **NOT** shifted indexes.

        :rtype: numpy.ndarray
        """

        distance_matrix = numpy.zeros((len(self), len(self)))
        all_b = numpy.vstack(list(a.position for a in self))
        for i, a in enumerate(self):
            t = numpy.linalg.norm(all_b - a.position, axis=1)
            distance_matrix[i] = t

        return distance_matrix

    def connectivities(self, threshold=.1):
        """
        Get connectivity atom by atom, using the VdW radii.

        :param threshold: threshold for acceptance (because bond are sometimes a little larger than the sum of VdWs)
        :type threshold: float
        :return: the connectivity, as a dictionary per atom index (as shifted_indexes)
        :rtype: dict
        """

        distances = self.distances()
        connectivities = {}

        for i in range(len(self)):
            connectivities[i + 1] = []

        radius = {}

        def get_radius(s):
            if s not in radius:
                radius[s] = mendeleev.element(s).covalent_radius_pyykko / 100  # radius are in pm !!

            return radius[s]

        for i, a1 in enumerate(self):
            covalent_radius1 = get_radius(a1.symbol)

            for j in numpy.where(distances[i, i + 1:] < 3.0)[0]:
                r_j = i + j + 1
                a2 = self.atom_list[r_j]
                covalent_radius2 = get_radius(a2.symbol)

                if distances[i, r_j] < (covalent_radius1 + covalent_radius2 + threshold):  # yup, there is a bond
                    connectivities[i + 1].append(r_j + 1)
                    connectivities[r_j + 1].append(i + 1)

        return connectivities

    def bonds(self, threshold=.1):
        """Get a list of the bond in the molecule.

        :param threshold:  threshold (see ``self. connectivites()``)
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

    def bounding_box(self, extra_space=0.0):
        """
        Get the bounding box corresponding to the molecule.

        :param extra_space: Add some space around the box
        :type extra_space: float
        :return: Bounding box
        :rtype: qcip_tools.bounding.AABoundingBox
        """

        origin = numpy.array([1e14, 1e14, 1e14])
        maximum = numpy.array([-1e14, -1e14, -1e14])
        for atm in self:
            coo = atm.position
            for i in [0, 1, 2]:
                if coo[i] < origin[i]:
                    origin[i] = coo[i]
                if coo[i] > maximum[i]:
                    maximum[i] = coo[i]

        return bounding.AABoundingBox([i - extra_space for i in origin], maximum=[i + extra_space for i in maximum])

    def formula(self, enhanced=False):
        """Get the formula, based on self.formula_holder and self. symbols_contained.
        Implementation note : logically, it can write the atoms in whatever order. But in chemistry,
        the usual way is to start in this order C,H,N,O,X ... Where X is the other(s) atom(s). So
        does this. "enhanced" use "X_{i}" to indicate indice i instead of just "Xi".
        """

        all_atoms = list(a.symbol for a in self)

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

    def linear(self, threshold=1e-5):
        """Compute if the molecule is linear or not. A molecule is linear if :math:`I_\\alpha \\approx 0` and
        :math:`I_\\beta \\approx I_\\gamma`, if :math:`I_\\alpha < I_\\beta \\approx I_\\gamma`.

        :param threshold: threshold for the moment of inertia
        :type threshold: float
        :rtype: bool
        """

        inertia_moments, _ = numpy.linalg.eigh(self.moments_of_inertia())
        if math.fabs(inertia_moments[0]) < threshold and math.fabs(inertia_moments[1] - inertia_moments[2]) < threshold:
            return True

        return False

    def square_of_spin_angular_moment(self):
        """Return :math:`\\hat{S}_z^2=S\\,(S+1)`, where :math:`S` is the spin value

        :rtype: float
        """
        m = (self.multiplicity - 1) / 2
        return m * (m + 1)


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


class MolecularSymmetryError(Exception):
    pass


class MolecularSymmetryFinder(symmetry.SymmetryFinder):
    """High-level symmetry finder on the molecule

    :param molecule: the molecule
    :type molecule: Molecule
    :param with_: used to differentiate atoms
    :type with_: str
    :param randomise_dummy: randomise dummy atoms so that they are all different
    :type randomise_dummy: bool
    :param tol: tolerance
    :type tol: float
    """

    def __init__(self, molecule, with_='Z', randomise_dummy=True, tol=1e-5):
        self.molecule = molecule
        self._prms = (with_, randomise_dummy)

        super().__init__(self.molecule.position_matrix(with_=with_, randomise_dummy=randomise_dummy), tol=tol)

    def _group_points(self):
        """Re generate points and group them
        """

        self.points = self.molecule.position_matrix(*self._prms)
        self.group_points_per_distance = symmetry.SymmetryFinder.group_points(self.points, self.decimals)

    def orient_molecule(self):
        """Find the symmetry of the molecule, then orient it
        (just translation, then rotation).

        :rtype: qcip_tools.symmetry.PointGroupDescription
        """

        description, com, rot = self.get_symmetry()

        self.molecule.translate_self(*[-i for i in com])
        nrot = numpy.eye(4)
        nrot[:3, :3] = rot
        self.molecule._apply_transformation_self(nrot)

        self._group_points()

        return description

    def symmetrise_molecule(self, lowers_infinite=-1):
        """After ```orient_molecule()``, and among the unique positions

        - Find the set of unique atoms: apply the symmetry elements and groups atoms together ;
        - Compute an average position by using the inverse of the symmetry elements on the position of the first atom ;
        - Set components that are almost zero (wrt ``self.tol``) to zero ;
        - Generate back the other positions using the symmetry elements on the average position.

        .. warning::

            + It needs to generate all the group elements first, so that may take time ;
            + This will therefore slightly change the geometry by moving the atoms.

        :return: the point group and a list of uniques atoms indexes
        :rtype: qcip_tools.symmetry.PointGroup, list
        """

        def pos_vec_in_vecs(vec, vecs):
            p = min((numpy.sum(numpy.abs(v - vec)), i) for i, v in enumerate(vecs))
            if p[0] > self.tol:
                raise MolecularSymmetryError('no position matches !')

            return p[1]

        description = self.orient_molecule()
        g = description.gen_point_group(lowers_infinite=lowers_infinite)

        uniques_atoms = []

        for lst in self.group_points_per_distance:
            r = self.points[lst, 1:]
            saw = {}

            for i in range(len(lst)):

                if i in saw:
                    continue

                # compute the set of unique positions:
                set_of_unique = set()
                transformations_between = {}

                ri = r[i]
                for e in g:
                    npoint = e.element.apply(ri)
                    index_2 = pos_vec_in_vecs(npoint, r)
                    if index_2 not in saw:
                        saw[index_2] = True
                        set_of_unique.add(index_2)
                        transformations_between[index_2] = e

                # compute average position
                average_position = \
                    sum(g.inverse(transformations_between[j]).element.apply(r[j]) for j in set_of_unique) / \
                    len(set_of_unique)

                # set "almost zero" to zero
                for j in range(3):
                    if numpy.abs(average_position[j]) < self.tol:
                        average_position[j] = .0

                # compute back position and set atoms
                for j in set_of_unique:
                    index = lst[j]
                    self.molecule[index].position = transformations_between[j].element.apply(average_position)

                # choose an unique atom among the set
                uniques_atoms.append(lst[min((sum(numpy.abs(r[j])), j) for j in set_of_unique)[1]])

        self._group_points()
        return g, uniques_atoms


def gen_molecule_from_unique_atoms(unique_atoms, group, tol=1e-3):
    """Generate a molecule out of the uniques atoms and the symmetry group

    :param unique_atoms: set of unique atoms (please use unique ones !!)
    :type unique_atoms: list
    :param group: point group
    :type group: qcip_tools.symmetry.PointGroup
    :rtype: Molecule
    """

    final_list = []
    for a in unique_atoms:
        positions = []
        p = a.position
        for e in group:
            np = e.element.apply(p)
            if not positions or not symmetry.SymmetryFinder.vec_in_vecs(np, positions, tol):
                positions.append(np)
                final_list.append(qcip_atom.Atom(atomic_number=a.atomic_number, position=np))

    return Molecule(atom_list=final_list)
