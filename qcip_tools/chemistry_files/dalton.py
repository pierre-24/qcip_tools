import tarfile
import math
import os
import collections
import io
import numpy
import re
import copy

from qcip_tools import molecule, atom, quantities, derivatives, derivatives_e, derivatives_g
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, ChemistryLogFile, \
    FormatError, WithIdentificationMixin, PropertyNotPresent


class InputFormatError(FormatError):
    pass


class MoleculeInput(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """Dalton mol input file.

    .. warning::

        Multiplicity seems to be given in the .dal file, so it may be wrong!!

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.title``: the title (``str``)
        + ``self.basis_set``: the basis set (``str``)

    """

    #: The identifier
    file_type = 'DALTON_MOL'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.title = ''
        self.basis_set = ''

    @classmethod
    def possible_file_extensions(cls):
        return ['mol']

    @classmethod
    def attempt_identification(cls, f):
        """A dalton molecule is quite unique ("charge=")
        """

        count = 0
        found_charge = 0

        for l in f.readlines():
            count += 1

            if 'charge=' in l.lower():
                found_charge += 1

            if count > 50:
                break

        return found_charge > 1

    @classmethod
    def from_molecule(cls, molecule, title='', basis_set='', *args, **kwargs):
        """Create a file from molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :param title: title of the run
        :type title: str
        :param basis_set: the basis set
        :type basis_set: str
        :rtype: qcip_tools.chemistry_files.dalton.MoleculeInput
        """

        obj = super().from_molecule(molecule, *args, **kwargs)
        obj.title = title
        obj.basis_set = basis_set
        return obj

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        lines = f.readlines()
        self.from_read = True

        if len(lines) < 6:
            raise InputFormatError('something wrong with dalton .mol: too short')

        self.basis_set = lines[1].strip()
        self.title = (lines[2] + '\n' + lines[3]).strip()

        info_line = lines[4].lower()

        in_angstrom = 'angstrom' in info_line
        atomic_number = .0

        charge = info_line.find('charge=')
        if charge != -1:
            next_space = info_line.find(' ', charge + len('charge='))
            self.molecule.charge = int(info_line[charge + len('charge='):next_space])
            self.molecule.multiplicity = 1 if self.molecule.charge % 2 == 0 else 1

        for line in lines[5:]:
            if 'charge=' in line.lower():
                atomic_number = int(float(line[7:line.lower().find('atoms')]))
                continue

            content = line.split()
            if len(content) != 4:
                continue

            self.molecule.insert(atom.Atom(
                atomic_number=atomic_number,
                position=[float(a) * (1 if in_angstrom else quantities.AuToAngstrom) for a in content[1:]])
            )

    def to_string(self, in_angstrom=True, nosym=False, group_atoms=False):
        """

        :param in_angstrom: gives the atomic coordinates in angstrom
        :type: in_angstrom: bool
        :param nosym: specify "nosymmetry"
        :type nosym: bool
        :param group_atoms: group of the same type together (the order may be lost)
        :type group_atoms: bool
        :rtype: str
        """
        r = 'BASIS\n{}\n'.format(self.basis_set)
        r += self.title

        if self.title.find('\n') == -1:
            r += '\n\n'

        r += 'Atomtypes={}{}{}{}\n'.format(
            len(self.molecule.symbols_contained) if group_atoms else len(self.molecule),
            ' Angstrom' if in_angstrom else '',
            ' Nosymmetry' if nosym else '',
            ' Charge={}'.format(self.molecule.charge) if self.molecule.charge != 0 else ''
        )

        if not group_atoms:
            for a in self.molecule:
                r += 'Charge={:.1f} Atoms={}\n'.format(a.atomic_number, 1)
                r += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(
                    a.symbol, *[p * (1 if in_angstrom else quantities.AuToAngstrom) for p in a.position])
        else:
            for symbol in self.molecule.symbols_contained:
                atms = self.molecule.atoms(symbol_in=[symbol])
                r += 'Charge={:.1f} Atoms={}\n'.format(atom.Definition[symbol][0], len(atms))

                for i in atms:
                    a = self.molecule[i]
                    r += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(
                        a.symbol, *[p * (1 if in_angstrom else quantities.AuToAngstrom) for p in a.position])

        return r

    def write(self, f, in_angstrom=True, nosym=False, group_atoms=False):
        """

        :param f: File
        :type f: file
        :param in_angstrom: gives the atomic coordinates in angstrom
        :type: in_angstrom: bool
        :param nosym: specify "nosymmetry"
        :type nosym: bool
        :param group_atoms: group of the same type together (the order may be lost)
        :type group_atoms: bool
        """

        f.write(self.to_string(in_angstrom, nosym, group_atoms))


class ArchiveOutput(ChemistryFile, WithMoleculeMixin, WithIdentificationMixin):
    """Archive output of Dalton. Contains lots of information for who can extract them.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.tar_file``: pipe to the tar file (``tarfile.TarFile``)

    """

    #: The identifier
    file_type = 'DALTON_ARCHIVE'
    requires_binary_mode = True

    def __init__(self):
        self.tar_file = None
        self.molecule = molecule.Molecule()

    @classmethod
    def possible_file_extensions(cls):
        return ['gz', 'tar.gz']

    @classmethod
    def attempt_identification(cls, f):
        """A dalton archive does contain a lot of files name "DALTON"
        """

        if 'b' not in f.mode:
            raise IOError('file {} must be open in binary mode!'.format(f.name))

        try:
            t = tarfile.open(fileobj=f)
        except tarfile.TarError:
            return False

        found_dalton = 0

        for n in t.getnames():
            if 'DALTON' in n:
                found_dalton += 1

        return found_dalton > 3

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        raise NotImplementedError('from_molecule')

    def read(self, f):
        """Expects ``f`` to be open on binary mode.

        .. note::

            In Dalton archive, the molecule is given in a old, weird format, called "molplt". Information about it is
            given in `this webpage <https://brettbode.github.io/wxmacmolplt/Manual_pages/MacMolPlt_Files.html>`_.

        :param f: File
        :type f: file
        """

        if 'b' not in f.mode:
            raise IOError('file {} must be open in binary mode!'.format(f.name))

        self.from_read = True
        self.tar_file = tarfile.open(f.name)

        # read molecule:
        mol_file = self.get_file('DALTON.MOL')

        lines = mol_file.readlines()

        info = lines[1].split()
        number_of_atoms = int(info[1])
        number_of_kinds = int(info[3])

        for i in range(number_of_atoms):
            a = lines[4 + number_of_kinds + i].split()
            self.molecule.insert(atom.Atom(symbol=str(a[0], 'utf-8'), position=[float(e) for e in a[1:]]))

        # correct masses if any
        fx = None
        try:
            fx = self.get_file('DALTON.NCA')
        except FileNotFoundError:
            pass

        if fx is not None:
            lines = fx.read().decode().splitlines()
            mode_line = 0
            for i, j in enumerate(lines[2:]):
                if 'MODE' in j:
                    mode_line = 2 + i
                    break
            masses = ' '.join(lines[2:mode_line]).split()
            if len(masses) != len(self.molecule):
                raise Exception('number of masses does not matches number of atoms')

            for i, mass in enumerate(masses):
                self.molecule[i].mass = float(mass)

    def get_file(self, name):
        """Get a ``io.BufferedReader`` from a filename, if exists. Note that you get ``bytes`` out of that kind of
        objects, so you need to ``decode()`` them.

        :param name: file name
        :type name: str
        :rtype: io.BufferedReader
        """

        if self.tar_file is None:
            raise IOError('no archive open')

        if name in self.tar_file.getnames():
            return self.tar_file.extractfile(self.tar_file.getmember(name))
        else:
            raise FileNotFoundError(name)


class WrongNumberOfData(ValueError):
    def __init__(self, n, p):
        super().__init__('wrong number of data ({}) for {}'.format(n, p))


@ArchiveOutput.define_property('computed_energies')
def dalton__archive_output__property__computed_energies(obj, *args, **kwargs):
    """Get the energies. Returns only the energy found on top of ``DALTON.PROP`` as the total energy (+labeled)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.ArchiveOutput
    :rtype: dict
    """

    try:
        f = obj.get_file('DALTON.PROP')
    except FileNotFoundError:
        raise PropertyNotPresent('computed_energies')

    energy_line = f.readline().decode()
    if 'ENERGY' not in energy_line:
        raise PropertyNotPresent('computed_energies')

    label = energy_line[13:23].strip()
    value = float(energy_line[23:46])

    return {'total': value, label: value}


@ArchiveOutput.define_property('electrical_derivatives')
def dalton__archive_output__property__electrical_derivatives(obj, *args, **kwargs):
    """Get electrical derivatives in Dalton archives. Returns a dictionary of dictionaries:

    .. code-block:: text

        + "F"
            + static : ElectricDipole
        + "FF":
            + static : PolarizabilityTensor
        + "FD":
            + 0.042 : PolarizabilityTensor
            + ...
        + "FDD":
            + static : FirstHyperpolarizabilityTensor
            + ...
        + ...

    .. note ::

        CC response calculations in Dalton are nice, since it allows form many order of response at the same time
        (thought it takes some memory) and it gives all the responses in one single file: DALTON.PROP.

        On the other hand, the ``**RESPONSE`` module of dalton does not allows for many order at the same time, and
        is f\*\*\*g inconsistent:

        + dipole moment is in DALTON.PROP
        + alpha is in DALTON.PROP, but only the non-zero components (and you need to assume permutation)
        + beta is in RESULTS.RSP, because why not. You need to assume Kleinman condition for SHG, and
          optical rectification is assumed to be the same as EOP, just permuted (and I'm not entirely sure
          of the validity of what I did there).
        + If you use the *patched* version of Dalton (see installation part), results are then available in
          DALTON.PROP (but you still need to assume permutation).
        + gamma is in DALTON.PROP (main source, but only non-zero components that contribute to :math:`\gamma_{||}`)
          and also in RESULTS.RSP.
          You need to use an hidden option, ``.GAMALL`` to get all components of the gamma tensor
          (otherwise **only the one contributing to :math:`\gamma_{||}` are computed**).
          The rest of the components are obtained from RESULT.RSP if ``.GAMALL`` is used.
          Note that the first frequency is non-null for dc-Kerr, while it is the last one with CC response module.
        + If you use the patched version, everything is in DALTON.PROP.

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.ArchiveOutput
    :rtype: dict
    """

    electrical_derivatives = {}
    translate_diplens = {'X': 0, 'Y': 1, 'Z': 2}

    # Check DALTON.PROP
    f = None

    try:
        f = obj.get_file('DALTON.PROP')
    except FileNotFoundError:
        pass

    beta_already_found = False
    missing_gamma_components = False

    if f is not None:
        lines = f.read().decode('utf-8').splitlines()
        splits = [0, 5, 8, 13, 23, 46, 55, 64, 73, 82, 105, 128, 151]
        CC_methods = ['CCS', 'CC2', 'CCSD', 'CC3']
        is_CC_calculation = False
        check_CC_calculation = False

        prop_sorted = {}
        for line in lines:
            info = []
            for index, split in enumerate(splits[1:]):
                info.append(line[splits[index]:split].strip())

            if not check_CC_calculation and info[5] != 'ENERGY':
                is_CC_calculation = info[3] in CC_methods

            derivative = int(info[2])
            if derivative not in prop_sorted:
                prop_sorted[derivative] = []
            prop_sorted[derivative].append(info)

        # dipole moment
        if 1 in prop_sorted:
            if len(prop_sorted[1]) != 3:
                raise WrongNumberOfData(len(prop_sorted[1]), 'dipole moment')

            dipole = numpy.zeros(3)

            for l in prop_sorted[1]:
                dipole[translate_diplens[l[5][0]]] = float(l[4])

            electrical_derivatives['F'] = {'static': derivatives_e.ElectricDipole(dipole=dipole)}

        # polarizability
        if 2 in prop_sorted:
            if is_CC_calculation and len(prop_sorted[2]) % 9 != 0:
                raise WrongNumberOfData(len(prop_sorted[2]), 'polarizability')

            num_of_tensors = int(len(prop_sorted[2]) / 9)
            found_frequencies = []
            data_per_frequencies = {}

            for l in prop_sorted[2]:
                freq = float(l[9])
                if freq not in found_frequencies:
                    found_frequencies.append(freq)
                    if is_CC_calculation and len(found_frequencies) > num_of_tensors:
                        raise WrongNumberOfData(num_of_tensors + 1, 'possible tensor for alpha')
                    data_per_frequencies[freq] = derivatives_e.PolarisabilityTensor(frequency=freq)

                c1, c2, v = translate_diplens[l[5][0]], translate_diplens[l[6][0]], float(l[4])
                data_per_frequencies[freq].components[c1, c2] = v

                if not is_CC_calculation and c1 != c2:
                    data_per_frequencies[freq].components[c2, c1] = v

            electrical_derivatives['dD'] = data_per_frequencies
            if .0 in electrical_derivatives['dD']:
                electrical_derivatives['FF'] = {'static': electrical_derivatives['dD'][.0]}
                electrical_derivatives['FF']['static'].input_fields = (0,)
                electrical_derivatives['FF']['static'].representation = derivatives.Derivative(from_representation='FF')
                del electrical_derivatives['dD'][.0]

        # hyperpolarizability
        if 3 in prop_sorted:
            if len(prop_sorted[3]) % 27 != 0:
                raise WrongNumberOfData(len(prop_sorted[3]), 'hyperpolarizability')

            beta_already_found = True
            num_of_tensors = int(len(prop_sorted[3]) / 27)
            found_frequencies = {'FFF': [], 'dDF': [], 'FDd': [], 'XDD': []}
            data = {}
            r_to_obj = {}

            for x in found_frequencies.keys():
                r_to_obj[x] = derivatives.Derivative(from_representation=x)

            for l in prop_sorted[3]:
                if len(l[9]) == 22:
                    freq_1 = float(l[9])
                    freq_2 = float(l[10])
                else:
                    freq_1 = float(l[9][:22])
                    freq_2 = float(l[9][22:45])

                if freq_1 == freq_2 and freq_1 == .0:
                    representation = 'FFF'
                    freq_1 = 'static'
                elif freq_1 == freq_2:
                    representation = 'XDD'
                elif freq_1 == -freq_2:
                    representation = 'FDd'
                elif freq_2 == .0:
                    representation = 'dDF'
                else:
                    raise Exception('unknown combination of field ({}, {})'.format(freq_1, freq_2))

                if freq_1 not in found_frequencies[representation]:
                    if representation not in data:
                        data[representation] = {}

                    found_frequencies[representation].append(freq_1)

                    if is_CC_calculation and sum(len(x) for x in found_frequencies.values()) > num_of_tensors:
                        raise WrongNumberOfData(num_of_tensors + 1, 'possible tensor for beta')

                    data[representation][freq_1] = derivatives_e.FirstHyperpolarisabilityTensor(
                        frequency=freq_1, input_fields=tuple(
                            derivatives_e.representation_to_field[x] for x in representation[1:]))

                c1, c2, c3 = translate_diplens[l[5][0]], translate_diplens[l[6][0]], translate_diplens[l[7][0]]
                v = float(l[4])
                data[representation][freq_1].components[c1, c2, c3] = v

                if not is_CC_calculation:
                    for k in r_to_obj[representation].inverse_smart_iterator((c1, c2, c3)):
                        data[representation][freq_1].components[k] = v

            electrical_derivatives.update(data)

        # second hyperpolarizability:
        if 4 in prop_sorted:
            if is_CC_calculation and len(prop_sorted[4]) % 81 != 0:
                raise WrongNumberOfData(len(prop_sorted[4]), 'hyperpolarizability')

            num_of_tensors = int(len(prop_sorted[4]) / 81)
            found_frequencies = {'FFFF': [], 'dFFD': [], 'dDFF': [], 'XDDF': [], 'dDDd': [], 'XDDD': []}
            data = {}

            for l in prop_sorted[4]:
                freq_1 = float(l[9])

                if len(l[10]) == 22:
                    freq_2 = float(l[10])
                    freq_3 = float(l[11])
                else:
                    freq_2 = float(l[10][:22])
                    freq_3 = float(l[10][22:])

                if freq_1 == freq_2 and freq_2 == freq_3 and freq_3 == .0:
                    representation = 'FFFF'
                    freq = 'static'
                elif freq_1 == freq_2 and freq_2 == .0:
                    representation = 'dFFD'
                    freq = freq_3
                elif freq_2 == freq_3 and freq_3 == .0:
                    representation = 'dDFF'
                    freq = freq_1
                elif freq_1 == freq_2 and freq_3 == .0:
                    representation = 'XDDF'
                    freq = freq_1
                elif freq_1 == freq_2 and freq_2 == -freq_3:
                    representation = 'dDDd'
                    freq = freq_1
                elif freq_1 == freq_2 and freq_2 == freq_3:
                    representation = 'XDDD'
                    freq = freq_1
                else:
                    raise Exception('unknown combination of field ({}, {}, {})'.format(freq_1, freq_2, freq_3))

                if freq not in found_frequencies[representation]:
                    if representation not in data:
                        data[representation] = {}

                    found_frequencies[representation].append(freq)

                    if sum(len(x) for x in found_frequencies.values()) > num_of_tensors:
                        if is_CC_calculation:
                            raise WrongNumberOfData(num_of_tensors + 1, 'possible tensor for gamma')
                        else:
                            missing_gamma_components = True

                    data[representation][freq] = derivatives_e.SecondHyperpolarizabilityTensor(
                        frequency=freq, input_fields=tuple(
                            derivatives_e.representation_to_field[x] for x in representation[1:]))

                data[representation][freq] \
                    .components[
                        translate_diplens[l[5][0]],
                        translate_diplens[l[6][0]],
                        translate_diplens[l[7][0]],
                        translate_diplens[l[8][0]]] \
                    = float(l[4])

            electrical_derivatives.update(data)

    # CHECK RESULTS.RSP (beta response at the HF/DFT level)
    f = None

    try:
        f = obj.get_file('RESULTS.RSP')
    except FileNotFoundError:
        pass

    if f is not None:
        lines = f.read().decode('utf-8').splitlines()
        data = {}

        if not beta_already_found:
            found_frequencies = {'FFF': [], 'dDF': [], 'FDd': [], 'XDD': []}
            splits = [0, 10, 20, 29, 40, 48, 53, 59, 73]
            is_optical_rectification_in_calculations = False
            is_EOP_in_calculations = False
            is_optical_rectification_use_later = []

            r_to_obj = {}  # to use inverse smart iterator
            for x in found_frequencies.keys():
                r_to_obj[x] = derivatives.Derivative(from_representation=x)

            for line in lines:
                if len(line) == 0:
                    break
                if 'Cubic response function' in line:
                    break

                info = []
                for index, split in enumerate(splits[1:]):
                    info.append(line[splits[index]:split].strip())

                freq_1 = float(info[1])
                freq_2 = float(info[3])
                components = tuple(translate_diplens[x] for x in re.split('\W', info[5]))
                value = float(info[7])

                if freq_1 == freq_2 and freq_1 == .0:
                    representation = 'FFF'
                    freq_1 = 'static'
                elif freq_1 == freq_2:
                    representation = 'XDD'
                elif freq_1 == -freq_2:  # will be treated later
                    representation = 'FDd'
                    is_optical_rectification_in_calculations = True
                    if is_EOP_in_calculations:
                        is_optical_rectification_use_later.append((freq_1, components, value))
                        continue
                elif freq_2 == .0:
                    representation = 'dDF'
                    is_EOP_in_calculations = True
                else:
                    raise Exception('unknown combination of field ({}, {})'.format(freq_1, freq_2))

                if freq_1 not in found_frequencies[representation]:
                    if representation not in data:
                        data[representation] = {}

                    found_frequencies[representation].append(freq_1)

                    data[representation][freq_1] = derivatives_e.FirstHyperpolarisabilityTensor(
                        frequency=freq_1, input_fields=tuple(
                            derivatives_e.representation_to_field[x] for x in representation[1:]))

                for components_ in r_to_obj[representation].inverse_smart_iterator(components):
                    data[representation][freq_1].components[components_] = value

            if is_optical_rectification_in_calculations and is_EOP_in_calculations:
                freqs = list(data['dDF'].keys())
                data['FDd'] = {}
                for freq in freqs:
                    data['FDd'][freq] = derivatives_e.FirstHyperpolarisabilityTensor(
                        frequency=freq, input_fields=(-1, 1))
                    for i in r_to_obj['dDF'].smart_iterator():
                        if data['dDF'][freq].components[i] == .0:
                            continue
                        for j in r_to_obj['dDF'].inverse_smart_iterator(i):
                            data['FDd'][freq].components[j[2], j[1], j[0]] = data['dDF'][freq].components[j]

                for freq, components_, value in is_optical_rectification_use_later:
                    data['FDd'][freq].components[components_] = value

        if missing_gamma_components:
            i = 0
            found_frequencies = {'FFFF': [], 'dDFF': [], 'XDDF': [], 'dDDd': [], 'XDDD': []}

            r_to_obj = {}
            for x in found_frequencies.keys():
                r_to_obj[x] = derivatives.Derivative(from_representation=x)

            while i < len(lines):
                if 'Cubic response function value in a.u. for' in lines[i]:
                    v = -float(lines[i + 6][25:])

                    freq_1, freq_2, freq_3 = tuple(float(lines[i + k][49:]) for k in range(2, 5))
                    if freq_1 == freq_2 and freq_2 == freq_3 and freq_3 == .0:
                        representation = 'FFFF'
                    elif freq_2 == freq_3 and freq_3 == .0:
                        representation = 'dDFF'
                    elif freq_1 == freq_2 and freq_3 == .0:
                        representation = 'XDDF'
                    elif freq_1 == freq_2 and freq_2 == -freq_3:
                        representation = 'dDDd'
                    elif freq_1 == freq_2 and freq_2 == freq_3:
                        representation = 'XDDD'
                    else:
                        raise Exception('unknown combination of field ({}, {}, {})'.format(freq_1, freq_2, freq_3))

                    components = tuple(translate_diplens[lines[i + k][36:37]] for k in range(1, 5))
                    i += 8

                    if not any(i == 1 for i in collections.Counter(components).values()):
                        continue

                    freq = 'static'
                    if representation != 'FFFF':
                        freq = next(f for f in electrical_derivatives[representation] if (f - freq_1) < 1e-5)

                    for components_ in r_to_obj[representation].inverse_smart_iterator(components):
                        electrical_derivatives[representation][freq].components[components_] = v
                else:
                    i += 1

        electrical_derivatives.update(data)

    if not electrical_derivatives:
        raise PropertyNotPresent('electrical_derivatives')

    return electrical_derivatives


@ArchiveOutput.define_property('geometrical_derivatives')
def dalton__archive_output__property__geometrical_derivatives(obj, *args, **kwargs):
    """Get geometrical derivatives in Dalton archives. Returns a dictionary of tensor

    .. note ::

        Only return geometrical hessian

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.ArchiveOutput
    :rtype: dict
    """

    geometrical_derivatives = {}

    # CHECK RESULTS.RSP (response at the HF/DFT level)
    f = None

    try:
        f = obj.get_file('DALTON.HES')
    except FileNotFoundError:
        pass

    if f is not None:
        lines = f.read().decode('utf-8').splitlines()

        n = int(lines[0].strip())
        if n != 3 * len(obj.molecule):
            raise Exception('hessian and molecule does not match')

        geometrical_derivatives['GG'] = derivatives_g.BaseGeometricalDerivativeTensor(
            representation='GG', spacial_dof=n, trans_plus_rot=5 if obj.molecule.linear() else 6)

        for i in range(n):
            for j in range(n):
                number = i * (n + 1) + j + 2
                geometrical_derivatives['GG'].components[i, j] = float(lines[number])

    if not geometrical_derivatives:
        raise PropertyNotPresent('geometrical_derivatives')

    return geometrical_derivatives


class OutputSection:
    """Remind when a section is entered and exited

    :param section: name of the section
    :type section: str
    :type line_start: int
    :type line_end: int
    """

    def __init__(self, section, line_start, line_end):

        self.section = section
        self.line_start = line_start
        self.line_end = line_end

    def __repr__(self):
        return 'Section {}: {}:{}'.format(self.section, self.line_start, self.line_end)


class OutputFormatError(FormatError):
    pass


class Output(ChemistryLogFile, WithMoleculeMixin, WithIdentificationMixin):
    """Output of Dalton.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.lines``: the lines of the file (``list`` of ``str``)
        + ``self.sections``: the different sections (``list`` of
          `OutputSection <#qcip_tools.chemistry_files.dalton.OutputSection>`_)

    """

    #: The identifier
    file_type = 'DALTON_LOG'
    chunk_title_variable = 'section'

    def __init__(self):
        self.molecule = molecule.Molecule()

    @classmethod
    def attempt_identification(cls, f):
        """A dalton output ... Does contains a lot of (european) countries and universities name at some point
        (and a little bit of "Dalton" as well)
        """

        count = 0
        found_countries = 0
        found_universities = 0
        found_dalton = 0
        countries = ['Norway', 'Denmark', 'Italy', 'Sweden', 'Germany']

        for l in f.readlines():
            count += 1

            if 'Dalton' in l or 'dalton' in l:
                found_dalton += 1

            if 60 < count < 150:
                if 'University' in l:
                    found_universities += 1
                for c in countries:
                    if c in l:
                        found_countries += 1
                        break
            if count > 150:
                break

        return found_countries > 20 and found_universities > 20 and found_dalton > 4

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        super().read(f)

        self.chunks.append(OutputSection('START', 0, -1))

        for index, line in enumerate(self.lines):
            if '.----------------' in line:
                if 'Starting in' in self.lines[index + 1]:
                    self.chunks[-1].line_end = index
                    section_end = self.lines[index + 1].rfind(')')
                    section_start = self.lines[index + 1].rfind('(')
                    if -1 in [section_end, section_start]:
                        raise Exception('no section found for line {}'.format(index))

                    self.chunks.append(
                        OutputSection(self.lines[index + 1][section_start + 1:section_end], index + 1, -1))

        # fetch molecule
        coordinates_line = self.search('Cartesian Coordinates (a.u.)', into='START')
        if coordinates_line != -1:
            number_of_atoms = int(math.floor(int(self.lines[coordinates_line + 3][-5:]) / 3))
            for i in range(number_of_atoms):
                content = self.lines[coordinates_line + 4 + i].split()
                self.molecule.insert(
                    atom.Atom(
                        symbol=content[0],
                        position=[a * quantities.AuToAngstrom for a in [
                            float(content[4]), float(content[7]), float(content[10])]]
                    ))
        else:  # if it is C1 molecule, it may not be reprinted
            coordinates_line = self.search('Content of the .mol file', into='START')
            if coordinates_line == -1:
                raise OutputFormatError('Unable to find molecule, not even in beginning!')
            line_end = self.search(
                'Output from DALTON general input processing', line_start=coordinates_line, into='START')

            mol_file = MoleculeInput()
            mol_file.read(io.StringIO(''.join(self.lines[coordinates_line + 3:line_end - 3])))

            self.molecule = copy.deepcopy(mol_file.molecule)

        mass_line = self.search('Isotopic Masses', into='START')
        for i in range(len(self.molecule)):
            line = self.lines[mass_line + i + 3].split()
            self.molecule[i].mass = float(line[-1])

        charge_line = self.search('@    Total charge of the molecule', into='SIRIUS')
        if charge_line == -1:
            raise OutputFormatError('cannot find charge of the molecule in SIRIUS')
        self.molecule.charge = int(self.lines[charge_line][-6:])
        self.molecule.multiplicity = int(self.lines[charge_line + 2][45:55])

        if self.molecule.number_of_electrons() != int(self.lines[charge_line - 2][-5:]):
            raise OutputFormatError(
                'Not the right number of electron, is the charge (found {}) wrong?'.format(self.molecule.charge))

    def get_inputs(self):
        """Fetch the inputs from the log file

        :rtype: tuple
        """

        line_dal = self.search('Content of the .dal input file', into='START')
        line_mol = self.search('Content of the .mol file', line_start=line_dal, into='START')
        line_end = self.search('Output from DALTON general input processing', line_start=line_mol, into='START')

        if -1 in [line_dal, line_mol, line_end]:
            raise OutputFormatError('inputs not found in START')

        dal_file = Input()
        mol_file = MoleculeInput()

        dal_file.read(io.StringIO(''.join(self.lines[line_dal + 2:line_mol - 2])))
        mol_file.read(io.StringIO(''.join(self.lines[line_mol + 3:line_end - 3])))

        return dal_file, mol_file


@Output.define_property('computed_energies')
def dalton__output__computed_energies(obj, *args, **kwargs):
    """Get the energies. Find the SCF/DFT energy, as well as the CCx one (which is given with MP2).

    .. note::

        + Not (yet?) CI, CASCI or MCSCF
        + Probably wrong on geometry optimization.

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.Output
    :rtype: dict
    """

    line_scf = obj.search('@    Final HF energy:', into='SIRIUS')
    if line_scf < 0:
        raise PropertyNotPresent('computed_energies')

    value = float(obj.lines[line_scf][-40:])
    energies = {'SCF/DFT': value, 'total': value}

    if obj.chunk_exists('CC'):
        line_ens = obj.search('Final results from the Coupled Cluster energy program', into='CC')

        if line_ens < 0:
            return energies

        energies['MP2'] = float(obj.lines[line_ens + 7][-25:])
        label = obj.lines[line_ens + 9][18:24].strip()
        energies[label] = float(obj.lines[line_ens + 9][-25:])
        energies['total'] = label

    return energies


@Output.define_property('geometrical_derivatives')
def dalton__output__get_geometrical_derivatives(obj, *args, **kwargs):
    """Get the cartesian gradient out of a dalton calculation output

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.Output
    :rtype: dict
    """

    geometrical_derivatives = {}

    dof = 3 * len(obj.molecule)
    trans_plus_rot = 5 if obj.molecule.linear() else 6

    # CC
    x = obj.search('Molecular gradient', into='CC')
    if x != -1:
        gradient = derivatives_g.BaseGeometricalDerivativeTensor(
            representation='G', spacial_dof=dof, trans_plus_rot=trans_plus_rot)

        for i in range(len(obj.molecule)):
            line = obj.lines[x + 3 + i].split()
            gradient.components[i * 3:(i + 1) * 3] = [float(a) for a in line[1:]]
        geometrical_derivatives['G'] = gradient

    # general !
    if obj.chunk_exists('ABACUS'):
        x = obj.search('Molecular gradient', into='ABACUS')
        if x != -1:
            gradient = derivatives_g.BaseGeometricalDerivativeTensor(
                representation='G', spacial_dof=dof, trans_plus_rot=trans_plus_rot)

            for i in range(len(obj.molecule)):
                line = obj.lines[x + 3 + i].split()
                gradient.components[i * 3:(i + 1) * 3] = [float(a) for a in line[1:]]
            geometrical_derivatives['G'] = gradient

        x = obj.search('Molecular Hessian (au)', into='ABACUS', line_start=x)
        if x != -1:
            hessian = derivatives_g.BaseGeometricalDerivativeTensor(
                representation='GG', spacial_dof=dof, trans_plus_rot=trans_plus_rot)

            spliting = [10, 22, 34, 46, 58, 70, 81]

            for i in range(dof):
                for j in range(dof):
                    if j > i:
                        continue

                    r = math.floor(j / 6)
                    s = math.floor(i / 3)
                    line = (dof + len(obj.molecule) + 7 - 4 * r) * r + i + s - 8 * r
                    hessian.components[i, j] = float(obj.lines[x + line + 5][spliting[j % 6]:spliting[j % 6 + 1]])
                    if i != j:
                        hessian.components[j, i] = hessian.components[i, j]

            geometrical_derivatives['GG'] = hessian

    if not geometrical_derivatives:
        raise PropertyNotPresent('geometrical_derivatives')

    return geometrical_derivatives


#: Name of the allowed modules in dalton (according to documentation)
ALLOWED_LEVEL_0_MODULES = ['**DALTON', '**INTEGRAL', '**WAVE F', '**START', '**EACH S', '**PROPERTIES', '**RESPONSE']


def check_module(name):
    """Check if the module name is allowed

    :param name: name of the module
    :type name: str
    :rtype: bool
    """

    n = '**' + name[:5]
    allowed = False
    for module in ALLOWED_LEVEL_0_MODULES:
        if module[:7] == n:
            allowed = True
            break

    return allowed


class InputModuleError(Exception):
    pass


class InputModule:
    """The dalton (sub)modules, with level indicating it they are module (level=0) or submodule (level=1)

    .. note::

        Since Dalton only interpret the 7 first characters of any input card or module (``*`` or ``.`` included),
        the storage key is reduced to that.

    :param level: level of the module (either 0 if principal or 1 if submodule)
    :type level: int
    :param name: name of the module
    :type name: str
    :param input_cards: The eventual input cards
    :type input_cards: collections.OrderedDict
    :param submodules: The eventual sub modules
    :type submodules: collections.OrderedDict
    """

    def __init__(self, level=1, name=None, input_cards=None, submodules=None):
        if level not in [0, 1]:
            raise InputModuleError('level should be zero or one')

        if level == 0 and name is not None:
            if not check_module(name):
                raise InputModuleError('not an allowed module: {}'.format(name))

        if level == 1 and submodules is not None:
            raise InputModuleError('subdmodule {} with subsubmodules'.format(name))

        self.name = name
        self.level = level
        self.input_cards = collections.OrderedDict() if input_cards is None else input_cards
        self.submodules = collections.OrderedDict() if submodules is None else submodules

    def set_submodule(self, submodule, force=False):
        """Set up a submodule

        :param submodule: module
        :type submodule: qcip_tools.chemistry_files.dalton.InputModule
        :param force: force replacing  the submodule if exists
        :type force: bool
        """

        if self.level != 0:
            raise InputModuleError('no subsubmodule allowed in {}!'.format(self.name))

        if not force and submodule.name[:6] in self.submodules:
            raise InputModuleError('submodule {} already exists in {}'.format(submodule.name, self.name))

        self.submodules[submodule.name[:6]] = submodule
        submodule.level = 1

    def set_input_card(self, input_card, force=True):
        """Set an input card

        :param input_card: module
        :type input_card: InputCard
        :param force: force replacing  the card if exists
        :type force: bool
        """

        name = input_card.name
        if name[0] == '.':
            name = name[1:]

        if not force and name[:6] in self.input_cards:
            raise InputModuleError('input card {} already exists in {}'.format(input_card.name, self.name))

        self.input_cards[name[:6]] = input_card

    def __contains__(self, item):
        if item[0] == '.' and item[1:7] in self.input_cards:
            return True

        if self.level == 0 and item[:6] in self.submodules:
            return True

        return False

    def __getitem__(self, item):
        if item[0] == '.' and item[1:7] in self.input_cards:
            return self.input_cards[item[1:7]]

        if self.level == 0 and item[:6] in self.submodules:
            return self.submodules[item[:6]]

        raise KeyError(item)

    def __setitem__(self, key, value):
        if type(value) not in [InputModule, InputCard]:
            raise TypeError(value)

        if type(value) is InputModule:
            if key[0] == '.':
                raise InputModuleError('module should not start with "."')

            if value.name is None:
                value.name = key

            if key[:6] != value.name[:6]:
                raise InputModuleError('key ({}) and name divergence ({})'.format(key[:6], value.name[:6]))

            self.set_submodule(value)

        else:
            if key[0] != '.':
                raise InputModuleError('input card should start with "."')

            if value.name is None:
                value.name = key[1:]

            name = value.name
            if name[0] == '.':
                name = name[1:]

            if key[1:7] != name[:6]:
                raise InputModuleError('key ({}) and name divergence ({})'.format(key[1:7], name[:6]))

            self.set_input_card(value)

    def __repr__(self):
        r = '*' + ('*' if self.level == 0 else '') + self.name + '\n'
        for i in self.input_cards.values():
            r += str(i)

        if self.level == 0:
            for i in self.submodules.values():
                r += str(i)

        return r

    def update(self, r):
        """Update a module with respect to another

        :param r: source module
        :type r: InputModule
        """

        # update submodule
        if self.level == 0:
            for k in r.submodules:
                if k in self.submodules:
                    self.submodules[k].update(r.submodules[k])
                else:
                    self.submodules[k] = copy.deepcopy(r.submodules[k])
        # update input card
        for k in r.input_cards:
            if k in self.input_cards:
                self.input_cards[k].parameters = r.input_cards[k].parameters.copy()
            else:
                self.input_cards[k] = r.input_cards[k]


class InputCard:
    """Dalton input card

    :param name: name of the module
    :type name: str
    :param parameters: the parameters
    :type parameters: list
    """

    def __init__(self, name=None, parameters=None):

        self.name = name

        if name is not None and self.name[0] == '.':
            self.name = name[1:]

        self.parameters = [] if parameters is None else parameters

    def __repr__(self):
        r = '.' + self.name + '\n'

        if self.parameters:
            for i in self.parameters:
                r += str(i) + '\n'

        return r


class Input(ChemistryFile, WithOutputMixin, WithIdentificationMixin):
    """Dalton dal input file.

    Do NOT contains a molecule!

    .. note::

        Since Dalton only interpret the 7 first characters of any module (``**`` included),
        the storage key is reduced to 5 characters.

    .. container:: class-members

        + ``self.module``: modules (``collections.OrderedDict`` of
          `InputModule <#qcip_tools.chemistry_files.dalton.InputModule>`_)

    """

    #: The identifier
    file_type = 'DALTON_DAL'

    def __init__(self):
        self.modules = collections.OrderedDict()

    @classmethod
    def possible_file_extensions(cls):
        return ['dal']

    @classmethod
    def attempt_identification(cls, f):
        """Very specific beginning (``**DALTON INPUT``) and end (``'*END OF INPUT``)
        """

        l0 = f.readline()
        if '**DALTON' == l0[:8]:

            with open(f.name, 'rb') as fh:
                # trick to get easily the last line, see https://stackoverflow.com/a/3346492
                # TODO: is it really necessary for a file this size ?
                fh.seek(-30, os.SEEK_END)
                last = fh.readlines()[-1].decode()

                if '*END OF' == last[:7]:
                    return True

        return False

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        lines = f.readlines()

        current_module = None
        current_submodule = None

        for index, line in enumerate(lines):  # comments
            if line[0] in ['!', '#']:
                continue
            if line[0] not in ['*', '.']:  # parameters will be read directly
                continue

            if line[0] == '*':
                if line[1] == '*':  # module
                    current_module = line[2:].strip()
                    self[current_module] = InputModule(level=0)
                    current_submodule = None

                else:  # submodule
                    if '*END OF' in line[:7]:
                        break

                    current_submodule = line[1:].strip()
                    self[current_module][current_submodule] = InputModule(level=1)

            elif line[0] == '.':
                current_input_card = line[1:].strip()
                parameters = []

                for line_param in lines[index + 1:]:
                    if line_param[0] in ['.', '*']:
                        break
                    elif line_param[0] in ['!', '#']:
                        continue
                    else:
                        parameters.append(line_param.strip())

                if current_submodule is None:
                    self[current_module]['.' + current_input_card] = InputCard(parameters=parameters)
                else:
                    self[current_module][current_submodule]['.' + current_input_card] = InputCard(parameters=parameters)

    def __contains__(self, item):

        if item[:5] in self.modules:
            return True

        return False

    def __getitem__(self, item):
        if item in self:
            return self.modules[item[:5]]

        raise KeyError(item)

    def __setitem__(self, key, value):
        if type(value) != InputModule:
            raise TypeError(value)

        if value.name is not None and key[:5] != value.name[:5]:
            raise InputModuleError('key ({}) and name ({}) divergence'.format(key[:5], value.name[:5]))

        elif value.name is None:

            if not check_module(key):
                raise InputModuleError('not an allowed module: {}'.format(key))

            value.name = key

        value.level = 0
        self.modules[key[:5]] = value

    def update(self, r):
        """Update the content from a given piece of string

        :param r: string
        :type r: str|Input
        """

        if type(r) is str:
            d = Input()
            try:
                d.read(io.StringIO(r))
            except (InputModuleError, InputFormatError) as e:
                raise e
        elif type(r) is Input:
            d = r
        else:
            raise TypeError('r should be either string or Input')

        for k in d.modules:
            if k in self.modules:
                self.modules[k].update(d.modules[k])
            else:
                self.modules[k] = copy.deepcopy(d.modules[k])

    def to_string(self):
        r = ''
        for i in self.modules.values():
            r += str(i)

        r += '*END OF\n'  # 7 characters is enough
        return r
