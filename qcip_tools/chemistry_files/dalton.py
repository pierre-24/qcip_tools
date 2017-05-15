import tarfile
import math

from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import ChemistryFile as qcip_ChemistryFile, apply_over_list


AuToAngstrom = 0.52917165


class MoleculeInput(qcip_ChemistryFile):
    """Dalton mol input file.

    .. warning::

        Multiplicity seems to be given in the .dal file, so it may be wrong!!

    **I/O class.**"""

    file_type = 'DALTON_MOL'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.title = ''
        self.basis_set = ''

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        lines = f.readlines()
        self.from_read = True

        if len(lines) < 6:
            raise Exception('something wrong with dalton .mol: too short')

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
                position=[float(a) * (1 if in_angstrom else AuToAngstrom) for a in content[1:]])
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
                    a.symbol, *[p * (1 if in_angstrom else AuToAngstrom) for p in a.position])
        else:
            for symbol in self.molecule.symbols_contained:
                atms = self.molecule.atoms(symbol_in=[symbol])
                r += 'Charge={:.1f} Atoms={}\n'.format(atom.Definition[symbol][0], len(atms))

                for i in atms:
                    a = self.molecule[i]
                    r += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(
                        a.symbol, *[p * (1 if in_angstrom else AuToAngstrom) for p in a.position])

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


class ArchiveOutput(qcip_ChemistryFile):
    """Archive output of Dalton. Contains lots of information for who can extract them.

    **Input only class.**"""

    file_type = 'DALTON_ARCHIVE'

    def __init__(self):
        self.tar_file = None
        self.molecule = molecule.Molecule()

    def read(self, f):
        """Expects ``f`` to be open on binary mode.

        .. note::

            In Dalton archive, the molecule is given in a old, weird format, called "molplt". Information about it is
            given in `this webpage <https://brettbode.github.io/wxmacmolplt/Manual_pages/MacMolPlt_Files.html>`_.

        :param f: File
        :type f: file
        """

        self.from_read = True
        self.tar_file = tarfile.open(fileobj=f)

        # read molecule:
        mol_file = self.get_file('DALTON.MOL')

        lines = mol_file.readlines()

        info = lines[1].split()
        number_of_atoms = int(info[1])
        number_of_kinds = int(info[3])

        for i in range(number_of_atoms):
            a = lines[4 + number_of_kinds + i].split()
            self.molecule.insert(atom.Atom(symbol=str(a[0], 'utf-8'), position=[float(e) for e in a[1:]]))

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


class Output(qcip_ChemistryFile):
    """Output of Dalton.

    **Input only class.**"""

    file_type = 'DALTON_LOG'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.lines = []
        self.sections = []

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.lines = f.readlines()
        self.sections.append(OutputSection('START', 0, -1))

        for index, line in enumerate(self.lines):
            if '.----------------' in line:
                if 'Starting in' in self.lines[index + 1]:
                    self.sections[-1].line_end = index
                    section_end = self.lines[index + 1].rfind(')')
                    section_start = self.lines[index + 1].rfind('(')
                    if -1 in [section_end, section_start]:
                        raise Exception('no section found for line {}'.format(index))

                    self.sections.append(
                        OutputSection(self.lines[index + 1][section_start + 1:section_end], index + 1, -1))

        # fetch molecule
        coordinates_line = self.search('Cartesian Coordinates (a.u.)', in_section='START')
        if coordinates_line == -1:
            raise Exception('cannot find molecule')

        number_of_atoms = math.floor(int(self.lines[coordinates_line + 3][-5:]) / 3)
        for i in range(number_of_atoms):
            content = self.lines[coordinates_line + 4 + i].split()
            self.molecule.insert(
                atom.Atom(
                    symbol=content[0],
                    position=[a * AuToAngstrom for a in [float(content[4]), float(content[7]), float(content[10])]]
                ))

        charge_line = self.search('@    Total charge of the molecule', in_section='SIRIUS')
        if charge_line == -1:
            raise Exception('cannot find charge of the molecule')
        self.molecule.charge = int(self.lines[charge_line][-6:])
        self.molecule.multiplicity = int(self.lines[charge_line + 2][45:55])

        if self.molecule.number_of_electrons() != int(self.lines[charge_line - 2][-5:]):
            raise Exception(
                'Not the right number of electron, is the charge (found {}) wrong?'.format(self.molecule.charge))

    def apply_function(self, func, line_start=0, line_end=None, in_section=None, **kwargs):
        """Apply ``apply_over_list()`` to the lines.

        :param lines: lines of the log
        :type lines: list
        :param func: function, for which the first parameter is the line, and the followings are the ``**kwargs``.
        :type func: callback
        :param line_start: starting index
        :type line_start: int
        :param line_end: end the search at some point
        :type line_end: int
        :param in_section: restrict the search to a given link
        :type in_section: str
        :type kwargs: dict
        :rtype: bool
        """

        if in_section is not None:

            for link_info in self.sections:
                if link_info.section == in_section:
                    r = apply_over_list(self.lines, func, link_info.line_start, link_info.line_end, **kwargs)
                    if r:
                        return True

            return False

        else:
            return apply_over_list(self.lines, func, line_start, line_end, **kwargs)

    def search(self, s, line_start=0, line_end=None, in_section=None):
        """Returns the line when the string is found in this line or -1 if nothing was found.

        :rtype: int
        """

        class FunctionScope:
            line_found = -1

        def find_string(line, current_index):
            if s in line:
                FunctionScope.line_found = current_index
                return True
            else:
                return False

        self.apply_function(find_string, line_start=line_start, line_end=line_end, in_section=in_section)
        return FunctionScope.line_found
