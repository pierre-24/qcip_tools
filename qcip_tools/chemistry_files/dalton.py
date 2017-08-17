import tarfile
import math
import os
import collections
import io

from qcip_tools import molecule, atom, quantities
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, ChemistryLogFile, FormatError


class InputFormatError(FormatError):
    pass


class MoleculeInput(ChemistryFile, WithOutputMixin, WithMoleculeMixin):
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
    def attempt_recognition(cls, f):
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


class ArchiveOutput(ChemistryFile, WithMoleculeMixin):
    """Archive output of Dalton. Contains lots of information for who can extract them.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.tar_file``: pipe to the tar file (``tarfile.TarFile``)

    """

    #: The identifier
    file_type = 'DALTON_ARCHIVE'

    def __init__(self):
        self.tar_file = None
        self.molecule = molecule.Molecule()

    @classmethod
    def possible_file_extensions(cls):
        return ['gz', 'tar.gz']

    @classmethod
    def attempt_recognition(cls, f):
        """A dalton archive does contain a lot of files name "DALTON"
        """

        try:
            t = tarfile.open(f.name)
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


class OutputFormatError(FormatError):
    pass


class Output(ChemistryLogFile, WithMoleculeMixin):
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
    def attempt_recognition(cls, f):
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
        if coordinates_line == -1:
            raise OutputFormatError('cannot find molecule in START')

        number_of_atoms = math.floor(int(self.lines[coordinates_line + 3][-5:]) / 3)
        for i in range(number_of_atoms):
            content = self.lines[coordinates_line + 4 + i].split()
            self.molecule.insert(
                atom.Atom(
                    symbol=content[0],
                    position=[a * quantities.AuToAngstrom for a in [
                        float(content[4]), float(content[7]), float(content[10])]]
                ))

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


class Input(ChemistryFile, WithOutputMixin):
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
    def attempt_recognition(cls, f):
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

    def to_string(self):
        r = ''
        for i in self.modules.values():
            r += str(i)

        r += '*END OF\n'  # 7 characters is enough
        return r
