from qcip_tools.molecule import Molecule
from qcip_tools.atom import Atom
from qcip_tools.chemistry_files import InputOutputChemistryFile


class InputModule():
    """GAMESS input module (ends with ``$END``).

    A "special" module is a multiline module (i.e. ``$DATA`` or ``$TDHFX``)
    """

    def __init__(self, name, options, special_module=False):
        self.name = name
        self.options = options
        self.special_module = special_module

    @staticmethod
    def from_string(s):
        """Create a module out of string

        :param s: the string
        :type s: str
        :rtype: qcip_tools.chemistry_files.gamess.InputModule
        """

        if s.count('$') != 2:
            raise Exception('Incorrect number of $ in {}'.format(s))

        lines = s.strip().split('\n')
        if lines[0][0] != '$':
            raise Exception('malformed input {} (not starting with $)'.format(s))
        if lines[-1][-4:].lower() != '$end':
            raise Exception('malformed input {} (not ending with $END)'.format(s))

        if len(lines) == 1:
            l = lines[0].split()
            return InputModule(l[0][1:].lower(), [a.split('=') for a in l[1:-1]])
        else:  # "special" (= multiline) module
            if len(lines[0].split()) != 1:
                raise Exception('not a correct input: {}'.format(lines[0]))
            return InputModule(
                lines[0][1:].strip().lower(),
                [a.split() for a in lines[1:-1]],
                special_module=True)

    def __repr__(self):
        """Note: shifts everything by one space"""

        r = ' $' + self.name.upper()

        if self.special_module:
            r += '\n'

        for option in self.options:
            if self.special_module:
                r += ' ' + ' '.join(option) + '\n'
            else:
                r += ' {}={}'.format(*option)

        r += ' $END'

        return r

    def __contains__(self, item):
        if self.special_module:
            raise Exception('cannot search keyword in special module')

        return item.lower() in [a[0].lower() for a in self.options]


class Input(InputOutputChemistryFile):
    """GAMESS (US) input file.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.modules``: the different modules (``dict`` of ``InputModule``,
          where the key is the name of the module **in lowercase**)
        + ``self.title``: the title of the run (from $DATA)

    .. warning::

        All module are stored in lowercase.

    .. warning::

        Does not implement Z-matrix molecule definition and other symmetries than C1.

    """

    #: The identifier
    file_type = 'GAMESS_INP'

    def __init__(self):
        self.molecule = Molecule()
        self.modules = {}
        self.title = ''

    @classmethod
    def possible_file_extensions(cls):
        return ['inp']

    @classmethod
    def attempt_recognition(cls, f):
        """A GAMESS input file contains too many "$"!!!"""

        num_dollar = 0
        num_end = 0

        for i in range(10):
            l = f.readline()
            num_dollar += l.count('$')
            num_end += l.lower().count('$end')

        return num_dollar > 4 and num_end > 2

    def __contains__(self, item):
        return item.lower() in self.modules

    def __getitem__(self, item):
        return self.modules[item.lower()]

    def __setitem__(self, key, value):
        if not isinstance(value, InputModule):
            raise TypeError(value)

        value.title = value.title.lower()

        if value.title != key.lower():
            raise Exception('cannot set {} while title is {}'.format(value.title, key.lower()))

        self.modules[value.title] = value

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        content = f.read()
        len_content = len(content)
        self.from_read = True
        prev_pos = 0

        while True:
            pos_begin = content.find('$', prev_pos + 1)

            if pos_begin < 0:
                break

            pos_end = content.find('$', pos_begin + 1)

            if pos_end < 0:
                raise Exception('module started but never ended')

            if pos_end + 4 > len_content:
                raise Exception('malformed end')

            if content[pos_end + 1:pos_end + 4].lower() != 'end':
                raise Exception('no end !')

            m = InputModule.from_string(content[pos_begin:pos_end + 4])
            self.modules[m.name] = m
            prev_pos = pos_end

        if len(self.modules) == 0:
            raise Exception('empty file')

        # get molecule:
        if 'data' not in self:
            raise Exception('no $DATA')

        data = self.modules['data'].options
        if len(data) < 3:
            raise Exception('$DATA is too small')
        if len(data[1]) != 1 or data[1][0].lower() != 'c1':
            raise NotImplementedError('Symmetry != c1')

        self.title = ' '.join(data[0])
        for atom_def in data[2:]:
            if len(atom_def) != 5:
                raise Exception('wrong atom definition: {}'.format(atom_def))

            try:
                atomic_number = int(atom_def[1][:-2])
            except ValueError:
                raise Exception('not a correct atomic number: {}'.format(atom_def[1]))

            self.molecule.insert(Atom(atomic_number=atomic_number, position=[float(a) for a in atom_def[2:]]))

    def to_string(self):
        r = ''

        # first, all modules except $DATA
        for module_ in self.modules:
            if module_ == 'data':
                continue

            r += str(self.modules[module_]) + '\n'

        # then, $DATA
        r += ' $DATA\n'
        r += ' '.join(self.title.split('\n')) + '\nC1\n'
        for a in self.molecule:
            r += '{:3}   {:>2}.0   {: 10.6f}   {: 10.6f}   {: 10.6f}\n'.format(a.symbol, a.atomic_number, *a.position)
        r += ' $END'

        return r
