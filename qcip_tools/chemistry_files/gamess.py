from qcip_tools.molecule import Molecule
from qcip_tools.atom import Atom
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, ChemistryLogFile, \
    FormatError, WithIdentificationMixin
from qcip_tools.quantities import AuToAngstrom


class InputFormatError(FormatError):
    pass


class InputModule(object):
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
            raise InputFormatError('Incorrect number of $ in {}'.format(s))

        lines = s.strip().split('\n')
        if lines[0][0] != '$':
            raise InputFormatError('malformed input {} (not starting with $)'.format(s))
        if lines[-1][-4:].lower() != '$end':
            raise InputFormatError('malformed input {} (not ending with $END)'.format(s))

        if len(lines[0].split()) != 1:
            full_module = ' '.join(x for x in lines if x.strip()[0] != '!')
            l = full_module.split()
            title = l[0][1:].lower()
            rest_of_module = [a.split('=') for a in l[1:-1]]
            return InputModule(title, rest_of_module)
        else:
            return InputModule(
                lines[0][1:].strip().lower(),
                [a.split() for a in lines[1:-1] if a.strip()[0] != '!'],
                special_module=True)

    def __repr__(self):
        """Note: shifts everything by one space"""

        r = ' $' + self.name.upper()

        if self.special_module:
            r += '\n'

        length_line = len(r)

        for option in self.options:
            if self.special_module:
                r += ' ' + ' '.join(option) + '\n'
            else:
                to_add = ' {}={}'.format(*option)
                if length_line + len(to_add) > 79:
                    r += '\n'
                    length_line = 0
                length_line += len(to_add)
                r += ' {}={}'.format(*option)

        r += ' $END'

        return r

    def __contains__(self, item):
        if self.special_module:
            raise Exception('cannot search keyword in special module')

        return item.lower() in [a[0].lower() for a in self.options]


class Input(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """GAMESS (US) input file.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.modules``: the different modules (``dict`` of
          `InputModule <#qcip_tools.chemistry_files.gamess.InputModule>`_,
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
    def attempt_identification(cls, f):
        """A GAMESS input file contains too many "$"!!!"""

        num_dollar = 0
        num_end = 0

        for i in range(10):
            l = f.readline()
            num_dollar += l.count('$')
            num_end += l.lower().count('$end')

        return num_dollar > 4 and num_end > 2

    @classmethod
    def from_molecule(cls, molecule, title='', modules=None, *args, **kwargs):
        """Create a file from molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :param title: title of the run
        :type title: str
        :param modules: running module
        :type modules: list
        :rtype: qcip_tools.chemistry_files.gamess.Input
        """

        obj = super().from_molecule(molecule, *args, **kwargs)
        obj.title = title
        obj.modules = {}

        for m in modules:
            if type(m) is str:
                m = InputModule.from_string(m)
            obj.modules[m.name] = m

        return obj

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
                raise InputFormatError('module started but never ended')

            if pos_end + 4 > len_content:
                raise InputFormatError('malformed end')

            if content[pos_end + 1:pos_end + 4].lower() != 'end':
                raise InputFormatError('no $END!')

            m = InputModule.from_string(content[pos_begin:pos_end + 4])
            self.modules[m.name] = m
            prev_pos = pos_end

        if len(self.modules) == 0:
            raise Exception('empty file')

        # get molecule:
        if 'data' not in self:
            raise InputFormatError('no $DATA')

        data = self.modules['data'].options
        if len(data) < 3:
            raise InputFormatError('$DATA is too small')
        if len(data[1]) != 1 or data[1][0].lower() != 'c1':
            raise NotImplementedError('Symmetry != c1')

        self.title = ' '.join(data[0])
        for atom_def in data[2:]:
            if len(atom_def) != 5:
                raise InputFormatError('wrong atom definition: {}'.format(atom_def))

            try:
                atomic_number = int(atom_def[1][:-2])
            except ValueError:
                raise InputFormatError('not a correct atomic number: {}'.format(atom_def[1]))

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


class OutputStep:
    """Remind when a "step" is entered and exited

    :param step_name: name of the step
    :type step_name: str
    :type line_start: int
    :type line_end: int
    """

    def __init__(self, step_name, line_start, line_end):

        self.step_name = step_name
        self.line_start = line_start
        self.line_end = line_end

    def __repr__(self):
        return 'Step {}: {}:{}'.format(self.step_name, self.line_start, self.line_end)


class OutputFormatError(FormatError):
    pass


class Output(ChemistryLogFile, WithMoleculeMixin, WithIdentificationMixin):
    """Output of GAMESS.

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.lines``: the lines of the file (``list`` of ``str``)
        + ``self.chunks``: the different sections (``list`` of
          `OutputStep <#qcip_tools.chemistry_files.gamess.OutputStep>`_)

    """

    #: The identifier
    file_type = 'GAMESS_LOG'
    chunk_title_variable = 'step_name'

    def __init__(self):
        self.molecule = Molecule()

    @classmethod
    def attempt_identification(cls, f):
        """Same trick as Dalton: check for universities and states name
        """

        count = 0
        found_states = 0
        found_universities = 0
        found_gamess = 0
        states = ['IOWA', 'MICHIGAN', 'CALIFORNIA', 'PENNSYLVANIA', 'NEBRASKA']

        for l in f.readlines():
            count += 1

            if 'GAMESS' in l:
                found_gamess += 1

            if 10 < count < 100:
                if 'UNIVERSITY' in l:
                    found_universities += 1
                for c in states:
                    if c in l:
                        found_states += 1
                        break
            if count > 100:
                break

        return found_states > 5 and found_universities > 10 and found_gamess > 2

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        super().read(f)

        self.chunks.append(OutputStep('X', 0, -1))

        for index, line in enumerate(self.lines):
            if 'EXECUTION OF GAMESS BEGUN' in line:
                self.chunks[-1].line_start = index
                continue
            elif '.....' in line and 'STEP CPU' in self.lines[index + 1]:
                title = line.replace('.', '')\
                    .replace('END OF ', '')\
                    .replace('DONE WITH ', '')\
                    .replace('DONE ', '')\
                    .strip()

                self.chunks[-1].line_end = index
                self.chunks[-1].step_name = title
                self.chunks.append(OutputStep('X', index + 1, -1))
                continue
            elif 'EXECUTION OF GAMESS TERMINATED' in line:
                self.chunks.pop()
                break

        # molecule
        line_coordinates = self.search('COORDINATES (BOHR)', into='SETTING UP THE RUN')

        if line_coordinates < 0:
            raise OutputFormatError('no molecule is found in SETTING UP THE RUN')

        for line in self.lines[line_coordinates + 2:]:
            atom_def = line.split()
            if len(atom_def) == 0:
                break

            elif len(atom_def) != 5:
                raise OutputFormatError('Wrong molecule definition: {}'.format(line))

            try:
                atomic_number = int(atom_def[1][:-2])
            except ValueError:
                raise OutputFormatError('not a correct atomic number: {} in {}'.format(atom_def[1], line))

            self.molecule.insert(
                Atom(atomic_number=atomic_number, position=[float(a) * AuToAngstrom for a in atom_def[2:]]))

        # charge and multiplicity
        num_of_electron_line = self.search(
            'NUMBER OF ELECTRONS', line_start=line_coordinates + len(self.molecule), into='SETTING UP THE RUN')

        if num_of_electron_line < 0:
            raise OutputFormatError('cannot find charge and multiplicity')

        number_of_electrons = int(self.lines[num_of_electron_line][-5:])
        self.molecule.charge = int(self.lines[num_of_electron_line + 1][-5:])
        self.molecule.multiplicity = int(self.lines[num_of_electron_line + 2][-5:])

        if number_of_electrons != self.molecule.number_of_electrons():
            raise OutputFormatError('not the same number of electron {} != {}'.format(
                number_of_electrons, self.molecule.number_of_electrons()))
