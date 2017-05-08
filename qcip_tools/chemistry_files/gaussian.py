import os
import math
import re

from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import ChemistryFile as qcip_ChemistryFile, apply_over_list

AuToAngstrom = 0.52917165


class Input(qcip_ChemistryFile):
    """Gaussian input file.

    **Input/Output class.**
    """

    file_type = 'GAUSSIAN_INPUT'

    def __init__(self):

        self.molecule = molecule.Molecule()

        self.raw_blocks = []
        self.other_blocks = []
        self.title = ''
        self.options = []
        self.options_dict = {}
        self.input_card = []

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        lines = f.readlines()

        temp = []
        for l in lines:
            if len(temp) > 0 and l.strip() == '':
                self.raw_blocks.append(temp)
                temp = []
            else:
                temp.append(l.strip())

        # add last
        if len(temp) > 0:
            self.raw_blocks.append(temp)
        if len(self.raw_blocks) < 3:
            raise Exception('The file contains less than 3 raw_blocks')

        # try to understand first block
        for line in self.raw_blocks[0]:
            if line[0] == '%':
                self.options.append(line[1:])
                opt = line[1:].split('=')
                self.options_dict[opt[0].lower()] = opt[1]
            else:
                self.input_card.append(line)

        # read title, charge and multiplicity:
        self.title = ' '.join(self.raw_blocks[1])
        charge_and_multiplicity = [int(i) for i in self.raw_blocks[2][0].strip().split()]
        self.molecule.charge = charge_and_multiplicity[0]
        self.molecule.multiplicity = charge_and_multiplicity[1]

        # read coordinates:
        coordinates = self.raw_blocks[2][1:]
        for atm in coordinates:
            s = atm.split()
            if len(s) != 4:
                raise Exception('Wrong atom definition "' + atm + '"')
            else:
                coords = [float(i) for i in s[1:]]
                atm_ = s[0]
                if atm_.isnumeric():
                    atm_ = atom.AtomicNumberToSymbol[int(atm_)]
                self.molecule.insert(atom.Atom(symbol=atm_, position=coords))

        # then store "other blocks"
        if len(self.raw_blocks) > 3:
            self.other_blocks = self.raw_blocks[3:]

    def to_string(self, chk='', original_chk=False):
        """

        :param chk: new chk name if needed
        :type chk: str
        :param original_chk: use original chk name
        :type original_chk: bool
        :rtype: str
        """

        rstr = ''

        # first, options:
        for key in self.options_dict:
            if not original_chk and key == 'chk':
                rstr += '%{}={}\n'.format(key, chk)
            else:
                rstr += '%{}={}\n'.format(key, self.options_dict[key])

        # then, input card and title
        rstr += '\n'.join(self.input_card) + '\n\n'
        rstr += self.title + '\n\n'

        # write charge and multiplicity
        rstr += str(self.molecule.charge) + ' ' + str(self.molecule.multiplicity) + '\n'

        # write coordinates
        rstr += self.molecule.output_atoms()
        rstr += '\n'

        # write other blocks :
        for i, block in enumerate(self.other_blocks):
            rstr += '\n'.join(block)
            rstr += '\n\n'
        return rstr

    def write(self, f, chk=None, original_chk=False):
        """

        :param f: File
        :type f: file
        """

        if chk is None:
            chk = os.path.basename('.'.join(f.name.split('.')[:-1]))

        f.write(self.to_string(chk=chk, original_chk=original_chk))


FCHK_AUTHORIZED_TYPES = ['I', 'R', 'C']
REGEX_FLOAT = re.compile(r'(?P<base>[0-9](\.[0-9]*)?)\-(?P<exp>[0-9]*)')


def transform_string_from_fchk(data, data_type):
    if data_type == 'I':
        return int(data)
    elif data_type == 'R':
        try:
            return float(data)
        except ValueError:
            if REGEX_FLOAT.match(data):
                return 0.0
    elif data_type == 'C':
        return data
    else:
        raise Exception('Unknown data type {}'.format(data_type))


class FCHKChunkInformation:
    """Where to find a given keyword in the FCHK file

    :param keyword: keyword related to this piece of information
    :type keyword: str
    :param data_type: type of this piece of information
    :type data_type: str
    :param data_length: length of this data
    :type data_length: int
    :param line_start: start in the source file
    :type line_start: int|object
    :param line_end: end in the source file
    :type line_end: int|object
    """

    def __init__(self, keyword, data_type, data_length, line_start, line_end):
        self.keyword = keyword
        self.data_type = data_type
        self.data_length = data_length
        self.line_start = line_start
        self.line_end = line_end


class FCHK(qcip_ChemistryFile):
    """A FCHK file. Based on the same principle as DataFile (split into chunks, interpret and store after).

    **Input only class.**
    """

    file_type = 'GAUSSIAN_FCHK'

    def __init__(self):
        self.molecule = molecule.Molecule()

        self.calculation_title = ''
        self.calculation_type = ''
        self.calculation_method = ''
        self.calculation_basis = ''

        self.chunks_information = {}
        self.chunks_parsed = {}
        self.lines = []

    def get(self, key):
        """Get the content from a given keyword

        :param key: tghe keyword
        :type key: str
        :rtype: list|int|float|str
        """
        content = []

        if key not in self.chunks_information:
            raise KeyError(key)

        if key not in self.chunks_parsed:
            chunk = self.chunks_information[key]

            if chunk.data_length != 1:
                matrix = self.lines[chunk.line_start:chunk.line_end]
                matrix = (''.join(matrix)).split()
                if len(matrix) == chunk.data_length:
                    for d in matrix:
                        content.append(transform_string_from_fchk(d, chunk.data_type))
                else:
                    raise Exception('Size is not the same as defined for {}'.format(key))
            else:
                content = self.lines[chunk.line_start][49:].strip()
                content = transform_string_from_fchk(content, chunk.data_type)

            self.chunks_parsed[key] = content
        else:
            content = self.chunks_parsed[key]

        return content

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        self.lines = f.readlines()

        # extract the first information :
        self.calculation_title = self.lines[0].strip()
        second_line = self.lines[1].split()
        self.calculation_type = second_line[0]
        self.calculation_method = second_line[1]
        self.calculation_basis = second_line[2]

        # now, crawl trough information
        for index, current_line in enumerate(self.lines[2:]):
            current_line_index = index + 2
            if current_line[0] != ' ':
                keyword = current_line[:40].strip()  # name is contained in the 40 first characters
                size = 1
                data_type = current_line[43]

                if data_type not in FCHK_AUTHORIZED_TYPES:
                    raise Exception('Type of {} ({}) is unknown'.format(keyword, data_type))

                is_matrix = current_line[47] == 'N'

                if is_matrix:
                    size = int(current_line[49:].strip())
                    values_per_line = 5
                    if data_type == 'I':
                        values_per_line = 6
                    num_lines = math.ceil(size / values_per_line)
                    line_start, line_end = current_line_index + 1, current_line_index + 1 + num_lines

                else:
                    line_start = line_end = current_line_index

                self.chunks_information[keyword] = FCHKChunkInformation(keyword, data_type, size, line_start, line_end)

        # finally, read the molecule
        nuclear_charges = self.get('Atomic numbers')
        atoms_coordinates = self.get('Current cartesian coordinates')
        weights = self.get('Real atomic weights')

        for index, c in enumerate(nuclear_charges):
            a = atom.Atom(
                atomic_number=c, position=[a * AuToAngstrom for a in atoms_coordinates[index * 3:index * 3 + 3]])
            a.mass = weights[index]
            self.molecule.insert(a)

        self.molecule.charge = self.get('Charge')
        self.molecule.multiplicity = self.get('Multiplicity')

        if self.get('Number of electrons') != self.molecule.number_of_electrons():
            raise ValueError('Number of electron does not match: {} != {}'.format(
                self.get('Number of electrons'), self.molecule.number_of_electrons()))

    def property(self, property_):
        """Rewritten to get access to keywords as well (which starts with uppercase letters,
        so cannot be variable name anyway).
        """

        if property_ in self.chunks_information:
            return self.get(property_)

        return super().property(property_)

    def __getitem__(self, item):
        """
        Implement access trough ``[]``
        `"""
        return self.get(item)

    def __contains__(self, item):
        """
        Implement test with ``in`` keyword
        """
        return item in self.chunks_information


class LinkCalled:
    """Remind when a link is called

    :param link: number of the link called
    :type link: int
    :type line_start: int
    :type line_end: int
    """

    def __init__(self, link, line_start, line_end):

        self.link = link
        self.line_start = line_start
        self.line_end = line_end

    def __repr__(self):
        return 'Link {}: {}:{}'.format(self.link, self.line_start, self.line_end)


class Output(qcip_ChemistryFile):
    """Log file of Gaussian. Contains a lot of informations, but not well printed or located.
    If possible, rely on the FCHK rather than on the LOG file.

    **Input only class.**

    .. note::

        ``self.geometry`` contains the input geometry, but it may be reoriented according to the symmetry if
        ``nosym`` was not provided. If it was not reoriented, ``self.input_orientation`` is set to ``True``.

    .. warning::

        Results are not guaranteed if the job ran without ``#P``.

    """

    file_type = 'GAUSSIAN_LOG'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.input_orientation = False

        self.links_called = []
        self.lines = []

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.lines = f.readlines()

        found_enter_link1 = False

        for index, line in enumerate(self.lines):
            if 'Entering Link 1' in line:
                self.links_called.append(LinkCalled(1, index, -1))
                found_enter_link1 = True
                break

        if not found_enter_link1:
            raise Exception('this does not seems to be a Gaussian output')

        for index, line in enumerate(self.lines):
            # List links (needs a job that ran with `#p`)
            if '(Enter ' in line:
                self.links_called[-1].line_end = index - 1

                link_num = line[-10:-6]
                if link_num[0] == 'l':
                    link_num = link_num[1:]
                current_link = int(link_num)
                self.links_called.append(LinkCalled(current_link, line_start=index, line_end=-1))

        # fetch molecule:
        charge_and_multiplicity_line = self.search('Charge = ', in_link=101)
        orientation_line = self.search('orientation:', in_link=202)

        if -1 in [charge_and_multiplicity_line, orientation_line]:
            raise Exception('Could not find geometry')

        self.input_orientation = self.lines[orientation_line][26] == 'I'

        for line in self.lines[orientation_line + 5:]:
            if '-------------' in line:
                break
            self.molecule.insert(atom.Atom(
                atomic_number=int(line[12:21].strip()),
                position=[float(a) for a in line[36:].split()]
            ))

        self.molecule.charge = int(self.lines[charge_and_multiplicity_line][9:13])
        self.molecule.multiplicity = int(self.lines[charge_and_multiplicity_line][27:])

    def link_called(self, link, all_times=False):
        """Is a link called in the file?

        :param link: the link to search
        :type link: int
        :param all_times: rater than stopping at the first time the link is found, go until the end and count
        :type all_times: bool
        :return: the number of times the link is found (max 1 if ``all_times`` is ``False``)
        :rtype: int
        """
        times = 0

        for link_info in self.links_called:
            if link_info.link == link:
                if not all_times:
                    return 1
                else:
                    times += 1

        return times

    def apply_function(self, func, line_start=0, line_end=None, in_link=None, **kwargs):
        """Apply ``apply_over_list()`` to the lines.

        :param lines: lines of the log
        :type lines: list
        :param func: function, for which the first parameter is the line, and the followings are the ``**kwargs``.
        :type func: callback
        :param line_start: starting index
        :type line_start: int
        :param line_end: end the search at some point
        :type line_end: int
        :param in_link: restrict the search to a given link
        :type in_link: int
        :type kwargs: dict
        :rtype: bool
        """

        if in_link is not None:
            if type(in_link) is not int:
                raise TypeError(in_link)

            for link_info in self.links_called:
                if link_info.link == in_link:
                    r = apply_over_list(self.lines, func, link_info.line_start, link_info.line_end, **kwargs)
                    if r:
                        return True

            return False

        else:
            return apply_over_list(self.lines, func, line_start, line_end, **kwargs)

    def search(self, s, line_start=0, line_end=None, in_link=None):
        """Returns the line when the string is found in this line or -1 if nothing was found.

        :rtype: int
        """

        class FunctionScope:
            line_found = -1

        def find_string(line, current_line_index):
            if s in line:
                FunctionScope.line_found = current_line_index
                return True
            else:
                return False

        self.apply_function(find_string, line_start=line_start, line_end=line_end, in_link=in_link)
        return FunctionScope.line_found
