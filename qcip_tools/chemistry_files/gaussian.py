import os
import math
import re

from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import File as qcip_File

AuToAngstrom = 0.52917165


class Input(qcip_File):
    """Gaussian input file. I/O class."""

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


class FCHK(qcip_File):
    """A FCHK file. Based on the same principle as DataFile (split into chunks, interpret and store after)
    """

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
