import os
import math
import re
import numpy
import copy

from qcip_tools import molecule, atom, quantities
from qcip_tools.chemistry_files import ChemistryFile as qcip_ChemistryFile, apply_over_list


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

    @classmethod
    def possible_file_extensions(cls):
        return ['com', 'inp']

    @classmethod
    def attempt_recognition(cls, f):
        """A gaussian input should contains '%' or '#' as the first character of the line of the "first block"
        """

        found_some_lines = 0
        found_input_line = False

        for l in f.readlines():
            if l.strip() == '':
                break
            if l[0] not in ['%', '#']:
                if not found_some_lines:
                    return False
            else:
                if l[0] == '#':
                    found_input_line = True
                found_some_lines += 1

        if not found_some_lines:
            return False

        if not found_input_line:
            return False

        # TODO: if ``found_some_lines`` is equal to 1, maybe an extra test ? (charge and multiplicity, for example)
        return True

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

        self.title = ''
        self.calculation_type = ''
        self.calculation_method = ''
        self.basis_set = ''

        self.chunks_information = {}
        self.chunks_parsed = {}
        self.lines = []

    @classmethod
    def possible_file_extensions(cls):
        return ['fchk']

    @classmethod
    def attempt_recognition(cls, f):
        """A gaussian fchk starts with two two line of info, then "Number of atoms"
        """

        i = 0

        for l in f.readlines():
            if i == 2:
                if 'Number of atoms' not in l:
                    return False
            if i == 3:
                if 'Info1-9' in l:
                    return True
                return False
            i += 1

        return False

    def get(self, key):
        """Get the content from a given keyword

        :param key: the keyword
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
        self.title = self.lines[0].strip()
        second_line = self.lines[1].split()
        self.calculation_type = second_line[0]
        self.calculation_method = second_line[1]
        self.basis_set = second_line[2]

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
                atomic_number=c,
                position=[a * quantities.AuToAngstrom for a in atoms_coordinates[index * 3:index * 3 + 3]])
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

    @classmethod
    def possible_file_extensions(cls):
        return ['log', 'out']

    @classmethod
    def attempt_recognition(cls, f):
        """A gaussian log ... Contains a lot of "gaussian" in the beginning (limit to the 150 first lines)"
        """

        count = 0
        num_of_gaussian = 0

        for l in f.readlines():
            if count > 100:
                break
            if 'Gaussian' in l or 'gaussian' in l:
                num_of_gaussian += 1
            count += 1

        return num_of_gaussian > 10

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

        def find_string(line, current_index):
            if s in line:
                FunctionScope.line_found = current_index
                return True
            else:
                return False

        self.apply_function(find_string, line_start=line_start, line_end=line_end, in_link=in_link)
        return FunctionScope.line_found


class ChargeTransferInformation:
    """Analyzed charge transfer.

    In it, ``charge`` is the amount of charge (in e) transfered from ``barycenter_m`` to ``barycenter_p``,
    linked by ``vector``, which lead to a certain charge transfer ``distance``.

    :param charge: charge
    :type charge: float
    :param barycenter_m: barycenter of the negative charges
    :type barycenter_m: numpy.array
    :param barycenter_p: barycenter of the positive charges
    :type barycenter_p: numpy.array
    """

    def __init__(self, charge, barycenter_m, barycenter_p):
        self.charge = charge
        self.barycenter_m = barycenter_m
        self.barycenter_p = barycenter_p
        self.vector = barycenter_m - barycenter_p
        self.distance = numpy.linalg.norm(self.vector)


class Cube(qcip_ChemistryFile):
    """Gaussian cube.

    Documentation on that format can be found `here <http://gaussian.com/cubegen/>`_.

    **Input/output class**."""

    file_type = 'GAUSSIAN_CUBE'

    def __init__(self):

        self.molecule = molecule.Molecule()
        self.cube_type = ''
        self.records_per_direction = [0, 0, 0]
        self.origin = numpy.zeros(3)
        self.increments = numpy.zeros(3)
        self.data_per_record = 1
        self.records = None
        self.title = ''
        self.subtitle = ''
        self.MOs = ''

    @classmethod
    def possible_file_extensions(cls):
        return ['cub', 'cube']

    @classmethod
    def attempt_recognition(cls, f):
        """A cube file contains, in line 3-6, quadruplets of numbers, the first one being an integer, then 3 floats.
        Then, after the different atoms, line should contains numbers in scientific notation.
        """

        count = 0
        number_of_atoms = 0

        for l in f.readlines():
            count += 1
            # line 3-6:
            if 7 > count > 2:
                c = l.split()
                if 6 < len(c) < 4:
                    return False

                if not c[0].isdigit():
                    return False

                for i in c[1:]:
                    try:
                        float(i)
                    except ValueError:
                        return False

                if count == 3:
                    number_of_atoms = abs(int(c[0]))
            # only numbers in scientific notation in records
            if count == 7 + number_of_atoms:
                c = l.split()

                for i in c:
                    if 'E' not in i:
                        return False
                return True

        return False

    def read(self, f):
        """

        .. warning::

            The charge/multiplicity may be wrong !

        .. note::

            Increments and origin are stored as found, in bohr.

        :param f: File
        :type f: file
        """

        self.from_read = True

        lines = f.readlines()
        additional_line = False
        self.cube_type = 'unknown'
        self.data_per_record = 1

        # titles:
        self.title = lines[0].strip()
        self.subtitle = lines[1].strip()

        self.MOs = []

        # parameters
        params = lines[2].split()
        num_of_atoms = int(params[0])
        if num_of_atoms < 0:
            additional_line = True
            num_of_atoms = - num_of_atoms
            self.cube_type = 'MO'

        self.origin = numpy.array([float(a) for a in params[1:4]])

        if len(params) > 4:
            self.data_per_record = int(params[4])

        # number of records and increments
        for i in range(0, 3):
            self.increments[i] = float(lines[3 + i].split()[i + 1])

        for i in range(0, 3):
            num = int(lines[i + 3][0:8].strip())
            self.records_per_direction[i] = num

        # atoms:
        self.molecule = molecule.Molecule()
        for i in range(0, num_of_atoms):
            atm = lines[i + 6].split()
            symbol = atom.AtomicNumberToSymbol[int(atm[0])]
            atom_ = atom.Atom(
                symbol=symbol, position=numpy.array([float(a) * quantities.AuToAngstrom for a in atm[2:5]]))
            self.molecule.insert(atom_)

        # records
        start_line = 6 + num_of_atoms

        if additional_line:
            # number of stuffs per line
            extra_line = lines[start_line].strip().split()
            self.data_per_record = int(extra_line[0])

            for i in extra_line[1:]:
                self.MOs.append(int(i))

            start_line += 1

        total_num_of_data = self.number_of_data()
        actual_num_of_data = 0

        self.records = numpy.zeros(total_num_of_data)

        for line in lines[start_line:]:
            raw = [float(a) for a in line.strip().split()]
            l = len(raw)
            self.records[actual_num_of_data:actual_num_of_data + l] = raw
            actual_num_of_data += l

        reshape = self.records_per_direction.copy()
        reshape.append(self.data_per_record)

        self.records = self.records.reshape(tuple(reshape))

        if actual_num_of_data != total_num_of_data:
            raise Exception(
                'num of record is not equal to total num of records ({}!={})'.format(
                    actual_num_of_data, total_num_of_data))

    def number_of_records(self):
        """Number of records in the file

        .. note::

            According to the documentation:
            "Cube files have one row per record (i.e., N1*N2 records each of length N3*NVal)" (where N1 and N2
            are the record per direction in the X and Y direction, while N3 is for the Z direction, the fastest
            running index).

        :rtype: int
        """

        return numpy.prod(self.records_per_direction[:2])

    def number_of_data(self):
        """Number of points in the file

        :rtype: int
        """
        return self.number_of_records() * self.records_per_direction[2] * self.data_per_record

    def dV(self):
        """Return the size of the smallest unit of volume in the cube, in angstrom**3.

        :rtype: float
        """

        return numpy.prod(self.increments) * quantities.AuToAngstrom ** 3

    def to_string(self):
        """

        :rtype: str
        """

        r = '{}\n{}\n'.format(self.title, self.subtitle)

        num_of_atoms = len(self.molecule)
        if self.cube_type == 'MO':
            num_of_atoms = - num_of_atoms

        # origin and increments
        r += '{:5d} {:11f} {:11f} {:11f} {:5d}\n'.format(
            num_of_atoms, self.origin[0], self.origin[1], self.origin[2], self.data_per_record)

        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[0], self.increments[0], 0.0, 0.0)
        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[1], 0.0, self.increments[1], 0.0)
        r += '{:5d} {:11f} {:11f} {:11f}\n'.format(self.records_per_direction[2], 0.0, 0.0, self.increments[2])

        # molecule
        for a in self.molecule:
            r += '{:5d} {:11f} {:11f} {:11f} {:11f}\n'.format(
                a.atomic_number, a.atomic_number, *(a.position / quantities.AuToAngstrom))

        if self.cube_type == 'MO':
            r += '    {}  {}\n'.format(len(self.MOs), '  '.join(str(m) for m in self.MOs))

        # records:
        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    for index_data in range(self.data_per_record):
                        if (z * self.data_per_record + index_data) % 6 == 0 and z != 0:
                            r += '\n'

                        r += ' {: .5E}'.format(self.records[x][y][z][index_data])

                r += '\n'

        return r

    def alike(self, other):
        """Compare two cubes to check if they have the same size

        :param other: the other cube
        :type other: Cube
        :rtype: bool
        """

        if self.number_of_data() != other.number_of_data():
            return False

        if not numpy.array_equal(self.increments, other.increments):
            return False

        if self.records.shape != other.records.shape:
            return False

        return True

    def __add__(self, other):
        """Sum the values of two cubes

        :param other: other cube
        :type other: Cube
        :rtype: Cube
        """

        if not self.alike(other):
            raise ArithmeticError('cannot sum with this cube')

        new_cube = copy.deepcopy(self)
        new_cube.records += other.records

        return new_cube

    def __sub__(self, other):
        """Subtract the values of two cubes

        :param other: other cube
        :type other: Cube
        :rtype: Cube
        """

        if not self.alike(other):
            raise ArithmeticError('cannot sub with this cube')

        new_cube = copy.deepcopy(self)
        new_cube.records -= other.records

        return new_cube

    def __pow__(self, power, modulo=None):
        """Power the records

        :param power: power
        :rtype power: float
        :param modulo: ?
        :rtype modulo: int
        :rtype: Cube
        """

        new_cube = copy.deepcopy(self)
        new_cube.records **= power

        return new_cube

    def slice(self, index=0):
        """Get a specific part of the data

        :param index: index of the record
        :type index: int
        :rtype: Cube
        """

        if index < 0 or index >= self.data_per_record:
            raise ValueError(index)

        new_cube = copy.deepcopy(self)
        new_shape = list(self.records.shape)
        new_shape[-1] = 1
        new_cube.records = self.records[:, :, :, index].reshape(new_shape)
        new_cube.data_per_record = 1

        if self.cube_type == 'MO':
            new_cube.MOs = [self.MOs[index]]

        return new_cube

    def positions(self):
        """Return an array with the position (in angstrom) of each point.

        :rtype: numpy.array
        """

        s = self.records_per_direction[:]
        s.append(3)
        positons = numpy.zeros(tuple(s))

        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    current_position = self.origin + numpy.array([x, y, z]) * self.increments
                    positons[x, y, z] = current_position * quantities.AuToAngstrom
        return positons

    def compute_charge_transfer(self, index=0):
        """
        Analyse CT by computing:

        - Transferred charge ;
        - Barycenter positions (the negative one being the starting point of the transfer) ;
        - Distance between the two barycenters.

        Made according to *J. Chem. Theory. Comput.* **7**, 2498 (2011).

        :param index: slice of the record to use
        :type index: int
        :rtype: ChargeTransferInformation
        """

        dV = self.dV()

        # positive and negative densities:
        rho_p = self.records.copy()
        rho_p[numpy.where(rho_p < 0.0)] = 0.0

        charge_transferred = rho_p.sum() * dV

        rho_m = self.records.copy()
        rho_m[numpy.where(rho_p > 0.0)] = 0.0

        # find barycenter positions (in angstrom)
        positions = self.positions()
        bar_p = numpy.zeros(3)
        bar_m = numpy.zeros(3)

        for i in range(3):
            bar_p[i] = numpy.sum(positions[:, :, :, i] * rho_p[:, :, :, index]) / charge_transferred * dV
            bar_m[i] = - numpy.sum(positions[:, :, :, i] * rho_m[:, :, :, index]) / charge_transferred * dV

        return ChargeTransferInformation(
            charge=charge_transferred / (quantities.AuToAngstrom ** 3),  # charge in |e|
            barycenter_m=bar_m,
            barycenter_p=bar_p)

    def sum_density_of_sets(self, bounding_sets, index=0):
        """Sum the density in a given space, defined by the bounding boxes

        :param bounding_sets: sets of bounding boxes
        :type bounding_sets: dict
        :param index: slice of the record to use
        :type index: int
        :return: dictionary with the sum of the density
        :rtype: dict
        """

        if self.data_per_record != 1:
            raise ValueError(self.data_per_record)

        if len(bounding_sets) == 0:
            return {}

        sums = {}

        for k in bounding_sets:
            sums[k] = .0

        dV = self.dV() / (quantities.AuToAngstrom ** 3)

        for x in range(self.records_per_direction[0]):
            for y in range(self.records_per_direction[1]):
                for z in range(self.records_per_direction[2]):
                    current_position = (self.origin + [x, y, z] * self.increments) * quantities.AuToAngstrom
                    current_value = self.records[x, y, z, index] * dV

                    for key, bounding_set in bounding_sets.items():
                        if current_position in bounding_set:
                            sums[key] += current_value

        return sums
