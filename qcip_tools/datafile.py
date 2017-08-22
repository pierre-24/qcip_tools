import math
import struct

#: Allowed data types in datafiles: ``I`` and ``R`` correspond to a **list** of integer or floats,
#: while ``S`` defines a list of characters, so a string.
ALLOWED_DATA_TYPES = ('I', 'R', 'S')


class InvalidDataFile(Exception):
    pass


def transform_string_from_type(data, data_type):
    """Transform a string to the good type

    :param data: a string
    :type data: str
    :param data_type: the type
    :type data_type: str
    :return: the data
    """
    if data_type == 'I':
        return int(data)
    elif data_type == 'R':
        return float(data)
    elif data_type == 'S':
        return data  # data is already a string ...
    else:
        raise DataFileTypeError(data_type)


class ChunkInformation:
    """Information on the location of a given chunk of information

    :param keyword: keyword related to this piece of information
    :type keyword: str
    :param data_type: type of this piece of information
    :type data_type: str
    :param data_length: length of this data
    :type data_length: int
    :param offset_start: start in the source file
    :type offset_start: int|object
    :param offset_end: end in the source file
    :type offset_end: int|object
    :param modified: indicate whether this chunk is was modifed during the execution
    :type modified: bool
    """

    def __init__(self, keyword, data_type, data_length, offset_start=-1, offset_end=-1, modified=False):
        if data_type not in ALLOWED_DATA_TYPES:
            raise DataFileTypeError(data_type)

        self.keyword = keyword
        self.data_type = data_type
        self.data_length = data_length
        self.offset_start = offset_start
        self.offset_end = offset_end
        self.modified = modified

    def __repr__(self):
        return '{} ({}{}): {}:{}'.format(
            self.keyword, self.data_type, self.data_length, self.offset_start, self.offset_end)


class DataFileTypeError(TypeError):
    def __init__(self, s):
        super().__init__('type {} is not in {}'.format(s, ', '.join(ALLOWED_DATA_TYPES)))


class DataFile:
    """
    Base class for data files.

    Implement a "lazy loading" behavior: ``read()`` should only slice the file in chunks, which are stored.
    When the value of the keyword is requested, the chunck is parsed with ``parse_chunk()``.

    To do so, ``read()`` should fill the variable ``self.raw`` with the content of the file, and parse it into
    chunks and filling ``self.chunks_information`` with information about that chunk (especially ``start_offset`` and
    ``end_offset``). Then, the value of ``self.raw`` and ``self.chunks_information`` can be used when ``parse_chunk()``
    is called (by the ``get()`` function, probably), which should triggers the parsing of the chunk, and should set
    the result in the corresponding ``self.chunks_parsed[keyword]``.

    :param filename: read file
    :type filename: str
    :param pipe: read a pipe
    :type pipe: file

    """

    def __init__(self):

        self.chunks_information = {}  #: Dict of ``ChunkInformation`` for **all** available keywords in the source file.
        self.chunks_parsed = {}  #: parsed chunks, so dict of stuffs
        self.raw = None  #: variable for storage of chunks for ``parse_chunk()``,  so "whatever floats you boat".

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

    def read(self, f):
        """Read in a pipe

        .. note::

            Must be implemented in child class.

        :param f: pipe
        :type f: file
        """
        raise NotImplementedError()

    def write(self, f):
        """Write in a pipe

        .. note::

            Must be implemented in child class.

        :param f: pipe
        :type f: file
        """
        raise NotImplementedError()

    def set(self, key, data_type, data, force_change_type=False):
        """Set a keyword to a given value

        :param key: the key
        :type key: str
        :param data_type: type of data
        :type data_type: str
        :param data: the data
        :type data: iterable
        :param force_change_type: do not check type of data before overwritting it
        :type force_change_type: bool
        """

        if data_type not in ALLOWED_DATA_TYPES:
            raise DataFileTypeError(data_type)

        if not force_change_type:
            if key in self.chunks_information and self.chunks_information[key].data_type != data_type:
                raise TypeError('new type {} != previous type {}'.format(
                    data_type, self.chunks_information[key].data_type))

        self.chunks_information[key] = ChunkInformation(key, data_type, len(data), modified=True)
        self.chunks_parsed[key] = data

    def get(self, key, data_type=None):
        """Get the value of a keyword

        :param key: the key
        :type key: str
        :param data_type: type of data
        :type data_type: str
        :rtype: iterable
        """

        if key in self.chunks_information:
            if key not in self.chunks_parsed:  # not yet parsed
                self.parse_chunk(keyword=key)

            if data_type is not None:
                if data_type != self.chunks_information[key].data_type:
                    raise TypeError('new type {} != previous type {}'.format(
                        data_type, self.chunks_information[key].data_type))

            return self.chunks_parsed[key]
        else:
            raise KeyError(key)

    def parse_chunk(self, keyword):
        """Parse a chunk

        .. note::

            Must be implemented in child class.

        :param keyword: keyword
        :type keyword: str
        """
        raise NotImplementedError()


class TextDataFile(DataFile):
    """
    Object that save data in the form of a text file, more or less readable by an human being
    """

    def write(self, f):
        for k in self.chunks_information:
            info = self.chunks_information[k]

            if info.modified:  # avoid interpreting data if not needed
                val = self.get(k)

                f.write('{}{} {}\n'.format(info.data_type, info.data_length, k))

                if info.data_type == 'S':
                    num_of_lines = math.ceil(len(val) / 79)
                    for i in range(num_of_lines):
                        f.write(' {}\n'.format(val[78 * i: 78 * (i + 1)]))
                else:
                    for index, val_ in enumerate(val):
                        if index != 0:
                            if info.data_type == 'R' and index % 3 == 0:
                                f.write('\n')
                            if info.data_type == 'I' and index % 5 == 0:
                                f.write('\n')

                        if info.data_type == 'R':
                            f.write(' {: .15e}'.format(val_))
                        if info.data_type == 'I':
                            f.write(' {: 15d}'.format(val_))
                    f.write('\n')
            else:
                f.write(''.join(self.raw[info.offset_start - 1:info.offset_end]))

    def read(self, f):

        self.raw = f.readlines()
        self.chunks_information = {}
        self.chunks_parsed = {}

        for index, line in enumerate(self.raw):
            if line.strip() == '':
                continue

            if line[0] != ' ':
                space = line.find(' ')
                data_type, data_length, data_keyword = line[:1], int(line[1:space]), line[space:].strip()

                line_start = line_end = index + 1
                if data_type == 'I':
                    line_end += math.ceil(data_length / 5)
                elif data_type == 'R':
                    line_end += math.ceil(data_length / 3)
                else:
                    line_end += math.ceil(data_length / 79)

                self.chunks_information[data_keyword] = \
                    ChunkInformation(data_keyword, data_type, data_length, line_start, line_end)

    def parse_chunk(self, keyword):
        if keyword not in self.chunks_information:
            raise KeyError(keyword)

        info = self.chunks_information[keyword]

        if info.data_type != 'S':

            matrix = self.raw[info.offset_start:info.offset_end]
            matrix = (''.join(matrix)).split()
            data = []

            if len(matrix) == info.data_length:
                for d in matrix:
                    data.append(transform_string_from_type(d, info.data_type))  # transform data in the right type
            else:
                raise InvalidDataFile('size is not the same as defined for ' + keyword)

            self.chunks_parsed[keyword] = data
        else:
            self.chunks_parsed[keyword] = ''.join([s.strip() for s in self.raw[info.offset_start:info.offset_end]])


TYPE_TO_BINARY_IDENTIFIER = {'I': 0, 'R': 1, 'S': 2}
BINARY_IDENTIFIER_TO_TYPE = dict((b, a) for a, b in TYPE_TO_BINARY_IDENTIFIER.items())
TYPE_TO_BINARY_SIZE = {'I': 8, 'R': 8, 'S': 1}
TYPE_TO_STRUCT_IDENTIFIER = {'I': 'q', 'R': 'd', 'S': 's'}


class BinaryDataFile(DataFile):
    """Binary data file. Same principle, but data stored in binary form: may be a bit faster, and there is no loss
    of information (other than the limit of the representation on computers).

    The file structure is

    + Header
    + Chunks table
    + Chunks

    The different part are detailed below. Note that the size are a consequence of the use of ``struct``
    (see `documentation <https://docs.python.org/3/library/struct.html#format-characters>`_).
    The file is little-endian (use ``<`` of the ``struct`` library).

    .. warning::

        String are encoded in UTF-8 before storage. So the ``data_length`` contains the size of the encoded string
        (so the number of bytes) rather than the number of letters (that ``len()`` gives).

    **File header** (``<iIII``)

    .. list-table::
        :header-rows: 1
        :widths: 30 25 20 25

        * - Name
          - Type
          - Size
          - Default value
        * - Magic number
          - Int
          - 4
          - ``0x17B1DAF1``
        * - Version
          - Unsigned int
          - 4
          - ``1``
        * - Number of chunks
          - Unsigned int
          - 4
          - ``0``
        * - Offset chunks start
          - Unsigned int
          - 4
          - ``16``

    **Chunk table** (For each chunk, ``<IIIII`` and then ``<XXs`` for the keyword, where ``XX`` is the length)

    .. list-table::
        :header-rows: 1
        :widths: 30 25 20 25

        * - Name
          - Type
          - Size
          - Default value
        * - Chunk type
          - Unsigned int
          - 4
          - ``0``
        * - Number of data
          - Unsigned int
          - 4
          - ``0``
        * - Offset start
          - Unsigned int
          - 4
          - ``16``
        * - Offset end
          - Unsigned int
          - 4
          - ``16``
        * - Keyword length
          - Unsigned int
          - 4
          - ``0``
        * - Keyword
          - String
          - 1 x ``length``
          - ``''``

    .. warning::

        Offset are always absolute positions with respect to the **beginning of the chunk section** !!

    For chunk type:

     +  ``0``: Integers (``I``). The long type (``q``) is used to pack these data.
     + ``1``: floats (``R``). The double type (``d``) is used to pack these data.
     + ``2``: string (``S``). The string (``s``) is used to pack these data, after encoding in UTF-8.

    **Chunks** (for each chunk)

    From offset start to offset end, data (size of the chunk is given by
    ``offset_end-offset_size``, and ``(offset_end-offset_size)/len(type)`` should give ``number_of_data``).
    """

    last_version = 1
    magic_number = 0x17B1DAF1

    def __init__(self, version=None):
        super().__init__()
        self.version = version if version is not None else self.last_version
        self.offset_chunk_start = 0

    def read(self, f):
        self.raw = f.read()
        self.chunks_information = {}
        self.chunks_parsed = {}

        # header
        magic_number, version, number_of_chunk, self.offset_chunk_start = struct.unpack('<iIII', self.raw[:16])

        if magic_number != self.magic_number:
            raise InvalidDataFile('wrong magic number 0x{:04X} != 0x{}'.format(magic_number, self.magic_number))
        if version > self.last_version:
            raise InvalidDataFile('Version {} unknown, should be <= {}'.format(version, self.last_version))

        # chunks:
        offset = 16
        for i in range(number_of_chunk):
            chunk_type, number_of_data, offset_start, offset_end, keyword_size = \
                struct.unpack('<IIIII', self.raw[offset:offset + 20])

            chunk_keyword = struct.unpack(
                '<{}s'.format(keyword_size), self.raw[offset + 20: offset + 20 + keyword_size])[0].decode('utf-8')

            if chunk_type not in BINARY_IDENTIFIER_TO_TYPE:
                raise InvalidDataFile('Invalid type {} near offset {} in chunks table'.format(chunk_type, offset))

            if offset_end - offset_start != number_of_data * TYPE_TO_BINARY_SIZE[BINARY_IDENTIFIER_TO_TYPE[chunk_type]]:
                raise InvalidDataFile(
                    'number of data does not match offset for {} in chunk table'.format(chunk_keyword))

            self.chunks_information[chunk_keyword] = ChunkInformation(
                chunk_keyword, BINARY_IDENTIFIER_TO_TYPE[chunk_type], number_of_data, offset_start, offset_end)

            offset += 20 + keyword_size

        if offset != self.offset_chunk_start:
            raise InvalidDataFile('{} informations in chunk table'.format(
                'missing' if offset < self.offset_chunk_start else 'extra'))

    def parse_chunk(self, keyword):
        if keyword not in self.chunks_information:
            raise KeyError(keyword)

        info = self.chunks_information[keyword]

        if info.data_type == 'S':
            self.chunks_parsed[keyword] = struct.unpack(
                '<{}s'.format(info.data_length),
                self.raw[info.offset_start + self.offset_chunk_start:info.offset_end + self.offset_chunk_start])[0]\
                .decode('utf-8')
        else:
            self.chunks_parsed[keyword] = list(struct.unpack(
                '<{}{}'.format(info.data_length, TYPE_TO_STRUCT_IDENTIFIER[info.data_type]),
                self.raw[info.offset_start + self.offset_chunk_start:info.offset_end + self.offset_chunk_start]))

    def set(self, key, data_type, data, force_change_type=False):
        """Override so that the length of a string correspond to the number of bytes in the encoded form"""

        r = super().set(key, data_type, data, force_change_type=force_change_type)

        if data_type == 'S':
            info = self.chunks_information[key]
            info.data_length = len(self.get(key).encode('utf-8'))

        return r

    def write(self, f):
        # first, create chunk table (to know its size) and chunks:
        chunk_table = b''
        chunks = b''
        offset = 0

        for k in self.chunks_information:
            info = self.chunks_information[k]

            offset_end = offset + info.data_length * TYPE_TO_BINARY_SIZE[info.data_type]
            keyword_in_binary = info.keyword.encode('utf-8')
            chunk_table += struct.pack(
                '<IIIII',
                TYPE_TO_BINARY_IDENTIFIER[info.data_type],
                info.data_length,
                offset,
                offset_end,
                len(keyword_in_binary))

            chunk_table += struct.pack('< {}s'.format(len(keyword_in_binary)), keyword_in_binary)

            if info.modified:  # avoid interpreting data if not needed
                val = self.get(k)
                if info.data_type == 'S':
                    encoded_string = val.encode('utf-8')
                    chunks += struct.pack(
                        '<{}{}'.format(info.data_length, TYPE_TO_STRUCT_IDENTIFIER[info.data_type]), encoded_string)
                else:
                    chunks += struct.pack(
                        '<{}{}'.format(info.data_length, TYPE_TO_STRUCT_IDENTIFIER[info.data_type]), *val)
            else:
                chunks += self.raw[
                    self.offset_chunk_start + info.offset_start:self.offset_chunk_start + info.offset_end]

            offset = offset_end

        # write:
        # header is 16 bytes long!
        f.write(struct.pack(
            '<iIII', self.magic_number, self.version, len(self.chunks_information), len(chunk_table) + 16))

        f.write(chunk_table)
        f.write(chunks)
