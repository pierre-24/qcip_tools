import math

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
        raise TypeError(data_type)


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
            raise ValueError(data_type)

        self.keyword = keyword
        self.data_type = data_type
        self.data_length = data_length
        self.offset_start = offset_start
        self.offset_end = offset_end
        self.modified = modified


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
        self.chunks_parsed = {}  #: parsed chunks, so list of stuffs
        self.raw = None  #: variable for storage of chunks for ``parse_chunk()``,  so "whatever floats you boat".
        self.data_source = None  #: source pipe, should not be used normally (memory access is faster than the others)

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
            raise TypeError(data_type)

        if not force_change_type:
            if key in self.chunks_information and self.chunks_information[key].data_type != data_type:
                raise TypeError(data_type)

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
                    raise TypeError(data_type)

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
        self.data_source = f

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
