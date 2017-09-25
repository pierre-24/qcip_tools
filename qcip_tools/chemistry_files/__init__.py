from qcip_tools.mixins import Dispatcher


class PropertyNotPresent(Exception):
    """Raised when a property is actually not in the file"""
    pass


class ChemistryFile(Dispatcher):
    """Purely abstract class that implement some basic methods that any child should implement if possible.

    Input methods:

    - ``read()``
    - ``property()`` (not yet implemented)

    Methods to help file recognition:

    - ``possible_file_extensions()`` (class method)
    - ``attempt_identification()`` (class method)

    """

    file_type = 'MISC'
    from_read = False
    use_binary = False

    def read(self, f):
        raise NotImplementedError

    @classmethod
    def define_property(cls, key):
        def wrapper(callback):
            return cls.set_callback(key, callback)
        return wrapper

    def property(self, property_, **kwargs):
        """

        :param property_: the property
        :type property_: str
        """

        callback = self.dispatch(property_)

        if callback is None:
            raise Exception('property {} not defined for {}'.format(property_, type(self)))

        return callback(self, **kwargs)

    def has_property(self, property_):
        """Checks whether a given property is available or not.

        :rtype: bool
        """

        return property_ in self.dispatcher


@ChemistryFile.define_property('file_type')
def chemistry_file__property__file_type(obj, **kwargs):
    """Get the file type trough ``file_type``"""
    return obj.file_type


@ChemistryFile.define_property('molecule')
def chemistry_file__property__molecule(obj, **kwargs):
    """Get the file geometry"""
    if issubclass(type(obj), WithMoleculeMixin):
        return obj.molecule
    else:
        raise PropertyNotPresent('geometry')


class WithIdentificationMixin(object):
    """Mixin for recogintion of the file by the *helpers*"""

    @classmethod
    def possible_file_extensions(cls):
        """Return the common extention of this kind of files

        :rtype: list
        """

        raise NotImplementedError

    @classmethod
    def attempt_identification(cls, f):
        """Attempt to identify the file as a possible file of this type

        :param f: file (in read mode)
        :type f: file
        :rtype: bool
        """

        raise NotImplementedError


class FormatError(Exception):
    """Raised when the format is different from what expected"""
    pass


class WithOutputMixin(object):
    """Mixin to add output methods.


    Output methods:

    - ``to_string()`` (preferable to override)
    - ``write()``
    """

    def to_string(self):
        raise NotImplementedError

    def write(self, f):
        f.write(self.to_string())

    def __repr__(self):
        return self.to_string()


class WithMoleculeMixin(object):
    """Mixin to add a molecule.
    """

    molecule = None

    def get_molecule(self):
        """Get the corresponding molecular geometry.

        :rtype: qcip_tools.molecule.Molecule
        """

        return self.molecule

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        """Create an object out of a molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :rtype: WithMoleculeMixin
        """

        obj = cls()
        obj.molecule = molecule
        obj.from_read = False
        return obj


def apply_over_list(lst, func, start=0, end=None, **kwargs):
    """
    Apply ``func()`` to a given element in list, and expects ``True`` if the iteration must stop, ``False`` otherwise.
    The prototype of ``func`` must be ``func(line,  current_index, **kwargs)``.

    :param lst: list over which the function is applied
    :type lst: list
    :param func: function, for which the first parameter is the index, and the followings are the ``**kwargs``.
    :type func: callback
    :param start: starting index
    :type start: int
    :param end: end the search at some point
    :type end: int
    :type kwargs: dict
    :return: True if it ``func()`` call for termination, False otherwise
    :rtype: bool
    """

    if end is None:
        end = len(lst)

    for index, line in enumerate(lst[start:end]):

        if func(line, current_index=index + start, **kwargs):
            return True

    return False


class ChemistryLogFile(ChemistryFile):
    """Mixin for log files, for which keeps lines in ``self.line``, allow the division of the file in chunks,
     and define search functions
     """

    lines = None
    chunks = None
    chunk_title_variable = 'chunk_title'

    @classmethod
    def possible_file_extensions(cls):
        """.log and .out are so common ...

        :rtype: list
        """

        return ['log', 'out']

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        raise NotImplementedError('from_molecule')

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.lines = f.readlines()
        self.from_read = True
        self.chunks = []

    def apply_function(self, func, line_start=0, line_end=None, into=None, **kwargs):
        """Apply ``apply_over_list()`` to the lines.

        :param lines: lines of the log
        :type lines: list
        :param func: function, for which the first parameter is the line, and the followings are the ``**kwargs``.
        :type func: callback
        :param line_start: starting index
        :type line_start: int
        :param line_end: end the search at some point
        :type line_end: int
        :param into: restrict the search to a given chunk
        :type into: str
        :type kwargs: dict
        :rtype: bool
        """

        if into is not None:
            for section_info in self.chunks:
                if getattr(section_info, self.chunk_title_variable) == into:
                    r = apply_over_list(self.lines, func, section_info.line_start, section_info.line_end, **kwargs)
                    if r:
                        return True

            return False

        else:
            return apply_over_list(self.lines, func, line_start, line_end, **kwargs)

    def search(self, s, line_start=0, line_end=None, into=None):
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

        self.apply_function(find_string, line_start=line_start, line_end=line_end, into=into)
        return FunctionScope.line_found

    def chunk_exists(self, title, all_times=False):
        """Is a chunk in the file?

        :param title: the chunk to search
        :type title: int|str
        :param all_times: rater than stopping at the first time the link is found, go until the end and count
        :type all_times: bool
        :return: the number of times the link is found (max 1 if ``all_times`` is ``False``)
        :rtype: int
        """
        times = 0

        for chunk_info in self.chunks:
            if getattr(chunk_info, self.chunk_title_variable) == title:
                if not all_times:
                    return 1
                else:
                    times += 1

        return times
