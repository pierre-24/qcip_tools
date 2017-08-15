from qcip_tools.mixins import Dispatcher


class InputChemistryFile(Dispatcher):
    """Purely abstract class that implement some basic methods that any child should implement if possible.

    Input methods:

    - ``read()``
    - ``get_molecule()``
    - ``property()`` (not yet implemented)

    Methods to help file recognition:

    - ``possible_file_extensions()`` (class method)
    - ``attempt_recognition()`` (class method)

    """

    molecule = None
    file_type = 'MISC'
    from_read = False

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

    def get_molecule(self):
        """Get the corresponding  molecular geometry. Raises ``NotImplementedError`` if ``self.molecule`` is ``None``.

        :rtype: qcip_tools.molecule.Molecule
        """

        molecule = self.molecule

        if molecule is None:
            raise NotImplementedError

        return molecule

    @classmethod
    def possible_file_extensions(cls):
        """Return the common extention of this kind of files

        :rtype: list
        """

        raise NotImplementedError

    @classmethod
    def attempt_recognition(cls, f):
        """Attempt to identify the file as a possible file of this type

        :param f: file (in read mode)
        :type f: file
        :rtype: bool
        """

        raise NotImplementedError


@InputChemistryFile.define_property('file_type')
def property__file_type(obj, **kwargs):
    """Get the file type trough ``file_type``"""
    return obj.file_type


class InputOutputChemistryFile(InputChemistryFile):
    """Purely abstract class that implement some basic methods that any child should implement if possible.


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
