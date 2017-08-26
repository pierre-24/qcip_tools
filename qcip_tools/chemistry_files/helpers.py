import inspect
from qcip_tools import chemistry_files
# Note: we need to explicitly import all modules, so that they are inspected!
from qcip_tools.chemistry_files import gamess, gaussian, dalton, xyz  # noqa


class ProbablyNotAChemistryFile(Exception):
    """Raised when the type of the file is not recognized"""
    pass


def identifiable_chemistry_file_objects():
    """Yield all objects with the ``WithIdentificationMixin``."""

    for module_name, module in inspect.getmembers(chemistry_files, inspect.ismodule):
        if module_name == 'helpers':
            continue

        for class_name, obj in inspect.getmembers(module, inspect.isclass):
            if issubclass(obj, chemistry_files.WithIdentificationMixin) and \
                    obj != chemistry_files.WithIdentificationMixin:
                yield obj


def open_chemistry_file(f, must_be=None, trust_extension=False):
    """Try to recognise a file in the different possibility, and return the (normally) correct class.

    .. note::

        Despite the fact that it uses the file extension to pre-sort the possibilities, the performances are bad.

    :param f: file
    :type f: file
    :param must_be: restrict the list of possibilities
    :type must_be: list
    :param trust_extension: do not run a second pass on every other possibility if the content does not match extension
    :type trust_extension: bool
    :rtype: qcip_tools.chemistry_files.ChemistryFile
    """

    def attempt_recognition(objs, exclusion=None):
        if not exclusion:
            exclusion = []

        for obj_ in objs:
            if obj_ in exclusion:
                continue

            f.seek(0, 0)

            try:
                if obj_.attempt_identification(f):
                    f.seek(0, 0)
                    o = obj_()
                    o.read(f)
                    return o
            except NotImplementedError:
                continue

        return None

    possible_matches = []
    file_extension = None

    if not must_be:
        must_be = identifiable_chemistry_file_objects()

    # first try is based on file extension:
    try:
        file_extension = f.name.split('.')[-1].lower()
    except:
        pass

    if file_extension:
        for obj in identifiable_chemistry_file_objects():
            if obj not in must_be:
                continue

            try:
                if file_extension in obj.possible_file_extensions():
                    possible_matches.append(obj)
            except NotImplementedError:
                continue

    if possible_matches:
        instance = attempt_recognition(possible_matches)

        if instance:
            return instance

    # then, if nothing matches, try again on every other possibilities
    if not trust_extension:
        instance = attempt_recognition(must_be, exclusion=possible_matches)

        if instance:
            return instance

    # at this point, we should raise an error if nothing matches
    raise ProbablyNotAChemistryFile(f)
