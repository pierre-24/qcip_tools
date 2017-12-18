import inspect
import argparse

from qcip_tools import chemistry_files
# Note: we need to explicitly import all modules, so that they are inspected!
from qcip_tools.chemistry_files import gamess, gaussian, dalton, xyz, chemistry_datafile  # noqa


class ProbablyNotAChemistryFile(Exception):
    """Raised when the type of the file is not recognized"""
    pass


def identifiable_chemistry_file_objects(identifier_must_be=None):
    """Yield all objects with the ``WithIdentificationMixin``.

    :param identifier_must_be: check if ``file_type`` is of the same type as the obj
    :type identifier_must_be: list of str
    """

    for module_name, module in inspect.getmembers(chemistry_files, inspect.ismodule):
        if module_name == 'helpers':
            continue

        for class_name, obj in inspect.getmembers(module, inspect.isclass):
            if issubclass(obj, chemistry_files.WithIdentificationMixin) and \
                    obj != chemistry_files.WithIdentificationMixin:
                if identifier_must_be is not None:
                    if obj.file_type in identifier_must_be:
                        yield obj
                else:
                    yield obj


def open_chemistry_file(f, must_be=None, trust_extension=False):
    """Try to recognise a file in the different possibility, and return the (normally) correct class.

    .. note::

        Despite the fact that it uses the file extension to pre-sort the possibilities, the performances are bad.

    .. note::

        Use the trick from https://stackoverflow.com/a/31123030 to serve a binary file if needed.

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
                f_in_good_mode = f if not obj_.requires_binary_mode else f.buffer.raw
                if obj_.attempt_identification(f_in_good_mode):
                    f_in_good_mode.seek(0, 0)
                    o = obj_()
                    o.read(f_in_good_mode)
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
    except KeyError:
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


def create_open_chemistry_file_action(must_be=None):
    """Create a custom argparse action.

    .. code-block:: python

        from qcip_tools.chemistry_files import helpers, gaussian

        parser = argparse.ArgumentParser()
        parser.add_argument('-a', action=helpers.create_open_chemistry_file_action())
        parser.add_argument('-b', action=helpers.create_open_chemistry_file_action(must_be=[gaussian.Input]))

    String (or list) is parsed, and identifier can be added with ``IDENTIFIER:/path/to/file``.
    Valid identifier are the values of ``file_type`` in the different ``ChemistryFile`` objects
    (a list can be found in the table on the beginning `of this section <chemistry_files.html#api-documentation>`_).

    :param must_be: restrict the list of possibilities
    :type must_be: list
    :rtype: argparse.Action
    """

    class CustomAction(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if type(values) is str:
                setattr(args, self.dest, self.open(values, parser))
            elif type(values) is list:
                files = []
                for value in values:
                    files.append(self.open(value, parser))
                setattr(args, self.dest, files)
            else:
                parser.error('Unknown type ?!?')

        def open(self, s, parser):
            x_must_be = must_be

            if ':' not in s:
                path = s
            else:
                r = s.split(':')
                path = ':'.join(r[1:])
                identifer = r[0]
                try:
                    obj = next(identifiable_chemistry_file_objects(identifier_must_be=[identifer]))
                except StopIteration:
                    parser.error('{}: iddentifier "{}" is not a valid identifier'.format(path, identifer))

                if x_must_be is not None and obj not in x_must_be:
                    parser.error('{}: identifier {} is not a correct possibility (must be {})'.format(
                        path, identifer, ', '.join(x.file_type for x in x_must_be)))
                x_must_be = [obj]

            try:
                with open(path) as f:
                    fx = open_chemistry_file(f, must_be=x_must_be)
                    return fx
            except ProbablyNotAChemistryFile:
                if must_be is not None:
                    parser.error('{}: must be {}'.format(path, ', '.join(x.file_type for x in must_be)))
                elif x_must_be is not None:
                    parser.error('{}: cannot recognize file as {}'.format(path, x_must_be[0].file_type))
                else:
                    parser.error('{}: not recognized as a chemistry file'.format(path))
            except Exception as e:
                parser.error('cannot open {}: {}'.format(s, e))

    return CustomAction
