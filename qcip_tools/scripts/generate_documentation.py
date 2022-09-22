#!/usr/bin/env python3
"""
Make the documentation of the different scripts, based on their content.
"""

import argparse
import os
import glob
import sys
import hashlib
import datetime

try:
    from importlib.machinery import SourceFileLoader
except ImportError:
    pass

try:
    import importlib.util
except ImportError:
    pass

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'
__longdoc__ = """
For each ``*.py`` is source directory, load it, looks for the autodoc variables (given below), then create a ``.rst``
file out of that in target directory.

+ ``__version__`` ;
+ ``__author___`` (comma separated list) ;
+ ``__maintainer__`` (which should match one of the author) ;
+ ``__email__`` ;
+ ``__status__`` (either "Prototype", "Development", or "Production") ;
+ ``__doc__`` ;
+ ``__longdoc__`` (see below).

``__longdoc__`` is not a real autodoc variable, but is added to generate the "rest" of the documentation (it will
go after the list of options in the documentation page, under the section "more information").
This variable should therefore contains reStructuredText code, which will later be processed by Sphinx.

Since the code load every script to fetch the variables, don't forget to "prevent" the execution with

.. code-block:: python

    if __name__ == '__main__':
        # (your program)

... and please avoid any computation outside this block.

The ``ArgumentParser`` in every script will be looked by name, then by function name. It will be used (if found) to
generate usage and list of option.

In the generated ``.rst`` (which have the same name as the python script), the first 3 lines looks like:

.. code-block:: text

    .. hash=d39b55382d773cd2bbe27ef23f58e431fe03c81c
    .. Generated: 25/08/17 17:36
    .. Do not edit!

The ``hash=`` line is used to check whether the documentation should be re-generated or not (the hash is based on
the script file), so deleting this line will force the generation.

Inspired by the code of
`autoprogram <https://bitbucket.org/birkenfeld/sphinx-contrib/src/7a7a570b8d07/autoprogram/?at=default>`_ and
the `sphinx documentation on domains <http://www.sphinx-doc.org/en/stable/domains.html#the-standard-domain>`_.
"""


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def hash_file(file_path):
    hash = hashlib.sha1()
    with open(file_path, 'rb') as f:
        while True:
            data = f.read(2**16)
            if not data:
                break
            hash.update(data)

    return hash.hexdigest()


def scan_arg_parser(arg_parser, fetch=''):
    """Go through the arg parser and yield options

    :param arg_parser: the arg parser
    :type arg_parser: argparse.ArgumentParser
    :param fetch: fetch type of argument to fetch (positional, required, optional)
    :type fetch: str
    :yield: tuple
    """

    if not fetch:
        fetch = 'positional'
    elif fetch not in ['positional', 'required', 'optional']:
        raise ValueError('fetch')

    for arg in arg_parser._actions:
        if arg.option_strings and fetch != 'positional':
            names = list(arg.option_strings)
            desc = arg.help or ''
            if (arg.required and fetch == 'required') or (not arg.required and fetch == 'optional'):
                yield (names, desc)
        elif not arg.option_strings and fetch == 'positional':
            names = [arg.metavar or arg.dest]
            desc = arg.help or ''
            yield (names, desc)


def get_arguments_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument(
        '-t', '--target',
        action='store',
        default='./docs/scripts/',
        help='Target directory',
        type=is_dir)
    parser.add_argument(
        '-s', '--source', action='store', default='./qcip_tools/scripts', help='Source directory', type=is_dir)
    parser.add_argument(
        '--suffix', action='store', default='', help='suffix for the command (i.e. ".py")')
    parser.add_argument(
        '-f', '--force', action='store_true', help='Skip hash check (force making)')
    parser.add_argument(
        '-d', '--dry-run', action='store_true', help='Check if file need change, but do not modify them')

    return parser


exclude_scripts = ['__init__.py', 'generate_documentation.py']
possible_argparse_parsers = ['arguments_parser', 'arg_parser', 'parser']
possible_argparse_parsers_func = ['get_arguments_parser', 'get_arg_parser']


def main():
    args = get_arguments_parser().parse_args()

    for script_path in glob.glob('{}/*.py'.format(args.source)):

        # find module name:
        script_name = os.path.basename(script_path)
        if script_name in exclude_scripts:
            continue
        module_name = script_name[:-3]

        # load module: (credits to https://stackoverflow.com/a/67692)
        if sys.version_info <= (3, 2):
            raise Exception('too old')
        elif sys.version_info < (3, 5):
            module_obj = SourceFileLoader(module_name, script_path).load_module()
        else:
            spec = importlib.util.spec_from_file_location(module_name, script_path)
            module_obj = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module_obj)

        # check whether we can skip the generation of the documentation for this file
        module_documentation_path = os.path.join(args.target, module_name + '.rst')
        skip_remake_documentation = False
        hash_source = hash_file(script_path)

        if not args.force:
            if os.path.isfile(module_documentation_path):
                with open(module_documentation_path) as f:
                    first_line = f.readline()
                    if first_line[:8] == '.. hash=':
                        hash_source_in_documentation = first_line[8:].strip()
                        if hash_source == hash_source_in_documentation:
                            skip_remake_documentation = True

        # do it!
        if not skip_remake_documentation:
            print(script_path)
        else:
            continue

        if args.dry_run:
            continue

        with open(module_documentation_path, 'w') as f:
            f.write('.. hash={}\n'.format(hash_source))
            f.write('.. Generated: {}\n'.format(datetime.datetime.strftime(datetime.datetime.now(), '%d/%m/%y %H:%M')))
            f.write('.. Do not edit!\n\n')

            module_name_in_code = '``{}{}``'.format(module_name, '{}'.format(args.suffix) if args.suffix else '')
            f.write('{}\n{}\n{}\n\n'.format(
                '=' * len(module_name_in_code), module_name_in_code, '=' * len(module_name_in_code)))

            # use autodoc variables:
            autodoc_variables = {
                '__author__': 'Unknown',
                '__version__': 'unknown',
                '__maintainer__': 'Unknown',
                '__email__': '',
                '__status__': 'Prototype',
                '__doc__': '',
                '__longdoc__': ''
            }

            for var in autodoc_variables:
                if var not in module_obj.__dict__:
                    print('   WARNING: no {}'.format(var))
                else:
                    autodoc_variables[var] = getattr(module_obj, var)

            f.write('By ')
            num_author = 0
            for author in autodoc_variables['__author__'].split(','):
                if num_author != 0:
                    f.write(', ')
                author = author.strip()

                is_maintainer = False
                if author == autodoc_variables['__maintainer__']:
                    is_maintainer = True

                f.write('{}{}{}'.format('**' if is_maintainer else '', author, '**' if is_maintainer else ''))

                if is_maintainer:
                    f.write('{}'.format(
                        ' (`{0} <{0}>`_)'.format(
                            autodoc_variables['__email__']) if autodoc_variables['__email__'] != '' else ''))
                num_author += 1
            f.write('.\n\n')

            f.write('**Version {}** ({}).\n\n'.format(
                autodoc_variables['__version__'], autodoc_variables['__status__']))

            f.write('Synopsis\n++++++++\n\n')

            f.write('{} - {}' .format(module_name_in_code, autodoc_variables['__doc__']))
            f.write('\n\n')

            # find and use argparser
            arg_parser = None
            for var in possible_argparse_parsers:
                if var in module_obj.__dict__:
                    if type(getattr(module_obj, var)) == argparse.ArgumentParser:
                        arg_parser = getattr(module_obj, var)
                        break

            if arg_parser is None:
                for func in possible_argparse_parsers_func:
                    if func in module_obj.__dict__:
                        o = getattr(module_obj, func)()
                        if type(o) == argparse.ArgumentParser:
                            arg_parser = o
                            break

            if arg_parser is not None:
                arg_parser.prog = '{}{}'.format(module_name, '{}'.format(args.suffix) if args.suffix else '')
                # usage:
                f.write('.. program:: {}\n\n'.format(module_name))
                f.write('.. code-block:: console\n\n')
                for line in arg_parser.format_usage().splitlines():
                    f.write('  {}\n'.format(line))
                f.write('\n\n')

                def format_option(names, desc):
                    r = '.. option:: {}\n\n'.format(', '.join(names))
                    if desc:
                        for line in desc.splitlines():
                            r += '  {}\n'.format(line)
                    return r + '\n'

                # arguments
                positional_arguments = ''
                for (names, desc) in scan_arg_parser(arg_parser, fetch='positional'):
                    positional_arguments += format_option(names, desc)

                if positional_arguments:
                    f.write('Positional arguments:\n\n')
                    f.write(positional_arguments)

                required_arguments = ''
                for (names, desc) in scan_arg_parser(arg_parser, fetch='required'):
                    required_arguments += format_option(names, desc)

                if required_arguments:
                    f.write('Required arguments:\n\n')
                    f.write(required_arguments)

                optional_arguments = ''
                for (names, desc) in scan_arg_parser(arg_parser, fetch='optional'):
                    optional_arguments += format_option(names, desc)

                if optional_arguments:
                    f.write('Optional arguments:\n\n')
                    f.write(optional_arguments)

            f.write('\n\n')

            if autodoc_variables['__longdoc__']:
                f.write('More information\n++++++++++++++++\n\n')
                f.write(autodoc_variables['__longdoc__'])


if __name__ == '__main__':
    main()
