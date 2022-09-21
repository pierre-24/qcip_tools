"""
Scripts to ease the manipulation of quantum chemistry results in Python 3
"""

import pathlib

# List of the scripts that are installed, without the .py extension. The name is used to give the command name.
# Please keep this list ordered alphabetically.
provide_scripts = [
    'check_chemistry_file',
    'to_xyz'
]


def make_console_scripts(package_name=pathlib.Path('qcip_tools/scripts')):
    """Function used to generate the ``console_scripts`` list for ``setup.py``

    :rtype: list
    """

    console_scripts = []
    for script in provide_scripts:
        path = package_name / (script + '.py')
        if not path.is_file():
            raise FileNotFoundError(path)

        console_scripts.append('{0} = {1}.{0}:main'.format(script, str(package_name).replace('/', '.')))

    return console_scripts
