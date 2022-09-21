"""
Scripts to ease the manipulation of quantum chemistry results in Python 3
"""

import pathlib

# List of the scripts that are installed, without the .py extension. The name is used to give the command name.
# Please keep this list ordered alphabetically.
provide_scripts = [
    'boltzmann_population',
    'check_chemistry_file',
#    'ct_analysis',  # noqa
#    'cube_radial_distribution',  # noqa
#    'electrical_derivatives',  # noqa
#    'excitations',  # noqa
#    'geometrical_derivatives',  # noqa
#    'gen_character_table',  # noqa
#    'gen_spectrum',  # noqa
#    'measure_mols',  # noqa
#    'symmetrise',  # noqa
#    'thermochemistry_analysis',  # noqa
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
