from setuptools import setup, find_packages
from os import path

import qcip_tools

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

with open(path.join(here, 'requirements/requirements-base.in')) as f:
    requirements = f.readlines()

with open(path.join(here, 'requirements/requirements.in')) as f:
    requirements_dev = f.readlines()[1:]

setup(
    name='qcip_tools',
    version=qcip_tools.__version__,

    # Description
    description=qcip_tools.__doc__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='website',

    project_urls={
        'Bug Reports': 'https://gitlab.unamur.be/chimie/lct/qcip_tools/issues',
        'Source': 'https://gitlab.unamur.be/chimie/lct/qcip_tools',
    },

    url='https://gitlab.unamur.be/chimie/lct/qcip_tools',
    author=qcip_tools.__author__,

    # Classifiers
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',

        # Specify the Python versions:
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    packages=find_packages(),
    python_requires='>=3.5',
    test_suite='tests',

    # requirements
    install_requires=requirements,

    extras_require={  # Optional
        'dev': requirements_dev,
    },
)
