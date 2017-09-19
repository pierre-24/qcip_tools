from setuptools import setup
from pip.download import PipSession
from pip.req import parse_requirements

import qcip_tools

session = PipSession()
pkgs = []
for pkg in parse_requirements('requirements.txt', session=session):
    if pkg.req:
        pkgs.append(str(pkg.req))

setup(
    name=qcip_tools.__name__,
    packages=['qcip_tools', 'qcip_tools.chemistry_files'],
    version=qcip_tools.__version__,
    author=qcip_tools.__author__,
    author_email=qcip_tools.__email__,
    description=qcip_tools.__doc__,
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'
    ],
    install_requires=pkgs,
    test_suite='tests'
)
