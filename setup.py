from setuptools import setup
from pip.download import PipSession
from pip.req import parse_requirements

session = PipSession()
pkgs = []
for pkg in parse_requirements('requirements.txt', session=session):
    if pkg.req:
        pkgs.append(str(pkg.req))

setup(
    name='qcip_tools',
    packages=['qcip_tools', 'qcip_tools.chemistry_files'],
    version='0.2',
    author='Pierre Beaujean',
    author_email='pierre.beaujean@unamur.be',
    description='Library to ease the manipulation of quantum chemistry results in Python 3',
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'
    ],
    install_requires=pkgs,
    test_suite='tests'
)
