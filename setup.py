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
    packages=['qcip_tools'],
    version='0.2',
    include_package_data=True,
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'
    ],
    install_requires=pkgs,
    test_suite='tests'
)
