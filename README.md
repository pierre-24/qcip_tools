# Quantum chemistry in python (QCiP)

## Installation

*Need to find out that part*.

## Installation for development

+ Clone it: git@git.pierrebeaujean.net:pierre/qcip_tools.git`.
+ Install scipy dependencies : ATLAS and Lapack, C and Fortran compilers, Python dev headers and cython : `sudo apt-get install libopenblas-dev libatlas-dev build-essential python-dev swig gfortran python-nose` (or something like this).
+ Create virtualenv and activate it: `virtualenv venv --python=python3 && source venv/bin/activate`.
+ Install dependencies : `pip install --upgrade -r requirements.txt` or `make install-dependencies` (for a developpement environement).
+ Finally, "install" the pakage: `pip install --editable .`
