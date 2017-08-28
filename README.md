# Quantum chemistry in python (QCiP) tools library

Library maintained by [Pierre Beaujean](pierre.beaujean@unamur.be) to ease the manipulation of quantum chemistry results in Python 3. Created in the frame of my PhD thesis in the [University of Namur](https://www.unamur.be).

Stable version: [`release-v0.2.2`](https://git.pierrebeaujean.net/pierre/qcip_tools/tree/release-v0.2.2) (August 25, 2017).

Previous releases:

+ [`release-v0.2.1`](https://git.pierrebeaujean.net/pierre/qcip_tools/tree/release-v0.2.1) (August 25, 2017).
+ [`release-v0.2`](https://git.pierrebeaujean.net/pierre/qcip_tools/tree/release-v0.2) (August, 24 2017).
+ [`release-v0.1`](https://git.pierrebeaujean.net/pierre/qcip_tools/tree/release-v0.1) (July, 11 2017).

## Purpose

+ General purpose objects (molecule, atom, basis set) ;
+ Handle and post-analyze results of calculations (currently, electrical and geometrical derivatives of the energy) ;
+ Retrieve data from quantum chemistry package and create input files (currently Gaussian, Dalton and GAMESS).

## Installation

With pip:

```bash
pip install git+ssh://git@gitlab.unamur.be/pierre.beaujean/qcip_tools.git@master
```

See the [installation page](./documentation/source/install.rst) for more information.

## API documentation

Build the documentation with `make doc`.

## Contributing

You can reports bugs and suggestions any time by email or using the bugtracker.

If you want to contribute to the code, see the [contribution page](./documentation/source/contributing.rst). 
Please note that the code is not actually developed on the git server of the University of Namur (which only contains the releases) but on a personal protected git server (with CI activated and properly configured). 
Feel free to ask access if needed.
