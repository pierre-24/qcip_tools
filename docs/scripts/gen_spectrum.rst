.. hash=08cd64948d5e4bc6a7893b3e95d5631190f686c4
.. Generated: 21/09/22 18:48
.. Do not edit!

================
``gen_spectrum``
================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``gen_spectrum`` - 
Generate a spectrum


.. program:: gen_spectrum

.. code-block:: console

  usage: gen_spectrum [-h] [-v] [-l LIMITS] [-e EACH] [-I] [-L LIMIT_IMPULSES]
                      [-m] [-n | -s SCALE]
                      [-E EXCLUDE | -i INCLUDE | -D DECONTAMINATE]
                      [infile] {uv,ir} ...


Positional arguments:

.. option:: infile

  source of the derivatives

.. option:: subparser_name

  source of the spectrum (`uv` or `ir`)

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -l, --limits

  Limit and units of the graph: `min:max:unit`

.. option:: -e, --each

  Interval (in unit of limits)

.. option:: -I, --impulses

  Add impulses (at the end)

.. option:: -L, --limit-impulses

  Only output impulses above that threshold

.. option:: -m, --maximums

  Add maximums (at the end, after the impulses if any)

.. option:: -n, --normalize

  Set maximum intensity to 1

.. option:: -s, --scale

  Scale the intensities

.. option:: -E, --exclude

  Exclude some peaks (list starts at 0)

.. option:: -i, --include

  Include only some peaks (list starts at 0)

.. option:: -D, --decontaminate

  Exclude peaks for which Δ<Ŝ²> is larger than this threshold (only relevant if <S²> is available)



More information
++++++++++++++++


Generate spectrum the point to plot a spectrum in a given window of energy (or wavelength).

The `y` values are in arbitrary units.

Currently, **only fetch UV intensities** and plot them, for

+ Gaussian (FCHK)
+ Dalton (archive)

To do, with the Hessian:

+ IR spectra (I need the derivatives of the dipole moment) ;
+ Raman spectra (I need the derivatives of the polarizability).

