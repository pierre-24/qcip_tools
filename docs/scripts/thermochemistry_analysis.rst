.. hash=5bdc9e07a785ff348129bf345393b1f9759f0fed
.. Generated: 21/09/22 22:41
.. Do not edit!

============================
``thermochemistry_analysis``
============================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``thermochemistry_analysis`` - 
Perform a thermochemistry analysis


.. program:: thermochemistry_analysis

.. code-block:: console

  usage: thermochemistry_analysis [-h] [-v] [-n SYMMETRY_NUMBER]
                                  [-t TEMPERATURE] [-p PRESSURE] [-s SCALE]
                                  [-x EXCLUDE] [-V] [-f FACTOR] [-g]
                                  [infile]


Positional arguments:

.. option:: infile

  source of the derivatives

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -n, --symmetry-number

  Symmetry number (for the rotation)

.. option:: -t, --temperature

  Temperature (in K)

.. option:: -p, --pressure

  Pressure (in Pa)

.. option:: -s, --scale

  Scaling factor for the vibrational frequencies

.. option:: -x, --exclude

  Exclude frequencies

.. option:: -V, --verbose

  Gives the detail of the different contributions

.. option:: -f, --factor

  multiply energies by a given factor

.. option:: -g, --guess-symmetry

  Guess the symmetry number



More information
++++++++++++++++


Try to fetch the hessian, and compute thermochemisty data out of that
Rely on the availability of ``geometrical_derivatives`` and ``computed_energies``.
