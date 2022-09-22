.. hash=e948d7d2a5d46805d9ae16c0a69faf582c261cce
.. Generated: 24/06/21 16:35
.. Do not edit!

========================
``boltzmann_population``
========================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``boltzmann_population`` - 
Compute the Boltzmann population of molecules


.. program:: boltzmann_population

.. code-block:: console

  usage: boltzmann_population [-h] [-v] [-t TEMPERATURE] [-p PRESSURE] [-c {E,H,G,U}] [-f FACTOR] [infiles ...]


Positional arguments:

.. option:: infiles

  sources

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -t, --temperature

  Temperature (in K)

.. option:: -p, --pressure

  Pressure (in Pa)

.. option:: -c, --criterion

  Criterion for the population

.. option:: -f, --factor

  multiply energies by a given factor



More information
++++++++++++++++


Compute the Boltzmann population of molecules, based on:

+ their total energy (``-c E``, default),
+ their internal energy (``-c U``),
+ their enthalpy (``-c H``),
+ their free Gibbs energy (``-c G``).

Checks if the molecular formula of all inputs matches
(but that's it, e.g. not that they have been obtained with the same method).
