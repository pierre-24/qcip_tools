.. hash=0b48528352ad8d9fa9c4c11b93f5147c8d7b0ee9
.. Generated: 21/09/22 18:09
.. Do not edit!

===============
``excitations``
===============

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``excitations`` - 
Extract the excited states information


.. program:: excitations

.. code-block:: console

  usage: excitations [-h] [-v] [-l LIMIT] [-L LIMIT_CSFS] [-E {au,eV,nm,cm-1}]
                     [infile]


Positional arguments:

.. option:: infile

  source of the derivatives

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -l, --limit

  Only print excitation above a given limit

.. option:: -L, --limit-CSFs

  Only print CSF above a limit

.. option:: -E, --unit

  Unit for the energy



More information
++++++++++++++++


You can

+ Change the unit of the energy
+ Limit the number of state that appears
