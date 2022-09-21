.. hash=101bfe10a5fb1c5e59b1487bfd762e49bc3c5da2
.. Generated: 21/09/22 17:25
.. Do not edit!

========================
``check_chemistry_file``
========================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``check_chemistry_file`` - 
checks whether qcip_tools identify a given file. If so, it can also check if it have a given property.


.. program:: check_chemistry_file

.. code-block:: console

  usage: check_chemistry_file [-h] [-v] [-p PROPERTY] [-T] [infile]


Positional arguments:

.. option:: infile

  file to check

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -p, --property

  List of atoms

.. option:: -T, --not-trust-extension

  Do not trust extension of the file (slower)



More information
++++++++++++++++


Prints 1 if it is the case, 0 if not.

.. note::

    Having this property does not mean that the property in question is indeed available in this given file,
    just that it *could* present this property.
