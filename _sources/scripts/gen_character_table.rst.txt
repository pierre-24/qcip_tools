.. hash=7babf638df201d609dd53311e4a945cc77f53935
.. Generated: 21/09/22 18:09
.. Do not edit!

=======================
``gen_character_table``
=======================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``gen_character_table`` - 
Generate character table for a given group


.. program:: gen_character_table

.. code-block:: console

  usage: gen_character_table [-h] [-v] group


Positional arguments:

.. option:: group

  group for which character table must be printed

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit



More information
++++++++++++++++


Generate the character table for a group.

Input is, for example, ``C2v`` (or ``C_2v``), ``D3d`` (or ``D_3d``), ``Td`` (or ``T_d``), ...

.. warning::

    + No warranty is given. This has been tested, but not on every group ;
    + The order of the class is (more than probably) not correct ;
    + Complex conjugate representations are not merged together, but are marked with a star.
    + Large groups may take time to generate.
    + Does not work for :math:`I_h` !
