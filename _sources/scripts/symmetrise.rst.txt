.. hash=d9c35f7c887b713be2df5f60058f5b9fd86823db
.. Generated: 21/09/22 18:48
.. Do not edit!

==============
``symmetrise``
==============

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``symmetrise`` - 
Detect point group, symmetrise molecule and outputs it as XYZ


.. program:: symmetrise

.. code-block:: console

  usage: symmetrise [-h] [-v] [-t TOLERANCE] [-s] [-u] [-S] [infile]


Positional arguments:

.. option:: infile

  source of the derivatives

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -t, --tolerance

  Tolerance threshold

.. option:: -s, --symmetrise

  Fully symmetrise the molecule (involve small displacements and the full generation of the group)

.. option:: -u, --uniques

  Get unique atoms (implies `-s`)

.. option:: -S, --same-output

  Output in the same format as the input



More information
++++++++++++++++


Detect point group and symmetrise molecule.

The tolerance option (``--tolerance``) depends of the quality of an eventual optimisation: ``1e-3`` (the default) is ok
in most cases, but may results in a larger group than it should in a few cases, so lower this threshold if any.

By default, the script only orient the molecule w.r.t. symmetry elements, so that it matches its group.
But ``--symmetrise`` and/or ``--uniques`` implies to generate the whole group, which may take time for
large ones (thus very symmetrical molecules, with :math:`\#\mathcal{G} \geq 60`). Then, it generates uniques atoms,
and eventually the whole molecule back (if only ``--symmetrise`` is set).

.. warning::

    ``--symmetrise`` and ``--uniques`` do not work for icosahedral molecules!!

