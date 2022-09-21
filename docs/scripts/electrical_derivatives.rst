.. hash=5dcf2ab22db36814212599c78974732d065c7a5a
.. Generated: 21/09/22 22:46
.. Do not edit!

==========================
``electrical_derivatives``
==========================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.3** (Development).

Synopsis
++++++++

``electrical_derivatives`` - 
Fetch the different derivatives of the energy with respect to electric field, and compute related quantities


.. program:: electrical_derivatives

.. code-block:: console

  usage: electrical_derivatives [-h] [-v] [infile]


Positional arguments:

.. option:: infile

  source of the derivatives

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit



More information
++++++++++++++++


Try to fetch the dipole moment and the frequency dependant (hyper)polarizabilities, then print the tensors, and
related quantities. Rely on the availability of ``electrical_derivatives``.

.. warning::

    In Dalton, by default, second hyperpolarizability with HF or DFT does not compute
    all components of the gamma tensor, but only the one that contribute to :math:`\gamma_{||}`.
    Use ``.GAMALL`` in ``*CUBIC`` to do so.
