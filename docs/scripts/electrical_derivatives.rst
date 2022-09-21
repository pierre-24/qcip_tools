.. hash=7d8ab69401afa169bec844428ceab8766d335c7e
.. Generated: 21/09/22 18:09
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

Currently implemented for:

+ Dalton archive output (CC and normal responses).
  Probably not working for explicit operator settings (please use ``.FREQUE`` along with things like ``.SHG`` or so),
  and other type of operators (different from ``DIPLEN``).
+ Gaussian FCHK, with sign correction (but **not the second hyperpolarizability tensors**).

To do:

+ GAMESS output (not easy, and don't forget the DFT one)
+ Gaussian LOG (all located in the same place)
+ Dalton LOG (not located in the same place)

.. warning::

    By default, second hyperpolarizability with HF or DFT does not compute all components of the gamma tensor, but only
    the one that contribute to :math:`\gamma_{||}`.
