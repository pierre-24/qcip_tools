.. hash=8be03ef9bb9400b35d2112ead47077efe99fb331
.. Generated: 21/09/22 22:41
.. Do not edit!

===========================
``geometrical_derivatives``
===========================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.2** (Development).

Synopsis
++++++++

``geometrical_derivatives`` - 
Fetch the different derivatives of the energy with respect to geometrical modifications


.. program:: geometrical_derivatives

.. code-block:: console

  usage: geometrical_derivatives [-h] [-v] [-N] [infile]


Positional arguments:

.. option:: infile

  source of the derivatives

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -N, --no-vibs

  Do not perform the vibrational analysis



More information
++++++++++++++++


Try to fetch the gradient and hessian. If the second one is found, compute vibrational frequencies.
Rely on the availability of ``geometrical_derivatives``.
