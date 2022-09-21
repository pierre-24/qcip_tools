.. hash=343dd826fbaea36b20cc66bbac8ee84c28500be0
.. Generated: 21/09/22 17:25
.. Do not edit!

============================
``cube_radial_distribution``
============================

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.1** (Development).

Synopsis
++++++++

``cube_radial_distribution`` - 
Radial distribution


.. program:: cube_radial_distribution

.. code-block:: console

  usage: cube_radial_distribution [-h] [-v] [-S] [-c CENTER] [--dr DR]
                                  [--n-polar N_POLAR]
                                  [--n-azimuthal N_AZIMUTHAL] [-d DATA]
                                  infile


Positional arguments:

.. option:: infile

  Density cube

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -S, --square

  square cube before

.. option:: -c, --center

  Center

.. option:: --dr

  Increment of radius

.. option:: --n-polar

  Number of subdivision for the integration over theta

.. option:: --n-azimuthal

  Number of subdivision for the integration over phi

.. option:: -d, --data

  Center



More information
++++++++++++++++


Report the radial distribution of a cube around a given center [by default :math:`(0,0,0)`].

The charge in a given region of the space, located by :math:`\mathbf{r}` and in an element of volume
:math:`d\mathbf{r}`, is given by

.. math::

    q(\mathbf{r}) = \rho(\mathbf{r})\,d\mathbf{r}.

Integration over whole space gives the number of particles, :math:`Q`.
In spherical coordinates, :math:`d\mathbf{r} = r^2\sin{\theta}\,dr\,d\theta\,d\phi`,
this integral becomes

.. math::

    Q =
    \int_0^{2\pi}\int_0^{\pi}\int_0^{\infty}
    \rho(r,\theta,\phi)\,r^2\,\sin{\theta}\,dr\,d\theta\,d\phi.

Thus, the radial distribution is given by

.. math::
    :label: tr

    \frac{dQ(r)}{dr} =
    r^2\,\int_0^{2\pi}\int_0^{\pi}
    \rho(r,\theta,\phi)\sin{\theta}\,d\theta\,d\phi,

Equation :eq:`tr` is obtained numerically by interpolation over the cube.
