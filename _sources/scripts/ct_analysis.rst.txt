.. hash=14725faf618fc39e28e85d3be4de42695f0a14b3
.. Generated: 21/09/22 17:25
.. Do not edit!

===============
``ct_analysis``
===============

By **Pierre Beaujean** (`pierre.beaujean@unamur.be <pierre.beaujean@unamur.be>`_).

**Version 0.2** (Development).

Synopsis
++++++++

``ct_analysis`` - 
Compute charge transfer (CT) quantities


.. program:: ct_analysis

.. code-block:: console

  usage: ct_analysis [-h] [-v] -g GROUND -e EXCITED [-S] [-o OUTPUT] [-D]
                     [-d DUMMY]


Required arguments:

.. option:: -g, --ground

  Ground state density

.. option:: -e, --excited

  Excited state density

Optional arguments:

.. option:: -h, --help

  show this help message and exit

.. option:: -v, --version

  show program's version number and exit

.. option:: -S, --square

  square cube before make the difference

.. option:: -o, --output

  difference cube

.. option:: -D, --output-with-diff

  Store all the densities (total, positive, negative)

.. option:: -d, --dummy

  Create an xyz file with dummy atoms in positions of both barycenters (negative then positive)



More information
++++++++++++++++


Theory
******

Based on an implementation of this theory by D. Jacquemin (see supporting informations of the corresponding paper).

The charge in a given region of the space, located by :math:`\mathbf{r}` and in an element of volume
:math:`d\mathbf{r}`, is given by

.. math::

    q(\mathbf{r}) = \int  \rho(\mathbf{r})\,d\mathbf{r}.

Charges at different point of the space (a "density") can be obtained by using the
``cubegen`` utility `provided by Gaussian <http://gaussian.com/cubegen/>`_. In particular,

+ ``cubegen 0 density=scf x.fchk out.cub`` permit to obtain the "density" of the ground state, and
+ ``cubegen 0 density=ci x.fchk out.cub`` permit to obtain the "density" of the excited state.

Note that you need to launch the Gaussian calculations with ``Density=(all)`` so that it stores the CI density
(if you use solvatation, please **make sure** to use ``TD=(NonEqSolv)``!).

Let :math:`\delta q(r)` be the difference between the charge of the excited and the ground state.
This quantity can be splitted into increasing areas and decreasing ones, :math:`\delta q_+(\mathbf{r})` and
:math:`\delta q_-(\mathbf{r})`, where

.. math::

    \delta q_+(\mathbf{r}) = \left\{
    \begin{array}{ll}
        \delta q(\mathbf{r}) & \text{if }\delta q(\mathbf{r}) > 0, \\
        0 & \text{otherwise.}
    \end{array}
    \right.

and,

.. math::

    \delta q_-(\mathbf{r}) = \left\{
    \begin{array}{ll}
        \delta q(\mathbf{r}) & \text{if }\delta q(\mathbf{r}) < 0, \\
        0 & \text{otherwise.}
    \end{array}
    \right.

One can therefore compute:

- The transferred charge between ground and excited state:

  .. math::

    q_{CT} = \frac{1}{2}\,\sum_{\mathbf{r}_i} \delta q_+(\mathbf{r}_i) - \delta q_-(\mathbf{r}_i).

  Note that the original papers states that

  .. math::

    q_{CT} =\sum_{\mathbf{r}_i} \delta q_+(\mathbf{r}_i) =-\sum_{\mathbf{r}_i} \delta q_-(\mathbf{r}_i),

  but the implementation of D. Jacquemin reports and uses the average.

- The starting and ending point of the charge transfer, :math:`\mathbf{r}_+` and :math:`\mathbf{r}_-`:

  .. math::

    \mathbf{r}_+ = \sum_{\mathbf{r}_i} \frac{\mathbf{r}_i\,q_+(\mathbf{r}_i)}{q_{CT}},

    \mathbf{r}_- = \sum_{\mathbf{r}_i} \frac{\mathbf{r}_i\,q_-(\mathbf{r}_i)}{q_{CT}}.

  Those are the barycenters of the positive and the negative densities.
  The vector between those two barycenter is the charge transfer vector, defined as

  .. math::

    \mathbf{v}_{CT} = \mathbf{r}_--\mathbf{r}_+.

  In particular, the charge transfer distance is the norm of this vector, :math:`d_{CT} = |\mathbf{v}_{CT}|`.
  Notice the usage of the so called *chemist convention*, where the dipole is defined from positive to negative
  positions.

+ The norm of variation of dipole moment between the ground and excited state:

  .. math::

    |\mu_{CT}| = q_{CT}\,d_{CT}.


Implementation
**************

Only works with gaussian cube.
Note that external programs may be abble to generate those as well
(if this is not a density but a probability, like with MO, squaring it gives the density, so use the ``-S`` option).

The program reports :math:`\mathbf{v}_{CT}` as well as :math:`q_{CT}`, :math:`d_{CT}`, and
:math:`|\mu_{CT}|`.

It allows to save the difference cube (for visualization) and an xyz file containing two dummy atoms (one for each
barycenter, first :math:`\mathbf{r}_-` and then :math:`\mathbf{r}_+`).

.. warning::

    Vector and :math:`d_{CT}` are given in Angstrom, :math:`q_{CT}` is in \|e\| (electron charge),
    and :math:`|\mu_{CT}|` is therefore in Angstrom \|e\|.

Source
******

+ T. Le Bahers *et al.* *J. Chem. Theory. Comput.* **7**, 2498 (2011)
  `10.1021/ct200308m <http://dx.doi.org/10.1021/ct200308m>`_.
+ D. Jacquemin *et al.* *Phys Chem Chem Phys.* **28**, 5383 (2012)
  `10.1039/c2cp40261k <http://dx.doi.org/10.1039/c2cp40261k>`_.
