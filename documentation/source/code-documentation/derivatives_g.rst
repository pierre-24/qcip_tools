====================================================================
Derivatives w.r.t. geometrical stuffs (``qcip_tools.derivatives_g``)
====================================================================

Tools to help the manipulation of derivatives of the energy with respect to cartesian (``G``) or normal mode (``N``) coordinates.

Theory
------

The displacement around the equilibrium position, :math:`\mathbf{x}_i=\mathbf{x}'_i-\mathbf{x}_{i,0}`, may be expressed as a Taylor series to describe the potential energy surface:

  .. math::

      V(\mathbf{x}) = V_0  +  \sum_i F_i\mathbf{x}_i + \frac{1}{2!}\,\sum_{ij}H_{ij}\mathbf{x}_i\mathbf{x}_j + \frac{1}{3!}\,\sum_{ij}F_{ijk}\mathbf{x}_i\mathbf{x}_j\mathbf{x}_k+ \ldots

where :math:`\mathbf{x}\in\mathbf{R}^{3N}`, :math:`F_i=\left.\frac{\partial V}{\partial\mathbf{x}_i}\right|_{\mathbf{x}=0}` is the gradient (so the forces, wich are close to 0 if the system is in equilibrium), :math:`H_{ij}=\left.\frac{\partial^2 V}{\partial\mathbf{x}_i\mathbf{x}_j}\right|_{\mathbf{x}=0}` is the hessian (force constant matrix) and :math:`F_{ijk}` is the cubic force constants matrix (and so on).
The sums runs over the :math:`3N` degrees of freedom of :math:`N` atoms.

Normal modes
============

By settings the Lagrange equations of the harmonic oscillator,

.. math::
    \begin{align}
    &T(\mathbf{x}) = \frac{1}{2} \sum_{ij} m_{ij}\,\dot{x}_i\dot{x}_j, \\
    &V(\mathbf{x}) = \frac{1}{2}\sum_{ij} H_{ij}\,x_i\,x_j,
    \end{align}

and :math:`L = T-V`. Taking advantage of the fact that the Hessian and the mass matrices are diagonals, :math:`H_{ij}=H_{ii}\delta_{ij}` and :math:`m_{ij}=m_i\delta_{ij}`, the :math:`3N` Lagrange equation are written

.. math::

    \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{x}_i}\right) - \left(\frac{\partial L}{\partial x_i}\right) = \sum_j m_{ij}\,\overset{\cdot\cdot}{x}_j + H_{ij}\,x_{j} = 0,

or, in a matrix form,

.. math::

    \mathbf{m}\,\overset{\cdot\cdot}{\mathbf{X}} + \mathbf{H}\,\mathbf{X} = 0.

One can assume, in the case of small displacements, a linear harmonic solution, :math:`\mathbf{X}=\mathbf{X}_0\,\cos{(\mathbf{\omega}\,t)}`, so that the previous equation reduce to

.. math::

    \mathbf{H}\,\mathbf{X}_0 = \mathbf{m}\,\mathbf{\omega}^2\,\mathbf{X}_0.

This is reduced to an eigenvalue problem by setting mass-weighted coordinates, :math:`\mathbf{Q}=\mathbf{m}^{\frac{1}{2}}\,\mathbf{X}_0`, which gives

.. math::

    (\mathbf{m}^{-\frac{1}{2}}\mathbf{H}\,\mathbf{m}^{-\frac{1}{2}})\,\mathbf{Q} = \mathbf{\omega}^2\,\mathbf{Q},

or,

.. math::
    :label: x1

    \mathbf{Q}^\dagger\,\mathbf{H}^m\,\mathbf{Q} = \omega^2

where :math:`\mathbf{H}^m=\mathbf{m}^{-\frac{1}{2}}\mathbf{H}\,\mathbf{m}^{-\frac{1}{2}}` is the mass-weighted Hessian.
The :math:`\mathbf{Q}`'s are therefore eigenvectors of the mass-weigted hessian, also called normal modes (which form a complete and orthogonal basis), while the :math:`\mathbf{\omega}`'s are the eigenvalues (vibrational frequencies).
There is 6 (or 5, if the system is linear) zeros eigenvalues for the translation and rotation modes, while the other normal modes describe the vibrations.


Computation of normal modes
===========================

The computation takes place in the  ``MassWeightedHessian`` class.

.. warning::

    Modes are not decontamined from translations and rotations.

1. Compute the mass weighted Hessian, :math:`\mathbf{H}^m`, as

    .. math::

        H^m_{i\alpha,j\beta} = \frac{1}{\sqrt{m_i\,m_j}}\,H_{i\alpha,j\beta},

  where :math:`i` and :math:`j` are the number of the atom, while :math:`\alpha` and :math:`\beta` ar the cartesian coordinates of the atoms.
  Masses, :math:`m_i` and :math:`m_j`, are in AMU.

2. Compute the :math:`3N` eigenvalues (:math:`\omega^2`) and egeinvectors (:math:`\mathbf{Q}`) of :math:`\mathbf{H}^m`. The explicit set of equations to solve, obtained from :eq:`x1`, is written

    .. math::
        :label: v1

        \omega_p^2\,\delta_{pq} = \sum_{i\alpha,j\beta} Q_{i\alpha,p}\,\frac{H_{i\alpha,j\beta}}{\sqrt{m_i\,m_j}}\,Q_{j\beta,q}.\label{eq:v1}

  Order them from the lowest to the largest eigenvalue.
  Eventually discard the 6 (or 5) first lowest ones.

3. Compute frequencies:

    .. math::

        \nu_p = \sqrt{\frac{\omega^2_p}{\mu}},

  where :math:`\mu` is the conversion factor from AMU to electron mass (about 1822, see `here <quantities.html#qcip_tools.quantities.AMUToElectronMass>`_) and :math:`p` is the normal mode. With the conversion, :math:`\nu_p` is therefore in atomic units.

4. Compute the carthesian displacements. From :eq:`v1`, they are defined by

    .. math::

        D_{p,i\alpha} = \frac{1}{\sqrt{m_i\,\mu}} Q_{i\alpha,p}.

  With the conversion, displacement are therefore in :math:`\text{Å}\,m_e^{-1/2}`, where :math:`m_e` is the electron mass (so atomic unit of mass).

.. note::

    The normal modes (``normal_modes``) and cartesian displacements (``displacements``) are stored in the form :math:`p,i\alpha`.

API documentation
-----------------

.. automodule:: qcip_tools.derivatives_g
    :members: