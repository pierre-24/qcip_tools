====================================================================
Derivatives w.r.t. geometrical stuffs (``qcip_tools.derivatives_g``)
====================================================================

Tools to help the manipulation of derivatives of the energy with respect to cartesian (``G``) or normal mode (``N``) coordinates.

Theory
------

The displacement around the equilibrium position, :math:`\mathbf{x}_i=\mathbf{x}'_i-\mathbf{x}_{i,0}`, may be expressed as a Taylor series to describe the potential energy surface:

  .. math::

      V(\mathbf{x}) = V_0  +  \sum_i F_i\mathbf{x}_i + \frac{1}{2!}\,\sum_{ij}H_{ij}\mathbf{x}_i\mathbf{x}_j +  + \frac{1}{3!}\,\sum_{ij}F_{ijk}\mathbf{x}_i\mathbf{x}_j\mathbf{x}_k+ \ldots

where :math:`\mathbf{x}\in\mathbf{R}^{3N}`, :math:`F_i=\left.\frac{\partial V}{\partial\mathbf{x}_i}\right|_{\mathbf{x}=0}` is the gradient (so the forces, wich are close to 0 if the system is in equilibrium), :math:`H_{ij}=\left.\frac{\partial^2 V}{\partial\mathbf{x}_i\mathbf{x}_j}\right|_{\mathbf{x}=0}` is the hessian (force constant matrix) and :math:`F_{ijk}` is the cubic force constants matrix (and so on).
The sums runs over the :math:`3N` degrees of freedom of :math:`N` atoms.

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

    \mathbf{M}\,\overset{\cdot\cdot}{\mathbf{X}} + \mathbf{H}\,\mathbf{X} = 0.

One can assume, in the case of small displacements, a linear harmonic solution, :math:`\mathbf{X}=\mathbf{X}_0\,\cos{(\mathbf{\omega}\,t)}`, so that the previous equation reduce to

.. math::

    \mathbf{H}\,\mathbf{X}_0 = \mathbf{\omega}^2\,\mathbf{m}\,\mathbf{X}_0.

This is reduced to an eigenvalue problem by setting mass-weighted coordinates, :math:`\mathbf{Q}=\mathbf{m}^{\frac{1}{2}}\,\mathbf{X}_0`, which gives

.. math::

    (\mathbf{m}^{-\frac{1}{2}}\mathbf{H}\,\mathbf{m}^{\frac{1}{2}})\,\mathbf{Q} = \mathbf{\omega}^2\,\mathbf{Q}.

The :math:`\mathbf{Q}`'s are therefore eigenvectors of the mass-weigted hessian, also called normal modes (which form a complete and orthogonal basis), while the :math:`\mathbf{\omega}`'s are the eigenvalues.
There is 6 (or 5, if the system is linear) zeros eigenvalues for the translation and rotation modes, while the other normal modes describe the vibrations.


.. note::

    *May be continued to explain the concept of displacements and the corresponding units*.

API documentation
-----------------

.. automodule:: qcip_tools.derivatives_g
    :members: