======================================================
Derivatives of the energy (``qcip_tools.derivatives``)
======================================================

Tools to help the representation of derivatives of the energy.

Theory
------

Supposing a function :math:`\mathbf{f}:\mathbb{R}^{n}\rightarrow \mathbb{R}^{m}`, a function which takes a vector :math:`\mathbf{x}\in\mathbb{R}^n` as input and produce a vector :math:`\mathbf{f}(\mathbf{x})\in\mathbb{R}^m`, one can write its `Maclaurin series expansion <https://en.wikipedia.org/wiki/Taylor_series>`_ as

.. math::
    :label: mclaurin

    \mathbf{M}(\mathbf{q}) = \sum_{r=0}^{r_{max}} \frac{1}{r!}\frac{\partial^r\mathbf{f}(0)}{\partial\mathbf{q}^r}\,\mathbf{q}^r,

where :math:`r_{max}=\infty` if the series is not truncated.

In this expansion, one recognise, for example, :math:`\frac{\partial\mathbf{f}(0)}{\partial\mathbf{q}}`, which is the jacobian of the function (and a linear application :math:`\mathbb{R}^n\rightarrow\mathbb{R}^m`) evaluated at :math:`\mathbf{q}=\mathbf{0}`, which correspond to a :math:`n\times m` matrix :math:`\alpha`, for which the elements :math:`i` and :math:`j` are given by

.. math::

    \alpha_{ij} = \left.\frac{\partial \mathbf{f}_i}{\partial \mathbf{q}_j}\right|_{\mathbf{q}=\mathbf{0}}

where :math:`\mathbf{f}_i` is the application that associate to :math:`\mathbf{x}\in\mathbb{R}^n` the :math:`i^\text{th}` coordinates of :math:`\mathbf{f}(\mathbf{x})\in\mathbb{R}^m`.
By facility, it will be referred as the term "tensor of order 2". For any :math:`r` in :eq:`mclaurin`, there is a corresponding tensor of order :math:`r+1` (corresponding to a multilinear application to :math:`\mathbb{R}^m` and wich can be flatten into a :math:`n^r\, m` dimension vector).

Due to the `Shwarz's theorem <https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives#Schwarz.27s_theorem>`_ (referred as "Kleinman symmetry" in the field of nonlinear optics), some component of those tensors are equals to each others.


Two kind of derivatives are considered for the moment:

+ Electrical, with respect to a static (``F``) or dynamic (``D``) electric field (though the first one is a special case of the second).
  The input space is :math:`\mathbf{f}\in\mathbb{R}^3`, and

  .. math::

      \begin{align}
        E(\mathbf{f}) = E_0  &+  \sum_i \mu_i\,\mathbf{f}_i+ \frac{1}{2!}\,\sum_{ij}\alpha_{ij}\mathbf{f}_i\mathbf{f}_j + \frac{1}{3!}\sum_{ijk}\beta_{ijk}\mathbf{f}_i\mathbf{f}_j\mathbf{f}_k \\
        &+ \frac{1}{4!}\,\sum_{ijkl}\gamma_{ijkl}\mathbf{f}_i\mathbf{f}_j\mathbf{f}_k\mathbf{f}_l + \ldots
      \end{align}

  where :math:`\mu` is the electric dipole moment, :math:`\alpha` is the polarizability, :math:`\beta` and :math:`\gamma` are the first and second hyperpolarizabilites.
  Note that in the case of dynamic electric fields, :math:`\mathbf{f}_i(t)=\mathbf{f}^0_i\,(e^{-i\omega t}+c.c.)`, where :math:`\omega` is the pulsation (:math:`\omega=2\pi\nu`).

+ Geometrical, with respect to cartesian (``G``) or normal mode (``N``) coordinates.
  The input space is :math:`\mathbf{q}\in\mathbf{R}^{3N}` where :math:`3N` is the number of degrees of freedom of :math:`N` atoms expressed in cartesian coordinates, and

  .. math::

      V(\mathbf{q}) = V_0  +  \sum_i F_i\mathbf{q}_i + \frac{1}{2!}\,\sum_{ij}H_{ij}\mathbf{q}_i\mathbf{q}_j +  + \frac{1}{3!}\,\sum_{ij}F_{ijk}\mathbf{q}_i\mathbf{q}_j\mathbf{q}_k+ \ldots

  where :math:`F_i=\left.\frac{\partial V}{\partial\mathbf{q}_i}\right|_{\mathbf{q}=0}` is the gradient (so the forces, wich are close to 0 if the system is in equilibrium), :math:`H_{ij}=\left.\frac{\partial^2 V}{\partial\mathbf{q}_i\mathbf{q}_j}\right|_{\mathbf{q}=0}` is the hessian (force constant matrix) and :math:`F_{ijk}` is the cubic force constants matrix (and so on).

All derivatives are written with respect to the energy (assuming a corresponding Maclaurrin expansion). For example, ``GG`` corresponds to the hessian expressed in cartesian coordinates (second order derivative of the energy), while ``F`` is the electric dipole moment.
Therefore, ``FG`` (or ``GF``) is the first-order derivative of the dipole moment with respect to the cartesian coordinates.

.. note::

    The name of a derivative is expressed as its electrical-derivative counterpart with respect to the geometrical parameters, though the opposite is valid (``GF`` is also the derivative of the gradient with respect to electric field).
    The first is just more used in the literature (for IR and Raman spectra, for example).

Implementation
--------------

*to be continued*

API documentation
-----------------

.. automodule:: qcip_tools.derivatives
    :members: