================================================================
Derivatives w.r.t. electric field (``qcip_tools.derivatives_e``)
================================================================

Tools to help the manipulation of derivatives of the energy with respect to a static (``F``) or dynamic (``D``) electric field (though the first one is a special case of the second).

Theory
------

Light can be described as an electromagnetic radiation, with time and space varying components.
Generally, only the phenomena associated with its electric components, because dominant, are considered.
If the space-dependent term is left constant, the electric part of the field reads

.. math::
    \mathbf{f}_i(\omega) = \mathbf{f}^0_i\,(e^{-i\omega t} + e^{i\omega t}),

where :math:`\mathbf{f}^0` is the amplitude of the the electric field, :math:`\omega` is the pulsation (equals to :math:`2\pi\nu`, with :math:`\nu` the frequency in :math:`s^{-1}`) and :math:`t` the time (in :math:`s`).
This situation refers to the *electric dipole approximation*.

A material is an ensemble of particles, negatively (electrons) and positively (nuclei) charged. An oscillating electric field interacts with all particles, though the motion of the electrons is the most significant part (because nuclei  have a smaller speed).
This leads to an oscillating induced dipole moment, or, at the macroscopic scale, a (change of) polarization. The total dipole of a molecule (sum of the induced and permanent dipole, :math:`\vec{\mu}_0`), can be described using a Taylor series. For the :math:`i`-components of the dipole moment,

.. math::

    \begin{align}
    \mu_i = \mu^0_i &+ \sum^{xyz}_j \alpha_{ij}(-\omega_{\sigma};\omega_1)\,\mathbf{f}_{j}(\omega_1) \nonumber\\
    &+ \frac{1}{2!} \sum^{xyz}_{jk} \beta_{ijk}(-\omega_{\sigma};\omega_1,\omega_2)\,\mathbf{f}_{j}(\omega_1)\,\mathbf{f}_{k}(\omega_2) \nonumber\\
    &+ \frac{1}{3!} \sum^{xyz}_{jkl} \gamma_{ijkl}(-\omega_{\sigma};\omega_1,\omega_2,\omega_3)\,\mathbf{f}_{j}(\omega_1)\,\mathbf{f}_{k}(\omega_2)\,\mathbf{f}_{l}(\omega_3)+ \ldots
    \end{align}

where :math:`\alpha` and :math:`\beta` are tensors of rank 2 and 3, called the polarizability and the first hyperpolarizability, respectively (generaly expressed in atomic units, *au*).
:math:`\omega_1`, :math:`\omega_2` and :math:`\omega_3` are the pulsations of the fields applied in the :math:`j`, :math:`k` and  :math:`l` directions, while :math:`\omega_\sigma = \sum_i \omega_i`.
For small electric fields, the polarizability term is dominant, but the hyperpolarizability terms becomes non-negligible when the field gets large.

Sources
-------

+ R.W. Boyd,  « Nonlinear Optics (3rd edition) ». *Elsevier*, 2008.
+ F\. Castet *et al*. *Acc. Chem. Res.* **46**, 2656 (2013).



API documentation
-----------------

.. automodule:: qcip_tools.derivatives_e
    :members: