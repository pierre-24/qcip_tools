====================================
Basis set (``qcip_tools.basis_set``)
====================================


Theory
------

Definition
**********

Quantum calculations are, most of the time, based on the LCAO (linear combination of atomic orbitals) principle,

.. math::

    \Psi_i = \sum_j^N c_{ij} \phi_{j},

where :math:`c_{ij}` are coefficient of the linear combination and :math:`\phi_j` are the atomic orbitals.
The atomic orbitals (or **basis functions**) are solution of the Hartree-Fock equation for the atom.
Slater type obritals (STO, :math:`\propto e^{-r}`) are sometimes used to describe the atomic orbitals, due to their similarity to hydrogen atom solutions, but Gaussian type orbitals (GTO) are more convenient from a computationnal point of view.
The general expression (in carthesian coordinates) of a GTO centered in :math:`\mathbf{R}` is

.. math::

    G_{lmn}^{\mathbf{R}, \alpha}(\mathbf{r}) = N_{lmn}^\alpha\,(r_x - R_x)^l\,(r_y-R_y)^m\,(r_z - R_z)^n\,e^{-\alpha\,|\mathbf{r}-\mathbf{R}|^2},

where :math:`\alpha` (the letter :math:`\zeta` is also used) is the exponent of the GTO (which controls the *width* of the orbital, small :math:`\alpha` giving large/diffuse functions), and :math:`N_{lmn}^\alpha` is a normalization factor.

To mimic STO's, a linear combination (a contraction) of GTO's (**primitives**) is used to form a **basis function**.

The sum of exponents :math:`L=l+m+n` (which are not quantum numbers) is used analogously to the angular quantum number for atom, to mark function as *s*-type (:math:`L=0`), *p*-type (:math:`L=1`), *d*-type  (:math:`L=2`), etc.
Primitives (and basis functions) are regrouped in shells (S, P, D, F ...), forming a collection of GTO's with the same :math:`L`: in each shell :math:`\mu`, the primitives are contracted in :math:`n_\mu` basis functions :math:`\mu\gamma` (for example, with *s*-type basis function, the usual atomic orbitals 1s, 2s, 3s ...).

The general expression for an atomic basis set centered in :math:`\mathbf{R}` is therefore:

.. math::
    :label: atomicbasisset

    \phi^{\mathbf{R}}(\mathbf{r}) = \sum_\mu^{\text{S, P, D, }\ldots}\sum_{(l,m,n)\in \mu} \sum_{\gamma}^{n_\mu} \underbrace{N_{\mu\gamma}\,\sum_k^{n_{\mu\gamma}} c_{\mu\gamma k}\,G_{lmn}^{\mathbf{R}, \alpha_{\mu\gamma k}}(\mathbf{r})}_{\begin{array}{c}\text{basis function}\,\mu\gamma\,\equiv\\\text{contraction of }n_{\mu\gamma}\text{ primitives}\end{array}},

where :math:`c_{\mu\gamma k}` is the contraction coefficient of the primitive :math:`k` in basis function :math:`\mu\gamma` in shell :math:`\mu`, to which correspond a exponent :math:`\alpha_{\mu\gamma k}`.
:math:`N_{\mu\gamma}` is a normalization constant for the whole basis function :math:`\gamma` in :math:`\mu`.

In the calculation, a LCAO coefficient is associated to each :math:`\mu\gamma` basis function.


Construction
************

Some basis sets are constucted using the same number of basis function in each shell. This number of contractions helps to label the basis set as *n*-uple :math:`\zeta`: simple :math:`\zeta` (for example the infamous STO-3G), double :math:`\zeta` (DZ),  triple :math:`\zeta` (TZ), quadruple :math:`\zeta` (QZ) and so all (5Z, 6Z, ...).

On the other hand, since valence orbitals are affected than the core (inner) orbitals by chemical processes, more basis functions (contractions) are sometimes used to describe these (which does not mean more primitives).
The corresponding basis sets (for example the Pople basis sets) are sometimes called *split-valence*.
These basis sets are usually structured in such a way that most diffuse primitives (small :math:`\alpha`) are left uncontracted.
For example, in the 6-31G basis set of Pople, a contraction of 6 primitives is used for core orbitals while 2 contractions (of 3 and 1 primitives, respectively) is used for valence orbitals. It is therefore referred as double :math:`\zeta`.

Basis sets are frequently augmented with others functions:

+ *Polarization functions*: functions with higher values of :math:`L` than those present in occupied atomic orbitals of the ground state of the corresponding atom (or, strictly speaking, of the noble gas belonging to the same row in the periodic table).
  They are important for correlated calculations as well as description of bonding.
  In Pople's notation, they are indicated in parentheses, for example 6-31G(d) include a *d*-type basis function.
+ *Diffuse functions*: primitives with very small (< 0.01 ?) exponent, contained in uncontracted basis functions, that therefore decay slowly with :math:`|\mathbf{r}-\mathbf{R}|^2`.
  They are important for the computation of properties, weak bonds (i.e. hydrogen bonds) and long distance interactions.
  In Pople's notation, they are indicated by `+`, for example 6+31G include a *p*-type diffuse basis function on heavy atoms.
  Note that "polarization-diffuse" functions are possible (in Dunning "aug" basis sets, for example), although not in Pople's basis sets.


.. note::

    **Shells**. Even thought it does not make sense strictly speaking (because of the orthogonality of primitives of different shells), the "SP" shell exists (in Pople basis sets), and define a shell when GTO share the same exponent but different contraction coefficients for the *s* and *p* orbitals.

    **Naming of basis sets**. Basis set are named based on arbitrary criterions (with the author name, generally) or more systematically (Pople or Dunning).
    One notation, used by Dalton and some authors in the literature is ``[uncontracted set|contracted set]``.
    For example ``[10s4p1d|3s2p1d]`` for second row atoms (6-31G* in the notation of Pople), which is constituted of 10 *s*-type primitives (contracted in 3 *s*-type basis functions), 4 *p*-type primitives (contracted in 2 *p*-type basis functions) and 1 *d*-type primitive (in one *d*-type polarization basis function).
    The representation of a atomic basis function in the code follows this notation.

    **Basis functions**. Therefore, when using a basis set defined by ``[10s4p1d|3s2p1d]`` for second row atoms, in means that there is, on an atomic point of view the [1s, 2s, 2p, 3s', 3p', 3d'] orbitals (where ' represents the added polarization orbitals).

    **Number of basis basis functions and primitives**. Don't forget that for each *p*-type "basis function", there is actually 3 basis function (:math:`p_x`, :math:`p_y` and :math:`p_z`) defined, as emphasized by the sum of over :math:`(l,m,n)\in\mu` in Eq. :eq:`atomicbasisset`, 5 for each *d*-type orbitals (or six if one consider the carthesian counterparts), and so all. The same goes for the number of primitives.
    For example, water with a 6-31+G(d) basis set (``[11s5p1d|4s3p1d]`` for oxygen, ``[4s|2s]`` for hydrogen) contains actually 18 basis functions for the oxygen (= :math:`4\times 1+3\times 3+1\times 5`) and 31 primitives (= :math:`11\times 1+5\times 3 + 1\times 5`), so in total there is 22 basis functions and 39 GTO's for the water molecule.

    **Output**. By convention, the primitives of a basis set are listed shell by shell, in decrasing order of number of primitive per basis function and starting from the highest coefficient.
    The code handle that order by sorting internal lists.


Sources
*******

+ http://www.ccl.net/cca/documents/basis-sets/img/basis-sets.shtml ;
+ http://vergil.chemistry.gatech.edu/courses/chem6485/pdf/basis-sets.pdf ;
+ http://www.helsinki.fi/kemia/fysikaalinen/opetus/jlk/luennot/Lecture5.pdf ;
+ https://en.wikipedia.org/wiki/Basis_set_(chemistry) ;

(last consultation on the 18th of August, 2017)

Implementation
--------------

A primitive is defined with the `Primitive <#qcip_tools.basis_set.Primitive>`_ object, and represent basically a contraction coefficient and and exponent.
Basis functions are defined with the `Function <#qcip_tools.basis_set.Function>`_ object, which contains the different primitives.
For a given atom, one defines a basis set with `AtomicBasisSet <#qcip_tools.basis_set.AtomicBasisSet>`_, in which you define the different basis functions for a given shell (S, P or SP, D, F ...).
Finally, a basis set is defined by the `BasisSet <#qcip_tools.basis_set.BasisSet>`_ object, which contains the different basis sets for each atoms.

.. code-block:: python

    from qcip_tools import basis_set

    # STO-3G for hydrogen or oxygen (with SP):
    # From ESML basis set exchange
    sto3g = basis_set.BasisSet('STO-3G')

    bf_H_S = basis_set.Function()
    bf_H_S.add_primitive(basis_set.Primitive(3.42525091, 0.15432897))  # exponent, then coefficient
    bf_H_S.add_primitive(basis_set.Primitive(0.62391373, 0.53532814))
    bf_H_S.add_primitive(basis_set.Primitive(0.16885540, 0.44463454))

    b_H = basis_set.AtomicBasisSet('H')  # you can use symbol or atomic number
    b_H.add_basis_function('s', bf_H_S)  # you can use upper or lowercase

    bf_O_S = basis_set.Function()
    bf_O_S.add_primitive(basis_set.Primitive(130.7093200, 0.15432897))
    bf_O_S.add_primitive(basis_set.Primitive(23.8088610, 0.53532814))
    bf_O_S.add_primitive(basis_set.Primitive(6.4436083, 0.44463454))

    bf_O_SP = basis_set.Function()
    bf_O_SP.add_primitive(basis_set.Primitive(5.0331513, -0.09996723, 0.15591627))  # "sp" primitives have two coefficients
    bf_O_SP.add_primitive(basis_set.Primitive(1.1695961, 0.39951283, 0.60768372))
    bf_O_SP.add_primitive(basis_set.Primitive(0.3803890, 0.70011547, 0.39195739))

    b_O = basis_set.AtomicBasisSet('O')
    b_O.add_basis_function('s', bf_O_S)
    b_O.add_basis_function('sp', bf_O_SP)

    sto3g.add_atomic_basis_set(b_O)  # add the atomic basis set for oxygen
    sto3g.add_atomic_basis_set(b_H)

API documentation
-----------------

.. automodule:: qcip_tools.basis_set
    :members:

.. automodule:: qcip_tools.basis_set_esml
    :members: