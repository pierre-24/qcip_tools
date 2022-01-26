==================================
Symmetry (``qcip_tools.symmetry``)
==================================

Handle symmetry.

Introduction
------------

What's a group anyway ?
=======================

A `group <https://en.wikipedia.org/wiki/Group_(mathematics)>`_ is a set :math:`G`, together with a binary operation :math:`\star` (group law),
form a group :math:`(G, \star)` if it satisfies 4 requirements:

- Closure :math:`\forall a, b \in \mathcal{G}: c=a\star b \land c\in\mathcal{G}` ;
- Associativity: :math:`\forall a, b, c \in \mathcal{G}: (a\star b)\star c = a\star (b\star c)`
- (unique) Identity: :math:`\exists! e\in\mathcal{G},\forall a\in\mathcal{G}: a\star e = a = e\star a`;
- (unique) Inverse: :math:`\forall a\in\mathcal{G}:\exists! a'\in\mathcal{G}: a\star a' = e`.

In the code, the binary operation is represented by ``BinaryOperation``, which works on a given codomain and uses a surjective function:

.. code-block:: python

    from qcip_tools import symmetry

    f = symmetry.BinaryOperation(
        symmetry.Set(range(4)),  # domain is {0, 1, 2, 3}
        lambda e: (e[0] + e[1]) % 4 # function is modulo 4
    )

    g = symmetry.Group(f)  # cyclic group of order 4

    e = g.identity()  # there is an identity element (here 0)
    inv1 = g.inverse(1)  # the inverse of "1" is "3"

    # and the inverse of 1 times 1 is identity
    assert symmetry.GroupElement(1, g) * inv1 == e

This implementation also compute the different classes (which includes, at least, an element and its inverse) and the Cayley table (multiplication table).

Point group
===========

By using symmetry operation as group element (and matrix multiplication as the binary operation), one can create a `point group <https://en.wikipedia.org/wiki/Finite_group>`_, which apply on finite objects (in oposition to space group, which aply to infinite objects).

Symmetry operation are of two types, basically:

+ Proper rotation of :math:`k\times \frac{2\pi}{n}` around an axis, noted :math:`C_n^k`.
+ Improper rotation of :math:`k\times \frac{2\pi}{n}` around an axis, noted :math:`S_n^k`
  (a rotation of :math:`\frac{2\pi}{n}` followed by a reflexion, repeated :math:`k` times).

From that, one can generate:

+ Identity, :math:`E=C_n^n` ;
+ Reflexion, :math:`\sigma=S_1` ;
+ Inversion, :math:`i = S_2^1`.

Those operation are generated trough function of the `Operation class <#qcip_tools.symmetry.Operation>`_ (internally, it uses quaternions).
One can generate simply any point group by using the methods of `PointGroup class <#qcip_tools.symmetry.PointGroup>`_.

.. code-block:: python

    C_3v = symmetry.PointGroup.C_nv(3)  # generate C3v

    self.assertEqual(C_3v.conjugacy_classes[0], {C_3v.e})  # first class is always identity
    self.assertEqual(  # second class contains the proper rotation
        {e.element for e in C_3v.conjugacy_classes[1]}, {symmetry.Operation.C(3), symmetry.Operation.C(3, 2)})

    # third class contains the reflexion planes (and, among other, the one that sets x)
    self.assertIn(symmetry.Operation.sigma(axis=numpy.array([0, 1., 0])), C_3v.conjugacy_classes[2])

Character table
===============

An useful tools when working with group theory is the `character table <https://en.wikipedia.org/wiki/Character_table>`_, which contains irreducible representations of character for each class of the group.

.. code-block:: python

    C_2v = symmetry.PointGroup.C_nv(2)
    t = C_2v.character_table()

    A1 = t.irreducible_representations[0]  # A1 (trivial) representation
    A2 = t.irreducible_representations[1]  # A2
    r = A1 + A2  # reducible representation contains A1 + A2
    assert A2 in A1 * A2 # direct product gives A2

One can also do ``print(t)`` to get the formated output:

.. code-block:: text

    C_2v (h=4)
    ------------------------------------------------------
           E           C_2         σ_v         σ_v
    ------------------------------------------------------
    A1     1.00        1.00        1.00        1.00
    A2     1.00        1.00       -1.00       -1.00
    B1     1.00       -1.00        1.00       -1.00
    B2     1.00       -1.00       -1.00        1.00
    ------------------------------------------------------

Discover symmetry
=================

Use `SymmetryFinder <#qcip_tools.symmetry.SymmetryFinder>`_.
In particular, the documentation of `find_symmetry() <#qcip_tools.symmetry.SymmetryFinder.find_symmetry>`_ describes the way the point group is found.

.. warning::

    This was not tested for all point groups. In particular, it detects icosahedral point groups correctly, but outputs a wrong orientation.

In the case of molecule, see `qcip_tools.molecule.MolecularSymmetryFinder <molecule.html#qcip_tools.molecule.MolecularSymmetryFinder>`_.

API documentation
-----------------

.. automodule:: qcip_tools.symmetry
    :undoc-members:
    :members: