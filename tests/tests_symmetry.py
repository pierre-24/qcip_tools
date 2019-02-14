import numpy
import random
# import unittest

import scipy.linalg as la

from tests import QcipToolsTestCase
from qcip_tools import symmetry


def old_simdiag(ops, evals=False):
    """Simulateous diagonalization of communting Hermitian matrices..

    Source: http://qutip.org/docs/latest/modules/qutip/simdiag.html
    """

    # print('---- old simdiag')
    tol = 1e-14
    start_flag = 0
    if not isinstance(ops, (list, numpy.ndarray)):
        ops = numpy.array([ops])
    num_ops = len(ops)
    for jj in range(num_ops):
        A = ops[jj]
        shape = A.shape
        if shape[0] != shape[1]:
            raise TypeError('Matricies must be square.')
        if start_flag == 0:
            s = shape[0]
        if s != shape[0]:
            raise TypeError('All matrices. must be the same shape')
        for kk in range(jj):
            B = ops[kk]
            if numpy.linalg.norm(A * B - B * A) > tol:
                raise TypeError('Matricies must commute.')

    A = ops[0]
    eigvals, eigvecs = numpy.linalg.eig(A)
    zipped = list(zip(-eigvals, range(len(eigvals))))
    zipped.sort()
    ds, perm = zip(*zipped)
    ds = -numpy.real(numpy.array(ds))
    perm = numpy.array(perm)
    eigvecs_array = numpy.array(
        [numpy.zeros((A.shape[0], 1), dtype=complex) for _ in range(A.shape[0])])

    for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
        eigvecs_array[kk][:, 0] = eigvecs[:, perm[kk]]

    # print(eigvecs_array)
    k = 0
    rng = numpy.arange(len(eigvals))
    while k < len(ds):
        # find degenerate eigenvalues, get indicies of degenerate eigvals
        inds = numpy.array(abs(ds - ds[k]) < max(tol, tol * abs(ds[k])))
        inds = rng[inds]
        if len(inds) > 1:  # if at least 2 eigvals are degenerate
            eigvecs_array[inds] = degen(
                tol, eigvecs_array[inds],
                numpy.array([ops[kk] for kk in range(1, num_ops)]))
        k = max(inds) + 1
    eigvals_out = numpy.zeros((num_ops, len(ds)), dtype=float)
    kets_out = numpy.array([(eigvecs_array[j] / la.norm(eigvecs_array[j]))
                            for j in range(len(ds))])
    if not evals:
        return kets_out
    else:
        for kk in range(num_ops):
            for j in range(len(ds)):
                eigvals_out[kk, j] = numpy.real(numpy.dot(
                    eigvecs_array[j].conj().T,
                    ops[kk].data * eigvecs_array[j]))
        return eigvals_out, kets_out


def degen(tol, in_vecs, ops):
    """
    Private function that finds eigen vals and vecs for degenerate matrices..
    """
    # print('--- old degen')
    n = len(ops)
    if n == 0:
        return in_vecs
    A = ops[0]
    vecs = numpy.column_stack(in_vecs)
    eigvals, eigvecs = numpy.linalg.eig(numpy.dot(vecs.conj().T, numpy.dot(A, vecs)))
    zipped = list(zip(-eigvals, range(len(eigvals))))
    zipped.sort()
    ds, perm = zip(*zipped)
    ds = -numpy.real(numpy.array(ds))
    perm = numpy.array(perm)
    vecsperm = numpy.zeros(eigvecs.shape, dtype=complex)
    for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
        vecsperm[:, kk] = eigvecs[:, perm[kk]]
    vecs_new = numpy.dot(vecs, vecsperm)
    vecs_out = numpy.array(
        [numpy.zeros((A.shape[0], 1), dtype=complex) for k in range(len(ds))])
    for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
        vecs_out[kk][:, 0] = vecs_new[:, kk]

    # print(vecs_out)
    k = 0
    rng = numpy.arange(len(ds))
    while k < len(ds):
        inds = numpy.array(abs(ds - ds[k]) < max(
            tol, tol * abs(ds[k])))  # get indicies of degenerate eigvals
        inds = rng[inds]
        if len(inds) > 1:  # if at least 2 eigvals are degenerate
            vecs_out[inds] = degen(tol, vecs_out[inds],
                                   numpy.array([ops[jj] for jj in range(1, n)]))
        k = max(inds) + 1
    return vecs_out


class SymmetryTestCase(QcipToolsTestCase):

    def test_closest_fraction(self):
        self.assertEqual(symmetry.closest_fraction(-.499999, max_denominator=1000), (-1, 2))
        self.assertEqual(symmetry.closest_fraction(.499999, max_denominator=1000), (1, 2))

    def test_group(self):
        f = symmetry.BinaryOperation(symmetry.Set(range(4)), lambda e: (e[0] + e[1]) % 4)
        self.assertTrue(f.check_surjectivity())
        self.assertTrue(f.check_associativity())

        g = symmetry.Group(f)  # just a cyclic group of order 4, then
        self.assertEqual(g.identity(), 0)
        self.assertEqual(g.inverse(0), 0)
        self.assertEqual(g.inverse(1), 3)
        self.assertEqual(g.inverse(2), 2)

        self.assertTrue(g.abelian)

        for i in range(4):
            self.assertEqual(g.inverse(i) ** -1, i)

        n3 = symmetry.GroupElement(3, g)
        self.assertEqual(n3, 3)
        self.assertEqual(n3 ** 1, 3)
        self.assertEqual(n3 ** 2, 2)
        self.assertEqual(n3 ** 3, 1)
        self.assertEqual(n3 ** -1, 1)
        self.assertEqual(n3 ** -2, 2)

    def test_symmetry_element(self):
        p = numpy.array([random.randint(-3, 3), random.randint(-3, 3), random.randint(-3, 3)])

        # basic elements
        E = symmetry.Operation.E()
        self.assertArraysAlmostEqual(E.apply(p), p)
        self.assertArraysAlmostEqual(numpy.dot(E.matrix_representation(), p), p)

        i = symmetry.Operation.i()
        self.assertArraysAlmostEqual(i.apply(p), -p)
        self.assertArraysAlmostEqual(numpy.dot(i.matrix_representation(), p), -p)

        Oz = symmetry.Operation.sigma()
        self.assertArraysAlmostEqual(Oz.apply(p), numpy.array([*p[:2], -p[2]]))
        self.assertArraysAlmostEqual(numpy.dot(Oz.matrix_representation(), p), numpy.array([*p[:2], -p[2]]))

        S52z = symmetry.Operation.S(5, 2, axis=numpy.array([1, 1., 0.]))
        self.assertArraysAlmostEqual(S52z.apply(p), numpy.dot(S52z.matrix_representation(), p))

        # test identity
        t = Oz * Oz
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)
        t = i * i
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that C_3 is cyclic
        C3z = symmetry.Operation.C(3)
        C23z = symmetry.Operation.C(3, 2)

        t = C3z * C3z
        self.assertEqual(t, C23z)
        self.assertEqual(t.get_description(), C23z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * C3z
        self.assertEqual(t, E)
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check that S_3 is cyclic
        S3z = symmetry.Operation.S(3)
        S53z = symmetry.Operation.S(3, 5)

        t = S3z * S3z
        self.assertEqual(t, C23z)
        self.assertEqual(t.get_description(), C23z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C23z.apply(p))

        t = t * S3z
        self.assertEqual(t, Oz)
        # TODO: self.assertEqual(t.get_description(), Oz.get_description()) → k = 2
        self.assertArraysAlmostEqual(t.apply(p), Oz.apply(p))

        t = t * S3z
        self.assertEqual(t, C3z)
        self.assertEqual(t.get_description(), C3z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C3z.apply(p))

        t = t * S3z
        self.assertEqual(t, S53z)
        self.assertEqual(t.get_description(), S53z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), S53z.apply(p))

        t = t * S3z
        self.assertEqual(t, E)
        # TODO: self.assertEqual(t.get_description(), E.get_description()) → all that matter is the symbol
        self.assertArraysAlmostEqual(t.apply(p), p)

        # Check composition: C_n(z) * O(z) = S_n(z)
        n = 5
        Cnz = symmetry.Operation.C(n)
        Snz = symmetry.Operation.S(n)
        t = Cnz * Oz
        self.assertEqual(t, Snz)
        self.assertEqual(t.get_description(), Snz.get_description())
        self.assertArraysAlmostEqual(t.apply(p), Snz.apply(p))

        # Check composition: C_2(y) * C_2(x) = C_2(z)
        C2z = symmetry.Operation.C(2, axis=numpy.array([0, 0, 1.]))
        C2y = symmetry.Operation.C(2, axis=numpy.array([0, 1., 0]))
        C2x = symmetry.Operation.C(2, axis=numpy.array([1., 0, 0]))
        t = C2y * C2x
        self.assertEqual(t, C2z)
        self.assertEqual(t.get_description(), C2z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(x) = C_2(z)
        Ox = symmetry.Operation.sigma(axis=numpy.array([1, 0., 0]))
        Oy = symmetry.Operation.sigma(axis=numpy.array([0, 1., 0]))
        t = Oy * Ox
        self.assertEqual(t, C2z)
        self.assertEqual(t.get_description(), C2z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C2z.apply(p))

        # Check composition: O(y) * O(xy') = C_4(z)
        Oxy = symmetry.Operation.sigma(axis=numpy.array([1., -1., 0]))
        C4z = symmetry.Operation.C(4)
        t = Oy * Oxy
        self.assertEqual(t, C4z)
        self.assertEqual(t.get_description(), C4z.get_description())
        self.assertArraysAlmostEqual(t.apply(p), C4z.apply(p))

    def test_point_groups(self):

        # checked against http://gernot-katzers-spice-pages.com/character_tables/index.html
        groups = [
            (symmetry.PointGroup.C_n(4), 4, True, 4),
            (symmetry.PointGroup.S_n(2), 2, True, 2),  # = C_i
            (symmetry.PointGroup.S_n(4), 4, True, 4),
            (symmetry.PointGroup.C_nv(2), 4, True, 4),
            (symmetry.PointGroup.C_nv(3), 6, False, 3),
            (symmetry.PointGroup.C_nv(4), 8, False, 5),
            (symmetry.PointGroup.C_nh(1), 2, True, 2),  # = C_s
            (symmetry.PointGroup.C_nh(2), 4, True, 4),  # != S_2
            (symmetry.PointGroup.C_nh(3), 6, True, 6),
            (symmetry.PointGroup.C_nh(4), 8, True, 8),  # != S_4
            (symmetry.PointGroup.D_n(4), 8, False, 5),
            (symmetry.PointGroup.D_nh(3), 12, False, 6),
            (symmetry.PointGroup.D_nh(4), 16, False, 10),
            (symmetry.PointGroup.D_nd(2), 8, False, 5),
            (symmetry.PointGroup.D_nd(3), 12, False, 6),
            (symmetry.PointGroup.T(), 12, False, 4),
            (symmetry.PointGroup.T_h(), 24, False, 8),
            (symmetry.PointGroup.T_d(), 24, False, 5),
            (symmetry.PointGroup.O(), 24, False, 5),
            (symmetry.PointGroup.O_h(), 48, False, 10),
            (symmetry.PointGroup.I(), 60, False, 5),  # ~2 seconds
            (symmetry.PointGroup.I_h(), 120, False, 10),  # ~ 5 seconds
        ]

        for g, n_elements, is_abelian, number_of_class in groups:
            self.assertEqual(len(g.G), n_elements)
            self.assertEqual(g.abelian, is_abelian)
            self.assertEqual(g.number_of_class, number_of_class)
            if len(g.G) <= 12:
                self.assertTrue(g.binary_operation.check_surjectivity())  # O(n²)
                self.assertTrue(g.binary_operation.check_associativity())  # O(n³)

    def test_conjugacy_classes(self):
        """Tests if the conjugacy classes comes in order"""

        C_3v = symmetry.PointGroup.C_nv(3)
        self.assertTrue(
            list(len(x) for x in C_3v.conjugacy_classes) == [1, 2, 3])  # ... Which gives 6 elements in total

        self.assertEqual(C_3v.conjugacy_classes[0], {C_3v.e})  # first class is always identity
        self.assertEqual(
            {e.element for e in C_3v.conjugacy_classes[1]}, {symmetry.Operation.C(3), symmetry.Operation.C(3, 2)})

    def test_character_table(self):
        """Try to generate the character table"""

        def simultaneous_diagonalization(matrices, tol=1e-14, in_vecs=None):
            """
            Simulateous diagonalization of (communting Hermitian) matrices.

            Seems to works by computing a *basis* for the first matrix, then by diagonalizing the different
            blocks (organized by the degeneracy of the eigenvalues).

            Original source: http://qutip.org/docs/latest/modules/qutip/simdiag.html.
            Also check
            https://ocw.mit.edu/courses/physics/8-05-quantum-physics-ii-fall-2013/lecture-notes/MIT8_05F13_Chap_05.pdf
            (which is the closest explanation I found of this implementation, which is not common).

            .. warning::

                For unknown reasons, fails for :math:`D_{4h}`, :math:`O_{h}`, :math:`T_{h}` and :math:`I_{h}`.

                + For :math:`D_{4h}` use ``del mat[1]`` ;
                + For :math:`T_{h}` use ``del mat[0]`` and ``del mat[3]`` (T representations need to be rescaled) ;
                + For :math:`O_{h}` use ``del mat[0]`` and ``del mat[4]`` (T representations need to be rescaled) :
                + For :math:`I_{h}` use ``del mat[2]`` (H representations need to be rescaled).

                That's the best I can do for the moment ;)

            :param matrices: list of matrices
            :type matrices: list
            :param tol: tolerance on the eigenvalue (to check if they are degenerate or not)
            :type tol: float
            :param in_vecs: previously computed eigenvectors
            :type in_vecs: numpy.ndarray
            """

            n = len(matrices)
            if n == 0:
                return in_vecs

            A = matrices[0]
            if in_vecs is not None:
                vecs = numpy.column_stack(in_vecs)
                eigvals, eigvecs = numpy.linalg.eig(numpy.dot(vecs.conj().T, numpy.dot(A, vecs)))
            else:
                eigvals, eigvecs = numpy.linalg.eig(A)

            perm = numpy.argsort(-eigvals)  # sort in reverse order (largest first)
            ds = eigvals[perm]

            if in_vecs is not None:
                vecsperm = numpy.zeros(eigvecs.shape, dtype=complex)
                for kk in range(len(perm)):
                    vecsperm[:, kk] = eigvecs[:, perm[kk]]
                vecs_new = numpy.dot(vecs, vecsperm)
                vecs_out = numpy.zeros((len(ds), A.shape[0]), dtype=complex)
                for kk in range(len(perm)):
                    vecs_out[kk] = vecs_new[:, kk]
            else:
                vecs_out = numpy.zeros(A.shape, dtype=complex)
                for kk in range(len(perm)):
                    vecs_out[kk] = eigvecs[:, perm[kk]]

            k = 0
            rng = numpy.arange(len(ds))

            while k < len(ds):
                indices = rng[numpy.where(abs(ds - ds[k]) < max(tol, tol * abs(ds[k])))]
                if len(indices) > 1:  # two (or more) degenerate eigenvalues
                    vecs_out[indices] = simultaneous_diagonalization(matrices[1:], tol, vecs_out[indices])
                k = max(indices) + 1

            return vecs_out

        def normalize_evec(evec):
            if abs(evec[0]) > 1e-14:
                evec /= evec[0]
            else:
                evec /= numpy.min(evec[numpy.where(abs(evec) > 1e-5)])

            return evec

        g = symmetry.PointGroup.T_d()

        class_inverses = numpy.zeros(g.number_of_class, dtype=int)
        class_sizes = numpy.zeros(g.number_of_class, dtype=int)

        matrices = []

        for i, c in enumerate(g.conjugacy_classes):
            e = next(iter(c))
            class_inverses[i] = g.to_conjugacy_class[g.inverse(e)]
            class_sizes[i] = len(c)

        for i in range(0, g.number_of_class):
            matrices.append(g.class_matrix(i))

        # del matrices[0]
        # del matrices[3]
        e = simultaneous_diagonalization(matrices)
        # ed = old_simdiag(matrices)

        final_eigenvectors = []
        for i in range(g.number_of_class):
            # final_eigenvectors.append(normalize_evec(ed[i].T[0]))
            final_eigenvectors.append(normalize_evec(e[i]))

        if len(final_eigenvectors) == g.number_of_class:
            print(g.conjugacy_classes)
            table = numpy.zeros((g.number_of_class, g.number_of_class), dtype=complex)
            for j in range(g.number_of_class):
                v = final_eigenvectors[j]
                degree = numpy.sqrt(g.order / numpy.einsum('i,i->', v, v[class_inverses] / class_sizes))
                # print('* d', degree)
                # print('* v', numpy.around(v, 3))
                table[j] = v * degree / class_sizes

            print(numpy.around(table, 2))

    def test_symmetry_finder(self):
        """Test if one is able to detect symmetry"""

        def s(points, tol=1e-5):
            d = symmetry.SymmetryFinder(points, tol).find_symmetry()[0]
            return d.symbol, d.order

        """O_h geometry:

           1 1
           |/
        1--3--1
          /|
         1 1
        """

        p = numpy.array([
            (3., 0, 0, 0),
            (1., 1, 0, 0),
            (1., -1, 0, 0),
            (1., 0, 1, 0),
            (1., 0, -1, 0),
            (1., 0, 0, 1),
            (1., 0, 0, -1),
        ])

        self.assertEqual(s(p), (symmetry.PointGroupType.octahedral_achiral, 0))
        self.assertEqual(s(p[:-1]), (symmetry.PointGroupType.pyramidal, 4))
        self.assertEqual(s(p[:-2]), (symmetry.PointGroupType.prismatic, 4))
        self.assertEqual(s(p[:-3]), (symmetry.PointGroupType.pyramidal, 2))
        self.assertEqual(s(p[:-4]), (symmetry.PointGroupType.prismatic, -1))
        self.assertEqual(s(p[:-5]), (symmetry.PointGroupType.pyramidal, -1))
