import numpy
import random
# import unittest

from tests import QcipToolsTestCase
from qcip_tools import symmetry


def simdiags(ops, tol=1e-14, in_vecs=None):
    """
    Private function that finds eigen vals and vecs for degenerate matrices..
    """
    n = len(ops)
    if n == 0:
        return in_vecs
    A = ops[0]

    if in_vecs is not None:
        vecs = numpy.column_stack(in_vecs)
        eigvals, eigvecs = numpy.linalg.eig(numpy.dot(vecs.conj().T, numpy.dot(A, vecs)))
    else:
        eigvals, eigvecs = numpy.linalg.eig(A)

    zipped = list(zip(-eigvals, range(len(eigvals))))
    zipped.sort()
    ds, perm = zip(*zipped)
    ds = -numpy.real(numpy.array(ds))

    print(eigvals[numpy.argsort(eigvals)][::-1], perm)
    perm = numpy.array(perm)  # numpy.argsort(eigvals)
    print(eigvals[perm])

    if in_vecs is not None:
        vecsperm = numpy.zeros(eigvecs.shape, dtype=complex)
        for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
            vecsperm[:, kk] = eigvecs[:, perm[kk]]
        vecs_new = numpy.dot(vecs, vecsperm)
        vecs_out = numpy.array(
            [numpy.zeros((A.shape[0], 1), dtype=complex) for k in range(len(ds))])
        for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
            vecs_out[kk][:, 0] = vecs_new[:, kk]
    else:
        vecs_out = numpy.array(
            [numpy.zeros((A.shape[0], 1), dtype=complex) for k in range(A.shape[0])])

        for kk in range(len(perm)):  # matrix with sorted eigvecs in columns
            vecs_out[kk][:, 0] = eigvecs[:, perm[kk]]

    k = 0
    rng = numpy.arange(len(ds))
    while k < len(ds):
        inds = numpy.array(abs(ds - ds[k]) < max(
            tol, tol * abs(ds[k])))  # get indices of degenerate eigvals
        inds = rng[inds]
        if len(inds) > 1:  # if at least 2 eigvals are degenerate
            vecs_out[inds] = simdiags(numpy.array([ops[jj] for jj in range(1, n)]), tol, vecs_out[inds])
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

    # @unittest.skip('not working, issues to find simultaneoous eigenvectors')
    def test_character_table(self):
        """Try to generate the character table"""

        def normalize_evec(evec):
            if abs(evec[0]) > 1e-5:
                evec /= evec[0]
            else:
                evec /= numpy.min(evec[numpy.where(abs(evec) > 1e-5)])

            return evec

        g = symmetry.PointGroup.D_nd(7)

        class_inverses = numpy.zeros(g.number_of_class, dtype=int)
        class_sizes = numpy.zeros(g.number_of_class, dtype=int)

        matrices = []

        for i, c in enumerate(g.conjugacy_classes):
            e = next(iter(c))
            class_inverses[i] = g.to_conjugacy_class[g.inverse(e)]
            class_sizes[i] = len(c)

        for i in range(0, g.number_of_class):
            matrices.append(g.class_matrix(i))

        e = simdiags(matrices, 1e-14)

        final_eigenvectors = []
        for i in range(g.number_of_class):
            final_eigenvectors.append(normalize_evec(e[i].T[0]))

        # print(final_eigenvectors)

        """for i in range(1, g.number_of_class):
            print('-----')
            final_eigenvectors = []
            if i not in matrices:
                matrices[i] = g.class_matrix(i)
            evals, evecs = numpy.linalg.eig(matrices[i])
            uniques, indexes, inverses, counts = numpy.unique(
                numpy.around(evals, 5), return_inverse=True, return_counts=True, return_index=True)

            evecs = normalize_evecs(evecs.T)
            for ic in range(len(counts)):
                count = counts[ic]
                pos_vecs = numpy.where(inverses == ic)
                if count == 1:  # save eigenvector
                    if abs(evecs[pos_vecs][0][0]) > 1e-5:
                        print('***', uniques[ic], 'saved one', evecs[pos_vecs][0])
                        final_eigenvectors.append(evecs[pos_vecs][0])
                else:
                    degenerated_evecs = evecs[pos_vecs]
                    print('***', uniques[ic], 'tries to extract', count)
                    for j in range(1, g.number_of_class):
                        if i == j:
                            continue
                        espace_mat = numpy.zeros((count, count))
                        if j not in matrices:
                            matrices[j] = g.class_matrix(j)
                        for k in range(count):
                            t = numpy.dot(matrices[j], degenerated_evecs[k])
                            for l in range(count):
                                espace_mat[k, l] = numpy.dot(t, degenerated_evecs[l])
                        # print(espace_mat)
                        espace_evals, espace_evecs = numpy.linalg.eig(espace_mat)
                        if not numpy.all(numpy.around(espace_evals, 5).astype(int) == numpy.around(espace_evals, 1)):
                            print(espace_evals, 'not integer evals, skipping')
                            continue
                        if len(numpy.unique(numpy.around(espace_evals, 5))) == count:
                            print('evals', espace_evals)
                            for k in range(count):
                                n_eigenvector = numpy.zeros(g.number_of_class)
                                for l in range(count):
                                    n_eigenvector += espace_evecs[k, l] * degenerated_evecs[l]

                                if abs(numpy.around(n_eigenvector[0], 5)) > 1e-5:
                                    print('extracted', normalize_evec(n_eigenvector))
                                    final_eigenvectors.append(normalize_evec(n_eigenvector))
                                else:
                                    print('got 0 as first component, it is not possible')
                            break

            # print(final_eigenvectors)"""

        if len(final_eigenvectors) == g.number_of_class:
            print('-------\n')
            print(g.conjugacy_classes)
            print(class_sizes)
            for j in range(g.number_of_class):
                v = final_eigenvectors[j]
                degree = numpy.around(numpy.real(
                    numpy.sqrt(g.order / numpy.einsum('i,i->', v, v[class_inverses] / class_sizes))), 1).astype(int)
                print(numpy.around(v * degree / class_sizes, 3))

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
