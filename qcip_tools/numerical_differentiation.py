import math
import numpy
import sys

from qcip_tools import assert_in_domain


def akf_(a, k, f):
    return (a ** k) * f


def sgn(a):
    if a < 0:
        return -1
    elif a > 1:
        return 1
    else:
        return 0


def ak_shifted(ratio, q):
    """
    Match the :math:`k` and the real :math:`a^k` value.

    For a given :math:`q`, returns

    .. math::
        \\left\\{
            \\begin{array}{ll}
                -a^{|q|-1} & \\text{if }q<0\\\\
                0 & \\text{if }k=0\\\\
                a^{q-1}&\\text{if }q>0
            \\end{array}\\right.

    :param q: q
    :type q: int
    :param ratio: ratio (a)
    :type ratio: float
    :rtype: float
    """
    if q == 0:
        return 0
    sgn = 1.
    if q < 0:
        sgn = -1.
    return sgn * (ratio ** (abs(q) - 1))


class Coefficients:
    """Class to store the coefficients for the derivative calculation.

        :param d: derivation order
        :type d: int
        :param p:  order of error
        :type p: int
        :param ratio: ratio (a)
        :type ratio: float
        :param method: derivation method, either forward (F), backward (B) or centered (C)
        :type method: str
    """

    def __init__(self, d, p, ratio=2., method='C'):

        if method not in ['C', 'F', 'B']:
            raise ValueError('{} not in CFB'.format(method))

        self.method = method
        self.i_min = 0

        if method == 'C':
            if (d + p) % 2 != 1:
                raise Exception('d+p should be odd for centered derivative')

            self.i_min = -(d + p - 1) / 2.
        if method == 'B':
            self.i_min = -(d + p - 1)

        self.derivation_order = d
        self.error_order = p

        self.ratio = ratio

        mat_solve = numpy.zeros((d + p, d + p))
        mat_res = numpy.zeros(d + p)
        self.mat_i = numpy.zeros(d + p)
        for n in range(d + p):
            for i_ in range(d + p):
                i = self.i_min + i_
                self.mat_i[i_] = i
                mat_solve[n, i_] = ak_shifted(ratio, i) ** n

            if n == d:
                mat_res[n] = 1

        self.mat_coefs = numpy.dot(numpy.linalg.inv(mat_solve), mat_res.transpose())

    @staticmethod
    def choose_p_for_centered(d, shift=0):
        """Choose the precision so that :math:`d+p-1` is even (for centered derivatives)

        :param d: order of the derivative
        :type d: int
        :param shift: increase the precision
        :type shift: int
        """
        if d % 2 == 0:
            return 1 + shift * 2
        else:
            return 2 + shift * 2

    def prefactor(self, k, h0):
        """return the :math:`\\frac{d!}{h^d}` prefactor, with :math:`h=a^kh_0`.

        :param h0: minimal value
        :type h0: float
        :param k: Minimal amplitude
        :type k: int
        :return: prefactor
        :rtype: float
        """
        return math.factorial(self.derivation_order) / akf_(self.ratio, k, h0) ** self.derivation_order


def compute_derivative_of_function(c, scalar_function, k, h0, input_space_dimension=1, **kwargs):
    """
    Compute the numerical derivative of a scalar, calling the appropriate results through ``scalar_function``,
    to which the ``**kwargs`` parameter is provided.

    :param c: list of tuple to describe the numerical derivative (Coefficient, coordinate)
    :type c: list
    :param scalar_function: a callback function to access the different quantities
    :type scalar_function: callback
    :param k: minimal amplitude
    :type k: int
    :param h0: minimal value
    :type h0: float
    :param input_space_dimension: dimension of the input vector space
    :type input_space_dimension: int
    :param kwargs: other arguments, transfered to ``scalar_function()``
    :type kwargs: dict
    :return: value of the derivative (as a scalar)
    :rtype: float
    """

    prev_coords = []

    for i in c:
        if i[1] in prev_coords:
            raise Exception('used coordinates `{}` twice!'.format(i[1]))

        prev_coords.append(i[1])

    coefficients = c[0][0].mat_coefs.copy()
    for s in c[1:]:
        ndc = s[0]
        coefficients = numpy.tensordot(coefficients, ndc.mat_coefs, axes=0)

    accum = .0
    it = numpy.nditer(coefficients, flags=['multi_index'])
    while not it.finished:
        if it[0] == .0:
            it.iternext()
            continue

        index = it.multi_index
        fields = [.0] * input_space_dimension

        for i, s in enumerate(c):
            init_val = s[0].mat_i[index[i]]
            fields[s[1]] = init_val + (k if init_val > 0 else -k if init_val < 0 else 0)

        # print(index, fields)
        accum += it[0] * scalar_function(fields, h0, **kwargs)

        it.iternext()

    for s in c:
        accum *= s[0].prefactor(k, h0)

    return accum


class RombergTriangle:
    """
    Do a Romberg triangle to lower (remove) the influence of higher-order derivatives.

    Note: Romberg triangle indices are n, then k (but printing require n to be column, so invert them !)

    :param vals: list of derivatives, the first one calculated with the lowest k, the last one the larger
    :type vals: list
    :param ratio: ratio (a)
    :type ratio: float
    :param r: derivatives to remove. r=1 for backward and forward approximation and 2 for centered approximation.
    :type r: int
    """

    def __init__(self, values, ratio=2., r=1):

        self.ratio = ratio
        self.side = len(values)
        self.romberg_triangle = numpy.zeros((self.side, self.side))
        self.romberg_triangle[:, 0] = values

        # compute Romberg triangle:
        for m in range(1, self.side):
            prev_iteration = self.romberg_triangle[:, m - 1]
            for index in range(self.side - m):
                mult_ = ratio ** (r * m)
                self.romberg_triangle[index, m] = \
                    (mult_ * prev_iteration[index] - prev_iteration[index + 1]) / (mult_ - 1)

    def amplitude_error(self, k=0, m=0):
        """Compute :math:`\\varepsilon_k(m) = H_{k+1,m} - H_{k,m}`

        :param k: amplitude
        :type k: int
        :param m: iteration
        :type m: int
        :rtype: float
        """

        assert_in_domain(m, 0, self.side - 2)
        assert_in_domain(k, 0, (self.side - m - 2))

        return self.romberg_triangle[k + 1, m] - self.romberg_triangle[k, m]

    def iteration_error(self, m=0, k=0):
        """Compute :math:`\\varepsilon_m(k) = H_{k,m+1} - H_{k,m}`

        :param k: amplitude
        :type k: int
        :param m: iteration
        :type m: int
        :rtype: float
        """

        assert_in_domain(m, 0, self.side - 2)
        assert_in_domain(k, 0, (self.side - m - 2))

        return self.romberg_triangle[k, m + 1] - self.romberg_triangle[k, m]

    def __repr__(self):
        return self.romberg_triangle_repr()

    def romberg_triangle_repr(self, with_decoration=False, k_min=0):
        """
        Print the Romberg triangle ...
        """
        r = ''

        if with_decoration:
            r += '-' * 6
            for m in range(self.side):
                r += '-' * 15
            r += '\n'

            r += ' ' * 6
            for m in range(self.side):
                r += '{:^15}'.format('m={}'.format(m))
            r += '\n'
            r += '-' * 6
            for m in range(self.side):
                r += '-' * 15
            r += '\n'

        for k in range(0, self.side):
            if with_decoration:
                r += '{:7}'.format('k={}'.format(k + k_min))
            for m in range(self.side):
                if with_decoration:
                    r += '  '
                if m < (self.side - k):
                    r += '{: 11.5e} '.format(self.romberg_triangle[k, m])
                elif with_decoration:
                    r += '{:12}'.format(' ')
            r += '\n'

        if with_decoration:
            r += '-' * 6
            for m in range(self.side):
                r += '-' * 15
            r += '\n'

        return r

    def find_best_value(self, threshold=1e-5, verbose=False, out=sys.stdout):
        """Find the "best value" in the Romberg triangle

        :param threshold: threshold for maximum iteration error
        :type threshold: float
        :param verbose: get an insight of how the value was chosen
        :type verbose: bool
        :return: a tuple of 3 information: position (tuple), value (float), iteration error (float)
        :rtype: tuple
        """

        if self.side == 1:  # if there is no Triangle, then the best value is the only value
            return (0, 0), self.romberg_triangle[0, 0], 0.0

        m = 0
        current_region = (0, self.side - 1)
        prev_region_error = .0

        if verbose:
            print('- starting with m=0', file=out)

        while True:
            stability_regions = [
                [current_region[0], current_region[0], self.amplitude_error(k=current_region[0], m=m)]]

            current_stability_region = 0
            prev_error = .0

            # dissect the column into stabiliy regions (and look for a value below the threshold, if any)
            for k in range(*current_region):
                amplitude_error = self.amplitude_error(k=k, m=m)
                if math.fabs(amplitude_error) < threshold:
                    if verbose:
                        print('→ value in ({}, {}) have an iteration error lower than threshold, stopping'.format(
                            k, m), file=out)

                    return (k, m), self.romberg_triangle[k, m], 0.0 if m < 1 else self.iteration_error(k=k, m=m)

                if k > current_region[0]:
                    if sgn(amplitude_error) != sgn(prev_error) or math.fabs(amplitude_error) < math.fabs(prev_error):
                        current_stability_region += 1
                        stability_regions.append([k, k, amplitude_error])
                    else:
                        stability_regions[current_stability_region][1] = k

                prev_error = amplitude_error

            stability_regions = sorted(stability_regions, key=lambda a: (a[1] - a[0], -math.fabs(a[2])), reverse=True)
            stability_region = stability_regions[0]

            if verbose:
                print('- stability region(s): {}'.format(
                    ', '.join(['[k={} to k={}, error={:.3e}]'.format(
                        r[0], r[1], r[2]) for r in stability_regions])), file=out)

                print('- select region between k={} and k={} as stability region'.format(
                    *stability_region[:-1]), file=out)

            if m > 0:
                if math.fabs(prev_region_error) < math.fabs(stability_region[2]):
                    if verbose:
                        print('→ amplitude error did not decrease, stopping', file=out)

                    return (stability_region[0], m), \
                        self.romberg_triangle[stability_region[0], m], \
                        self.iteration_error(k=stability_region[0], m=m - 1)

            prev_region_error = stability_region[2]

            if stability_region[1] - stability_region[0] < 2:
                if verbose:
                    if len(stability_regions) == 1:
                        print('→ only one value in the remaining stability region, stopping', file=out)
                    else:
                        print('→ column with amplitude errors of alternating sign, stopping', file=out)

                return (stability_region[0], m + 1), \
                    self.romberg_triangle[stability_region[0], m + 1], \
                    self.iteration_error(k=stability_region[0], m=m)

            m += 1
            current_region = (stability_region[0], stability_region[1])

            if verbose:
                print('- set m={} and continue analysis on this column'.format(m), file=out)
