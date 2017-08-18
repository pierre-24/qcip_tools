import math

from qcip_tools import atom as qcip_atom, assert_in_domain, ValueOutsideDomain

#: Shell to total angular momentum
SHELL_TO_TAM = {
    'SP': -1,
    'S': 0,
    'P': 1,
    'D': 2,
    'F': 3,
    'G': 4,
    'H': 5,
    'I': 6
}

#: total angular momentum to shell
TAM_TO_SHELL = dict((b, a) for a, b in SHELL_TO_TAM.items())


class NotAValidAngularMomentum(ValueError):
    pass


class Primitive:
    """Gaussian type orbital (**primitive** gaussian function)

    :param contraction_coefficient: contraction coefficient
    :type contraction_coefficient: float
    :param p_coefficient: contraction coefficient for the p function (SP orbital)
    :type p_coefficient: float
    :param exponent: exponent of the gaussian (also called "zeta")
    :type exponent: float
    """

    def __init__(self, exponent=1.0, contraction_coefficient=1.0, p_coefficient=.0):
        self.exponent = exponent
        self.contraction_coefficient = contraction_coefficient
        self.p_coefficient = p_coefficient

    def __add__(self, other):
        if type(other) is not Primitive:
            raise TypeError(other)

        c = Function()
        c.add_primitive(self)
        c.add_primitive(other)

        return c


class Function:
    """Contracted Gaussian type orbital (basis function)
    """

    def __init__(self):
        self.primitives = []

    def __len__(self):
        return len(self.primitives)

    def add_primitive(self, primitive):
        """Add a GTO to the contraction

        :param primitive: the GTO
        :type primitive: Primitive
        """

        self.primitives.append(primitive)

        if len(self.primitives) > 1:
            self.primitives.sort(key=lambda x: x.exponent, reverse=True)

    def is_diffuse(self, exponent_threshold=.1):
        """Test if the function correspond to a diffuse one, based on the following criterions:

        + number of primitive is one
        + this primitive have a contraction coefficient of 1
        + this primitive have a exponent lower than ``exponent_threshold``

        :param exponent_threshold: threshold for the exponent of the primitive
        :type exponent_threshold: float
        :rtype: bool
        """

        if len(self) != 1:
            return False

        p = self.primitives[0]
        if p.contraction_coefficient != 1.:
            return False

        if p.p_coefficient != .0 and p.p_coefficient != 1.:
            return False

        if p.exponent > exponent_threshold:
            return False

        return True


class AtomicBasisSet:
    """Set of basis functions for a given atom.

    :param atomic_number: the atomic number
    :type atomic_number: int
    :param symbol: the atom symbol (starting with an upercase letter)
    :type symbol: str
    """

    def __init__(self, symbol=None, atomic_number=None):
        self.atom = qcip_atom.Atom(symbol=symbol, atomic_number=atomic_number)
        self.basis_functions_per_shell = {}

    @classmethod
    def to_key(cls, key):
        if type(key) is str:
            u = key.upper()
            if u not in SHELL_TO_TAM:
                raise NotAValidAngularMomentum(u)
            return u
        elif type(key) is int:
            if key not in TAM_TO_SHELL:
                NotAValidAngularMomentum(key)
            return TAM_TO_SHELL[key]
        else:
            raise TypeError(key)

    def __contains__(self, item):
        return AtomicBasisSet.to_key(item) in self.basis_functions_per_shell

    def __getitem__(self, item):
        return self.basis_functions_per_shell[AtomicBasisSet.to_key(item)]

    def add_basis_function(self, shell, basis_function):
        """

        :param shell: the shell
        :type shell: str
        :param basis_function: the basis function
        :type basis_function: Function
        """

        if len(basis_function) < 1:
            raise Exception('Empty basis function')

        k = AtomicBasisSet.to_key(shell)
        k_sp = AtomicBasisSet.to_key('SP')
        k_p = AtomicBasisSet.to_key('P')

        if k == k_sp:
            if k_p in self:
                raise Exception('SP basis function with P function already defined')
        if k == k_p:
            if k_sp in self:
                raise Exception('P basis function with SP functions already defined')

        if k not in self.basis_functions_per_shell:
            self.basis_functions_per_shell[k] = []

        self.basis_functions_per_shell[k].append(basis_function)

        if len(self.basis_functions_per_shell[k]) > 1:
            self.basis_functions_per_shell[k].sort(
                key=lambda x: len(x) * 100 + math.log10(x.primitives[0].exponent), reverse=True)

    def contracted_representation(self):
        r = ''
        start = 0

        if 'SP' in self:
            start = 2
            len_s = len(self['S'])
            len_sp = len(self['SP'])
            r += '{}s{}p'.format(len_s + len_sp, len_sp)

        for i in range(start, len(TAM_TO_SHELL) - 1):
            if i in self:
                shell = TAM_TO_SHELL[i]
                r += '{}{}'.format(len(self[shell]), shell.lower())
        return r

    def full_representation(self):
        r = ''
        start = 0

        if 'SP' in self:
            start = 2
            len_s = sum(len(a) for a in self['S'])
            len_sp = sum(len(a) for a in self['SP'])
            r += '{}s{}p'.format(len_s + len_sp, len_sp)

        for i in range(start, len(TAM_TO_SHELL) - 1):
            if i in self:
                shell = TAM_TO_SHELL[i]
                r += '{}{}'.format(sum(len(a) for a in self[shell]), shell.lower())
        return r

    def __repr__(self):
        return '{} [{}|{}]'.format(self.atom.symbol, self.full_representation(), self.contracted_representation())


class BasisSet:
    """Basis set

    :param nickname: name of the basis set
    :type nickname: str
    """

    def __init__(self, nickname=''):
        self.atomic_basis_sets = {}
        self.nickname = nickname

    @classmethod
    def to_key(cls, key):
        if type(key) is str:
            if key not in qcip_atom.Definition:
                raise KeyError(key)
            return qcip_atom.Definition[key][0]
        elif type(key) is int:
            try:
                assert_in_domain(key, 1, 92, 'atomic number')
                return key
            except ValueOutsideDomain as e:
                raise KeyError(e)
        else:
            raise KeyError(key)

    def __contains__(self, item):
        return BasisSet.to_key(item) in self.atomic_basis_sets

    def __getitem__(self, item):
        return self.atomic_basis_sets[BasisSet.to_key(item)]

    def __setitem__(self, key, value):
        if type(value) is not AtomicBasisSet:
            raise TypeError(value)

        t_key = BasisSet.to_key(key)
        e_key = BasisSet.to_key(value.atom.atomic_number)
        if e_key != t_key:
            raise ValueError('key and atom mismatch, Z {} != {}'.format(t_key, e_key))

        self.atomic_basis_sets[t_key] = value

    def add_atomic_basis_set(self, atomic_basis_set, force_replace=True):
        """

        :param atomic_basis_set: the atomic basis set
        :type atomic_basis_set: AtomicBasisSet
        :param force_replace: replace if existing (otherwise, raise an exception if already exists)
        """

        e_key = BasisSet.to_key(atomic_basis_set.atom.atomic_number)

        if not force_replace and e_key in self:
            raise Exception('atomic basis set for atom {} already defined'.format(atomic_basis_set.atom.symbol))

        self[e_key] = atomic_basis_set
