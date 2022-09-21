import csv

import numpy

from qcip_tools import chemistry_files, derivatives_e, derivatives_g, derivatives


class CSVTensor(chemistry_files.ChemistryFile, chemistry_files.WithIdentificationMixin):
    """Special CSV input. Assume that the first cell contains ``CSV_TENSOR_IN``, and that the first row contains
    different information:

    + In second column, the derivative type (i.e. ``FFF``, ``GG``, ...) ;
    + In third column: spacial DOF, if needed ;
    + In fourth column: frequency, if needed ;
    + In fifth column: number of (excited) states, if needed.
    + For the latter rows, values of the first column does not matter.

    Then, the rest correspond to the values of the last dimension of the tensor (F/D/X = 3,
    G = spacial DOF).
    The number of row is equal to the product of the dimensions of the tensor (except the last one).

    New tensor starts with ``CSV_TENSOR_IN``.

    For example,

    .. code:: text

        CSV_TENSOR_IN,G,6,,,
        1,1,1,1,1,1
        CSV_TENSOR_IN,GG,6
        1,0,0,0,0,0
        0,2,0,0,0,0
        0,0,3,0,0,0
        0,0,0,4,0,0
        0,0,0,0,5,0
        0,0,0,0,0,6

    """

    file_type = 'CSV_IN'
    keyword = 'CSV_TENSOR_IN'

    def __init__(self):
        self.tensors = {}

    @classmethod
    def possible_file_extensions(cls):
        return ['csv']

    @classmethod
    def attempt_identification(cls, f):
        """
        Check for "CSV_TENSOR_IN"

        :rtype: bool
        """

        x = f.read(len(CSVTensor.keyword))
        return x[:len(CSVTensor.keyword)] == CSVTensor.keyword

    @staticmethod
    def create_empty_tensor(deriv, spacial_dof=None, frequency=None, nstates=None):

        if deriv == 'F':
            return derivatives_e.ElectricDipole()
        if deriv in ['FF', 'dD']:
            return derivatives_e.PolarisabilityTensor(frequency=frequency)
        if deriv in ['FFF', 'dDF', 'XDD']:
            return derivatives_e.FirstHyperpolarisabilityTensor(frequency=frequency)
        if deriv in ['FFFF', 'dFFD', 'dDFF', 'XDDF', 'dDDd', 'XDDD']:
            return derivatives_e.SecondHyperpolarizabilityTensor(frequency=frequency)
        if deriv in ['G', 'GG']:
            if spacial_dof is None:
                raise Exception('spacial DOF should be given for {}'.format(deriv))
            return derivatives_g.BaseGeometricalDerivativeTensor(representation=deriv, spacial_dof=int(spacial_dof))
        else:
            return derivatives.Tensor(
                deriv,
                spacial_dof=int(spacial_dof) if spacial_dof not in [None, ''] else None,
                frequency=frequency if frequency not in [None, ''] else None,
                nstates=int(nstates) if nstates not in [None, ''] else None
            )

    def read(self, f, reader_kwargs={'delimiter': ','}):
        reader = csv.reader(f, **reader_kwargs)

        vstack = []
        t_info = []
        t_all = []
        first = True

        for r in reader:
            if r[0] == 'CSV_TENSOR_IN':
                if first:
                    first = False
                else:
                    t_all.append((t_info, vstack))

                t_info = r[1:]
                vstack = []
                continue

            try:
                vstack.append(list(float(a) for a in r if a != ''))
            except ValueError:
                raise Exception('float point conversion error at line `{}`'.format(','.join(r)))

        t_all.append((t_info, vstack))  # add last read

        for t_info, vstack in t_all:
            t = CSVTensor.create_empty_tensor(*t_info[:4])
            r = t.representation
            rs = str(r)

            if len(vstack) != r.dimension() / r.shape()[-1]:
                raise Exception('first dimensions mismatch for {}: expected {} lines, got {}'.format(
                    rs, r.dimension() / r.shape()[-1], len(vstack)))

            t.components[:] = numpy.vstack(vstack).reshape(r.shape())

            if derivatives.is_electrical(r) and not derivatives.is_excitation(r):
                if rs not in self.tensors:
                    self.tensors[rs] = {}
                self.tensors[rs][t.frequency] = t
            else:
                self.tensors[rs] = t

    @staticmethod
    def output_tensor(tensor: derivatives.Tensor):
        info = [
            'CSV_TENSOR_IN',
            tensor.representation.representation(),
            tensor.spacial_dof if tensor.spacial_dof is not None else '',
            tensor.frequency if tensor.frequency is not None else '',
            tensor.representation.nstates if tensor.representation.nstates is not None else ''
        ]

        out = [info]
        shape = tensor.components.shape

        for line in tensor.components.reshape((int(numpy.product(shape[:-1])), shape[-1])):
            out.append(list(line))

        return out


@CSVTensor.define_property('geometrical_derivatives')
def csvtensor__property__geometrical_derivatives(obj, *args, **kwargs):
    """Get the geometrical derivatives, if available

    :param obj: object
    :type obj: CSVTensor
    :rtype: dict
    """

    d = {}

    for t in obj.tensors:
        if not derivatives.is_electrical(t) and derivatives.is_geometrical(t):
            d[t] = obj.tensors[t]

    return d


@CSVTensor.define_property('electrical_derivatives')
def csvtensor__property__electrical_derivatives(obj, *args, **kwargs):
    """Get the geometrical derivatives, if available

    :param obj: object
    :type obj: CSVTensor
    :rtype: dict
    """

    d = {}

    for t in obj.tensors:
        if derivatives.is_electrical(t) and not derivatives.is_geometrical(t):
            d[t] = obj.tensors[t]

    return d


@CSVTensor.define_property('excitations')
def csvtensor__property__excitations(obj, *args, **kwargs):
    """

    :param obj: object
    :type obj: CSVTensor
    :rtype: dict
    """

    d = {}

    for t in obj.tensors:
        if t[0] in derivatives.SYMBOL_EXCITATIONS:
            d[t] = obj.tensors[t]

    return d
