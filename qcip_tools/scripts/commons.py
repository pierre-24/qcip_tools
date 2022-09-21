import sys
import numpy
import math
import csv

from qcip_tools import chemistry_files, quantities, derivatives, derivatives_g, derivatives_e
from qcip_tools.chemistry_files import gamess, gaussian, helpers


def exit_failure(msg, status=1):
    """Write a message in stderr and exits

    :param msg: the msg
    :type msg: str
    :param status: exit status (!=0)
    :type status: int
    """

    sys.stderr.write(msg)
    sys.stderr.write('\n')
    return sys.exit(status)


# STDA UV/VIS output
class StdaDATOutput(chemistry_files.ChemistryFile, chemistry_files.WithIdentificationMixin):
    """STDA "tda.dat" output
    """

    file_type = 'STDA_DAT'

    allowed_keywords = ['UV', 'VELO', 'NM', 'WIDTH', 'SHIFT', 'LFAKTOR', 'RFAKTOR', 'MMASS']
    keywords_with_values = ['WIDTH', 'SHIFT', 'LFAKTOR', 'RFAKTOR', 'MMASS']

    def __init__(self):
        self.num_excitations = 0
        self.excitation_energies = None
        self.transitions = None

        self.keywords = {}

    @classmethod
    def possible_file_extensions(cls):
        return ['dat']

    @classmethod
    def attempt_identification(cls, f):
        """Looks for "RFAKTOR" and "LFAKTOR"
        """

        i = 0
        found_lfaktor = False
        found_lfaktor_l = 0
        found_rfaktor = False
        found_rfaktor_l = 0

        for line in f.readlines():
            if i > 10:
                break

            if 'LFAKTOR' in line:
                found_lfaktor = True
                found_lfaktor_l = i
            elif 'RFAKTOR' in line:
                found_rfaktor = True
                found_rfaktor_l = i
                if found_lfaktor:
                    break

            i += 1

        return found_rfaktor and found_lfaktor and (found_rfaktor_l - found_lfaktor_l) == 2

    def read(self, f):
        lines = f.readlines()
        i = 0
        expects_value = False
        current_keyword = ''

        for line in lines:
            if 'DATXY' in line:
                break

            if expects_value:
                self.keywords[current_keyword] = float(line.strip())
                expects_value = False
            else:
                k = line.strip()
                if k not in self.allowed_keywords:
                    raise Exception('not allowed: {}'.format(k))
                self.keywords[k] = None
                current_keyword = k
                expects_value = k in self.keywords_with_values

            i += 1

        energies = [.0]
        trs = [[.0, .0, .0]]
        eV_to_Ha = quantities.convert(quantities.ureg.eV, quantities.ureg.hartree)

        nstates = 0
        for line in lines[i + 1:]:
            inf = line.split()
            e = float(inf[1]) * eV_to_Ha
            energies.append(e)
            trs.append([.0, .0, math.sqrt(1.5 * float(inf[2]) / e)])
            nstates += 1

        self.excitation_energies = derivatives.Tensor('!', components=energies, nstates=nstates + 1)
        self.transitions = derivatives.Tensor('!F', components=numpy.vstack(trs), nstates=nstates + 1)
        self.num_excitations = nstates


@StdaDATOutput.define_property('excitations')
def stda__out__property__excitations(obj, *args, **kwargs):
    """

    :param obj: object
    :type obj: StdaDATOutput
    :rtype: dict
    """

    return {
        '!': obj.excitation_energies,
        '!F': obj.transitions
    }


helpers.EXTRA_CHEMISTRY_FILES.append(StdaDATOutput)


# GAMESS output
@gamess.Output.define_property('computed_energies')
def gamess__log__property__computed_energies(obj, *args, **kwargs):
    """Get the energies (actually only the HF one)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gamess.Output
    :rtype: dict
    """
    found = obj.search('TOTAL ENERGY', into='PROPERTY EVALUATION')
    if found < 0:
        raise chemistry_files.PropertyNotPresent('computed_energies')

    e = float(obj.lines[found][-20:].strip())
    return {'total': e, 'SCF': e}


@gamess.Output.define_property('geometrical_derivatives')
def gamess__log__property__geometrical_derivatives(obj, *args, **kwargs):
    """Get the geometrical derivatives, if available

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gamess.Output
    :rtype: dict
    """

    if not obj.chunk_exists('NORMAL COORDINATE ANALYSIS'):
        raise chemistry_files.PropertyNotPresent('geometrical_derivatives')

    n = len(obj.molecule)
    spacial_dof = 3 * n
    trans_plus_rot_dof = 5 if obj.molecule.linear() else 6
    geometrical_derivatives = {}

    found = obj.search('ENERGY GRADIENT', into='NORMAL COORDINATE ANALYSIS')
    if found > 0:
        gradient = numpy.zeros((spacial_dof))
        for i in range(n):
            gradient[3 * i:3 * (i + 1)] = [float(a) for a in obj.lines[found + 4 + i].split()[-3:]]

        geometrical_derivatives['G'] = derivatives_g.BaseGeometricalDerivativeTensor(
            spacial_dof=spacial_dof, trans_plus_rot=trans_plus_rot_dof, representation='G', components=gradient)

    found = obj.search('CARTESIAN FORCE CONSTANT MATRIX', line_start=found, into='NORMAL COORDINATE ANALYSIS')
    if found > 0:
        hessian = numpy.zeros((spacial_dof, spacial_dof))
        num_of_pack = int(math.ceil(n / 2))
        shift = found + 6
        for i in range(num_of_pack):
            num_lines = 3 * (n - 2 * i)
            jinf = 1
            for j in range(num_lines):
                jatom = 2 * i
                iatom = j // 3 + jatom
                coord = j % 3

                for jc in range(jinf):
                    e = float(obj.lines[shift + j][20 + 9 * jc: 20 + 9 * jc + 9].strip())
                    hessian[3 * iatom + coord, 3 * jatom + jc] = e

                if jinf != 6:
                    jinf += 1

            shift += num_lines + 4

        hessian_ = hessian.T + hessian
        numpy.fill_diagonal(hessian_, numpy.diag(hessian))

        geometrical_derivatives['GG'] = derivatives_g.BaseGeometricalDerivativeTensor(
            spacial_dof=spacial_dof, trans_plus_rot=trans_plus_rot_dof, representation='GG', components=hessian_)

    return geometrical_derivatives


# CSV Tensor file
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


helpers.EXTRA_CHEMISTRY_FILES.append(CSVTensor)


# gaussian output
@gaussian.Output.define_property('electrical_derivatives')
def gaussian__property__electrical_derivatives(obj, *args, **kwargs):
    """

    :param obj: object
    :type obj: gaussian.Output
    :rtype: dict
    """

    line_start = obj.search(' Property number 1', into=801)
    if line_start < 0:
        return {}

    offset = 0
    data = {}
    while True:
        line = obj.lines[line_start + offset]
        if 'Property' not in line:
            break
        elif 'Alpha' in line:
            frequency = float(line[-10:-2])
            tensor = derivatives_e.PolarisabilityTensor(frequency=frequency if frequency != .0 else 'static')
            for i in range(3):
                cpts = obj.lines[line_start + offset + 2 + i].split()
                tensor.components[i, :] = [float(cpts[j + 1].replace('D', 'e')) for j in range(3)]
            if frequency == .0:
                data.update({'FF': {'static': tensor}})
            else:
                if 'dD' not in data:
                    data['dD'] = {}
                data['dD'].update({frequency: tensor})
            offset += 5
        elif 'Beta' in line:
            frequency = float(line[-10:-2])
            is_SHG = 'w,w,-2w' in line
            tensor = derivatives_e.FirstHyperpolarisabilityTensor(
                input_fields=(0, 0) if frequency == .0 else ((1, 1) if is_SHG else (0, 1)),
                frequency=frequency if frequency != .0 else 'static')

            compressed_tensor = []
            for i in range(18):
                compressed_tensor.append(float(obj.lines[line_start + offset + 2 + i][-15:].replace('D', 'e')))
            if is_SHG:
                tensor.components = gaussian.beta_SHG_from_fchk(compressed_tensor)
                if frequency == .0:
                    data.update({'FFF': {'static': tensor}})
                else:
                    if 'XDD' not in data:
                        data['XDD'] = {}
                    data['XDD'].update({frequency: tensor})
            else:
                tensor.components = gaussian.beta_EOP_from_fchk(compressed_tensor)
                if frequency != .0:
                    if 'dDF' not in data:
                        data['dDF'] = {}
                    data['dDF'].update({frequency: tensor})

            offset += 20
        else:  # probably Gamma
            break

    return data
