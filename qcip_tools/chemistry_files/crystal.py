from qcip_tools.chemistry_files import ChemistryLogFile, WithIdentificationMixin, PropertyNotPresent, \
    WithMoleculeMixin, FormatError
from qcip_tools import molecule, derivatives_e, atom
import numpy
import qcip_tools.math as qmath


class Output(ChemistryLogFile, WithIdentificationMixin, WithMoleculeMixin):
    """Output of Crystal

    .. container:: class-members

        + ``self.molecule``: the primitive cell (``qcip_tools.molecule.Molecule``)
        + ``self.lattice_vectors``, ``self.lattice_abc`` and ``self.lattice_albega``: cell vectors, abc lengths
          (in angstrom) and alpha-beta-gamma (in degree).
        + ``self.lines``: the lines of the file (``list`` of ``str``)

    """

    #: The identifier
    file_type = 'CRYSTAL_LOG'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.lattice_vectors = numpy.eye(3)

        self.lattice_abc = numpy.ones(3)
        self.lattice_albega = numpy.pi / 2 * numpy.ones(3)

    @classmethod
    def attempt_identification(cls, f):
        """A Crystal output, like DALTON, contains a bit of countries, universities and sometimes its own name.
        """

        count = 0
        found_universities = 0
        found_crystal = 0
        found_countries = 0
        countries = ['ITALY', 'FRANCE', 'GERMANY', 'SPAIN', 'UK', 'USA']

        for line in f.readlines():
            count += 1

            if 'CRYSTAL' in line:
                found_crystal += 1
            elif 'UNIVERSITY' in line or 'UNIVERSITA' in line:
                found_universities += 1
                for c in countries:
                    if c in line:
                        found_countries += 1
                        break
            if count > 100:
                break

        return found_countries > 5 and found_universities > 5 and found_crystal > 2

    def read(self, f):
        """There does not seem to be any obvious sectioning"""
        super().read(f)

        line_primitive = self.search('CARTESIAN COORDINATES - PRIMITIVE CELL')
        if line_primitive < 0:
            raise FormatError('not primitive cell')

        # lattice vectors
        l_vecs = [[float(a) for a in line.split()] for line in self.lines[line_primitive - 5:line_primitive - 2]]
        self.lattice_vectors = numpy.array(l_vecs).reshape((3, 3))

        self.lattice_abc = numpy.linalg.norm(self.lattice_vectors, axis=1)
        self.lattice_albega = numpy.array([
            qmath.angle(self.lattice_vectors[1], numpy.zeros(3), self.lattice_vectors[2]),
            qmath.angle(self.lattice_vectors[0], numpy.zeros(3), self.lattice_vectors[2]),
            qmath.angle(self.lattice_vectors[0], numpy.zeros(3), self.lattice_vectors[1]),
        ])

        # geometry of the primitive cell
        for line in self.lines[line_primitive + 4:]:
            if line == '\n':
                break

            atom_def = line.split()

            self.molecule.insert(
                atom.Atom(symbol=atom_def[2][0] + atom_def[2][1:].lower(), position=[float(a) for a in atom_def[3:]]))


@Output.define_property('electrical_derivatives')
def crystal__output__property__electrical_derivatives(obj, *args, **kwargs):
    """Get electrical derivatives in Dalton archives. Returns a dictionary of dictionaries:

    .. code-block:: text

        + "F"
            + static : ElectricDipole
        + "FF":
            + static : PolarizabilityTensor
        + "FD":
            + 0.042 : PolarizabilityTensor
            + ...
        + "FDD":
            + static : FirstHyperpolarizabilityTensor
            + ...
        + ...

    .. note::

        Only works with the CPKS module

    :param obj: object
    :type obj: qcip_tools.chemistry_files.crystal.Output
    :rtype: dict
    """

    electrical_derivatives = {}
    translate_c = {'X': 0, 'Y': 1, 'Z': 2}

    # look for the CPKS module
    line_cphf = obj.search('COUPLED-PERTURBED KOHN-SHAM CALCULATION (CPKS)')
    if line_cphf < 0:
        raise PropertyNotPresent('electrical_derivatives')

    # if so, gather info
    for i, l in enumerate(obj.lines[line_cphf:]):
        if ' POLARIZABILITY TENSOR' in l or 'POLARIZABILITY (ALPHA)' in l:  # static or dynamic polarizability
            is_slab = 'TENSOR' in l
            freq = obj.lines[line_cphf + i + (4 if not is_slab else 3)][2:8].strip()  # TODO: additional digits?
            tensor = derivatives_e.PolarisabilityTensor(
                frequency=derivatives_e.convert_frequency_from_string(freq + 'eV'))

            for j in range(6):
                info = obj.lines[line_cphf + i + (4 if not is_slab else 3) + j].split()
                if len(info) == 0:
                    break

                component = tuple(translate_c[k] for k in info[1])
                tensor.components[component] = tensor.components[component[1], component[0]] = float(info[2])

            if tensor.frequency == 0:
                electrical_derivatives['FF'] = {'static': tensor}
            else:
                if 'dD' not in electrical_derivatives:
                    electrical_derivatives['dD'] = {}
                electrical_derivatives['dD'][tensor.frequency] = tensor

        elif 'FIRST HYPERPOLARIZABILITY' in l:
            is_slab = 'TENSOR' in l
            tensor = derivatives_e.FirstHyperpolarisabilityTensor(input_fields=(0, 0))

            for j in range(10):
                line = obj.lines[line_cphf + i + (6 if not is_slab else 5) + j]
                if '***' in line:
                    break

                info = line.split()
                component = tuple(translate_c[k] for k in info[0])

                for c in tensor.representation.inverse_smart_iterator(component):
                    tensor.components[c] = float(info[1])

            electrical_derivatives['FFF'] = {'static': tensor}

        elif 'SHG Effect beta(-2omega;' in l or 'Pockels Effect beta(' in l:
            is_pockel = 'Pockels' in l

            freq = l[-12:-1]
            tensor = derivatives_e.FirstHyperpolarisabilityTensor(
                input_fields=(0, 1) if is_pockel else (1, 1),
                frequency=derivatives_e.convert_frequency_from_string(freq))

            line_beta = obj.search('BETA', line_start=line_cphf + i)
            for j in range(27):
                info = obj.lines[line_beta + 3 + j].split()
                component = tuple(translate_c[k] for k in info[0])
                # NOTE: it seems that the component is reversed!!
                tensor.components[component] = -float(info[1])

            if is_pockel:
                if 'dDF' not in electrical_derivatives:
                    electrical_derivatives['dDF'] = {}
                electrical_derivatives['dDF'][tensor.frequency] = tensor
            else:
                if 'XDD' not in electrical_derivatives:
                    electrical_derivatives['XDD'] = {}
                electrical_derivatives['XDD'][tensor.frequency] = tensor

    return electrical_derivatives
