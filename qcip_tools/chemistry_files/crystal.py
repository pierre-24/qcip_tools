from qcip_tools.chemistry_files import ChemistryLogFile, WithIdentificationMixin, PropertyNotPresent, WithMoleculeMixin
from qcip_tools import molecule, derivatives_e


class Output(ChemistryLogFile, WithIdentificationMixin, WithMoleculeMixin):
    """Output of Crystal

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.lines``: the lines of the file (``list`` of ``str``)

    """

    #: The identifier
    file_type = 'CRYSTAL_LOG'

    def __init__(self):
        self.molecule = molecule.Molecule()

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

        print(found_countries, found_universities, found_crystal)

        return found_countries > 5 and found_universities > 5 and found_crystal > 2

    def read(self, f):
        """There does not seem to be any obvious sections"""
        super().read(f)


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
        if 'POLARIZABILITY (ALPHA)' in l:  # static or dynamic polarizability
            freq = obj.lines[line_cphf + i + 4][2:8]  # TODO: maybe somewhere else with additional digits?
            tensor = derivatives_e.PolarisabilityTensor(
                frequency=derivatives_e.convert_frequency_from_string(freq + 'eV'))

            for j in range(6):
                info = obj.lines[line_cphf + i + 4 + j].split()
                component = tuple(translate_c[k] for k in info[1])
                tensor.components[component] = tensor.components[component[1], component[0]] = float(info[2])

            if tensor.frequency == 0:
                electrical_derivatives['FF'] = {'static': tensor}
            else:
                if 'dD' not in electrical_derivatives:
                    electrical_derivatives['dD'] = {}
                electrical_derivatives['dD'][tensor.frequency] = tensor

        elif 'FIRST HYPERPOLARIZABILITY (BETA)' in l:
            tensor = derivatives_e.FirstHyperpolarisabilityTensor(input_fields=(0, 0))

            for j in range(10):
                info = obj.lines[line_cphf + i + 6 + j].split()
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
                tensor.components[component] = float(info[1])

            if is_pockel:
                if 'dDF' not in electrical_derivatives:
                    electrical_derivatives['dDF'] = {}
                electrical_derivatives['dDF'][tensor.frequency] = tensor
            else:
                if 'XDD' not in electrical_derivatives:
                    electrical_derivatives['XDD'] = {}
                electrical_derivatives['XDD'][tensor.frequency] = tensor

    return electrical_derivatives
