import math

import numpy

from qcip_tools import chemistry_files, quantities, derivatives


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
