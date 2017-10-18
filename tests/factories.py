import numpy

from qcip_tools import derivatives_e


class FakeElectricDipole(derivatives_e.ElectricDipole):
    def __init__(self, factor=1., **kwargs):
        super().__init__(**kwargs)
        self.components = numpy.array([0., 1. * factor, 2. * factor])


class FakePolarizabilityTensor(derivatives_e.PolarisabilityTensor):
    def __init__(self, isotropy_factor=1., anisotropy_factor=1., **kwargs):
        super().__init__(**kwargs)

        dummy = .0

        for i in self.representation.smart_iterator():
            dummy += 1
            val = dummy * (isotropy_factor if i[0] == i[1] else anisotropy_factor)
            for j in self.representation.inverse_smart_iterator(i):
                self.components[j] = val


class FakeFirstHyperpolarizabilityTensor(derivatives_e.FirstHyperpolarisabilityTensor):
    def __init__(self, dipolar_factor=1., octupolar_factor=1., **kwargs):
        super().__init__(**kwargs)

        dummy = .0

        for i in self.representation.smart_iterator():
            dummy += 1
            val = dummy

            if i[0] != i[1] and i[1] != i[2] and i[0] != i[2]:
                val *= octupolar_factor
            else:
                val *= dipolar_factor

            for j in self.representation.inverse_smart_iterator(i):
                self.components[j] = val
