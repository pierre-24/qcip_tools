from pint import UnitRegistry
from scipy import constants

#: Unit registry. Include definitions for:
#:
#:  - ``bohr``
#:  - ``wavenumber``.
ureg = UnitRegistry()

# Define extra units
# NOTE: most of it is already there :
# - https://github.com/hgrecco/pint/blob/master/pint/constants_en.txt
# - https://github.com/hgrecco/pint/blob/master/pint/default_en.txt
ureg.define('bohr_radius = {} m = bohr'.format(constants.value('Bohr radius')))
ureg.define('wavenumber = 100 * planck_constant * speed_of_light / meter = cm-1')

#: Shortcut for quantities
Q_ = ureg.Quantity


class ValueOutsideDomain(ValueError):
    """Raised when a value is larger or lower that given boundaries"""
    def __init__(self, val, min_, max_):
        super().__init__('{} outside [{};{}]'.format(val, min_, max_))


def assert_in_domain(val, min_, max_):
    if val > max_ or val < min_:
        raise ValueOutsideDomain(val, min_, max_)
