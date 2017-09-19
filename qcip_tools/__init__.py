"""
Quantum Chemistry In Python (QCIP) tools package.
"""

from pint import UnitRegistry
from scipy import constants

__name__ = 'qcip_tools'
__version__ = '0.3a'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

#: Unit registry. Include definitions for:
#:
#:  - ``bohr``
#:  - ``wavenumber``.
ureg = UnitRegistry()

# Define extra units
# NOTE: most of it is already there :
# - https://github.com/hgrecco/pint/blob/master/pint/constants_en.txt
# -
ureg.define('bohr_radius = {} m = bohr'.format(constants.value('Bohr radius')))
ureg.define('wavenumber = 100 * planck_constant * speed_of_light / meter = cm-1')

#: Shortcut for quantities
Q_ = ureg.Quantity


class ValueOutsideDomain(ValueError):
    """Raised when a value is larger or lower that given boundaries"""
    def __init__(self, val, min_, max_, help_=''):
        super().__init__('{}{} outside [{};{}]'.format('' if not help_ else help_ + ' ', val, min_, max_))


def assert_in_domain(val, min_, max_, help_=''):
    if val > max_ or val < min_:
        raise ValueOutsideDomain(val, min_, max_, help_='')
