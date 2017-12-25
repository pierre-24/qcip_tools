"""
Quantum Chemistry In Python (QCIP) tools package.
"""

__name__ = 'qcip_tools'
__version__ = '0.4.3'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


class ValueOutsideDomain(ValueError):
    """Raised when a value is larger or lower that given boundaries"""
    def __init__(self, val, min_, max_, help_=''):
        super().__init__('{}{} outside [{};{}]'.format('' if not help_ else help_ + ' ', val, min_, max_))


def assert_in_domain(val, min_, max_, help_=''):
    if val > max_ or val < min_:
        raise ValueOutsideDomain(val, min_, max_, help_='')
