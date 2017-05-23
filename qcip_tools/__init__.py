from pint import UnitRegistry
from scipy import constants

ureg = UnitRegistry()
Q_ = ureg.Quantity

# Define extra units
# NOTE: most of it is already there :
# - https://github.com/hgrecco/pint/blob/master/pint/constants_en.txt
# - https://github.com/hgrecco/pint/blob/master/pint/default_en.txt

ureg.define('bohr_radius = {} m = bohr'.format(constants.value('Bohr radius')))
ureg.define('wavenumber = 100 * planck_constant * speed_of_light / meter = cm-1')
