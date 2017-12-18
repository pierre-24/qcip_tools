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
# -
ureg.define('bohr_radius = {} m = bohr'.format(constants.value('Bohr radius')))
ureg.define('wavenumber = 100 * planck_constant * speed_of_light / meter = cm-1')
ureg.define('atomic_unit_of_time = {} s'.format(constants.value('atomic unit of time')))

#: Shortcut for quantities
Q_ = ureg.Quantity


def convert(from_, to_, value=1.0):
    """Conversion factor from one unit to the other

    :param from_: base unit
    :type from_: pint.unit.Unit
    :param to_: converted unit
    :type to_: pint.unit.Unit
    :param value: value to convert
    :type value: float
    """

    return (value * from_).to(to_).magnitude


# widely used conversion factors
#: Convert bohr to angstrom
AuToAngstrom = convert(ureg.bohr, ureg.angstrom)
#: Convert atomic mass units to electron mass
AMUToElectronMass = convert(ureg.atomic_mass_unit, ureg.electron_mass)
