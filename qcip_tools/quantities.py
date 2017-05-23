from . import ureg


def convert(from_, to_):
    """Conversion factor from one unit to the other

    :param from_: base unit
    :type from_: pint.unit.Unit
    :param to_: converted unit
    :type to_: pint.unit.Unit
    """

    return (1.0 * from_).to(to_).magnitude

# widely used conversion factors
#: Convert bohr to angstrom
AuToAngstrom = convert(ureg.bohr, ureg.angstrom)
#: Convert atomic mass units to electron mass
AMUToElectronMass = convert(ureg.atomic_mass_unit, ureg.electron_mass)
