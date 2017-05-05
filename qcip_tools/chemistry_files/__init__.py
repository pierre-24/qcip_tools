

class File:
    """Purely abstract class that implemnent some basic methods that any child should implement if possible.

    Input methods:

    - ``read()``
    - ``property()``
    - ``get_molecule()``

    Output methods:

    - ``to_string()``
    - ``write()``

    """

    molecule = None
    file_type = 'MISC'
    from_read = False

    def read(self, f):
        raise NotImplementedError

    def write(self, f):
        f.write(self.to_string())

    def to_string(self):
        raise NotImplementedError

    def property(self, property_):
        """

        :param property_: the property
        :type property_: str
        """

        return getattr(self, property_)

    def get_molecule(self):
        """Get the corresponding  molecular geometry. Raises ``NotImplementedError`` if ``self.molecule`` is ``None``.

        :rtype: qcip_tools.molecule.Molecule
        """

        molecule = self.molecule

        if molecule is None:
            raise NotImplementedError

        return molecule
