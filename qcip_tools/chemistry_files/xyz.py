from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, FormatError, \
    WithIdentificationMixin


class XYZFormatError(FormatError):
    pass


class File(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    """The (in)famous XYZ file

    .. container:: class-members

        + ``self.molecule``: the molecule (``qcip_tools.molecule.Molecule``)
        + ``self.title``: title of the molecule (``str``)
    """

    #: The identifier
    file_type = 'XYZ'

    def __init__(self):

        self.molecule = molecule.Molecule()
        self.title = ''

    @classmethod
    def possible_file_extensions(cls):
        return ['xyz']

    @classmethod
    def attempt_identification(cls, f):
        """An XYZ is kinda simple: it starts with a number, then contains series of 4 things (the 3 last being floats)
        """

        try:
            number = int(f.readline())
        except ValueError:
            return False

        f.readline()

        for i in range(3):
            if i >= number:
                break

            info = f.readline().split()
            if len(info) != 4:
                return False

            for j in range(3):
                try:
                    float(info[j + 1])
                except ValueError:
                    return False

        return True

    @classmethod
    def from_molecule(cls, molecule, title='', *args, **kwargs):
        """Create a file from molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :param title: title of the run
        :type title: str
        :rtype: qcip_tools.chemistry_files.xyz.File
        """

        obj = super().from_molecule(molecule, *args, **kwargs)
        obj.title = title
        return obj

    def read(self, f, trust_number_of_atoms=True):

        self.from_read = True

        try:
            number_of_atoms = int(f.readline())
        except ValueError:
            raise XYZFormatError('XYZ must start with number of atoms')

        self.title = f.readline().strip()
        count = 0

        while True:
            info = f.readline().split()

            if len(info) == 0:
                break

            if len(info) != 4:
                raise XYZFormatError('wrong number of data in {}'.format(info))

            position = [float(a) for a in info[1:]]
            if info[0].isnumeric():
                self.molecule.insert(atom.Atom(atomic_number=int(info[0]), position=position))
            else:
                self.molecule.insert(atom.Atom(symbol=info[0], position=position))

            count += 1

        if trust_number_of_atoms and count != number_of_atoms:
            raise XYZFormatError('the actual number of atom ({}) is no coherent (!={})'.format(count, number_of_atoms))

    def to_string(self):
        r = str(len(self.molecule)) + '\n'
        r += self.title.replace('\n', '') + '\n'

        for a in self.molecule:
            r += '{}  {: .8f} {: .8f} {: .8f}\n'.format(a.symbol, *a.position)

        return r
