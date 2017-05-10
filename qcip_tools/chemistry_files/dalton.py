import tarfile

from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import ChemistryFile as qcip_ChemistryFile


AuToAngstrom = 0.52917165


class MoleculeInput(qcip_ChemistryFile):
    """Dalton mol input file.

    **I/O class.**"""

    file_type = 'DALTON_MOL'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.title = ''
        self.basis_set = ''

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        lines = f.readlines()
        self.from_read = True

        if len(lines) < 6:
            raise Exception('something wrong with dalton .mol')

        self.basis_set = lines[1].strip()
        self.title = (lines[2] + '\n' + lines[3]).strip()

        in_angstrom = 'angstrom' in lines[4].lower()
        atomic_number = .0

        for line in lines[5:]:
            if 'charge=' in line.lower():
                atomic_number = int(float(line[7:line.lower().find('atoms')]))
                continue

            content = line.split()
            if len(content) != 4:
                continue

            self.molecule.insert(atom.Atom(
                atomic_number=atomic_number,
                position=[float(a) * (1 if in_angstrom else AuToAngstrom) for a in content[1:]])
            )

    def to_string(self, in_angstrom=True, nosym=False, group_atoms=False):
        """

        :param in_angstrom: gives the atomic coordinates in angstrom
        :type: in_angstrom: bool
        :param nosym: specify "nosymmetry"
        :type nosym: bool
        :param group_atoms: group of the same type together (the order may be lost)
        :type group_atoms: bool
        :rtype: str
        """
        r = 'BASIS\n{}\n'.format(self.basis_set)
        r += self.title

        if self.title.find('\n') == -1:
            r += '\n\n'

        if not group_atoms:
            r += 'Atomtypes={}{}{}\n'.format(
                len(self.molecule), ' Angstrom' if in_angstrom else '', ' Nosymmetry' if nosym else '')

            for a in self.molecule:
                r += 'Charge={:.1f} Atoms={}\n'.format(a.atomic_number, 1)
                r += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(
                    a.symbol, *[p * (1 if in_angstrom else AuToAngstrom) for p in a.position])
        else:
            r += 'Atomtypes={}{}{}\n'.format(
                len(self.molecule.symbols_contained),
                ' Angstrom' if in_angstrom else '',
                ' Nosymmetry' if nosym else ''
            )

            for symbol in self.molecule.symbols_contained:
                atms = self.molecule.atoms(symbol_in=[symbol])
                r += 'Charge={:.1f} Atoms={}\n'.format(atom.Definition[symbol][0], len(atms))

                for i in atms:
                    a = self.molecule[i]
                    r += '{:3} {:16.9f} {:16.9f} {:16.9f}\n'.format(
                        a.symbol, *[p * (1 if in_angstrom else AuToAngstrom) for p in a.position])

        return r

    def write(self, f, in_angstrom=True, nosym=False, group_atoms=False):
        """

        :param f: File
        :type f: file
        :param in_angstrom: gives the atomic coordinates in angstrom
        :type: in_angstrom: bool
        :param nosym: specify "nosymmetry"
        :type nosym: bool
        :param group_atoms: group of the same type together (the order may be lost)
        :type group_atoms: bool
        """

        f.write(self.to_string(in_angstrom, nosym, group_atoms))


class ArchiveOutput(qcip_ChemistryFile):
    """Archive output of Dalton. Contains lots of information for who can extract them.

    **Input only class.**"""

    def __init__(self):
        self.tar_file = None
        self.molecule = molecule.Molecule()

    def read(self, f):
        """Expects ``f`` to be open on binary mode.

        .. note::

            In Dalton, the molecule is given in a old, weird format, called "molplt". Informations about it is
            given in https://brettbode.github.io/wxmacmolplt/Manual_pages/MacMolPlt_Files.html ...

        :param f: File
        :type f: file
        """

        self.from_read = True
        self.tar_file = tarfile.open(fileobj=f)

        # read molecule:
        mol_file = self.get_file('DALTON.MOL')

        lines = mol_file.readlines()

        info = lines[1].split()
        number_of_atoms = int(info[1])
        number_of_kinds = int(info[3])

        for i in range(number_of_atoms):
            a = lines[4 + number_of_kinds + i].split()
            self.molecule.insert(atom.Atom(symbol=str(a[0], 'utf-8'), position=[float(e) for e in a[1:]]))

    def get_file(self, name):
        """Get a ``io.BufferedReader`` from a filename, if exists. Note that you get ``bytes`` out of that kind of
        objects, so you need to ``decode()`` them.

        :param name: file name
        :type name: str
        :rtype: io.BufferedReader
        """

        if self.tar_file is None:
            raise IOError('no archive open')

        if name in self.tar_file.getnames():
            return self.tar_file.extractfile(self.tar_file.getmember(name))
        else:
            raise FileNotFoundError(name)
