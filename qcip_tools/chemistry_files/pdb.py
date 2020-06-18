from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin


class File(ChemistryFile, WithOutputMixin, WithMoleculeMixin, WithIdentificationMixin):
    #: The identifier
    file_type = 'PDB_FILE'

    def __init__(self):
        self.molecule = molecule.Molecule()
        self.TERs = []

    @classmethod
    def possible_file_extensions(cls):
        return ['pdb']

    @classmethod
    def attempt_identification(cls, f):
        """A PDB file has a lot of specific keywords
        """

        count = 0
        found = 0
        keywords = ['HEADER', 'TITLE', 'COMPND', 'REMARK', 'MODEL', 'ATOM', 'KEYWDS', 'AUTHOR']

        for line in f.readlines():
            if count > 50:
                break

            count += 1
            if line[:6].strip() in keywords:
                found += 1

        return found > 30

    @classmethod
    def from_molecule(cls, molecule, *args, **kwargs):
        """Create a file from molecule

        :param molecule: the molecule
        :type molecule: qcip_tools.molecule.Molecule
        :rtype: qcip_tools.chemistry_files.pdb.File
        """

        obj = super().from_molecule(molecule, *args, **kwargs)
        return obj

    def read(self, f):
        """

        Format is described in
        - http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        - http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
        - http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#TER

        :param f: File
        :type f: file
        """

        lines = f.readlines()
        self.from_read = True

        i = 0

        for line in lines:
            if line[:6].strip() in ['ATOM', 'HETATM']:
                symbol = line[76:78].strip()
                atm = atom.Atom(symbol=symbol)
                atm.position = [float(x) for x in [line[30:38], line[38:46], line[46:54]]]

                atm.extra['pdb_name'] = line[12:16]
                atm.extra['pdb_altLoc'] = line[16]
                atm.extra['pdb_resName'] = line[17:20]
                atm.extra['pdb_chainId'] = line[21]
                atm.extra['pdb_resSeq'] = line[22:26]
                atm.extra['pdb_iCode'] = line[26]
                atm.extra['pdb_occupancy'] = float(line[54:60])
                atm.extra['pdb_tempFactor'] = float(line[60:66])
                atm.extra['pdb_charge'] = line[78:80]

                self.molecule.insert(atm)
                i += 1

            elif line[:3] == 'TER':
                self.TERs.append(i)

    def to_string(self):
        """

        :rtype: str
        """

        r = ''
        i = 1

        def fx(at, n, alt=''):
            return at.extra[n] if n in at.extra else alt

        for a in self.molecule:

            if i - 1 in self.TERs:
                r += 'TER\n'

            r += '{:6}{:5} {:4}{:1}{:4}{}{:3}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           {:2}{:2}\n'.format(
                'ATOM',
                i,
                fx(a, 'pdb_name', 'UKN' + str(i)),
                fx(a, 'pdb_altLoc', ''),
                fx(a, 'pdb_resName', 'UKR' + str(i)),
                fx(a, 'pdb_chainId', 1),
                fx(a, 'pdb_resSeq', i),
                fx(a, 'pdb_iCode', ''),
                a.position[0],
                a.position[1],
                a.position[2],
                fx(a, 'pdb_occupancy', 1.),
                fx(a, 'pdb_tempFactor', 0.0),
                a.symbol,
                fx(a, 'pdb_charge', '')
            )
            i += 1

        r += 'END'

        return r
