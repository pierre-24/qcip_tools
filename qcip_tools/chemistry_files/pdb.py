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

        for l in f.readlines():
            if count > 50:
                break

            count += 1
            if l[:6].strip() in keywords:
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

        for l in lines:
            if l[:6].strip() in ['ATOM', 'HETATM']:
                symbol = l[77:79].strip()
                atm = atom.Atom(symbol=symbol)
                atm.position = [float(x) for x in [l[30:38], l[38:46], l[46:54]]]

                atm.extra['pdb_name'] = l[13:17]
                atm.extra['pdb_resName'] = l[17:20]
                atm.extra['pdb_chainId'] = l[21]
                atm.extra['pdb_resSeq'] = l[22:26]
                atm.extra['pdb_iCode'] = l[26]
                atm.extra['pdb_occupancy'] = float(l[54:60])
                atm.extra['pdb_tempFactor'] = float(l[60:66])
                atm.extra['pdb_charge'] = l[78:80]

                self.molecule.insert(atm)
                i += 1

            elif l[:3] == 'TER':
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

            r += '{:6}{:5}  {:4}{:4}{}{:3}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.3f}{:6.3f}           {:2}{:2}\n'.format(
                'ATOM',
                i,
                fx(a, 'pdb_name', 'UKN' + str(i)),
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
