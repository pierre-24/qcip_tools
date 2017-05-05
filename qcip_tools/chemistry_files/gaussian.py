import os

from qcip_tools import molecule, atom
from qcip_tools.chemistry_files import File as qcip_File


class Input(qcip_File):
    """Gaussian input file. I/O class."""

    file_type = 'GAUSSIAN_INPUT'

    def __init__(self):

        self.molecule = molecule.Molecule()

        self.raw_blocks = []
        self.other_blocks = []
        self.title = ''
        self.options = []
        self.options_dict = {}
        self.input_card = []

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        self.from_read = True
        lines = f.readlines()

        temp = []
        for l in lines:
            if len(temp) > 0 and l.strip() == '':
                self.raw_blocks.append(temp)
                temp = []
            else:
                temp.append(l.strip())

        # add last
        if len(temp) > 0:
            self.raw_blocks.append(temp)
        if len(self.raw_blocks) < 3:
            raise Exception('The file contains less than 3 raw_blocks')

        # try to understand first block
        for line in self.raw_blocks[0]:
            if line[0] == '%':
                self.options.append(line[1:])
                opt = line[1:].split('=')
                self.options_dict[opt[0].lower()] = opt[1]
            else:
                self.input_card.append(line)

        # read title, charge and multiplicity:
        self.title = ' '.join(self.raw_blocks[1])
        charge_and_multiplicity = [int(i) for i in self.raw_blocks[2][0].strip().split()]
        self.molecule.charge = charge_and_multiplicity[0]
        self.molecule.multiplicity = charge_and_multiplicity[1]

        # read coordinates:
        coordinates = self.raw_blocks[2][1:]
        for atm in coordinates:
            s = atm.split()
            if len(s) != 4:
                raise Exception('Wrong atom definition "' + atm + '"')
            else:
                coords = [float(i) for i in s[1:]]
                atm_ = s[0]
                if atm_.isnumeric():
                    atm_ = atom.AtomicNumberToSymbol[int(atm_)]
                self.molecule.insert(atom.Atom(symbol=atm_, position=coords))

        # then store "other blocks"
        if len(self.raw_blocks) > 3:
            self.other_blocks = self.raw_blocks[3:]

    def to_string(self, chk='', original_chk=False):
        """

        :param chk: new chk name if needed
        :type chk: str
        :param original_chk: use original chk name
        :type original_chk: bool
        :rtype: str
        """

        rstr = ''

        # first, options:
        for key in self.options_dict:
            if not original_chk and key == 'chk':
                rstr += '%{}={}\n'.format(key, chk)
            else:
                rstr += '%{}={}\n'.format(key, self.options_dict[key])

        # then, input card and title
        rstr += '\n'.join(self.input_card) + '\n\n'
        rstr += self.title + '\n\n'

        # write charge and multiplicity
        rstr += str(self.molecule.charge) + ' ' + str(self.molecule.multiplicity) + '\n'

        # write coordinates
        rstr += self.molecule.output_atoms()
        rstr += '\n'

        # write other blocks :
        for i, block in enumerate(self.other_blocks):
            rstr += '\n'.join(block)
            rstr += '\n\n'
        return rstr

    def write(self, f, chk=None, original_chk=False):
        """

        :param f: File
        :type f: file
        """

        if chk is None:
            chk = os.path.basename('.'.join(f.name.split('.')[:-1]))

        f.write(self.to_string(chk=chk, original_chk=original_chk))
