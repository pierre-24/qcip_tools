import unittest
import tempfile
import shutil
import os

from qcip_tools.chemistry_files import gaussian


class GaussianTestCase(unittest.TestCase):
    """Gaussian stuffs"""

    random_input_content = \
        '%mem=3gb\n'\
        '%chk=calculation_0002\n'\
        '%nprocshared=4\n'\
        '#p MP4 6-311+G(d) nosym field=read\n'\
        'scf=(Conver=11,NoVarAcc,MaxCyc=600,vshift=1000) IOP(9/6=600,9/9=11)\n'\
        '\n'\
        'finite field calculation (-1, 0, 0)\n'\
        '\n'\
        '-1 2\n'\
        'O	  0.00000000   0.00000000   0.00000000\n'\
        'H	  0.75698805   0.00000000   0.58632839\n'\
        'H	 -0.75698805   0.00000000   0.58632839\n'\
        '\n'\
        '-0.0004000000	 0.0000000000	 0.0000000000\n'\
        ''  # the blank line at the end is important ;)

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.temp_dir, 'input.com')

        with open(self.input_file, 'w') as f:
            f.write(self.random_input_content)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_input_files(self):
        """Test the behavior of gaussian input file class"""

        text = 'some random text'

        # test with existing input
        gi1 = gaussian.Input()

        with open(self.input_file) as f:
            gi1.read(f)

        self.assertIsNotNone(gi1.molecule)
        self.assertEqual(len(gi1.molecule), 3)
        self.assertEqual(gi1.molecule.multiplicity, 2)
        self.assertEqual(gi1.molecule.charge, -1)
        self.assertEqual(gi1.molecule.number_of_electrons(), 11)  # = 2 (H and H) + 8 (O) + 1 (charge)

        self.assertEqual(len(gi1.options_dict), 3)
        self.assertIn('mem', gi1.options_dict)
        self.assertIn('nprocshared', gi1.options_dict)
        self.assertIn('chk', gi1.options_dict)

        self.assertEqual(len(gi1.other_blocks), 1)  # one block for the electric field

        # write it
        other_file = os.path.join(self.temp_dir, 'other.com')
        gi1.title = text

        with open(other_file, 'w') as f:
            gi1.write(f)

        self.assertTrue(os.path.exists(other_file))

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk=other\n' in content)  # chk has changed based on the file name
            self.assertTrue('%mem=3gb' in content)
            self.assertTrue('%nprocshared=4' in content)
            self.assertTrue(text in content)
            self.assertTrue('{} {}'.format(gi1.molecule.charge, gi1.molecule.multiplicity) in content)
            self.assertTrue(gi1.molecule.output_atoms() in content)
            self.assertTrue(gi1.other_blocks[0][0] in content)

        # write it with a custom chk
        custom_chk = 'another_one'
        with open(other_file, 'w') as f:
            gi1.write(f, chk=custom_chk)

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk={}\n'.format(custom_chk) in content)  # chk is custom

        # write with initial chk
        with open(other_file, 'w') as f:
            gi1.write(f, original_chk=True)

        with open(other_file) as f:
            content = f.read()
            self.assertTrue('%chk=calculation_0002\n' in content)  # chk remains the same as the original
