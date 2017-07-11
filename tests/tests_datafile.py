import os
import tempfile
import math

from tests import QcipToolsTestCase
from qcip_tools import datafile


class DataFileTestCase(QcipToolsTestCase):

    def setUp(self):
        self.temporary_file = tempfile.mkstemp()[1]

    def tearDown(self):
        os.unlink(self.temporary_file)

    def test_data_file(self):
        """Test the base class (DataFile)"""

        some_integers = [1, 2, 3, 4, 5, 6, 7]

        # setting and getting a value
        f = datafile.DataFile()
        f.set('integers', 'I', some_integers)

        self.assertIn('integers', f.chunks_information)
        self.assertTrue(f.chunks_information['integers'].modified)
        self.assertEqual(f.chunks_information['integers'].data_type, 'I')
        self.assertIn('integers', f.chunks_parsed)
        self.assertEqual(f.chunks_parsed['integers'], some_integers)

        self.assertTrue('integers' in f)
        self.assertFalse('whatever' in f)

        self.assertEqual(f.get('integers'), some_integers)
        self.assertEqual(f['integers'], some_integers)

        # fool the storage system to prove it is used:
        f.chunks_parsed['integers'] = some_integers[:-1]
        self.assertEqual(f.get('integers'), some_integers[:-1])

        # fire exceptions !
        with self.assertRaises(TypeError):  # non-existing type
            f.set('whatever', 'X', 'test')
        with self.assertRaises(TypeError):  # wrong type
            f.set('integers', 'R', [.1, .2, .3])
        with self.assertRaises(TypeError):  # wrong type
            f.get('integers', 'R')
        with self.assertRaises(TypeError):  # non-existing type
            datafile.transform_string_from_type('test', 'X')
        with self.assertRaises(KeyError):  # non-existing variable
            _ = f['whatever']  # noqa
        with self.assertRaises(KeyError):  # non-existing variable
            _ = f.get('whatever')  # noqa

        # force change variable, even if type is wrong
        some_floats = [.1, .2, .3]
        f.set('integers', 'R', some_floats, force_change_type=True)
        self.assertEqual(f.chunks_information['integers'].data_type, 'R')
        self.assertEqual(f.chunks_information['integers'].data_length, len(some_floats))

    def test_text_data_file(self):
        """Test the behavior of TextDataFile"""

        some_integers = [1, -20002, 3, -405, 574, -6, -700000000000]
        some_floats = [.1, -.2e4, -.3e-3, .4e5, -.5, .6, -.7e-2, math.pi]
        some_text = 'This is a' + ' very' * 20 + ' long text!'

        # saving stuffs:
        f = datafile.TextDataFile()
        f.set('integers', 'I', some_integers)
        f.set('floats', 'R', some_floats)
        f.set('text', 'S', some_text)

        self.assertEqual(os.path.getsize(self.temporary_file), 0)

        with open(self.temporary_file, 'w') as fx:
            f.write(fx)

        self.assertNotEqual(os.path.getsize(self.temporary_file), 0)

        file_ = open(self.temporary_file)
        file_content = file_.read()
        file_.close()

        self.assertTrue('I' + str(len(some_integers)) + ' integers' in file_content)
        self.assertTrue('R' + str(len(some_floats)) + ' floats' in file_content)
        self.assertTrue('S' + str(len(some_text)) + ' text' in file_content)

        # reading stuffs:
        with open(self.temporary_file) as fx:
            g = datafile.TextDataFile()
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 3)

        self.assertIn('integers', g)
        self.assertIn('floats', g)
        self.assertIn('text', g)

        self.assertEqual(g.chunks_information['integers'].data_length, len(some_integers))
        self.assertEqual(g.chunks_information['integers'].data_type, 'I')
        self.assertFalse(g.chunks_information['integers'].modified)

        self.assertEqual(g.chunks_information['floats'].data_length, len(some_floats))
        self.assertEqual(g.chunks_information['floats'].data_type, 'R')
        self.assertFalse(g.chunks_information['floats'].modified)

        self.assertEqual(g.chunks_information['text'].data_length, len(some_text))
        self.assertEqual(g.chunks_information['text'].data_type, 'S')
        self.assertFalse(g.chunks_information['text'].modified)

        self.assertEqual(len(g.chunks_parsed), 0)

        self.assertEqual(g['text'], some_text)
        self.assertEqual(g['integers'], some_integers)
        self.assertEqual(g['floats'], some_floats)

        self.assertEqual(len(g.chunks_parsed), 3)

        # open, modify and save:
        with open(self.temporary_file) as fx:
            g = datafile.TextDataFile()
            g.read(fx)

        some_floats.append(math.e)
        g.set('floats', 'R', some_floats)
        self.assertTrue(g.chunks_information['floats'].modified)
        self.assertEqual(g.chunks_information['floats'].data_length, len(some_floats))

        self.assertFalse(g.chunks_information['text'].modified)
        self.assertFalse(g.chunks_information['integers'].modified)

        with open(self.temporary_file, 'w') as fx:
            g.write(fx)

        file_ = open(self.temporary_file)
        file_content = file_.read()
        file_.close()

        self.assertTrue('I' + str(len(some_integers)) + ' integers' in file_content)
        self.assertTrue('R' + str(len(some_floats)) + ' floats' in file_content)
        self.assertTrue('S' + str(len(some_text)) + ' text' in file_content)

        # and reopen:
        with open(self.temporary_file) as fx:
            g = datafile.TextDataFile()
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 3)

        self.assertEqual(g['text'], some_text)
        self.assertEqual(g['integers'], some_integers)
        self.assertEqual(g['floats'], some_floats)  # ok, good!
