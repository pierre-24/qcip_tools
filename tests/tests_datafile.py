import os
import tempfile
import math
import numpy

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

    def test_binary_data_file(self):
        """Test the binary version. Special care is taken with UTF-8."""

        some_integers = [1, -20002, 3, -405, 574, -6, -700000000000]
        some_floats = [.1, -.2e4, -.3e-3, .4e5, -.5, .6, -.7e-2, math.pi]
        some_text = 'This is a' + ' very' * 20 + ' long text!'

        utf8_stuff = '_ħħłł'

        # saving stuffs:
        f = datafile.BinaryDataFile()
        f.set('integers', 'I', some_integers)
        f.set('floats', 'R', some_floats)
        f.set('text', 'S', some_text)
        f.set('text' + utf8_stuff, 'S', some_text + utf8_stuff)

        self.assertEqual(os.path.getsize(self.temporary_file), 0)

        with open(self.temporary_file, 'wb') as fx:
            f.write(fx)

        self.assertNotEqual(os.path.getsize(self.temporary_file), 0)

        # reading stuffs:
        g = datafile.BinaryDataFile()

        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 4)

        self.assertIn('integers', g)
        self.assertIn('floats', g)
        self.assertIn('text', g)
        self.assertIn('text' + utf8_stuff, g)

        self.assertEqual(g.chunks_information['integers'].data_length, len(some_integers))
        self.assertEqual(g.chunks_information['integers'].data_type, 'I')
        self.assertFalse(g.chunks_information['integers'].modified)

        self.assertEqual(g.chunks_information['floats'].data_length, len(some_floats))
        self.assertEqual(g.chunks_information['floats'].data_type, 'R')
        self.assertFalse(g.chunks_information['floats'].modified)

        self.assertEqual(g.chunks_information['text'].data_length, len(some_text.encode('utf-8')))
        self.assertEqual(g.chunks_information['text'].data_type, 'S')
        self.assertFalse(g.chunks_information['text'].modified)

        self.assertEqual(
            g.chunks_information['text' + utf8_stuff].data_length, len((some_text + utf8_stuff).encode('utf-8')))

        self.assertEqual(g.chunks_information['text' + utf8_stuff].data_type, 'S')
        self.assertFalse(g.chunks_information['text' + utf8_stuff].modified)

        self.assertEqual(len(g.chunks_parsed), 0)

        self.assertEqual(g['text'], some_text)
        self.assertEqual(g['integers'], some_integers)
        self.assertEqual(g['floats'], some_floats)
        self.assertEqual(g['text' + utf8_stuff], some_text + utf8_stuff)

        self.assertEqual(len(g.chunks_parsed), 4)

        # open, modify and save:
        g = datafile.BinaryDataFile()
        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        some_floats.append(math.e)
        g.set('floats', 'R', some_floats)
        self.assertTrue(g.chunks_information['floats'].modified)
        self.assertEqual(g.chunks_information['floats'].data_length, len(some_floats))

        self.assertFalse(g.chunks_information['text'].modified)
        self.assertFalse(g.chunks_information['integers'].modified)
        self.assertFalse(g.chunks_information['text' + utf8_stuff].modified)

        with open(self.temporary_file, 'wb') as fx:
            g.write(fx)

        # and reopen:
        g = datafile.BinaryDataFile()
        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 4)

        self.assertEqual(g['text'], some_text)
        self.assertEqual(g['text' + utf8_stuff], some_text + utf8_stuff)
        self.assertEqual(g['integers'], some_integers)
        self.assertEqual(g['floats'], some_floats)  # ok, good!

    def test_array_type(self):
        """A type added a bit after"""

        array = numpy.eye(3)

        # TextDataFile:
        ft = datafile.TextDataFile()
        ft.set('array', 'A', array)

        self.assertEqual(ft.chunks_information['array'].data_length, 9)
        self.assertEqual(ft.chunks_information['array'].shape, (3, 3))

        with open(self.temporary_file, 'w') as fx:
            ft.write(fx)

        file_ = open(self.temporary_file)
        file_content = file_.read()
        file_.close()
        self.assertTrue('A3x3 array' in file_content)

        g = datafile.TextDataFile()
        with open(self.temporary_file) as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 1)

        self.assertEqual(ft.chunks_information['array'].data_length, 9)
        self.assertEqual(ft.chunks_information['array'].shape, (3, 3))
        self.assertArrayAlmostEqual(g['array'], array)  # ok.

        g = datafile.TextDataFile()  # open-read-close, then reopen
        with open(self.temporary_file) as fx:
            g.read(fx)

        g.set('test', 'I', [1, 2])

        with open(self.temporary_file, 'w') as fx:
            g.write(fx)

        g = datafile.TextDataFile()
        with open(self.temporary_file) as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 2)
        self.assertArrayAlmostEqual(g['test'], [1, 2])

        self.assertEqual(ft.chunks_information['array'].data_length, 9)
        self.assertEqual(ft.chunks_information['array'].shape, (3, 3))
        self.assertArrayAlmostEqual(g['array'], array)  # ok.

        # BinaryDataFile:
        fb = datafile.BinaryDataFile()
        fb.set('array', 'A', array)

        self.assertEqual(fb.chunks_information['array'].data_length, 9)
        self.assertEqual(fb.chunks_information['array'].shape, (3, 3))

        with open(self.temporary_file, 'wb') as fx:
            fb.write(fx)

        g = datafile.BinaryDataFile()
        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 1)

        self.assertEqual(g.chunks_information['array'].data_length, 9)
        self.assertEqual(g.chunks_information['array'].shape, (3, 3))
        self.assertArrayAlmostEqual(g['array'], array)  # ok.

        g = datafile.BinaryDataFile()  # open-read-close, then reopen
        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        g.set('test', 'I', [1, 2])

        with open(self.temporary_file, 'wb') as fx:
            g.write(fx)

        g = datafile.BinaryDataFile()
        with open(self.temporary_file, 'rb') as fx:
            g.read(fx)

        self.assertEqual(len(g.chunks_information), 2)
        self.assertArrayAlmostEqual(g['test'], [1, 2])

        self.assertEqual(g.chunks_information['array'].data_length, 9)
        self.assertEqual(g.chunks_information['array'].shape, (3, 3))
        self.assertArrayAlmostEqual(g['array'], array)  # ok.
