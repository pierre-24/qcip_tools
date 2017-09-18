===================================
Data File (``qcip_tools.datafile``)
===================================

Implement a key-to-value storage, with a "lazy-loading" behavior for some types of data (lists, basically, but see below).

+ ``read()`` should only slice the file in chunks, for which the information (via `ChunkInformation <#qcip_tools.datafile.ChunkInformation>`_) are stored.
+ When the value of the keyword is requested, the chunck is parsed with ``parse_chunk()`` if not present in the dict of previously parsed chunks.
+ ``get()`` and ``set()`` are used to access and modify data.
+ ``write()`` is finally used to write the file accordingly.

This basic mechanic is set in place by `DataFile <#qcip_tools.datafile.DataFile>`_, from which there is two derived classes:

+ `TextDataFile <#qcip_tools.datafile.TextDataFile>`_: save data in the form of a text file, more or less readable by an human being (and a little bit of data loss in the case of floatting numbers) ;
+ `BinaryDataFile <#qcip_tools.datafile.BinaryDataFile>`_: save data in the form of a binary file (details in the docstring of the class).

Four data types exists:

+ String (``S``): a string ;
+ Integer list (``I``): list of integers ;
+ Float list (``R``): list of floats ;
+ Numpy arrays (``A``): array of floats (the shape of the array is also stored).

It is therefore easy to store and retrieve some data:

.. code-block:: python

    import numpy
    from qcip_tools import datafile

    some_integers = [1, -20002, 3, -405, 574, -6, -700000000000]
    some_floats = [.1, -.2e4, -.3e-3, .4e5, -.5, .6, -.7e-2, math.pi]
    some_text = 'This is a' + ' very' * 20 + ' long text!'
    some_array = numpy.eye(3)

    # saving:
    f = datafile.TextDataFile()  # or datafile.BinaryDataFile()
    f.set('integers', 'I', some_integers)
    f.set('floats', 'R', some_floats)
    f.set('text', 'S', some_text)
    f.set('array', 'A', some_array)

    with open(file_path, 'w') as fx:
        f.write(fx)

    # reading:
    g = datafile.TextDataFile()
    with open(file_path) as fx:
        g.read(fx)

    # >>> print(g['integers'])
    # [1, -20002, 3, -405, 574, -6, -700000000000]

API documentation
-----------------

.. automodule:: qcip_tools.datafile
    :members: