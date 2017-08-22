===================================
Data File (``qcip_tools.datafile``)
===================================

Implement a key-to-value storage, with a "lazy-loading" behavior.

+ ``read()`` should only slice the file in chunks, for which the information (via `ChunkInformation <#qcip_tools.datafile.ChunkInformation>`_) are stored.
+ When the value of the keyword is requested, the chunck is parsed with ``parse_chunk()`` if not present in the dict of previously parsed chunks.
+ ``get()`` and ``set()`` are used to access and modify data.
+ ``write()`` is finally used to write the file accordingly.

This basic mechanic is set in place by `DataFile <#qcip_tools.datafile.DataFile>`_, from which there is two derived classes:

+ `TextDataFile <#qcip_tools.datafile.TextDataFile>`_: save data in the form of a text file, more or less readable by an human being ;
+ `BinaryDataFile <#qcip_tools.datafile.BinaryDataFile>`_: save data in the form of a binary file (details in the docstring of the class).

API documentation
-----------------

.. automodule:: qcip_tools.datafile
    :members: