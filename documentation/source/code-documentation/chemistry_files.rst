================================================
Chemistry files (``qcip_tools.chemistry_files``)
================================================

Read/write files commonly used in quantum chemistry.

Concepts
--------

Human readable text files are widely used in the field of quantum chemistry, but there is large syntaxic differences from one to another (when there IS a syntax).
This is especially true for input files, which follows rigorous syntax, while output files are more or less incoherent. This module implement access to different formats.

All classes derive from `ChemistryFile <chemistry-files/base.html#qcip_tools.chemistry_files.ChemistryFile>`_.
You need to override the ``read()`` method to read from a file.

Thanks to the `dispatcher pattern <mixins.html#qcip_tools.mixins.Dispatcher>`_, derivation from this class allow to define taylor-made accessors: use the ``@Class.define_property(property_name)`` decorator to give the ability to use this callback function with ``property(property_name[, optional kwargs ...])``.
Callbacks should have the form ``callback_func(obj, **kwargs)``, where ``obj`` is the current object. In the code, those callbacks have the form ``property__*(obj, **kwargs)``. One can check if a property is available with ``has_property()`` (which does not mean that the file actually contains this property, but that it can be carry out).

There is some *mixins* and derived classes to add functionnalities:

+ `WithMoleculeMixin <chemistry-files/base.html#qcip_tools.chemistry_files.WithMoleculeMixin>`_ adds ``self.molecule``, which is accessible trough ``get_molecule()``.
  It also provides a ``Class.from_molecule(molecule, *args, **kwargs)`` class method to ease creation of a file from scratch (if defined in subclass).
+ `WithIdentificationMixin <chemistry-files/base.html#qcip_tools.chemistry_files.WithIdentificationMixin>`_ adds ``Class.possible_file_extensions()`` and ``Class.attempt_identification()`` to allow file recognition by `helpers <chemistry-files/helpers.html#helpers-helpers>`_.
    * ``possible_file_extensions()`` helps to avoid some possibilities if the extension is trusted. The helper expects a list of possible extentions (without the dot).
    * ``attempt_identification()`` actually make the recognition of the file: given ``f`` (the file open in read mode), the helper function expects ``True`` or ``False`` if it is a file corresponding to a given objects.
      This recognition is usually based on patterns or words (like the name of the program) that frequently appears.
+ `WithOutputMixin <chemistry-files/base.html#qcip_tools.chemistry_files.WithOutputMixin>`_, which gives the ability to ``write()`` in a new file.
  Note that it is easier to override ``to_string()`` in this case.
+ `ChemistryLogFile <chemistry-files/base.html#qcip_tools.chemistry_files.ChemistryLogFile>`_ (derived class) keeps the whole file in memory (in ``self.lines``) and allow to divide the file in ``chunks`` to ease their manipulation.
  The ``chunk_exists()`` function check if a given chunk is present in the file.
  The mixins defines a ``apply_function()`` method to apply a given function (with ``apply_over_list()``, `see here <chemistry-files/base.html#qcip_tools.chemistry_files.apply_over_list>`_) to the lines (to carry out information for example).
  Based on that, it also defines a ``search()`` function, which returns the line where the information was found.


.. note::

    If one would like to construct the chemistry files from scratch (without using ``read()`` on an existing file) or get additional information from it, the different class members are given in the documentation.
    Also check if ``from_molecule()`` is defined.

API documentation
-----------------

.. toctree::
    :maxdepth: 2

    chemistry-files/base.rst
    chemistry-files/helpers.rst
    chemistry-files/xyz.rst
    chemistry-files/gaussian.rst
    chemistry-files/dalton.rst
    chemistry-files/gamess.rst
    chemistry-files/chemistry_datafile.rst

More precisely,

.. list-table::
    :header-rows: 1
    :widths: 35 45 10 10 10

    * - Class
      - Description
      - Identifier
      - Input
      - Ouput
    * - `xyz.File <./chemistry-files/xyz.html#qcip_tools.chemistry_files.xyz.File>`_
      - XYZ file (.xyz)
      - ``XYZ``
      - Yes
      - Yes
    * - `gaussian.BasisSet <./chemistry-files/gaussian.html#qcip_tools.chemistry_files.gaussian.BasisSet>`_
      - Gaussian basis set (.gbs ?)
      - ``GAUSSIAN_BS``
      - Yes
      - Yes
    * - `gaussian.Cube <./chemistry-files/gaussian.html#qcip_tools.chemistry_files.gaussian.Cube>`_
      - `Gaussian Cube file <http://gaussian.com/cubegen/>`_ (.cub)
      - ``GAUSSIAN_CUBE``
      - Yes
      - Yes
    * - `gaussian.FCHK <./chemistry-files/gaussian.html#qcip_tools.chemistry_files.gaussian.FCHK>`_
      - `Gaussian formated checkpoint file <http://gaussian.com/formchk/>`_ (.fchk)
      - ``GAUSSIAN_FCHK``
      - Yes
      - **No**
    * - `gaussian.Input <./chemistry-files/gaussian.html#qcip_tools.chemistry_files.gaussian.Input>`_
      - Gaussian input file (.com, .gau or .inp)
      - ``GAUSSIAN_INP``
      - Yes
      - Yes
    * - `gaussian.Output <./chemistry-files/gaussian.html#qcip_tools.chemistry_files.gaussian.Output>`_
      - Gaussian output (.log or .out)
      - ``GAUSSIAN_LOG``
      - Yes
      - **No**
    * - `dalton.ArchiveOutput <./chemistry-files/dalton.html#qcip_tools.chemistry_files.dalton.ArchiveOutput>`_
      - Dalton archive (.tar.gz)
      -  ``DALTON_ARCHIVE``
      - Yes
      - **No**
    * - `dalton.Input <./chemistry-files/dalton.html#qcip_tools.chemistry_files.dalton.Input>`_
      - Dalton input file (.dal)
      - ``DALTON_DAL``
      - Yes
      - Yes
    * - `dalton.MoleculeInput <./chemistry-files/dalton.html#qcip_tools.chemistry_files.dalton.MoleculeInput>`_
      - Dalton molecule file (.mol)
      - ``DALTON_MOL``
      - Yes
      - Yes
    * - `dalton.Output <./chemistry-files/dalton.html#qcip_tools.chemistry_files.dalton.Output>`_
      - Dalton output file (.out)
      - ``DALTON_LOG``
      - Yes
      - **No**
    * - `gamess.Input <./chemistry-files/gamess.html#qcip_tools.chemistry_files.gamess.Input>`_
      - GAMESS input (.inp)
      - ``GAMESS_INP``
      - Yes
      - Yes
    * - `gamess.Output <./chemistry-files/gamess.html#qcip_tools.chemistry_files.gamess.Output>`_
      - GAMESS output (.log or .out)
      - ``GAMESS_LOG``
      - Yes
      - **No**
    * - `chemistry_datafile.ChemistryDataFile <./chemistry-files/chemistry_datafile.html#qcip_tools.chemistry_files.chemistry_datafile.ChemistryDataFile>`_
      - QCIP chemistry data file (.chdf)
      - ``QCIP_CDF``
      - Yes
      - Yes