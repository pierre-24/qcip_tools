================================================
Chemistry files (``qcip_tools.chemistry_files``)
================================================

Read/write files commonly used in quantum chemistry.

Concepts
--------

Human readable text files are widely used in the field of quantum chemistry, but there is large syntaxic differences from one to another (when there IS a syntax).
This is especially true for input files, which follows rigorous syntax, while output files are more or less incoherent. This module implement access to different formats.

+ All classes derive from `InputChemistryFile <#qcip_tools.chemistry_files.InputChemistryFile>`_. Some function to override:

  - ``read()``: read from a file. One should also (if it makes sense) fill ``self.molecule``, which is accessible trough ``get_molecule()`` ;
  - ``possible_file_extensions()`` and ``attempt_recognition()`` (both class methods): helps file recognition.

  Thanks to the dispatcher pattern, derivation from this class allow to define taylor-made accessors: use the ``@Class.define_property(property_name)`` decorator to give the ability to use this callback function with ``property(property_name[, optional kwargs ...])``.
  Callbacks should have the following form: ``callback_func(obj, **kwargs)``, where ``obj`` is the current object. In the code, those callbacks have the form ``property__*(obj, **kwargs)``.

+ Some of them (typically input files of programs) also derive from `InputOutputChemistryFile <#qcip_tools.chemistry_files.InputOutputChemistryFile>`_, which gives the ability to ``write()`` in a new file.
  Note that it is easier to override ``to_string()`` in this case.

API documentation
-----------------

.. list-table::
    :header-rows: 1
    :widths: 35 55 10 10

    * - Class
      - Description
      - Input
      - Ouput
    * - `xyz.File <#qcip_tools.chemistry_files.xyz.File>`_
      - XYZ file (.xyz)
      - Yes
      - Yes
    * - `gaussian.Cube <#qcip_tools.chemistry_files.gaussian.Cube>`_
      - `Gaussian Cube file <http://gaussian.com/cubegen/>`_ (.cub)
      - Yes
      - Yes
    * - `gaussian.FCHK <#qcip_tools.chemistry_files.gaussian.FCHK>`_
      - `Gaussian formated checkpoint file <http://gaussian.com/formchk/>`_ (.fchk)
      - Yes
      - **No**
    * - `gaussian.Input <#qcip_tools.chemistry_files.gaussian.Input>`_
      - Gaussian input file (.com, .gau or .inp)
      - Yes
      - Yes
    * - `gaussian.Output <#qcip_tools.chemistry_files.gaussian.Output>`_
      - Gaussian output (.log or .out)
      - Yes
      - **No**
    * - `dalton.ArchiveOutput <#qcip_tools.chemistry_files.dalton.ArchiveOutput>`_
      - Dalton archive (.tar.gz)
      - Yes
      - **No**
    * - `dalton.Input <#qcip_tools.chemistry_files.dalton.Input>`_
      - Dalton input file (.dal)
      - Yes
      - Yes
    * - `dalton.MoleculeInput <#qcip_tools.chemistry_files.dalton.MoleculeInput>`_
      - Dalton molecule file (.mol)
      - Yes
      - Yes
    * - `dalton.Output <#qcip_tools.chemistry_files.dalton.Output>`_
      - Dalton output file (.out)
      - Yes
      - **No**
    * - `gamess.Input <#qcip_tools.chemistry_files.gamess.Input>`_
      - GAMESS input (.inp)
      - Yes
      - Yes

.. note::

    If one would like to construct the chemistry files from scratch (without using ``read()`` on an existing file) or get additional information from it, the different class members are given in the documentation.


Base
====

.. automodule:: qcip_tools.chemistry_files
    :members:


Helpers (``helpers``)
=====================

Based on ``possible_file_extensions()`` and ``attempt_recognition()``, if defined, the ``helpers`` submodule helps to open any chemistry files.

.. automodule:: qcip_tools.chemistry_files.helpers
    :members:

XYZ files (``xyz``)
===================

.. automodule:: qcip_tools.chemistry_files.xyz
    :members:

Gaussian files (``gaussian``)
=============================

.. automodule:: qcip_tools.chemistry_files.gaussian
    :members:

Dalton files (``dalton``)
=========================

.. automodule:: qcip_tools.chemistry_files.dalton
    :members:

GAMESS files (``gamess``)
=========================

.. automodule:: qcip_tools.chemistry_files.gamess
    :members: