================================================
Chemistry files (``qcip_tools.chemistry_files``)
================================================

Read/write files commonly used in quantum chemistry.

Concepts
--------

Human readable text files are widely used in the field of quantum chemistry, but there is large syntaxic differences from one to another (when there IS a syntax).
This is especially true for input files, which follows rigorous syntax, while output files are more or less incoherent. This module implement access to different formats.

All classes derive from `ChemistryFile <#qcip_tools.chemistry_files.ChemistryFile>`_. Some functions to override:

- ``read()``: read from a file.
- ``Class.possible_file_extensions()`` and ``Class.attempt_recognition()`` (both class methods): helps file recognition by `helpers <#helpers-helpers>`_.

Thanks to the `dispatcher pattern <mixins.html#qcip_tools.mixins.Dispatcher>`_, derivation from this class allow to define taylor-made accessors: use the ``@Class.define_property(property_name)`` decorator to give the ability to use this callback function with ``property(property_name[, optional kwargs ...])``.
Callbacks should have the following form: ``callback_func(obj, **kwargs)``, where ``obj`` is the current object. In the code, those callbacks have the form ``property__*(obj, **kwargs)``. One can check if a property is available with ``has_property()``.

There is some *mixins* to add functionnalities:

+ `WithMolecule <#qcip_tools.chemistry_files.WithOutput>`_ adds ``self.molecule``, which is accessible trough ``get_molecule()``.
  It also provides a ``Class.from_molecule(molecule, *args, **kwargs)`` class method to ease creation of a file from scratch (if defined in subclass).
+ `WithOutput <#qcip_tools.chemistry_files.WithOutput>`_, which gives the ability to ``write()`` in a new file.
  Note that it is easier to override ``to_string()`` in this case.

Specific design remarks:

+ Input file classes (``Input*``) are I/O and generally intepret data as much as they can.
+ Log file classes (``Output*``), keep the whole file in memory (trough ``self.lines``) and are generally divided in portions to ease their manipulation.
  They define a ``search()`` function, which returns the line where the information was found.
  They also define a ``apply_function()`` method to apply a given function (trough `apply_over_list() <#qcip_tools.chemistry_files.apply_over_list>`_) to the lines (to carry out information for example).

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
    * - `gamess.Output <#qcip_tools.chemistry_files.gamess.Output>`_
      - GAMESS output (.log or .out)
      - Yes
      - **No**


.. note::

    If one would like to construct the chemistry files from scratch (without using ``read()`` on an existing file) or get additional information from it, the different class members are given in the documentation.
    Also check if ``from_molecule()`` is defined.


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