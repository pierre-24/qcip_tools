============
Contributing
============
You first need to `install <./install.html>`_ if you wan to contribute to the code.

Design rules
------------

+ The code is written in Python 3, and follows the (in)famous `PEP-8 <http://legacy.python.org/dev/peps/pep-0008/>`_. You can check it by running ``make lint``, which launch the ``flake`` utility.
+ Codes and comments are written in english.
+ The code is documented using docstrings and Sphinx. The docstrings must contains the basic description of the function, as well as a description of the paramters (with the ``:type`` instruction, please).
+ The code is tested. You can launch the test series by using ``make test``.
  Every functionality should be provided with at least one unit test.
  Every script should be provided with at least one unit test.
  You may need test files to do so, but try to make them small (say, don't use d-aug-cc-pVDZ while STO-3G could do the job).
+ The package is documented. You can generate this documentation by using ``make doc``. Non-basic stuffs should be explained in this documentation. Don't forget to cite some articles or website if needed.
+ Before reimplementing something, please consider if there is no library that already exists to do the job.

Workflow
--------

Adapted from the (in)famous `Git flow <http://nvie.com/posts/a-successful-git-branching-model/>`_.

+ Development is mad in ``dev`` branch, while ``master`` contains the production version (and is protected from edition).
+ Functionality are added through merge request (MR) in the ``dev`` branch. Do not work in ``dev`` directly, but create a new branch (``git checkout -b my_branch origin/dev``).
+ Theses merge requests should be unitary, and include unit test(s) and documentation if needed. The test suite must succeed for the merge request to be accepted.
+ At some (random) points, ``dev`` will be merged by the maintainer into ``master`` to create a new version, with a tag of the form ``release-vXX``.

.. note::

    Since ``qcip_tools`` now rely on `pip-tools <https://github.com/jazzband/pip-tools>`_, the workflow is currently the following :

    1. Normal installation use ``pip-sync && pip install -e .`` (``make init``)
    2. To update the dependencies from upstream, ``pip-sync``  (``make sync``).
    3. To update the ``requirements.txt`` (and thus the actual version of the dependencies), a **specific** merge request is done, with the result of ``pip-compile requirements/requirements.in -o requirements.txt`` (followed by ``make sync`` on the dev's machine).

A note about units
------------------

Unless mentioned in the docstring, atomic units (Hartree) are used. The two exceptions are:

+ Angstrom for distances, since it is widely used in day-to-day chemistry, and widely used by quantum chemistry packages (Gaussian and GAMESS, to name a few).
+ Atomic units for masses, for the same kind of reason.
+ Kelvin for temperature (there is no conversion to atomic unit for the temperature, I think)
+ Pascal for the pressure (force atomic unit could be used to define pressure atomic unit, but it is a bit overthinking for thermochemistry)

For unit conversions, the `Pint library <http://pint.readthedocs.io>`_ is used, which already define a whole bunch of `units <https://github.com/hgrecco/pint/blob/master/pint/default_en.txt>`_ and `constants <https://github.com/hgrecco/pint/blob/master/pint/constants_en.txt>`_.
Also, take a look in the `code documentation for quantities <./code-documentation/quantities.html>`_ for extra units and conversion factors.