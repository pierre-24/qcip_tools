============
Contributing
============

Please note that the code is not actually developed on the git server of the University of Namur (which only contains the releases) but on a personal protected git server (with CI activated and properly configured).
Feel free to ask access if needed.

You first need to `install <./install.html>`_ if you wan to contribute to the code.

Design rules
------------

+ The code is written in Python 3, and follows the (in)famous `PEP-8 <http://legacy.python.org/dev/peps/pep-0008/>`_. You can check it by running ``make lint``, which launch the ``flake`` utility.
+ Codes and comments are written in english.
+ The code is documented using docstrings. The docstrings must contains the basic description of the function, as well as a description of the paramters (with the ``:type`` instruction, please).
+ The code is tested. You can launch the test series by using ``make test``. Every functionality should be provided with at least one unit test. Even though *coverage* is not monitored (yet?), it should be ok.
+ The package is documented. You can generate this documentation by using ``make doc``. Non-basic stuffs should be explained in this documentation. Don't forget to cite some articles or website if needed.
+ Before reimplementing something, please consider if there is no library that already exists to do the job.

Workflow
--------

Adapted from the (in)famous `Git flow <http://nvie.com/posts/a-successful-git-branching-model/>`_.

+ Development is mad in ``dev`` branch, while ``master`` contains the production version.
+ Functionalities are added through merge request (MR) in the ``dev`` branch.
+ Theses merge requests should be unitary, and include unit test(s) and documentation if needed. The test suite must succeed for the merge request to be accepted.
+ At some (random) points, ``dev`` will be merged into ``master`` to create a new version.

Licence
-------

This code belong to me, `Pierre Beaujean <pierre.beaujean@unamur.be>`_, and to the `University of Namur <https://www.unamur.be>`_ since it is developed and used in the frame of my PhD thesis.

A note about units
------------------

Unless mentioned in the docstring, atomic units (Hartree) are used. The two exceptions are:

+ Angstrom for distances, since it is widely used in day-to-day chemistry, and widely used by quantum chemistry packages (Gaussian and GAMESS, to name a few).
+ Atomic units for masses, for the same kind of reason.

For unit conversions, the `Pint library <http://pint.readthedocs.io>`_ is used, which already define a whole bunch of `units <https://github.com/hgrecco/pint/blob/master/pint/default_en.txt>`_ and `constants <https://github.com/hgrecco/pint/blob/master/pint/constants_en.txt>`_.
Also, take a look in the `code documentation for quantities <./code-documentation/quantities.html>`_ for extra units and conversion factors.