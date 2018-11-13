=====================
Installing qcip_tools
=====================

Normal install
--------------

To just use ``qcip_tools`` in your Python projects, simply use pip:

.. code-block:: bash

    pip3 install git+ssh://git@gitlab.unamur.be/chimie/lct/qcip_tools.git@release-vXX

Change ``@release-vXX`` at the end of the line to fetch a given version (listed `in the README <https://gitlab.unamur.be/chimie/lct/qcip_tools/blob/master/README.md>`_).

Note that ``--user`` allow you to install the package without being superuser (see `here <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_).
On the other hand, you can install it in a *virtualenv* (see below).


(*optional*) Patching Dalton
----------------------------

You can download a patch here: :download:`dalton.patch` to outputs responses functions in a better place (``DALTON.PROP`` in the archive), with more digits.
This is generally not important, except if you are looking for accuracy (for numerical differentiation, for example).

To apply this patch you need to recompile Dalton (so you need git, a fortran compiler, and eventual BLAS/LAPACK/MKL).
The following  commands allow you to `get the sources <https://gitlab.com/dalton/dalton>`_, and check if everything is ok:

.. code-block:: bash

  git clone --recursive https://gitlab.com/dalton/dalton.git
  cd dalton/
  git checkout -b special_version origin/release/2016
  # use "origin/master" if you want the latest (maybe buggy) version
  git apply --check /path/to/dalton.patch

If it is the case (no ``error`` in the output of the previous command), you can apply it and compile a new version of Dalton:

.. code-block:: bash

  git am --signoff < /path/to/dalton.patch
  #  You may get "trailing whitespace" error, but it's ok
  ./setup --prefix=/path/to/destination # [--int64]
  cd build
  make
  make install


Installation for contributors
-----------------------------

To contribute to the project, you need to clone the repository:

+ Clone it: ``git clone git@git.pierrebeaujean.net:pierre/qcip_tools.git``.
+ Install scipy dependencies : ``sudo apt-get install libopenblas-dev libatlas-dev build-essential libpng12-dev libfreetype6-dev libpython3.*-dev`` (or something like this).
+ Install pipenv: ``pip3 install pipenv``
+ Install virtualenv and dependencies: ``make init``.

You can launch the tests series with ``make test``