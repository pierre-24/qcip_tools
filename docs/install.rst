=====================
Installing qcip_tools
=====================

Normal install
--------------

To just use ``qcip_tools`` in your Python projects, simply use pip:

.. code-block:: bash

    pip3 install --user git+https://github.com/pierre-24/qcip_tools.git@vXX

Change ``@vXX`` at the end of the line to fetch a given version (listed `in the README <https://github.com/pierre-24/qcip_tools#readme>`_).

Note that ``--user`` allow you to install the package without being superuser (see `here <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_).
You will probably need to add ``$HOME/.local/bin`` to ``$PATH`` for this to work:

.. code-block:: bash

  echo 'PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc

On the other hand, you can install it in a *virtualenv* (see below).


(*optional*) Patching Dalton
----------------------------

.. warning::

    Not checked with recent (> 2018) version of Dalton.

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

To contribute to the project,

+ `Fork it <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
+ Clone your fork: ``git clone git@github.com:<USERNAME>/qcip_tools.git``.
+ Go in it: ``cd qcip_tools``
+ Install pip-tools: ``pip3 install pip-tools``
+ Install virtualenv ``python3 -m venv venv; source venv/bin/activate``
+ Install dependencies: ``make init``.
+ Add upstream: ``git remote add upstream https://github.com/pierre-24/qcip_tools.git``
+ Don't forget to create a separate branch to implement your changes: ``git checkout -b my_branch upstream/dev``.

See `the contribution part <contributing.html>`_.

You can launch the tests series with ``make test``