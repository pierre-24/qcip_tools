=====================
Installing qcip_tools
=====================

+ To just use ``qcip_tools`` in your Python projects, simply use pip:

  .. code-block:: bash

    pip install git+ssh://git@gitlab.unamur.be/pierre.beaujean/qcip_tools.git@release-vXX

  Change ``@release-vXX`` at the end of the line to fetch a given version (listed `in the README <https://gitlab.unamur.be/pierre.beaujean/qcip_tools/blob/master/README.md>`_).

+ To contribute to the project, you need to clone the repository:

  + Clone it: ``git clone git@git.pierrebeaujean.net:pierre/qcip_tools.git``.
  + Install scipy dependencies : ``sudo apt-get install libopenblas-dev libatlas-dev build-essential libpng12-dev libfreetype6-dev libpython3.*-dev`` (or something like this).
  + Create virtualenv and activate it: ``virtualenv venv --python=python3 && source venv/bin/activate``.
  + Install dependencies : ``make install-dependencies``.
  + Finally, "install" the pakage: ``pip install --editable .``

  You can launch the tests series with ``python3 setup.py test``