=====================
Installing qcip_tools
=====================

For a development environment,

+ Clone it: ``git clone git@git.pierrebeaujean.net:pierre/qcip_tools.git``.
+ Install scipy dependencies : ``sudo apt-get install libopenblas-dev libatlas-dev build-essential libpng12-dev libfreetype6-dev`` (or something like this).
+ Create virtualenv and activate it: ``virtualenv venv --python=python3 && source venv/bin/activate``.
+ Install dependencies : ``pip install --upgrade -r requirements.txt`` or ``make install-dependencies`` (for a development environment).
+ Finally, "install" the pakage: ``pip install --editable .``

You can launch the tests series with ``python3 setup.py test-back``