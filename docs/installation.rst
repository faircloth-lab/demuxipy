.. _installation:

*************
Installation
*************

demuxipy_ has several dependencies.  You need to install both of:

- numpy_
- seqtools_

Installing numpy
================

You first need to install numpy_ (>= 1.5).  Barring platform-specific
options, you can accomplish this by running:

.. code-block:: bash

    easy_install numpy

Sometimes on OSX, this is problematic.  There is a `binary installer
<http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1-py2.6-python.org-macosx10.3.dmg/download>`_
that you can use in this case.

Installing seqtools
===================

You also need to The easiest option to get up and running quickly is to
use easy_install:

.. code-block:: bash

    easy_install seqtools

Pip should also work:

.. code-block:: bash

    pip install seqtools

And, if none of those work, you can also install seqtools from the tarball:

.. code-block:: bash

    wget http://pypi.python.org/packages/source/s/seqtools/seqtools-0.5.tar.gz
    tar -xzvf seqtools-0.5.tar.gz
    cd seqtools-*
    python setup.py install

Installing `demuxipy`
=====================

Finally, you need to install demuxipy, itself.  Similar to the above,
you can use easy_install:

.. code-block:: bash

    easy_install demuxipy

Pip should also work

.. code-block:: bash

    pip install demuxipy

And, if none of those work, you can also install demuxipy from the tarball:

.. code-block:: bash

    wget localhost (link not up yet)
    tar -xzvf seqtools-0.5.tar.gz
    cd seqtools-*
    python setup.py install

.. _numpy: http://numpy.scipy.org
.. _seqtools:  https://github.com/faircloth-lab/seqtools/
.. _demuxipy:  http://github.com/faircloth-lab/demuxipy/
