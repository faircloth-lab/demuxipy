Introduction
============

demuxipy_ is a program for the parsing, error-correcting, and tracking
of DNA reads identified using molecular sequence tags.  Specifically
demuxipy_ is for demultiplexing ("demuxing") *hierarchically tagged*
sequence reads - meaning those identified using multiple sequence tags
as plate and well markers, respectively.  demuxipy_ differs from 
previous software/methodological approaches because:

- it demultiplexes `hierarchical reads`_.  Hierarchical tagging vastly
  expands the number of sequence pools that may be mixed during any
  single next-generation sequencing run

- it corrects sequencing errors within tags using fuzzy string matching
  along with sequence tags of Levenshtein distance >= X to identify and
  correct sequencing errors, which are somewhat common on certain
  next-generation sequencing platforms.  The Levenshtein_ distance
  differs from other implementations (i,e.  Hamming_ distance) in that
  the distance represents the number of insertions, deletions and
  substitutions needed to get from one sequence of characters to
  another.  For additional information see Suggested Readings, below.

- it organizes sequence read data, by tag or other metadata, in a
  relational database sqlite_ to ease downstream processing by sequence
  group(s)

- it can take advantage of multiprocessing_ for the parallel parsing and
  error correcting of sequence reads to reduce overall processing time

- it allows the user to specify linkers and hierarchical tag combinations
  in a flexible and easily-edited configuration file.  Once the
  appropriate sequence tags are added to the file, it is only a matter
  of providing the combinations used within a particular run.

- it intelligently creates regular expression groups based on tag
  combinations so that only those combinations within a given run are
  search for.

- it can search for potential concatemers within sequence reads

Dependencies
------------

- python_    (>=2.7)
- numpy_     (>=1.3)
- seqtools_  (>=0.5)

Installation
------------

demuxipy_ requires numpy_ and seqtools_, as mentioned above.  After
installing these dependencies, to install demuxipy:

- from source::

    tar -xzvf ~/your/download/location/demuxipy-*.tar.gz
    python setup.py install

- using easy_install::

    easy_install demuxipy

- using pip::

    pip install demuxipy

Tests
-----

To run the unit tests associated with demuxipy_, you have several
options.  If you have installed python-nose_, then you can do the following:

>>> import demuxipy
>>> demuxipy.test()


Alternatively, you can run:::

    python setup.py test

Or, you can run:::

    python demuxipy/tests/run.py

If you would like to run the tests for seqtools_, please see the README_
file for that library.

Sequence tagging
----------------

We have written a software package to help you design edit-distance
sequence tags named edittag_.  edittag_ is available from github_ and
pypi_.

Running demuxi.py
-----------------

To run demuxipy, prepare a valid configuration file and run:::

    demuxi.py my_configuration_file.conf

Using your data
---------------

After you have demultiplexed your data, you probably want to be able to
use it, rather than just keep it within the database.  You have several
options.

- we provide helper scripts to pull data out of the database and write
  those data to appropriate output files

- the database stores an object for each read containing the read itself
  and the read's quality, and tag data.  You can write a query to grab
  the objects your looking for, unpickle those, and use them directly in
  whatever code/pipeline you'd like.

Frequently asked questions
--------------------------

Please see the FAQ_.

Suggested reading
-----------------

If you'd like to know more about the history of sequence tagging,
sequence tags, error correction, and error-correcting codes, the
following should be of interest:

- `Meyer M, Stenzel U, Myles S, Prüfer K, Hofreiter M (2007) Targeted
  high-throughput sequencing of tagged nucleic acid samples.  Nucleic
  Acids Research 35(15):e97 <http://dx.doi.org/10.1093/nar/gkm566>`_

- `Hamady M, Walker JJ, Harris JK, Gold NJ, Knight R (2008)
  Error-correcting barcoded primers for pyrosequencing hundreds of
  samples in multiplex.  Nature Methods 5 (3):235-237
  <http://dx.doi.org/10.1038/nmeth.1184>`_

- `Binladen J, Gilbert MT, Bollback JP, Panitz F, Bendixen C, Nielsen R,
  Willerslev E (2009) The use of coded PCR primers enables
  high-throughput sequencing of multiple homolog amplification products
  by 454 parallel sequencing.  BMC Bioinformatics
  10:362<http://dx.doi.org/10.1371/journal.pone.0000197>`_

- `Levenshtein VI (1966). Binary codes capable of correcting deletions,
  insertions, and reversals. Soviet Physics Doklady 10:707–10
  <http://sascha.geekheim.de/wp-content/uploads/2006/04/levenshtein.pdf>`_

- `Hamming RW (1950) Error detecting and error correcting codes. Bell
  System Technical Journal 26 (2):147–160
  <http://www.caip.rutgers.edu/~bushnell/dsdwebsite/hamming.pdf>`_


.. _multiprocessing: http://en.wikipedia.org/wiki/Multiprocessing
.. _demuxipy: http://github.com/faircloth-lab/demuxipy
.. _Levenshtein: http://en.wikipedia.org/wiki/Levenshtein_distance
.. _Hamming: http://en.wikipedia.org/wiki/Hamming_distance
.. _edittag: http://github.com/faircloth-lab/edittag
.. _github: http://github.com/faircloth-lab/edittag
.. _pypi: http://pypi.python.org/pypi/edittag
.. _seqtools: http://github.com/faircloth-lab/seqtools/
.. _README: http://github.com/faircloth-lab/seqtools/README.rst
.. _python-nose: http://code.google.com/p/python-nose/
.. _FAQ: https://github.com/faircloth-lab/demuxipy/wiki/faq
.. _sqlite:  http://www.sqlite.org/
