************
Introduction
************

Purpose
=======

*demuxipy* is a program for demultiplexing or "demuxing" DNA sequencing
reads that are identified using sequence tags, which are also known by
other names, including "multiplex identification tags (MIDs)",
"barcodes", "barcoded adapters", etc.  Specifically, demuxipy_ is a
program for demultiplexing `hierarchically tagged`_ DNA sequencing
reads, although the program is flexible with respect to it's processing
regimes.

Features
========

Although similar programs exist for this purpose, including commercial
(Roche, Illumina) and non-commercial (e.g., split_libraries_fastq.py_)
exist, demuxipy_ differs from the other software/methodological approaches 
because:

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
  another.  For additional information see :ref:`Suggested Readings`, below.

- it organizes sequence read data, by tag or other metadata, in a
  relational database (sqlite_) to ease downstream processing, allowing you to 
  query or work by the resulting sequence "group(s)"

- it can take advantage of multiprocessing_ for the parallel parsing and
  error correcting of sequence reads to reduce overall processing time

- it allows the user to specify linkers and hierarchical tag combinations
  in a flexible and easily-edited configuration file.  Once the
  appropriate sequence tags are added to the file, it is only a matter
  of providing the combinations used within a particular run.

- it intelligently creates regular expression groups based on tag
  combinations so that only those combinations within a given run are
  searched for, reducing processing time.

- it can search for potential concatemers within sequence reads

Availability
============

We provide several methods of installing eddittag, see
:ref:`installation` for additional details:

Dependencies
============

* `python   (>= 2.7.x) <http://www.python.org>`_
* `numpy    (>= 1.5.x) <http://numpy.scipy.org>`_
* `seqtools (>= 0.5) <https://github.com/faircloth-lab/seqtools/>`_


.. _demuxipy: https://github.com/faircloth-lab/demuxipy/
.. _split_libraries_fastq.py: http://qiime.org/scripts/split_libraries_fastq.html
.. _`hierarchically tagged`: http://127.0.0.1/
.. _`hierarchical reads`: http://127.0.0.1/
.. _Levenshtein: http://en.wikipedia.org/wiki/Levenshtein_distance
.. _Hamming: http://en.wikipedia.org/wiki/Hamming_distance
.. _multiprocessing: http://en.wikipedia.org/wiki/Multiprocessing
.. _sqlite: http://www.sqlite.org/

