.. _setting-up-a-configuration-file:

*******************************
Setting up a configuration file
*******************************

demuxipy_ has lots of options, so I thought it would be easiest to setup
the options the program needs in a configuration file.  Below, I'll walk
you though the various options.  Note that the example configuration
file (demuxipy-example.conf) and the test configuration file 
(demuxipy/tests/test-data/demuxi-test.conf) are heavily annotated,
giving most of the text presented below.  The configuration file follows
the conventions of the Python_ `ConfigParser module`_.

Below, I'll describe the options used in each portion of the file by
using the ``[section headers]`` of the file as section headers below.

[Multiprocessing]
=================

.. code-block:: python

    [Multiprocessing]
    MULTIPROCESSING     = True
    PROCESSORS          = Auto

The multiprocessing section allows you so specify whether (or not) you
wish to run multiprocessing on the data.  By "multiprocessing", I mean
that the data will be split into chunks, those chunks will be passed to
separate physical or virtual CPUs on your computer, your computer will
process those chunks in parallel, and demuxipy_ will capture the
resulting data and store it in a database for you.

One downside to multiprocessing is that not everyone is using a
computer with multiple physical processors or multiple virtual cores (or
both).  That said if you are working with massively parallel sequencing
data, chances are you have reasonable processing power at your
fingertips.  If not, you can always `rent some`_.

The other downside to multiprocessing is that for small jobs, the
overhead involved with sending data to multiple cores may add more work
to be done that you actually need to do.  This means that small jobs may
run faster on a single computing core.  To turn off multiprocessing:

.. code-block:: python

    [Multiprocessing]
    MULTIPROCESSING     = False

Generally speaking, you will want to turn on multiprocessing

.. code-block:: python

    MULTIPROCESSING     = True

When you do, you can choose how many cores to use

.. code-block:: python

    PROCESSORS          = 2

or let the program choose automatically:

.. code-block:: python

    PROCESSORS          = Auto

When you choose automatically, the program will use n - 2 cores to
process your data and one core to put data into the database.  One core
will be left to the system.

When you choose to specify a number of cores, 1 of that number will be
reserved for putting data into the database, and the rest will be used
to process the data.

[Database]
==========

.. code-block:: python

    [Database]
    DATABASE            = demuxipy-results.sqlite

The database is where you specify the name of the database in which
your results will be stored. The database engine is SQLITE.

You may alter the database engine by writing your own `demuxi/db.conf`.

[Sequence]
==========

This section is where you tell demuxipy where your input sequence files
are located.  These files can be separate fasta and quality files:

.. code-block:: python

    [Sequence]
    fasta               = 'path/to/my/file.fasta'
    quality             = 'path/to/my/file.quality'

Or, you can pass the name of a single fastq file:

.. code-block:: python

    [Sequence]
    fastq               = 'path/to/my/file.fastq'

[Quality]
=========

If you would like to trim your reads **before** looking for sequence tags,
set:

.. code-block:: python

    [Quality]
    QualTrim                = True
    MinQualScore            = 10

Else, you need to set:

.. code-block:: python

    [Quality]
    QualTrim                = False
    MinQualScore            = 10

[MID (Outer) Tags]
==================

.. Note::

    Before you use these options, understand the way that demuxipy
    searches for Outer and Inner tags, given different options.  See
    :ref:`hierarchical-tagging`.

.. code-block:: python

    MidTrim                 = True
    MidGap                  = 5
    MidFuzzyMatching        = True
    MidAllowedErrors        = 1

If you have read the section on :ref:`hierarchical-tagging`, then the
concept of inner- and outer-tags should make sense to you.  Because we
started all this using 454 sequencing, we sometimes refer to the
"outer" tags as the MID (for multiplex identification) tags.

To specify that you wish to identify and trim the outer (MID) tag and to
provide the amount of sequence (in bp) prior to the MID/Outer tag that you are
willing to tolerate, set:

.. code-block:: python

    MidTrim                 = True
    MidGap                  = 5

To use fuzzy-matching+error correction to compensate for synthesis and
sequencing errors, set:

.. code-block:: python

    MidFuzzyMatching        = True
    MidAllowedErrors        = 1

This will turn on fuzzy matching and allow the program to correct
sequence tags having a single error.  The number of allowed errors
in Outer/MID tags can be independent of the number of allowed errors in
the Linker/Inner tags.

.. warning::

    DO NOT set these options unless you are using
    Levenshtein-distance-based (AKA edit-distance) sequence tags -
    otherwise the error correction may return erroneous results.
    
    - Sequence tags generated by edittag_ are Levenshtein-based
    - The MID and RL-MID tags from Roche, Inc. are Levenshtein-based
    
    - The multiplex indexes available from Illumina Inc. **ARE NOT**.

If you do not turn on fuzzy matching, then only tags matching the
expected sequence **perfectly** will be matched.


[Linker (Inner) Tags]
=====================

.. Note::

    Before you use these options, understand the way that demuxipy
    searches for Outer and Inner tags, given different options.  See
    :ref:`hierarchical-tagging`.

.. code-block:: python

    [Linker (Inner) Tags]
    LinkerTrim              = True
    LinkerGap               = 5
    LinkerFuzzyMatching     = True
    LinkerAllowedErrors     = 1

If you have read the section on :ref:`hierarchical-tagging`, then the
concept of inner- and outer-tags should make sense to you. We sometimes
refer to the "inner" tags as the Linker (for multiplex identification)
tags.

To specify that you wish to identify and trim the Inner tag and to
provide the amount of sequence (in bp) that you will tolerate between
the Outer tag and the Inner tag, set:

.. code-block:: python

    LinkerTrim                 = True
    LinkerGap                  = 5

To use fuzzy-matching+error correction to compensate for synthesis and
sequencing errors in the inner tags, set:

.. code-block:: python

    LinkerFuzzyMatching     = True
    LinkerAllowedErrors     = 1

This will turn on fuzzy matching and allow the program to correct
sequence tags having a single error.  The number of allowed errors
in Inner tags can be independent of the number of allowed errors in
the Outer tags.

.. warning::

    DO NOT set these options unless you are using
    Levenshtein-distance-based (AKA edit-distance) sequence tags -
    otherwise the error correction may return erroneous results.

    - Sequence tags generated by edittag_ are Levenshtein-based

    - The MID and RL-MID sequences from Roche, Inc. are Levenshtein-based

    - The multiplex indexes available from Illumina Inc. **ARE NOT**.

[Concatemers]
=============

If you are interested in searching your sequences for potential
concatemers, turn these options on.  Realize that when you search for
concatemers, demuxipy_ will only look within your sequence for all
possible Inner tags.  Also realize that searching for concatemers is
**very slow**.

.. code-block:: python

    [Concatemers]
    ConcatemerChecking      = True
    ConcatemerFuzzyMatching = False
    ConcatemerAllowedErrors = 1

Generally speaking, I would not use fuzzy matching to search for
concatemers.  However, if you like, you can turn on that option.  Be
aware that the higher the number of allowed errors, the more likely tou
are to match something that is not a true concatemer.  Using the fuzzy
matching for concatemers makes things **even slower**.

[Search]
========

Here, you give the section heading of the Inner/Outer tag grouping for
you have used to identify your reads.  If you only used Outer/MID tags
to identify reads, then you will set this parameter to:

.. code-block:: python

    [Search]
    SearchFor               = MidGroups

And define a :ref:`MidGroups` section, as below.  If you only used
Inner tags to identify reads, then you will set this parameter to:

.. code-block:: python

    [Search]
    SearchFor               = LinkerGroups

And define a :ref:`LinkerGroups` section, as below.  Finally, if you
used both Outer and Inner tags to identify your reads, then you will set
this parameter to:

.. code-block:: python

    [Search]
    SearchFor               = MidLinkerGroups

And define a :ref:`MidLinkerGroups` section, as below.

.. _MidGroups:

[MidGroups]
===========

If you have used only Outer/MID tags to identify your reads, define
which MID tag name identifies which group of reads.  The MID tag name
should correspond to the `Name:Sequence` defined in :ref:`Mids`.

This approach allows you to define the tags in the :ref:`Mids` section
once (many are predefined), and just pick which ones you want to use in
a given run in this section.i

.. code-block:: python

    [MidGroups]
    MID13 = crabs
    MID14 = lobsters
    MID15 = shrimp

This means that reads found starting with the sequence corresponding to
MID13 (as defined in :ref:`Mids`) will be denoted as belonging to the
group of reads derived from "crabs".  Similarly, reads starting with the
sequence corresponding to MID14 (as defined in :ref:`Mids`) will be
denoted as belonging to the group of reads derived from "lobsters", etc.

.. _LinkerGroups:

[LinkerGroups]
==============

If you have used only Inner/Linker tags to identify your reads, define
which Linker tag name identifies which group of reads.  The Linker tag name
should correspond to the `Name:Sequence` defined in :ref:`Linkers`.

This approach allows you to define the tags in the :ref:`Linkers` section
once (many are predefined), and just pick which ones you want to use in
a given run in this section.

.. code-block:: python

    [LinkerGroups]
    SimpleX1 = walrus
    SimpleX2 = whale
    SimpleX3 = porpoise

This means that reads found starting with the sequence corresponding to
SimpleX1 (as defined in :ref:`Linkers`) will be denoted as belonging to the
group of reads derived from "walrus".  Similarly, reads starting with the
sequence corresponding to SimpleX2 (as defined in :ref:`Linkers`) will be
denoted as belonging to the group of reads derived from "whale", etc.

.. _MidLinkerGroups:

[MidLinkerGroups]
=================

If you have used `Hierarchical tags`_ to identify your reads, define
which combination of Outer/Inner tag identifies which group of reads.  
The combination used should include an Outer tag as defined in
:ref:`Mids` and an inner tag as defined in the :ref:`Linkers`
section.

.. code-block:: python

    [MidLinkerGroups]
    MID15, SimpleX1 = cat
    MID15, SimpleX2 = unicorn
    MID16, SimpleX1 = pony

This means that reads that have an Outer tag sequence identical to MID15
and an inner tag sequence identical to SimpleX1 will be denoted as
belonging to the group of reads derived from "cat".  Similarly, reads
having an outer tag of MID15 and an Inner tag of SimpleX2 will be
denoted as belonging to the group named "unicorn".

.. _Linkers:

[Linkers]
=========

Here, define the name and sequence of those sequence tags used as Inner
tags.  This list may include tags you did not use during the run being
analyze, so you can keep altag of your Inner tag sequences within this
file, and select those you actually used as described above.

.. code-block:: python

    [Linkers]
    SimpleX1    = CGTCGTGCGGAATC
    SimpleX2    = GCTGCTGGCGAATC
    SimpleX3    = CGTGCTGCGGAACT
    IDX1        = AACC
    IDX2        = ACAG
    IDX3        = AGGA

.. _Mids:

[Mids]
======

Here, define the name and sequence of those sequence tags used as Outer
tags.  This list may include tags you did not use during the run being
analyze, so you can keep all of your Outer tag sequences within this
file, and select those you actually used as described above.


.. code-block:: python

    [Mids]
    RLMID1    = ACACGACGAC
    RLMID2    = ACACGTAGTA
    RLMID3    = ACACTACTCG
    MID1      = ACGAGTGCGT
    MID10     = TCTCTATGCG
    MID100    = TAGACTGCAC


.. _Python:  http://www.python.org
.. _`ConfigParser module`: http://docs.python.org/library/configparser.html
.. _edittag:  https://github.com/faircloth-lab/edittag/
