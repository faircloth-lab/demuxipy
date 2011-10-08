.. _running:

*******************
Running `demuxipy`
*******************

demuxipy_ is a program that is run from the command line (i,e.
"terminal").  Because there are many options, demuxipy_ requires that you setup a
configuration file containing the options that will drive the program
and then pass that file to the program on the command line.  To setup a
configuration file, see :ref:`setting-up-a-configuration-file`.

After installing demuxipy_ and generating a configuration file, 
running demuxipy_ is simply a matter of:

.. code-block:: bash

    python demuxi.py path/to/my_configuration_file.conf

After you start a run, you will see something like the following:

.. code-block:: none

    % demuxi.py demuxipy/tests/test-data/demuxi-test.conf

    ###############################################################
    #                       demuxi.py                             #
    #                                                             #
    # Demultiplexing of hierarchically tagged, massively parallel #
    # DNA sequence reads                                          #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2011 Brant C. Faircloth                  #
    #                                                             #
    #                                                             #
    # Ecology and Evolutionary Biology                            #
    # 621 Charles E. Young Drive                                  #
    # University of California, Los Angeles, 90095, USA           #
    ###############################################################


    Started:  Fri Oct 07, 2011  14:28:17
    Parsing reads into groups of 9 reads
    Starting 2 workers

    .........%.........%.........%...

    Ended: Fri Oct 07, 2011  14:28:18 (run time 0.011 minutes)


Which indicates progress of the program.  Generally speaking, one dot in
the above progress indicator is output for every 1000 sequences
processed, and one percent symbol is output every 10,000 sequences.

Once the program runs, you should proceed to :ref:`getting_data`.

.. _demuxipy: https://github.com/faircloth-lab/demuxipy/
