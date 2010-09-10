# Running linker-py #

There are several things you will likely want to do before you run the program:

1. you need to install the dependencies of the code and get MySQL setup with both a database to hold your results and a user that will connect to the database to write the results of this program (the user needs to be able to WRITE, CREATE, and DROP tables).  You also need to install MySQLdb.

1. you will need to fill out the linker-py.conf file appropriately.  this means that if you are using hierarchical tags, you need to list which MID+SimpleX go together to give what combination (e.g. MID16, SimpleX2 = unicorns).  If you used RLMIDs, then use those in the combinations instead (e.g. RLMID1, SimpleX2 = unicorns) Only list those combinations you used in your run, otherwise, things will go slower (the program will look for all additional combinations that are not there).  You SHOULD NOT have to input the actual MID, RLMID, or SimpleX sequences, provided they are the same as those provided at the bottom of the file ("Sequence Tag Section").

2. you may want to change MinQualScore (in linker-py.conf) from 10 to something higher - this is the min score to trim off the ends, which happens before we search for tags2. you may want to change MinQualScore (in linker-py.conf) from 10 to something higher - this is the min score to trim off the ends, which happens before we search for tags.

3. you will want to add your database parameters to linker-py.conf and make sure that it's setup for multiprocessing.you will want to add your database parameters to the file and make sure that it's setup for multiprocessing.

4. once you do all of that, you probably want to install pylevenshtein (link on the project page above), and then run helpers/levenshtein.py on your configuration file, like so:

        python helpers/levenshtein.py --configuration=linker-py.conf \
        --section=MidLinkerGroups --verbose

    this will let you know the minimum edit distance in your linker combinations, which determines how many errors you can correct.  you need a minimum edit distance of 3 to correct 1 error; 5 to correct 2 errors; 7 to correct 3 errors, etc.  the general formula is:

        (edit_distance - 1)/2 = number_of_correctable_errors

    where n is the number of errors you want to correct.  You *cannot* round up a fractional result.  Once you know that, you can edit the AllowedErrors section of the `linker-py.conf` file to reflect how many error you can correct, so it should read like:

        AllowedErrors       = 1

1. once all thats done, you can run the program:

        python linker.py --configuration=linker-py.conf

    this should begin running without errors, placing your reads into a table named `sequence` in whatever database you told it to use.     

1. you can look at the data in this database with the following query (which leaves off the BLOB column):

    select id, name, mid, mid_seq, mid_match, mid_method, linker, linker_seq,
    linker_match, linker_method, cluster, concat_seq, concat_match,
    concat_method, n_count, untrimmed_len, trimmed_len, seq_trimmed from
    sequence;

1. You can get the data out and into some fasta files with:

        python helpers/get_parsed_reads.py --configuration=linker-py.conf \
        --all --trimmed

    which will create fasta and qual files by whatever groups your hierarchical tags go into, but i haven't extensively tested this on other machines, so it may not work so smoothly.