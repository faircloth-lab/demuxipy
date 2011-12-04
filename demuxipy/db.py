"""
File: db.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 October 2011 15:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Database creation and entry functions for demuxi.py 

"""

import os
import sys
import sqlite3
try:
    import cPickle as pickle
except:
    print "Using pickle instead of cPickle"
    import pickle

def create_db_and_new_tables(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_keys = ON")
    try:
        cur.execute('''CREATE TABLE tags (
                id integer PRIMARY KEY AUTOINCREMENT,
                name text,
                mid text,
                mid_seq text,
                mid_match text,
                mid_method text,
                linker text,
                linker_seq text,
                linker_match text,
                linker_method text,
                cluster text,
                concat_seq text,
                concat_match text,
                concat_method text
            )'''
        )
        cur.execute('''CREATE TABLE sequence (
                id INTEGER,
                untrimmed_len integer,
                trimmed_len integer,
                n_count integer,
                seq_trimmed text,
                record blob,
                FOREIGN KEY(id) REFERENCES tags(id) DEFERRABLE INITIALLY
                DEFERRED
            )'''
        )
        cur.execute("CREATE INDEX idx_sequence_cluster on tags(cluster)")
    except sqlite3.OperationalError, e:
        #pdb.set_trace()
        if "already exists" in e[0]:
            answer = raw_input("\nDatabase already exists.  Overwrite [Y/n]? ")
            #pdb.set_trace()
            if answer == "Y" or answer == "YES":
                os.remove(db_name)
                conn, cur = create_db_and_new_tables(db_name)
            else:
                sys.exit()
        else:
            raise sqlite3.OperationalError, e
    return conn, cur

def insert_record_to_db(cur, tagged):
    name = tagged.read.identifier.split(' ')[0].lstrip('>')
    cur.execute('''INSERT INTO tags (
            name,
            mid,
            mid_seq,
            mid_match,
            mid_method,
            linker,
            linker_seq,
            linker_match,
            linker_method,
            cluster,
            concat_seq,
            concat_match,
            concat_method
        ) 
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', 
        (
            name,
            tagged.mid_name,
            tagged.mid,
            tagged.seq_match,
            tagged.m_type,
            tagged.l_name,
            tagged.l_tag,
            tagged.l_seq_match,
            tagged.l_m_type,
            tagged.l_critter,
            tagged.concat_tag,
            tagged.concat_seq_match,
            tagged.concat_m_type
        )
    )
    key = cur.lastrowid
    # pick the actual sequence
    sequence_pickle = pickle.dumps(tagged.read,pickle.HIGHEST_PROTOCOL)
    cur.execute('''INSERT INTO sequence (
            id,
            untrimmed_len,
            trimmed_len,
            n_count,
            seq_trimmed,
            record
        )
        VALUES (?,?,?,?,?,?)''',
        (
            key,
            tagged.original_len,
            len(tagged.read.sequence),
            tagged.read.sequence.lower().count('n'),
            tagged.read.sequence,
            sqlite3.Binary(sequence_pickle)
        )
    )
