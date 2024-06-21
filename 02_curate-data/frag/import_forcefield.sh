#!/bin/bash

set -xe

# drop the existing forcefields table and import a new one from hpc3 in a file
# named `ff.sql`

sqlite3 store.sqlite 'DROP TABLE forcefields'
sqlite3 store.sqlite < ff.sql
