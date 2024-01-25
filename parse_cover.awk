#!/usr/bin/awk -f

# parse output from check_coverage.py looking for molecules with few matching
# molecules

$1 ~ /^t/ && $3 < 15 { print $1 }
