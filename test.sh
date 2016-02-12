#!/usr/bin/env bash
PYTHONPATH=. ptipython
#sort -m <(grep "^chr1\b.*+" examples/data/control.bed |
#    cut -f 1-3 |
#    sort -u -k2,3n |
#    perl -a -ne '$F[1]+=75; # shift fragments
#    $F[1]=$F[1]-($F[1] % 200); # turn exact start into bin
#    print "@F[0,1]\n"')
#    <(grep "^chr1\b.*-" examples/data/control.bed |
#    cut -f 1-3 |
#    sort -u -k2,3n |
#    perl -a -ne '$F[2]-=75; # shift fragments
#    $F[2]=$F[2]-($F[2] % 200); # turn exact start into bin
#    print "@F[0,2]\n"') >deleteme.txt

# sort -m <(grep "^${2}\b.*+" $1  | cut -f 2,3 | sort -u -k2,3n | awk '{t = $2+75; t2 = t - (t % 200); print t2}') <(grep "^${2}\b.*-" $1  | cut -f 2,3 | sort -u -k2,3n | awk '{t = $2-75; t2 = t - (t % 200); print t2}') | uniq -c | sed -e 's/^[ ]*//' > ${3}
