#!/bin/bash

if [ -f time.txt ]
then
    rm time.txt
fi

top=20
for i in {1..20}
do
    { time ./Analyzer OutTree.root test.root > /dev/null; } 2>> time.txt
done


sed -r -ne 's|real\s+0m([[:digit:]].[[:digit:]]{3})s| \1|p' <time.txt | awk '{for(i=0; i<NF; i++){ total += $i; sofar +=1;} if(NR == '$top') print "real: ", total/sofar;}'

sed -r -ne 's|user\s+0m([[:digit:]].[[:digit:]]{3})s| \1|p' <time.txt | awk '{for(i=0; i<NF; i++){ total += $i; sofar +=1;} if(NR == '$top') print "user: ", total/sofar;}'

sed -r -ne 's|sys\s+0m([[:digit:]].[[:digit:]]{3})s| \1|p' <time.txt | awk '{for(i=0; i<NF; i++){ total += $i; sofar +=1;} if(NR == '$top') print "sys:  ", total/sofar;}'