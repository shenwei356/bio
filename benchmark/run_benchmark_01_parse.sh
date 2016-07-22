#!/bin/sh

echo Test: FASTA/Q Parsing
echo Output sequences of all apps are not wrapped to fixed length.


echo == seqtk
for f in dataset_*.f{a,q}; do    
    echo data: $f;
    
    echo read file once by cat
    cat $f > t; /bin/rm t; # warm up

    memusg -t -H seqtk seq $f > t;
    
    /bin/rm t;
done



echo == seqkit
for f in dataset_*.f{a,q}; do    
    echo data: $f;
    
    echo read file once by cat
    cat $f > t; /bin/rm t; # warm up

    memusg -t -H seqkit seq $f -w 0 > t;
    
    /bin/rm t;
done
