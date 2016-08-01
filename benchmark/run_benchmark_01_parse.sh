#!/bin/sh

echo Test: FASTA/Q Parsing
echo Output sequences of all apps are not wrapped to fixed length.


echo -en "\n============================================\n";

echo == seqtk
for f in dataset_*.f{a,q}; do
    echo -en "\n------------------------------------\n";
    echo data: $f;
    
    echo read file once by cat
    cat $f > t; /bin/rm t; # warm up

    memusg -t -H seqtk seq $f > $f.seqtk.fa;
    
    md5sum $f.seqtk.fa;
    seqkit stat $f.seqtk.fa;
    
    /bin/rm $f.seqtk.fa;
done


echo -en "\n============================================\n";

echo == seqkit
for f in dataset_*.f{a,q}; do
    echo -en "\n------------------------------------\n";
    echo data: $f;
    
    echo read file once by cat
    cat $f > t; /bin/rm t; # warm up

    memusg -t -H seqkit seq $f -w 0 > $f.seqkit.fa;    
    
    md5sum $f.seqkit.fa;
    seqkit stat $f.seqkit.fa;
    
    /bin/rm $f.seqkit.fa;
done
