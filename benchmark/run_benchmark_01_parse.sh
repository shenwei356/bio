#!/bin/sh

echo Test: FASTA/Q Parsing
echo Output sequences of all apps are not wrapped to fixed length.


echo -en "\n============================================\n";

for f in dataset_*.f{a,q}; do

    echo read file once with cat
    cat $f > /dev/null;
    
    
    echo -en "\n------------------------------------\n";    
    
    echo == seqkit
    echo data: $f;
    
    memusg -t -H seqkit seq $f -w 0 > $f.seqkit.fa;    
    
    md5sum $f.seqkit.fa;    
    /bin/rm $f.seqkit.fa;
    
    
    echo -en "\n------------------------------------\n";  
    
    echo == seqtk
    echo data: $f;

    memusg -t -H seqtk seq $f > $f.seqtk.fa;
    
    md5sum $f.seqtk.fa;    
    /bin/rm $f.seqtk.fa;        
done

