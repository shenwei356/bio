#!/bin/sh

echo Test: Plain text
echo Output sequences of all apps are not wrapped to fixed length.


echo -en "\n============================================\n";

for f in dataset_*.f{a,q}; do

    echo read file once with cat
    cat $f > /dev/null;
    
    
    echo -en "\n------------------------------------\n";    
    
    echo == seqkit
    echo data: $f;
    
    memusg -t -H seqkit seq $f -w 0 > $f.seqkit.fx;    
    
    md5sum $f.seqkit.fx;    
    /bin/rm $f.seqkit.fx;
    
    
    
    echo -en "\n------------------------------------\n";    
    
    echo == seqkit_t1
    echo data: $f;
    
    memusg -t -H seqkit seq $f -w 0 -j 1 > $f.seqkit.fx;    
    
    md5sum $f.seqkit.fx;    
    /bin/rm $f.seqkit.fx;
    
    
    
    echo -en "\n------------------------------------\n";  
    
    echo == seqtk
    echo data: $f;

    memusg -t -H seqtk seq $f > $f.seqtk.fx;
    
    md5sum $f.seqtk.fx;    
    /bin/rm $f.seqtk.fx;        
done

