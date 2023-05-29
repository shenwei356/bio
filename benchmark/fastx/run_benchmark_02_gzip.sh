#!/bin/sh

echo Test: Gzip compressed
echo Output sequences of all apps are not wrapped to fixed length.


echo -en "\n============================================\n";

for f in dataset_*.f{a,q}.gz; do

    echo read file once with cat
    cat $f > /dev/null;
    
    
    echo -en "\n------------------------------------\n";    
    
    echo == seqkit
    echo data: $f;
    
    memusg -t -H seqkit seq $f -w 0 -o $f.seqkit.gz --compress-level 6;    
    
    pigz -cd $f.seqkit.gz | md5sum;
    /bin/rm $f.seqkit.gz;
    
    
    echo -en "\n------------------------------------\n";    
    
    echo == seqkit_t1
    echo data: $f;
    
    memusg -t -H seqkit seq $f -w 0 -j 1 -o $f.seqkit.gz --compress-level 6;    
    
    pigz -cd $f.seqkit.gz | md5sum;
    /bin/rm $f.seqkit.gz;
    
    
    echo -en "\n------------------------------------\n";  
    
    echo == seqtk+gzip
    echo data: $f;

    memusg -t -H seqtk seq $f | gzip -c > $f.seqtk.gz;
    
    pigz -cd $f.seqtk.gz | md5sum
    /bin/rm $f.seqtk.gz;  
    
        
    echo -en "\n------------------------------------\n";  
    
    echo == seqtk+pigz
    echo data: $f;

    memusg -t -H seqtk seq $f | pigz -p 4 -c > $f.seqtk.gz;
    
    pigz -cd $f.seqtk.gz | md5sum
    /bin/rm $f.seqtk.gz;  
done

