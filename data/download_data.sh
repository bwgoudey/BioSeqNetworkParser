#!/bin/bash

#wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
#zcat < ./uniprot_sprot.dat.gz | bgzip -c > ./uniprot_sprot.dat.bgz

for i in `seq 0 5 300`; do 
    wget https://ftp.ncbi.nlm.nih.gov/genbank/gbbct${i}.seq.gz 
done
