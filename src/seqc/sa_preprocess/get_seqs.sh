#!/bin/bash

wget ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
wget ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

gunzip *.fa.gz
cat *.fa > Mus_musculus.GRCm38.cdna_all_plus_ncrna.fa

rm Mus_musculus.GRCm38.ncrna.fa
rm Mus_musculus.GRCm38.cdna.all.fa

