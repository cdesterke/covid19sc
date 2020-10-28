#!/bin/bash

mkdir GSE145926

cd GSE145926

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145926/suppl/GSE145926_RAW.tar
tar -xvf  GSE145926_RAW.tar
rm GSE145926_RAW.tar
rm *.gz


