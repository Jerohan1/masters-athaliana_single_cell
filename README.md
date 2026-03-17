# masters-athaliana_single_cell

## Data

*PRJEB77115 scRNA col0hybrids*
```
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR140/ERR14083819/col0xmany.sorted.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR140/ERR14083819/col0xmany.sorted.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR140/ERR14083820/col0xdb1.sorted.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR140/ERR14083820/col0xdb1.sorted.bam
```

*GCA_028009825.2_Col-CC reference genome*
```
wget -r -np -nH --cut-dirs=6 \
  -R "index.html*" \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/009/825/GCA_028009825.2_Col-CC/
```

*Organellar sequence and gene features*
```
https://www.ncbi.nlm.nih.gov/nuccore/NC_037304.1?report=genbank #NC_037304.1, mitochondrion, send to > complete record > fasta and gff3
https://www.ncbi.nlm.nih.gov/nuccore/7525012 #NC_000932.1, chloroplast, send to > complete record > fasta and gff3
```

*doi-10.17617-3.aeojbl assemblies for parent2*
```
https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.AEOJBL
```

*Schneeberger paper files*
```
git clone https://github.com/schneebergerlab/snrna_eqtl_mapping.git # pipeline files and plotting notebooks
wget https://zenodo.org/api/records/14864054/files-archive # supplementary data
```

*Cellranger barcode whitelist for STAR*
```
wget https://github.com/f0t1h/3M-february-2018/blob/master/3M-february-2018.txt.gz
```