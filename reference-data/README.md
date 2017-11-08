# Preparing reference data

## Preparing reference genome sequence data

We need reference genome of Homo sapiens. We can easily download it from [ensemble FTP](http://www.ensembl.org/info/data/ftp/index.html). If you are confused about selecting appropriate reference genomes (since there are some variants, such as `dna_sm`, `dna_rm` or `primary_assembly`, `toplevel`), there is an execellent [blog post](http://genomespot.blogspot.kr/2015/06/mapping-ngs-data-which-genome-version.html) which might help you. We are going to select *soft masked and primary assembly* reference genome here.

```shell
mkdir GRCh38_rel90

wget \
ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
-P GRCh38_rel90
```

## Preparing reference transcriptome sequence data

For some alignment-free transcript quantification tools, we need to prepare reference cDNA/transcriptome sequence data. Reference cDNA sequence data can also be downloaded easily from [ensembl FTP](http://www.ensembl.org/info/data/ftp/index.html). For our tutorial, *cdna.all* sequence file will do.

```shell
mkdir GRCh38_rel90  # only if it doesn't exist
wget \
ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-P GRCh38_rel90
```