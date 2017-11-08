# Preparing test data
Gene expression profiles of human primary monocytes from six individuals are available at Gene Expression Omnibus ([GSE80095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80095)). One of these data ([SRR3352080]) is selected for use.

There are a few ways to obtain this data. In this tutorial, we are going to use Europian Nucleotide Archive (ENA). Fortunately, ENA FTP is well-managed so that we can programmatically download sequence data. `scripts/fastq-dump-ena.py` allows us to download short read sequence data conveniently. 

**Usage**

```
./scripts/fastq-dump-ena.py ACCESSION [-p] [-o output]
```
**Options**
`-p` : Add this option if the data to be downloaded is paired-end sequences.
`-o output`: Specify output directory