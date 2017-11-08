# kallisto Tutorial

# Installation
The latest version of kallisto can be installed with bioconda:

```shell
conda install kallisto
```

Executables and source codes are available at [official kallisto project page](https://pachterlab.github.io/kallisto/download).

At 2017.11.08, the latest version of kallisto is v0.43.1, which can be downloaded via command below (in Linux OS):

```shell
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz -O kallisto_linux-v0.43.1.tar.gz
```
Extract executables and change into the directory.

```shell
tar xvf kallisto_linux-v0.43.1.tar.gz; cd kallisto_linux-v0.43.1
```

Official installation of kallisto is over. It will be convenient to add environment variable `$KALLISTO` which stores the path to the directory containing kallisto binary file. Then,

```shell
echo $KALLISTO
```
should give path to kallisto directory.

# Running kallisto

## Preparing reference transcriptome
To run kallisto, you don't need reference whole genome sequence, but *reference cDNA/transcriptome*. For more information on preparing reference transcriptome, refer to [reference-data](../reference-data) part.

## kallisto index
We need to index reference transcriptome for efficient quantification of transcripts. This can be done by running `kallisto index`

**Usage**

```shell
kallisto index [arguments] FASTA-files
```

- Required argument:

`-i, --index=STRING` : Filename for the kallisto index to be constructed

- Optional argument:

`-k, --kmer-size=INT` : k-mer (odd) length (default: 31, max value: 31)

`--make-unique` : Replace repeated target names with unique names

For example, run kallisto index with command below:

```shell
mkdir index

$KALLISTO/kallisto index -i index/Homo_sapiens.GRCh38.cdna.all.kallistoindex \
../reference-data/GRCh38_rel90/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

This will create `Homo_sapiens.GRCh38.cdna.all.kallistoindex` file under the directory `index/`.
## kallisto quant
Now it is time to quantify your RNA-seq data. Of course you should have your RNA-seq result data, if not, please refer to [test-data](../test-data) part.