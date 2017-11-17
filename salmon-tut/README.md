# Salmon Tutorial

# Installation
You can find the latest binary releases of Salmon [here](https://github.com/COMBINE-lab/salmon/releases). Current (2017.11.17) latest version of Salmon is `Salmon-0.8.2`, which can be downloaded by the command below.

```shell
wget \
https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz \
-O Salmon-0.8.2_linux_x86_64.tar.gz
```

Decompress it and add the directory to `$SALMON` environment variable for further convenience by modifying your `.bashrc` file.

```shell
tar xvf Salmon-0.8.2_linux_x86_64.tar.gz
```

Echoing `$SALMON` should give the directory containing salmon binary file.

```shell
echo $SALMON
```

and test your installation. Running

```shell
$SALMON/salmon -h
```

should give

```shell
Salmon v0.8.2

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     cite  Show salmon citation information
     index Create a salmon index
     quant Quantify a sample
     swim  Perform super-secret operation
```

if it is successfully installed.

# Running Salmon

## Preparing reference transcriptome
To run Salmon, you don't need reference whole genome sequence, but *reference cDNA/transcriptome*. For more information on preparing reference transcriptome, refer to [reference-data](../reference-data) part.

## salmon index
Creating Salmon index which facilitates the quasi-mapping process of Salmon only need to be done once before quantifying your RNA-seq data. 

**Usage**

```shell
salmon index -t <TRANSCRIPT_FASTA> -i <INDEX>
```

- Required argument

  - `-t, --transcripts`: Transcript fasta file.

  - `-i, --index`: Salmon index file name.

- Optional argument

  - `-k, --kmerlen`: The size of k-mers that should be used for the quasi-index

  - `-p, --threads`: Number of threads to use.

For example, `salmon index` can be executed as follows, and will be done in few minutes:

```shell
mkdir index 

$SALMON/salmon index \
-t ../reference-data/Homo_sapiens.GRCh38.cdna.all.fa.gz  \
-i index/Homo_sapiens.GRCh38.cdna.all.fa.salmonindex
```

## salmon quant
Quantification is simple.

```shell
salmon quant -i <INDEX> -l <LIBRARY_TYPE> \
-1 <READ_1> -2 <READ_2> -o <OUTPUT>
```

- Required argument

  - `-i, --index`: Salmon index file name.

  - `-1, --mates1`: File containing the #1 mates

  - `-2, --mates2`: File containing the #2 mates

  - `-o, --output`: Output quantification file.

- Optional argument

  - `-p, --threads`: Number of threads to use.

For example, `salmon index` can be executed as follows, and will be done in few minutes:

```shell
mkdir -p results/SRR1313090

$SALMON/salmon quant \
-i index/Homo_sapiens.GRCh38.cdna.all.fa.salmonindex \
-l A \
-1 ../test-data/data/SRR1313090/large_1.fastq \
-2 ../test-data/data/SRR1313090/large_2.fastq \
-o results/SRR1313090/large_quant \
-p 8
```