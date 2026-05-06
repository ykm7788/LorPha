# LorPha

LorPha is a pixi based one-stop phageome analysis pipeline designed for biologists. LorPha can simultaneously extract both viral and prokaryotic sequences at contig level. LorPha can also parallelly analysize multiple samples and automatically generate vOTU-abundance table for downstream analysis.

## Usage

### Step 1 Installation

#### 1.1 Install pixi

As LorPhais developed based on pixi, please refer to pixi to first make sure pixi is installed.&#x20;

[https://pixi.prefix.dev/](https://pixi.prefix.dev/ "https://pixi.prefix.dev/")

#### 1.2 Install LorPha

```bash 
git clone https://github.com/ykm7788/LorPha.git
```


#### 1.3 Install databases

Overall, LorPha requires 9 databases. Users can easily download and uncompress (7zip format) our pre-configured database file from netdisk.&#x20;

1. Quark:
2. Aliyun:
3. Baidu:

Alternatively, Users can manually download these databases from official sources. Please create a folder to contain all databases. Each subfolder within it corresponds to one database, and the names (case sensitive) of these subfolders are listed below:&#x20;

1. `checkv`
   [https://bitbucket.org/berkeleylab/CheckV/src/master/#markdown-header-checkv-database](https://bitbucket.org/berkeleylab/CheckV/src/master/#markdown-header-checkv-database "https://bitbucket.org/berkeleylab/CheckV/src/master/#markdown-header-checkv-database")
2. `dbcan`
   [https://pro.unl.edu/dbCAN2/browse\_download.php](https://pro.unl.edu/dbCAN2/browse_download.php "https://pro.unl.edu/dbCAN2/browse_download.php")
   run\_dbCAN\_database and V14 are required
3. `defensefinder`
   [https://github.com/mdmparis/defense-finder](https://github.com/mdmparis/defense-finder "https://github.com/mdmparis/defense-finder")
   Models based on MacSyFinder and CasFinder should be installed.
4. `genomad`
   [https://github.com/apcamargo/genomad/](https://github.com/apcamargo/genomad/ "https://github.com/apcamargo/genomad/")
5. `gtdb`
   [https://github.com/soedinglab/mmseqs2/wiki#downloading-databases](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases "https://github.com/soedinglab/mmseqs2/wiki#downloading-databases")
   Please use mmseq2 to install gtdb taxonomy database and make sure the database is named as 'GTDB' (case sensitive). Current version: r232.
6. `phabox`
   [https://github.com/KennthShang/PhaBOX](https://github.com/KennthShang/PhaBOX "https://github.com/KennthShang/PhaBOX")
   Follow phabox team to download the latest database file. Current version: V2\_2.
7. `PhaStyle`
   [https://huggingface.co/neuralbioinfo/prokbert-mini-long-phage](https://huggingface.co/neuralbioinfo/prokbert-mini-long-phage "https://huggingface.co/neuralbioinfo/prokbert-mini-long-phage")
   By default, model prokbert-mini-long-phage was used.
8. `vibrant`
   [https://github.com/AnantharamanLab/VIBRANT](https://github.com/AnantharamanLab/VIBRANT "https://github.com/AnantharamanLab/VIBRANT")
9. `vhost`&#x20;
   [https://www.genome.jp/ftp/db/virushostdb/](https://www.genome.jp/ftp/db/virushostdb/ "https://www.genome.jp/ftp/db/virushostdb/")
   We used the data from Virus-host DB to generate validated virus-host links. To accelerate the matching speed, LexicMap was applied to replace traditional blast. Therefore, users need to create a LexicMap-style database using `virushostdb.formatted.genomic.fna` and attach a file connecting geneID and virus-host relationship (find `virus.seqID.txt` in our netdisk). Users can also create their customized databases by generating new LexicMap-style database and modifying relationship file.&#x20;

### Step 2 Configure input sample file and parameter file

#### 2.1 sample file

Before runing LorPha, users need to manually edit a sample file containing the path to sample data. The sample file is a tab-separated file with four columns. Each row represents one sample that you want to analysize. All samples in this file will be parallelly processed and their results will be merged together. Please do not modify the column name of sample file. The `sample.txt` file should like this:

| SampleID | Contig                  | Cleandata\\\_R1           | Cleandata\\\_R2           |
| -------- | ----------------------- | ------------------------- | ------------------------- |
| sample1  | /path/sample1.contig.fa | /path/sample1\\\_R1.fq.gz | /path/sample1\\\_R2.fq.gz |
| sample2  | /path/sample2.contig.fa | /path/sample2\\\_R1.fq.gz | /path/sample2\\\_R2.fq.gz |
| sample3  | /path/sample3.contig.fa | /path/sample3\\\_R1.fq.gz | /path/sample3\\\_R2.fq.gz |

SampleID: Name of sample.

Contig: The absolute path to contig sequence (Should be an uncompressed fasta file).

Cleandata\_R1: The absolute path to paired cleandata R1 (Can be compressed).

Cleandata\_R2: The absolute path to paired cleandata R2 (Can be compressed).

#### 2.2 parameter file

In addition, users need to provide a parameter file containing the parameter of each software. The parameter file is a tab-separated file with two columns. First column indicates the software and second column shows the exact command-line parameter. Our `parameter.txt` gives a default setting while uses can also adjust parameters following the protocal of each software (phabox, genomad, lexicmap, vclust, dbcan, mmseqs, coverm). Please do not modify the column name and first column of this file.

| Software            | Parameter                                                                                                                                      |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| phabox              | \\--len 1500                                                                                                                                   |
| genomad             | \\--lenient-taxonomy --full-ictv-lineage                                                                                                       |
| lexicmap            | \\--align-min-match-pident 80 --min-qcov-per-hsp 0 --min-qcov-per-genome 50 --align-min-match-len 1500 --top-n-genomes 1                       |
| vclust prefilter    | \\--min-ident 0.95                                                                                                                             |
| vclust cluster      | \\--algorithm cd-hit --metric ani --ani 0.95 --qcov 0.85 --out-repr                                                                            |
| dbcan               | --e\_value\_ threshold 0.00001 --evalue\_cutoff 0.00001                                                                                        |
| mmseqs easy-cluster | \\--min-seq-id 0.8 -c 0.8 --cov-mode 5 --remove-tmp-files 1 --write-lookup 1                                                                   |
| mmseqs taxonomy     | \\--tax-lineage 1 -e 0.00001                                                                                                                   |
| coverm              | --min-read-aligned-percent 50 --min-read-percent-identity 90 --min-covered-fraction 0 --contig-end-exclusion 0 -m rpkm count reads\_per\_ base |

### Step 3 Run LorPha

#### 3.1 (Optional) Change to the LorPha fold

```bash 
cd /path/LorPha
```


#### 3.2 Setup LorPha jobs

```bash 
 pixi run lorpha_setup -r /AbsolutePath/rawdata.txt -p /AbsolutePath/parameter.txt -t core_num -d /AbsolutePath/database -o /path/result
```


With this command, a `job.sh` file should be created that contains all jobs LorPha will run.&#x20;

#### 3.3 Run

```bash 
sh job.sh
```


We recommend users to use `nohup `or `tmux `to run this script as normally it takes a long time to fully process all jobs. Moreover, users can also manually edit `job.sh` to selectively run certain jobs.
