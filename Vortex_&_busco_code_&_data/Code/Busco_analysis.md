```bash
wd: /hpc-home/ryanc/scratch/BUSCO/
source lmod-6.1
ml git python
git clone https://gitlab.com/ezlab/busco

python setup.py install

wget https://busco.ezlab.org/datasets/aves_odb9.tar.gz


```



**<u>Config.ini file</u>**

```
## BUSCO specific configuration
## It overrides default values in code and dataset cfg, and is overridden by arguments in command line
## Uncomment lines with single # when appropriate
[busco]
## Input file
# in = ./sample_data/target.fa
## Run name, used in output files and folder
# out = SAMPLE
## Where to store the output directory
# out_path = ./sample_data
## Path to the BUSCO dataset
# lineage_path = ./sample_data/example
## Which mode to run (genome / protein / transcriptome)
# mode = genome
## How many threads to use for multithreaded steps
# cpu = 1
## Domain for augustus retraining, eukaryota or prokaryota
# domain = eukaryota # do not change this unless you know exaclty why !!!
## Force rewrite if files already exist (True/False)
# force = False
## Restart mode (True/False)
# restart = False
## Blast e-value
# evalue = 1e-3
## Species to use with augustus, for old datasets only
# species = fly
## Augustus extra parameters
# augustus_parameters = '' # nothing here, use single quotes, like this: '--param1=1 --param2=2'
## Tmp folder
# tmp_path = ./tmp/
## How many candidate regions (contigs, scaffolds) to consider for each BUSCO
# limit = 3
## Augustus long mode for retraining (True/False)
# long = False
## Quiet mode (True/False)
# quiet = False
## Debug logs (True/False), it needs Quiet to be False
# debug = True
## tar gzip output files (True/False)
# gzip = False

[tblastn]
## path to tblastn
path = /tgac/software/testing/blast/2.2.30/x86_64/bin/

[makeblastdb]
## path to makeblastdb
path = /tgac/software/testing/blast/2.2.30/x86_64/bin/

[augustus]
## path to augustus
path = /tgac/software/testing/augustus/3.2.1/x86_64/bin/

[etraining]
## path to augustus etraining
path = /tgac/software/testing/augustus/3.2.1/x86_64/bin/

## path to augustus perl scripts, redeclare it for each new script
[gff2gbSmallDNA.pl]
path = /tgac/software/testing/augustus/3.2.1/x86_64/scripts/
[new_species.pl]
path = /tgac/software/testing/augustus/3.2.1/x86_64/scripts/
[optimize_augustus.pl]
path = /tgac/software/testing/augustus/3.2.1/x86_64/scripts/

[hmmsearch]
## path to HMMsearch executable
path = /tgac/software/testing/hmmer/3.1b2/x86_64/bin/

[Rscript]
## path to Rscript, if you wish to use the plot tool
path = /tgac/software/testing/R/3.3.1/x86_64/bin/
```



**<u>Script to run</u>**

```bash
source lmod-6.1
ml python
run_BUSCO.py -o pp -i pink_pigeon.prep.scafSeq.fasta -l aves_odb9 -c 4 -m genome
```



**<u>Method</u>**

An estimate of the completeness of the assembled genome was calculated using BUSCO v. 3.1.0 (Robert et al. 2017) and the aves databes (aves_odb9). A total of 4915 BUSCO groups were searched and of those 4639 (94.38%) were identified in the genome assembly. There were 4586  (93%)  BUSCO groups identified as single copy and complete, 168 (3.42%) as fragmented and 108 (2.20%) were missing. 



**<u>Results</u>**

| Type of  BUSCO                      | Percentage of total BUSCOs | Number of BUSCOs |
| ----------------------------------- | -------------------------- | ---------------- |
| Complete BUSCOs (C)                 | 94.38453713                | 4639             |
| Complete and single-copy BUSCOs (S) | 93.30620549                | 4586             |
| Complete and duplicated BUSCOs (D)  | 1.078331638                | 53               |
| Fragmented BUSCOs (F)               | 3.418107833                | 168              |
| Missing BUSCOs (M)                  | 2.197355036                | 108              |

A total of 4915 BUSCO groups were searched and of those 93% were single copy and complete.
