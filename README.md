###This is the repository of code used for alignment and variant detection for APL sequencing project.

#### Parameters used for alignment and variant detection are hardcoded into the script. May have to edit accordingly.

#### Alignment:
Assumes bwa, samblaster, sambaba and samtools are under path.
Edit paths for faidx indexed reference genome, gatk jar file, known indels (from gatk [bundle](ftp://ftp.broadinstitute.org/bundle/2.8/hg19/)) and exome baits manually. (From line 15 to 26)

```bash
$ fqtobam.py sample_id sample_name sample_library fq1.gz fq2.gz
```
#### Variant calling using VarScan2:
Assuems samtools and [bam-readcount](https://github.com/genome/bam-readcount) are installed under path.
Edit paths for faidx indexed reference genome, varscan2 jar file, [fpFilter](https://github.com/ckandoth/variant-filter) perl script maunally. (From line 15 to 18) 

```bash
$ bash run_varScan.sh normal.bam tumor.bam basename
```