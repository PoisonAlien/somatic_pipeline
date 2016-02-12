#!/bin/bash
#
# Author: Anand Mayakonda 
# Runs varScan2 somatic on given tumor and normal pair and outpts snp as well as indel files followed by filtering for false positives.
#
# Usage: bash run_varScan.sh normal.bam tumor.bam basename

NORMALBAM=$1
TUMORBAM=$2
bn=$3

echo -e "$NORMALBAM\t$TUMORBAM\t$bn"

#Set path to faidx indexed reference genome, varscan jar file, fpFilter perl script and annovar.
REFFILE="/home/csipk/NGS/ref_genomes/hg19.fa"
VARSCAN="/home/csipk/NGS/varscan/VarScan.v2.4.0.jar"
fpFilter="/home/csipk/NGS/varscan/fpfilter.pl"
annovar="/home/csipk/NGS/annovar/table_annovar.pl"

NOW=$(date)
echo -e "Started run\t$NOW" 

samtools mpileup -B -f $REFFILE -q 15 -L 10000 -d 10000 $NORMALBAM $TUMORBAM |
java -Xmx16g -d64 -jar $VARSCAN somatic -mpileup $bn --min-coverage-normal 10 --min-coverage-tumor 14 --min-var-freq 0.02 --strand-filter 1

echo -e "processing raw output \n"

java -jar $VARSCAN processSomatic $bn.snp --min-tumor-freq 0.07 --max-normal-freq 0.02 --p-value 0.05

NOW=$(date)
echo -e "DONE!\n${NOW}"


#### Filter processSomatic output (high confidence snps) for potential false positives.

#Make an output directory for output.
mkdir -p $bn"_fpFilter"

#prepare file for bam-readcount
sed 1d $bn.snp.Somatic.hc | awk '{print $1"\t"$2"\t"$2}' > ./$bn"_fpFilter"/$bn.var 

#Run bam-readcount to get variant statistics.
bam-readcount -q15 -w1 -b15 -l ./$bn"_fpFilter"/$bn.var -f $ref $bam > ./$bn"_fpFilter"/$bn.readCounts

#prepare --var-file file for fpfilter
cat $var | sed '/^#/d' | sed 1d | cut -f 1-4 > ./$bn"_fpFilter"/$bn.var 

#Run fpfilter script.
perl $fpFilter --var-file ./$bn"_fpFilter"/$bn.var --readcount-file ./$bn"_fpFilter"/$bn.readCounts --min-depth 14 --min-read-pos 0.1 --min-strandedness 0.01 --max-mmqs-diff 90 --min-var-count 4 --min-var-frac 0.07 --min-var-dist-3 0.1 --max-var-mmqs 100 --output-file ./$bn"_fpFilter"/$bn.fpfilter

#Keep only variants tagged as PASS by fpfilter.
grep 'PASS' ./$bn"_fpFilter"/$bn.fpfilter | awk '{OFS="\t" ; print $1,$2,$2,$3,$4,$5,$6,$7}' > ./$bn"_fpFilter"/$bn.pass.annovar

#Run annovar on passed variants
#perl  $annovar ./$bn"_fpFilter"/$bn.pass.annovar /home/csipk/NGS/annovar/humandb/ -buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2014oct_all,snp135,snp135NonFlagged,snp131,clinvar_20150330,cosmic70,ljb2_all -operation g,r,r,f,f,f,f,f,f,f,f -nastring NA --otherinfo --remove --outfile ./$bn"_fpFilter"/$bn
