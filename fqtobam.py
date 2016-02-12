#!/usr/bin/python

__author__ = 'Anand Mayakonda'

import argparse, datetime, subprocess, os, shutil


############## little function to make a read groups ##############
def rg_maker(rg_id, rg_sm, rg_lb, rg_pl, rg_pu, rg_dt):
    rg = "@RG\\tID:" + rg_id + "\\tSM:" + rg_sm + "\\tLB:" + rg_lb + "\\tPL:" + rg_pl + "\\tPU:" + rg_pu + "\\tDT:" + rg_dt
    return rg


############# source files for reference genome and gatk bundle kit ###############
ref = "/home/csipk/NGS/ref_genomes/hg19.fa"
gatk = "/home/csipk/NGS/GenomeAnalysisTK.jar"

# one can obtain these from gatk ftp (knonw as gatk bundle)
# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
millsIndels = "/home/csipk/NGS/ref_genomes/Mills_and_1000G_gold_standard.indels.hg19.vcf"
KGIndels = "/home/csipk/NGS/ref_genomes/1000G_phase1.indels.hg19.vcf"
dbSNP138 = "/home/csipk/NGS/ref_genomes/dbSNP_138.vcf"
# exome bait file in list format
# cat baits.bed | awk '{print $1":"$2"-"$3}' > exome_baits.list

regions = "/home/csipk/NGS/ref_genomes/exome_baits/agilent_v5_v3_merged_baits.list"  # exome sure select 50mb intervals
########################################################################################


parser = argparse.ArgumentParser(
    description="aligns given two sequences using bwa mem and processes them according to GATK best practices.")

# Positional arguments
parser.add_argument('sample_id', help='rg_id')
parser.add_argument('sample_name', help='rg_sm')
parser.add_argument('sample_library', help='rg_lb')
parser.add_argument('fq1', help='forward reads')
parser.add_argument('fq2', help='reverse reads')

# Optional arguments
parser.add_argument('-t', '--threads', help='threads to use [4]', default=4)
parser.add_argument('-pl', '--platform', help='read group platform [Illumina]', default='Illumina')
parser.add_argument('-pu', '--platform_unit', help='read group platform unit [Unknown]', default='Unknown')
parser.add_argument('-lb', '--library_strategy', help='strategy used for making library [WXS]', default='WXS',
                    choices=['WXS', 'WGS'], dest="lb")

# samblaster arguments
parser.add_argument('-rd', '--remove_dups', help='if specified, removes PCR duplicates instead marking them', dest='rd',
                    action='store_true')

# just align (no gatk best practices)
parser.add_argument('-ja', '--just_align',
                    help='if specified, bam files will not be processed according to gatk best practices', dest='gatk',
                    action='store_false')

# align and do indel realignment(no gatk best practices)
parser.add_argument('-indel', '--indel_realign',
                    help='if specified, alignment and indel realignment step is performed - BQSR ignored', dest='indel',
                    action='store_false')

args = parser.parse_args()

# date in iso8601 format
dt = str(datetime.date.today().isoformat())

# make read group
rg = rg_maker(args.sample_id, args.sample_name, args.sample_library, args.platform, args.platform_unit, dt)

# create an output directory
opdir = args.sample_id + "_log"

if not os.path.exists(opdir):
    os.mkdir(opdir)

print("Forward reads\t%s" % args.fq1)
print("Reverse reads\t%s" % args.fq2)
print("Using %s threads\n" % args.threads)

#print(args.gatk)
#print(args.indel)

if args.rd:
    print('[%s]\tAligning and Removing duplicates \n' % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

    alignCommand = 'bwa mem -M -v 1 -t %s -R \'%s\' %s %s %s 2>%s/%s | samblaster -r 2>%s/%s | samtools view -h -@ %s -bS - | sambamba sort -t %s -o %s /dev/stdin' % (
        args.threads, rg, ref, args.fq1, args.fq2, opdir, args.sample_id + "_bwa.log",
        opdir, args.sample_id + "_samblaster.log", args.threads, args.threads, args.sample_id + ".bam")

else:
    print('[%s]\tAligning and Marking duplicates. \n' % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

    alignCommand = 'bwa mem -M -v 1 -t %s -R \'%s\' %s %s %s 2>%s/%s | samblaster 2>%s/%s | samtools view -h -@ %s -bS - | sambamba sort -t %s -o %s /dev/stdin' % (
        args.threads, rg, ref, args.fq1, args.fq2, opdir, args.sample_id + "_bwa.log",
        opdir, args.sample_id + "_samblaster.log", args.threads, args.threads, args.sample_id + ".bam")


os.system(alignCommand)
#print(alignCommand)

cat = "cat " + opdir + "/" + args.sample_id + "_samblaster.log"
os.system(cat)

if args.gatk:

    print('[%s]\tCreating targets for indel realignment.\n' % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

    if args.lb == 'WXS':
        targetCreatorCommand = 'java -d64 -Xmx12g -jar %s -T RealignerTargetCreator -R %s -L %s -ip 50 -I %s -o %s -nt %s -known %s -known %s 2>%s/%s' % (
            gatk, ref, regions, args.sample_id + ".bam", args.sample_id + ".intervals", args.threads, millsIndels,
            KGIndels,
            opdir, args.sample_id + ".indel.log")

    else:
        targetCreatorCommand = 'java -d64 -Xmx12g -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s -known %s -known %s 2>%s/%s' % (
            gatk, ref, args.sample_id + ".bam", args.sample_id + ".intervals", args.threads, millsIndels, KGIndels,
            opdir, args.sample_id + ".indel.log")

    os.system(targetCreatorCommand)
    #print(targetCreatorCommand)

    print('[%s]\tPerforming  indel realignment.\n' % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

    if args.lb == 'WXS':
        indelAlignerCommand = 'java -d64 -Xmx12g -jar %s -T IndelRealigner -R %s -L %s -ip 50 -I %s -targetIntervals %s -known %s -known %s -o %s 2>%s/%s' % (
            gatk, ref, regions, args.sample_id + ".bam",
            args.sample_id + ".intervals", millsIndels, KGIndels, args.sample_id + "_processed.bam", opdir,
            args.sample_id + ".indel2.log")

    else:
        indelAlignerCommand = 'java -d64 -Xmx12g -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -known %s -known %s -o %s 2>%s/%s' % (
            gatk, ref, args.sample_id + ".bam",
            args.sample_id + ".intervals", millsIndels, KGIndels, args.sample_id + "_processed.bam", opdir,
            args.sample_id + ".indel2.log")


    os.system(indelAlignerCommand)
    #print(indelAlignerCommand)

    if args.indel:

        print("[%s]\tPerforming  base quality recalibration.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

        if args.lb == 'WXS':
            bqsrCommand = 'java -Xmx12g -d64 -jar %s -T BaseRecalibrator -L %s -ip 50 -I %s -R %s -knownSites %s -knownSites %s -knownSites %s -o %s  2>%s/%s' % (
                gatk, regions, args.sample_id + "_processed.bam", ref, KGIndels, millsIndels, dbSNP138,
                args.sample_id + "_recal.table", opdir, args.sample_id + '.BQSR.log')

        else:
            bqsrCommand = 'java -Xmx12g -d64 -jar %s -T BaseRecalibrator -I %s -R %s -knownSites %s -knownSites %s -knownSites %s -o %s  2>%s/%s' % (
                gatk, args.sample_id + "_processed.bam", ref, KGIndels, millsIndels, dbSNP138,
                args.sample_id + "_recal.table", opdir, args.sample_id + '.BQSR.log')


        os.system(bqsrCommand)
        #print(bqsrCommand)

        print("[%s]\tPrinting final bam file.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

        printReads = 'java -d64 -Xmx12g -jar %s -T PrintReads -R %s -I %s -nct %s -BQSR %s -o %s 2>%s/%s' % (
            gatk, ref, args.sample_id + "_processed.bam", args.threads, args.sample_id + "_recal.table",
            args.sample_id + "_recal.bam", opdir, args.sample_id + '.BQSR2.log')

        os.system(printReads)
        #print(printReads)

        print("[%s]\tRemoving intermediate files and cleaning.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

        os.remove(args.sample_id + "_processed.bam")
        os.remove(args.sample_id + "_processed.bai")

        shutil.move(args.sample_id + ".intervals", opdir+"/"+args.sample_id + ".intervals")
        shutil.move(args.sample_id + "_recal.table", opdir+"/"+args.sample_id + "_recal.table")

        shutil.move(args.sample_id+"_recal.bam", args.sample_id+".bam")
        shutil.move(args.sample_id+"_recal.bai", args.sample_id+".bam.bai")

        print("[%s]\tFinished.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

    else:
        print("[%s]\tRemoving intermediate files and cleaning.\n" % datetime.datetime.now().strftime(
            "%d-%m-%Y %H:%M"))

        shutil.move(args.sample_id + ".intervals", opdir+"/"+args.sample_id + ".intervals")
        shutil.move(args.sample_id + "_recal.table", opdir+"/"+args.sample_id + "_recal.table")

        shutil.move(args.sample_id+"_processed.bam", args.sample_id+".bam")
        shutil.move(args.sample_id+"_processed.bai", args.sample_id+".bam.bai")

        print("[%s]\tFinished.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))

else:
    print("[%s]\tFinished.\n" % datetime.datetime.now().strftime("%d-%m-%Y %H:%M"))
