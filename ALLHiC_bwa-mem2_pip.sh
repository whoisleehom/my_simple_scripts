#!/bin/bash
usage()
{
	echo "    Usage: `basename $0` -r reference -1 R1.fq -2 R2.fq -k group_count [-e enzyme] [-t threads] [-b bin_size]"
	echo "          -r: reference genome"
	echo "          -1: Lib_R1.fq.gz"
	echo "          -2: Lib_R2.fq.gz"
	echo "          -k: group_count"
	echo "          -e: enzyme_sites (HindIII: AAGCTT; MboI: GATC), default: HindIII"
	echo "          -t: threads, default: 10"
	echo "          -b: bin_size for hic heatmap, can be divided with comma, default: 500k"
	exit 0
}
### edit your bwa-mem2 path ###
export BWAMEM2='/share/home/stu_gaoyuan/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2'
### get options
while getopts ':r:1:2:k:e:t:b:' OPT; do
	case $OPT in
		r)
			ref="$OPTARG";;
		1)
			R1="$OPTARG";;
		2)
			R2="$OPTARG";;
		e)
			enzyme="$OPTARG";;
		k)
			group_count="$OPTARG";;
		t)
			threads="$OPTARG";;
		b)
			bin_size="$OPTARG";;
		?)
			usage;;
	esac
done
### check required variants
if [ -z $ref ] || [ -z $R1 ] || [ -z $R2 ] || [ -z $group_count ]; then
	usage
fi
### set default values while optional variants were not set
if [ -z $threads ]; then
	threads=10
fi
if [ -z $bin_size ]; then
	bin_size=500k
fi
if [ -z $enzyme ]; then
	enzyme=AAGCTT
fi
enzyme=`echo $enzyme | tr '[a-z]' '[A-Z]'`
if [ $enzyme = HINDIII ]; then
	enzyme=AAGCTT
fi
if [ $enzyme = MBOI ]; then
	enzyme=GATC
fi
### step_1 link required files
ln -s ${ref} ./seq.fasta
ln -s ${R1} ./Lib_R1.fastq.gz
ln -s ${R2} ./Lib_R2.fastq.gz
### step_2 1st round of mapping
$BWAMEM2 index seq.fasta
samtools faidx seq.fasta
$BWAMEM2 mem -M -t $threads -K 10000000 seq.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz | samtools sort --reference seq.fasta -O bam -o sorted.bam -@ $threads -
### step_3 uncorrected assembly
mkdir uncorrected_assembly;cd uncorrected_assembly
ln -s ../sorted.bam
ln -s ../seq.fasta* ./ 
ALLHiC_partition -r seq.fasta -e $enzyme -k $group_count -b sorted.bam
for ((I=1;I<=$group_count;I++));do j=$(printf "%02s" $I | tr " " "0");mv sorted.counts_GATC.${group_count}g${I}.txt group${j}.txt;echo "allhic optimize group${j}.txt sorted.clm" >> cmd.list;done
ParaFly -c cmd.list -CPU $threads
while [[ -e FailedCommands ]]; do rm cmd.list; rm cmd.list.completed; mv FailedCommands cmd.list; ParaFly -c cmd.list -CPU $threads; done
rm cmd.list; rm cmd.list.completed
### step_4 first build and plot
ALLHiC_build seq.fasta
perl ~/scripts/getFaLen.pl -i groups.asm.fasta -o len.txt
grep group len.txt > chrn.list
ALLHiC_plot sorted.bam groups.agp chrn.list $bin_size pdf
cd ..
### step_5 correct contig
samtools index sorted.bam
ALLHiC_corrector -m sorted.bam -r seq.fasta -o seq.HiCcorrected.fasta -t $threads
### step_6 2nd round of mapping
$BWAMEM2 index seq.HiCcorrected.fasta
samtools faidx seq.HiCcorrected.fasta
$BWAMEM2 mem -M -t $threads -K 10000000 seq.HiCcorrected.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz | samtools sort --reference seq.HiCcorrected.fasta -O bam -o sample.bwa_mem.bam -@ $threads -
### step_7 assembly
mkdir assembly;cd assembly
ln -s ../seq.HiCcorrected.fasta* ./
ln -s ../sample.bwa_mem.bam
ALLHiC_partition -r seq.HiCcorrected.fasta -e $enzyme -k $group_count -b sample.bwa_mem.bam
for ((I=1;I<=$group_count;I++));do j=$(printf "%02s" $I | tr " " "0");mv sample.bwa_mem.counts_GATC.${group_count}g${I}.txt group${j}.txt;echo "allhic optimize group${j}.txt sample.bwa_mem.clm" >> cmd.list;done
ParaFly -c cmd.list -CPU $threads
while [[ -e FailedCommands ]]; do rm cmd.list; rm cmd.list.completed; mv FailedCommands cmd.list; ParaFly -c cmd.list -CPU $threads; done
rm cmd.list; rm cmd.list.completed
### step_8 second build and plot
ALLHiC_build seq.HiCcorrected.fasta
perl ~/scripts/getFaLen.pl -i groups.asm.fasta -o len.txt
grep group len.txt > chrn.list
ALLHiC_plot sample.bwa_mem.bam groups.agp chrn.list $bin_size pdf
cd ..
