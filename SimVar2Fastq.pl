use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw/$Bin/;

my ($fa,$reg,$samtools,$sample_name,$InsP1,$InsP2,$InsLen,$vaf,$minCovLen,$reads_num,$read_len,$frag_len,$outdir);

GetOptions(
    "fa:s" => \$fa,                     # Need
    "r:s" => \$reg,                     # Need <format is: chr1:1000-2000 or 1:1000-2000>
    "samtools:s" => \$samtools,         # Need
    "name:s" => \$sample_name,          # Need
    "InsP1:i" => \$InsP1,               # Need
    "InsP2:i" => \$InsP2,               # Need
    "InsLen:i" => \$InsLen,             # Need
    "vaf:s" => \$vaf,                   # default: 0.2 (means 20% mutant allel)
    "minlen:i" => \$minCovLen,          # default: 5
    "N:i" => \$reads_num,               # default: 5000
    "rl:i" => \$read_len,               # default: 150
    "fl:i" => \$frag_len,               # default: 260
    "od:s" => \$outdir,                 # Need
    ) or die "unknown args\n";


# check args
if (not defined $fa || not defined $reg || not defined $sample_name || not defined $InsP1 || not defined $InsP2 || not defined $InsLen || not defined $outdir){
    die "please check your args\n";
}

# default value
$reads_num = 5000 if (not defined $reads_num);
$read_len = 150 if (not defined $read_len);
$frag_len = 260 if (not defined $frag_len);

if (not defined $fa){
    $fa = "/data1/database/b37/human_g1k_v37.fasta";
}

if (not defined $samtools){
    $samtools = "/home/fulongfei/miniconda3/bin/samtools";
}

if (!-d $outdir){
    `mkdir -p $outdir`;
}

if (not defined $vaf){
    $vaf = 0.2;
}

my $runsh = "$outdir/$sample_name\.SimVar2Fastq.sh";
open O, ">$runsh" or die;
print O "cd $outdir\n";
print O "perl $Bin/bin/makeInsSeq.pl -fa $fa -reg $reg -samtools $samtools -insP1 $InsP1 -insP2 $InsP2 -insLen $InsLen -od $outdir\n\n";

print O "perl $Bin/bin/SimIns.pl -in $outdir/ins.log -n $sample_name -vaf $vaf -num $reads_num -rlen $read_len -flen $frag_len -minlen $minCovLen -od $outdir\n\n";

# aln
my $r1 = "$outdir/$sample_name\.R1.fastq";
my $r2 = "$outdir/$sample_name\.R2.fastq";

my $cmd = "bwa mem -M -R \"\@RG\\tID\:$sample_name\\tSM\:$sample_name\\tPL\:illumina\" $fa $r1 $r2 \| $samtools view -b -o $sample_name\.bam -";
print O "$cmd\n";

$cmd = "$samtools sort $sample_name\.bam -o $sample_name\.sort.Ins\.$InsLen\.bam";
print O "$cmd\n";

$cmd = "$samtools index $sample_name\.sort.Ins\.$InsLen\.bam";
print O "$cmd\n";

$cmd = "rm $sample_name\.bam";
print O "$cmd\n";

close O;





