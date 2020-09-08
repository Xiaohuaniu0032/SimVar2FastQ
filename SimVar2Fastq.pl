use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw/$Bin/;

my ($fa,$reg,$samtools,$sample_name,$InsP1,$InsP2,$InsLen,$reads_num,$read_len,$frag_len,$outdir);

GetOptions(
    "fa:s" => \$fa,                     # Need
    "r:s" => \$reg,                     # Need <format is: chr1:1000-2000 or 1:1000-2000>
    "samtools:s" => \$samtools,         # Need
    "name:s" => \$sample_name,          # Need
    "InsP1:i" => \$InsP1,               # Need
    "InsP2:i" => \$InsP2,               # Need
    "InsLen:i" => \$InsLen,             # Need
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

if (!-d $outdir){
    `mkdir -p $outdir`;
}

my $runsh = "$outdir/$sample_name\.SimVar2Fastq.sh";
open O, ">$runsh" or die;
print O "cd $outdir\n";

# s1. extract target region DNA sequence
my @pos = split /\:/, $reg;
my $chr = $pos[0];
my @reg = split /\-/, $pos[1];
my $s = $reg[0];
my $e = $reg[1];
my $outfile = "$outdir/$chr\.$s\.$e\.fa";
print O "$samtools faidx $fa $reg >$outfile\n";



# s2. get new DNA sequence
my $in_name = basename $outfile;
my $od = dirname $outfile;
$outfile = "$od/$in_name\.Ins\.$InsP1\-$InsP2\.Len\.$InsLen";
print O "perl $Bin/bin/Intro_Ins.pl $outdir/$chr\.$s\.$e\.fa $InsP1 $InsP2 $InsLen $outfile\n";

# s3. simulate N reads
print O "perl $Bin/bin/SimIns.pl -in $outfile -n $sample_name -od $outdir\n";

# s4. aln
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





