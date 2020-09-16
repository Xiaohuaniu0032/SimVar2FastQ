use strict;
use warnings;
use Getopt::Long;

my ($fa,$region,$samtools,$dupStart,$dupLen,$repTimes,$outdir);

GetOptions(
    "fa:s" => \$fa,                    # Optional
    "reg:s" => \$region,               # NEED
    "samtools:s" => \$samtools,        # Optional
    "dupStart:i" => \$dupStart,        # NEED
    "dupLen:i" => \$dupLen,            # NEED
    "repTimes:i" => \$repTimes,        # Optional
    "outdir:s" => \$outdir,            # NEED
    ) or die;


# default value
if (not defined $fa){
    $fa = "/data1/database/b37/human_g1k_v37.fasta";
}

if (not defined $samtools){
    $samtools = "/home/fulongfei/miniconda3/bin/samtools";
}

if (not defined $repTimes){
    $repTimes = 1;
}

# check region
# PASS

# extract seq
my $outfile = "$outdir/tandemDup.log";
open O, ">$outfile" or die;
my $dupEndPos = $dupStart + $dupLen - 1;
print O "\# $region\.TandemDup\.$dupLen\.DupPos\.$dupStart\-$dupEndPos\n";

my $cmd = "$samtools faidx $fa $region >$outdir/origin.seq";
system($cmd) == 0 or die "samtools faidx failed\n";

# read origin seq
my $origin_seq = "";
open IN, "$outdir/origin.seq" or die;
<IN>;
while (<IN>){
    chomp;
    $origin_seq = $origin_seq.uc($_);
}
close IN;

`rm $outdir/origin.seq`;

print O "origin seq:$origin_seq\n";

# get dup seq
my $chr = (split /\:/, $region)[0];
my $dup_region = "$chr\:$dupStart\-$dupEndPos";
$cmd = "$samtools faidx $fa $dup_region >$outdir/dup.seq";
system($cmd) == 0 or die "samtools faidx failed\n";

# read dup seq
my $dup_seq = "";
open IN, "$outdir/dup.seq" or die;
<IN>;
while (<IN>){
    chomp;
    $dup_seq = $dup_seq.uc($_);
}
close IN;

$dup_seq = $dup_seq x $repTimes;

`rm $outdir/dup.seq`;

print O "dup seq:$dup_seq\n";


# make sure that the start pos of dup is among the specified region <-r>


# make new seq
$chr = (split /\:/, $region)[0];
my $reg = (split /\:/, $region)[1];
my @reg = split /\-/, $reg;
my $start_pos = $reg[0];
my $end_pos = $reg[1];

my $left_len = $dupEndPos - $start_pos + 1;
my $left_seq = substr($origin_seq,0,$left_len);
my $right_seq = substr($origin_seq,$left_len);
my $new_seq = $left_seq.$dup_seq.$right_seq;

print O "new seq:$new_seq\n";
close O;


