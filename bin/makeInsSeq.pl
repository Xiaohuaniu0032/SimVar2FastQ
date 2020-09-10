use strict;
use warnings;
use Getopt::Long;

my ($fa,$region,$samtools,$insP1,$insP2,$insLen,$outdir);

GetOptions(
    "fa:s" => \$fa,                    # Optional
    "reg:s" => \$region,               # NEED
    "samtools:s" => \$samtools,        # Optional
    "insP1:i" => \$insP1,              # NEED
    "insP2:i" => \$insP2,              # NEED
    "insLen:i" => \$insLen,            # NEED
    "od:s" => \$outdir,                # NEED
    ) or die "unknown args\n";

# default value
if (not defined $fa){
    $fa = "/data1/database/b37/human_g1k_v37.fasta";
}

if (not defined $samtools){
    $samtools = "/home/fulongfei/miniconda3/bin/samtools";
}

# extract seq
my $outfile = "$outdir/ins.log";

open O, ">$outfile" or die;
print O "\# $region\.Ins\.$insLen\.InsPos\.$insP1\-$insP2\n";

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

# make new seq
my @InsSeq;
my @base = qw/A T C G/;
for (my $i=1;$i<=$insLen;$i++){
    my $int = int(rand(4)); # [0,3]
    my $base = $base[$int];
    push @InsSeq, $base;
}

my $InsSeq = "";
for (@InsSeq){
    $InsSeq = $InsSeq.$_;
}

print O "ins seq:$InsSeq\n";

my $pos = (split /\:/, $region)[1];
my @pos = split /\-/, $pos;
my $seq_start = $pos[0];
my $seq_end = $pos[1];

# check is Ins pos is included by -region
die "insP2 - insP1 must == 1\n" if ($insP2 - $insP1 != 1);

if ($insP1 < $seq_start || $insP2 > $seq_end){
    die "ins pos is out of extracted region\n";
}

my $left_len = $insP1 - $seq_start + 1;
my $left_seq = substr($origin_seq,0,$left_len);
my $right_seq = substr($origin_seq,$left_len);
my $new_seq = $left_seq.$InsSeq.$right_seq;

print O "new seq:$new_seq\n";
close O;



