use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($infile,$InsP1,$InsP2,$InsLen,$outfile) = @ARGV;

my $in_name = basename $infile;
#my $od = dirname $infile;
#my $outfile = "$od/$in_name\.Ins\.$InsP1\-$InsP2\.Len\.$InsLen";


my $seq = "";
open IN, "$infile" or die;
<IN>;
while (<IN>){
    chomp;
    $seq = $seq.uc($_);
}
close IN;

print "origin seq is: $seq\n";

my $seq_start = (split /\./, $in_name)[1];
my $seq_end = (split /\./, $in_name)[2];
my $seq_len = length($seq);

# check if InsP1/InsP2 is included by seq_start/seq_end
die "error InsP1 and InsP2 args\n "if ($InsP2 - $InsP1 != 1);

if ($InsP1 < $seq_start || $InsP2 > $seq_end){
    die "$InsP1 and $InsP2 is out of origin seq's genome range ($seq_start to $seq_end)\n";
}

my @InsSeq;
my @base = qw/A T C G/;
for (my $i=1;$i<=$InsLen;$i++){
    my $int = int(rand(4)); # [0,3]
    my $base = $base[$int];
    push @InsSeq, $base;
}

my $InsSeq = "";
for (@InsSeq){
    $InsSeq = $InsSeq.$_;
}

my $len = length($InsSeq);
print "ins seq is $InsSeq (len:$len)\n";

my $left_len = $InsP1 - $seq_start + 1;
my $left_seq = substr($seq,0,$left_len);
my $right_seq = substr($seq,$left_len);
my $newseq = $left_seq.$InsSeq.$right_seq;

print "new seq is: $newseq\n";

open O, ">$outfile" or die;
print O "$newseq\n";
close O;



