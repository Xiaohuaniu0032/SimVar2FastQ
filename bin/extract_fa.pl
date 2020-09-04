use strict;
use warnings;
use Getopt::Long;


my ($fa,$pos,$samtools,$outdir) = @ARGV;
# pos format: chr1:1-10 or 1:1-10

my @pos = split /\:/, $pos;
my $chr = $pos[0];
my @reg = split /\-/, $pos[1];
my $s = $reg[0];
my $e = $reg[1];

my $outfa = "$outdir\/$chr\.$s\.$e\.fa";

my $cmd = "$samtools faidx $fa $pos >$outfa";
system($cmd) == 0 or die "failed to extract fa sequence\n";

