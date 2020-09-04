use strict;
use warnings;
use Getopt::Long;

my ($infile,$name,$reads_num,$read_len,$frag_len,$outdir);

GetOptions(
    "in:s" => \$infile,           # Need
    "n:s" => \$name,              # Need
    "num:i" => \$reads_num,       # Default: 5000
    "rlen:i" => \$read_len,       # Default: 150
    "flen:i" => \$frag_len,       # Default: 260
    "od:s" => \$outdir,           # Need
    ) or die "unknown args\n";

# default value
if (not defined $read_len){
    $read_len = 150;
}

if (not defined $frag_len){
    $frag_len = 260;
}

if (not defined $reads_num){
    $reads_num = 5000;
}

# read_len must <= frag_len
die "read len must <= frag len" if ($read_len > $frag_len);

my $seq = "";
open IN, "$infile" or die;
while (<IN>){
    chomp;
    $seq = $seq.uc($_);
}
close IN;


my $len = length($seq);

my $r1 = "$outdir/$name\.R1.fastq";
my $r2 = "$outdir/$name\.R2.fastq";

open O1, ">$r1" or die;
open O2, ">$r2" or die;

for (my $i=1;$i<=$reads_num;$i++){
    # get sub-frag
    my $rand_pos = int(rand($len-$frag_len+1));
    my $sub_seq = substr($seq,$rand_pos,$frag_len);

    my $r1 = substr($sub_seq,0,$read_len);
    my $rc_sub_seq = &rc_seq($sub_seq); # reverse complement
    my $r2 = substr($rc_sub_seq,0,$read_len);

    my $base_Q = "F" x $read_len; # illumina 1.8+, Phred+33, F means Q38

    print O1 "\@seq$i R1\n";
    print O1 "$r1\n";
    print O1 "+\n";
    print O1 "$base_Q\n";

    print O2 "\@seq$i R2\n";
    print O2 "$r2\n";
    print O2 "+\n";
    print O2 "$base_Q\n";
}

close O1;
close O2;

sub rc_seq{
    my ($seq) = @_;
    my $rev = reverse($seq);
    my $rev_c = "";
    my @rev = split //, $rev;
    for my $base (@rev){
        my $new_base;
        if ($base eq "A"){
            $new_base = "T";
        }elsif ($base eq "T"){
            $new_base = "A";
        }elsif ($base eq "C"){
            $new_base = "G";
        }else{
            $new_base = "C";
        }
        $rev_c = $rev_c.$new_base;
    }

    return($rev_c);
}





