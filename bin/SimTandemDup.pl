use strict;
use warnings;
use Getopt::Long;

my ($sim_log_file,$name,$mutFrac,$reads_num,$read_len,$frag_len,$outdir);

GetOptions(
    "in:s" => \$sim_log_file,                # NEED
    "n:s" => \$name,                         # NEED
    "vaf:f" => \$mutFrac,                    # Default: 0.2 (means 20% mutant allel)
    "num:i" => \$reads_num,                  # Default: 5000
    "rlen:i" => \$read_len,                  # Default: 150
    "flen:i" => \$frag_len,                  # Default: 260
    "outdir:s" => \$outdir,                  # NEED
    ) or die "unknown args\n";

# default value
if (not defined $mutFrac){
    $mutFrac = 0.2;
}

if (not defined $reads_num){
    $reads_num = 5000;
}

if (not defined $read_len){
    $read_len = 150;
}

if (not defined $frag_len){
    $frag_len = 260;
}

my $outfile = "$outdir/sim.log";
open O, ">$outfile" or die;

# read origin seq and dup info
my ($dup_info,$seq,$new_seq,$dup_seq);
open IN, "$sim_log_file" or die;
while (<IN>){
    chomp;
    if (/^\#/){
        s/^\#\s//;
        $dup_info = $_;
    }

    if (/origin seq/){
        $seq = (split /\:/, $_)[1];
    }

    if (/new seq/){
        $new_seq = (split /\:/, $_)[1];
    }

    if (/dup seq/){
        $dup_seq = (split /\:/, $_)[1];
    }
}
close IN;

print O "orgin seq is: $seq\n";
print O "dup seq is: $dup_seq\n";
print O "new seq is: $new_seq\n";
close O;


# make fastq

my $alt_frag_n = int($reads_num * $mutFrac);
my $ref_frag_n = $reads_num - $alt_frag_n;

my $r1 = "$outdir/$name\.R1.fastq";
my $r2 = "$outdir/$name\.R2.fastq";

open O1, ">$r1" or die;
open O2, ">$r2" or die;

# make fq for ref

my $ref_len = length($seq);
for (my $i=1;$i<=$ref_frag_n;$i++){
    my $rand_pos = int(rand($ref_len-$frag_len+1)); # start point to get random DNA fragment
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

my $start = $ref_frag_n + 1;
my $alt_len = length($new_seq);
for (my $i=$start;$i<=$reads_num;$i++){
    my $rand_pos = int(rand($alt_len-$frag_len+1)); # start point to get random DNA fragment
    my $sub_seq = substr($new_seq,$rand_pos,$frag_len);

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


