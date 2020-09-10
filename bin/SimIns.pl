use strict;
use warnings;
use Getopt::Long;

my ($infile,$name,$mutFrac,$reads_num,$read_len,$frag_len,$frag_cover_len,$outdir);

GetOptions(
    "in:s" => \$infile,                         # Need
    "n:s" => \$name,                            # Need
    "vaf:s" => \$mutFrac,                       # Default: 0.2 (means 20% mutant allel)
    "num:i" => \$reads_num,                     # Default: 5000
    "rlen:i" => \$read_len,                     # Default: 150
    "flen:i" => \$frag_len,                     # Default: 260
    "minlen:i" => \$frag_cover_len,             # Default: 5
    "od:s" => \$outdir,                         # Need
    ) or die "unknown args\n";

# default value
if (not defined $mutFrac){
    $mutFrac = 0.2; # 20% alt allel
}

if (not defined $read_len){
    $read_len = 150;
}

if (not defined $frag_len){
    $frag_len = 260;
}

if (not defined $reads_num){
    $reads_num = 5000;
}

if (not defined $frag_cover_len){
    $frag_cover_len = 5;
}


# read_len must <= frag_len
die "read len must <= frag len" if ($read_len > $frag_len);

my $outfile = "$outdir/sim.log";
open O, ">$outfile" or die;

# read origin seq and ins info
my ($ins_info, $seq, $new_seq, $ins_seq);
open IN, "$infile" or die;
while (<IN>){
    chomp;
    if (/^\#/){
        s/^\#\s//;
        $ins_info = $_;
    }
    if (/origin seq/){
        $seq = (split /\:/, $_)[1];
    }
    if (/new seq/){
        $new_seq = (split /\:/, $_)[1];
    }
    if (/ins seq/){
        $ins_seq = (split /\:/, $_)[1];
    }
}
close IN;

print O "orgin seq is: $seq\n";
print O "ins seq is: $ins_seq\n";
print O "new seq is: $new_seq\n";

# 13:28608000-28608500.Ins.30.InsPos.28608250-28608251
my $reg = (split /\./, $ins_info)[0]; # 13:28608000-28608500
my $ins_pos = (split /\./, $ins_info)[-1]; # 28608250-28608251
my $ins_len = (split /\./, $ins_info)[2]; # 30
my @ins_pos = split /\-/, $ins_pos;
my $ins_left_pos = $ins_pos[0];
my $ins_right_pos = $ins_pos[1];


# origin seq's start and end pos
my $reg_pos = (split /\:/, $reg)[1];
my @reg_pos = split /\-/, $reg_pos;
my $seq_start = $reg_pos[0];
my $seq_end = $reg_pos[1];


# make new ref seq
my $left_len = $ins_left_pos - $seq_start + 1;
my $left_seq = substr($seq,0,$left_len);
my $max_use_len = $frag_len - $frag_cover_len;
my $left_seq_r = reverse($left_seq);
my $left_seq_r_cut = substr($left_seq_r,0,$max_use_len);
my $new_ref_left = reverse($left_seq_r_cut);

my $right_seq = substr($seq,$left_len);
my $new_ref_right = substr($right_seq,0,$max_use_len);

my $new_ref_seq = $new_ref_left.$new_ref_right;

print O "new ref seq is: $new_ref_seq\n";

# make new alt seq
my $new_alt_seq = $new_ref_left.$ins_seq.$new_ref_right;
print O "new alt seq is: $new_alt_seq\n";
close O;


# make fastq
my $r1 = "$outdir/$name\.R1.fastq";
my $r2 = "$outdir/$name\.R2.fastq";

open O1, ">$r1" or die;
open O2, ">$r2" or die;

# ref/alt fragment number
my $alt_frag_n = int($reads_num * $mutFrac);
my $ref_frag_n = $reads_num - $alt_frag_n;


# make fq for ref
my $new_ref_seq_len = length($new_ref_seq);

for (my $i=1;$i<=$ref_frag_n;$i++){
    my $rand_pos = int(rand($new_ref_seq_len-$frag_len+1)); # start point to get random DNA fragment
    my $sub_seq = substr($new_ref_seq,$rand_pos,$frag_len);

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


# make fq for alt
my $new_alt_seq_len = length($new_alt_seq);

my $start = $ref_frag_n + 1;
for (my $i=$start;$i<=$reads_num;$i++){
    my $rand_pos = int(rand($new_alt_seq_len-$frag_len+1)); # start point to get random DNA fragment
    my $sub_seq = substr($new_alt_seq,$rand_pos,$frag_len);

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




