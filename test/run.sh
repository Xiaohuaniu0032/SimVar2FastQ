perl ../SimVar2Fastq.pl -fa /data1/database/b37/human_g1k_v37.fasta -r 13:28608000-28609000 -samtools /home/fulongfei/miniconda3/bin/samtools -name test -InsP1 28608500 -InsP2 28608501 -InsLen 30 -vaf 0.2 -minlen 5 -N 1000 -rl 150 -fl 260 -od $PWD
