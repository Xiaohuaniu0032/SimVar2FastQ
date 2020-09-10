cd /home/fulongfei/workdir/git_repo/SimVar2Fastq/test
perl /data1/workdir/fulongfei/git_repo/SimVar2Fastq/bin/makeInsSeq.pl -fa /data1/database/b37/human_g1k_v37.fasta -reg 13:28608000-28609000 -samtools /home/fulongfei/miniconda3/bin/samtools -insP1 28608500 -insP2 28608501 -insLen 30 -od /home/fulongfei/workdir/git_repo/SimVar2Fastq/test

perl /data1/workdir/fulongfei/git_repo/SimVar2Fastq/bin/SimIns.pl -in /home/fulongfei/workdir/git_repo/SimVar2Fastq/test/ins.log -n test -vaf 0.2 -num 1000 -rlen 150 -flen 260 -minlen 5 -od /home/fulongfei/workdir/git_repo/SimVar2Fastq/test

bwa mem -M -R "@RG\tID:test\tSM:test\tPL:illumina" /data1/database/b37/human_g1k_v37.fasta /home/fulongfei/workdir/git_repo/SimVar2Fastq/test/test.R1.fastq /home/fulongfei/workdir/git_repo/SimVar2Fastq/test/test.R2.fastq | /home/fulongfei/miniconda3/bin/samtools view -b -o test.bam -
/home/fulongfei/miniconda3/bin/samtools sort test.bam -o test.sort.Ins.30.bam
/home/fulongfei/miniconda3/bin/samtools index test.sort.Ins.30.bam
rm test.bam
