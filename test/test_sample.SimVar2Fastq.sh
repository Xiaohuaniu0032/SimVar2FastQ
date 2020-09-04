cd /home/fulongfei/workdir/git_repo/SimVar2FastQ/test
/home/fulongfei/miniconda3/bin/samtools faidx /data1/database/b37/human_g1k_v37.fasta 13:28608000-28608500 >/home/fulongfei/workdir/git_repo/SimVar2FastQ/test/13.28608000.28608500.fa
perl /data1/workdir/fulongfei/git_repo/SimVar2FastQ/bin/Intro_Ins.pl /home/fulongfei/workdir/git_repo/SimVar2FastQ/test/13.28608000.28608500.fa 28608250 28608251 30 /home/fulongfei/workdir/git_repo/SimVar2FastQ/test/13.28608000.28608500.fa.Ins.28608250-28608251.Len.30
perl /data1/workdir/fulongfei/git_repo/SimVar2FastQ/bin/SimIns.pl -in /home/fulongfei/workdir/git_repo/SimVar2FastQ/test/13.28608000.28608500.fa.Ins.28608250-28608251.Len.30 -n test_sample -od /home/fulongfei/workdir/git_repo/SimVar2FastQ/test
bwa mem -R "@RG\tID:test_sample\tSM:test_sample\tPL:illumina" /data1/database/b37/human_g1k_v37.fasta /home/fulongfei/workdir/git_repo/SimVar2FastQ/test/test_sample.R1.fastq /home/fulongfei/workdir/git_repo/SimVar2FastQ/test/test_sample.R2.fastq | /home/fulongfei/miniconda3/bin/samtools view -b -o test_sample.bam -
/home/fulongfei/miniconda3/bin/samtools sort test_sample.bam -o test_sample.sort.bam
/home/fulongfei/miniconda3/bin/samtools index test_sample.sort.bam
