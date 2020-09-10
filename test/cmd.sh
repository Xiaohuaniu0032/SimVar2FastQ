samtools view test.sort.Ins.30.bam|awk '$6~/I/'|awk '{print $1}' >seq.ins
samtools view test.sort.Ins.30.bam|awk '$6~/S/'|awk '{print $1}' >seq.softclip

cat seq.ins seq.softclip|sort|uniq -c|wc -l

