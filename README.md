# SimVar2Fastq
Simulate SNV/InDel at a specific genome position to make Fastq files

## Usage
`perl SimVar2Fastq.pl -fa /path/to/*.fasta -r 13:28608000-28608500 -samtools /path/to/samtools -name test_sample -InsP1 28608250 -InsP2 28608251 -InsLen 30 -od /path/to/outdir`

> at present, only support Ins type. SNV and Del modules have not been implemented. 

### Parameter specification
`-fa`: the fasta file

`-r`: the region that you will get it DNA sequence. for example: -r chr1:100-200 or -r 1:100-200

`-samtools`: the path to the binary samtools

`-name`: sample name

`-InsP1`: the left pos where you want to insert a DNA sequence

`-InsP2`: the right pos where you wang to insert a DNA sequence

`-InsLen`: the length of insertion

`-N`: how many reads you want to gererate. dafault: 5000

`-rl`: read length. default: 150

`-fl`: DNA fragment length. default: 260

`-od`: output dir

## Testing
1. `sh /path/SimVar2Fastq/test/run.sh` 
2. `sh test_sample.SimVar2Fastq.sh`

after running step2, you will get the following files:

```
-rw-rw-r-- 1 fulongfei fulongfei     532 9月   4 16:02 13.28608000.28608500.fa
-rw-rw-r-- 1 fulongfei fulongfei     532 9月   4 16:02 13.28608000.28608500.fa.Ins.28608250-28608251.Len.30
-rw-rw-r-- 1 fulongfei fulongfei     212 9月   4 16:02 run.sh
-rw-rw-r-- 1 fulongfei fulongfei  271062 9月   4 16:02 test_sample.bam
-rw-rw-r-- 1 fulongfei fulongfei 1578893 9月   4 16:02 test_sample.R1.fastq
-rw-rw-r-- 1 fulongfei fulongfei 1578893 9月   4 16:02 test_sample.R2.fastq
-rw-rw-r-- 1 fulongfei fulongfei    1219 9月   4 16:02 test_sample.SimVar2Fastq.sh
-rw-rw-r-- 1 fulongfei fulongfei  101251 9月   4 16:02 test_sample.sort.bam
-rw-rw-r-- 1 fulongfei fulongfei   14728 9月   4 16:02 test_sample.sort.bam.bai
```

### File content specification
1. `run.sh` is the main shell
2. `test_sample.SimVar2Fastq.sh` is the generated file after runing `run.sh`
3. `13.28608000.28608500.fa` this file contain the DNA sequence you want to extract (-r arg)
4. `13.28608000.28608500.fa.Ins.28608250-28608251.Len.30` this file contain the new DNA sequence (add the inserted seq into the original seq). `Ins` means this file is a Insertion type, `28608250-28608251` means you want to insert a seq into this position, `Len` stand for length, `30` means the inserted seq is 30bp
5. `test_sample.R1.fastq` and `test_sample.R2.fastq` are the generated PE fastq files, based on above new DNA sequence
6. `test_sample.bam` is the origin unsorted bam
7. `test_sample.sort.bam` is the sorted bam


### Notes
* the -r format should same as your fasta file. if your fasta file's chr naming is `chr1`, then you should specify `-r chr1:100-200`
* `-InsP1` `-InsP2` must in the `-r` region. and `-InsP1` + `1` must equal `-InsP2`
* read length must <= fragment length

