# SimVar2Fastq
Simulate SNV/InDel at a specific genome position to make Fastq files

## Usage
`perl SimVar2Fastq.pl -fa /path/to/*.fasta -r 13:28608000-28609000 -samtools /path/to/samtools -name test_sample -InsP1 28608500 -InsP2 28608501 -InsLen 30 -vaf 0.2 -minlen 5 -N 1000 -rl 150 -fl 260 -od /path/to/outdir`

> at present, only support Ins type. SNV and Del modules have not been implemented. 

### Parameter specification
`-fa`: the fasta file

`-r`: the region that you will get it DNA sequence. for example: -r chr1:100-200 or -r 1:100-200

`-samtools`: the path to the binary samtools

`-name`: sample name

`-InsP1`: the left pos where you want to insert a DNA sequence

`-InsP2`: the right pos where you want to insert a DNA sequence

`-InsLen`: the length of insertion

`-vaf`: mutant allel fraction

`-minlen`: the least number of bases that the fragment must overlap with the ins position. default: 5(bp)

`-N`: how many reads you want to gererate. dafault: 5000

`-rl`: read length. default: 150

`-fl`: DNA fragment length. default: 260

`-od`: output dir

## Testing
1. `sh /path/SimVar2Fastq/test/run.sh` 
2. `sh test_sample.SimVar2Fastq.sh`

after running step2, you will get the following files:

```
-rw-r--r--  1 lffu  staff   2.1K  9 10 15:12 ins.log
-rw-r--r--  1 lffu  staff   248B  9 10 15:12 run.sh
-rw-r--r--  1 lffu  staff   3.1K  9 10 15:12 sim.log
-rw-r--r--  1 lffu  staff   308K  9 10 15:12 test.R1.fastq
-rw-r--r--  1 lffu  staff   308K  9 10 15:12 test.R2.fastq
-rw-r--r--  1 lffu  staff   1.0K  9 10 15:12 test.SimVar2Fastq.sh
-rw-r--r--  1 lffu  staff    27K  9 10 15:12 test.sort.Ins.30.bam
-rw-r--r--  1 lffu  staff    14K  9 10 15:12 test.sort.Ins.30.bam.bai
```

### File content specification
1. `run.sh` is the main shell
2. `test_sample.SimVar2Fastq.sh` is the generated file after runing `run.sh`
3. `ins.log`: contain insertion info. this file is created by `/SimVar2Fastq/bin/makeInsSeq.pl` and was used by `/SimVar2Fastq/bin/SimIns.pl`
4. `sim.log`: contain insertion info and new ref/alt seq. this file is created by `/SimVar2Fastq/bin/SimIns.pl` and was used to check the results.
5. `test.R1.fastq`: simulated read1 fastq file
6. 	`test.R2.fastq`: simulated read2 fastq file
7. `test.sort.Ins.30.bam`: aligned bam file
8. `test.sort.Ins.30.bam.bai`: indexed bam file


### Notes
* the -r format should same as your fasta file. if your fasta file's chr naming is `chr1`, then you should specify `-r chr1:100-200`
* `-InsP1` `-InsP2` must in the `-r` region. and `-InsP1` + `1` must equal `-InsP2`
* read length must <= fragment length

