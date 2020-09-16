# SimVar2Fastq
Simulate SNV/InDel at a specific genome position to make Fastq files

## Usage
### for tandem dup
`python3 ../SimVar2Fastq.py TD -n test -reg 13:28608000-28609000 -dupStart 28608500 -dupLen 50 -outdir $PWD`

you can see detail arguments by `python3 SimVar2Fastq.py TD -h`


```
lffu@bogon:~/git_repo/SimVar2Fastq/test$ python3 ../SimVar2Fastq.py TD -h
usage: Simulate SNV/InDel/TandemDup into Fastq files TD [-h] [-fa FASTA] [-n NAME] [-reg REGION] [-samtools SAMTOOLS]
                                                        [-dupStart DUPSTART] [-dupLen DUPLEN] [-repTimes REPTIMES]
                                                        [-vaf VAF] [-num NUM] [-rlen RLEN] [-flen FLEN]
                                                        [-outdir OUTDIR]

optional arguments:
  -h, --help          show this help message and exit
  -fa FASTA           fasta file
  -n NAME             sample name
  -reg REGION         genome region, format: <chr:start-end>
  -samtools SAMTOOLS  samtools bin
  -dupStart DUPSTART  dup start pos
  -dupLen DUPLEN      dup length
  -repTimes REPTIMES  how many times this DUP exists
  -vaf VAF            mut allele percentages
  -num NUM            how many reads
  -rlen RLEN          read len
  -flen FLEN          fragment len
  -outdir OUTDIR      outdir to store files
```


### for insertion
you can see detail arguments by `python3 SimVar2Fastq.py Ins -h`

```
lffu@bogon:~/git_repo/SimVar2Fastq/test$ python3 ../SimVar2Fastq.py Ins -h
usage: Simulate SNV/InDel/TandemDup into Fastq files Ins [-h] [-fa FASTA] [-n NAME] [-reg REGION] [-samtools SAMTOOLS]
                                                         [-insP1 INSP1] [-insP2 INSP2] [-insLen INSLEN] [-vaf VAF]
                                                         [-num NUM] [-rlen RLEN] [-flen FLEN] [-outdir OUTDIR]

optional arguments:
  -h, --help          show this help message and exit
  -fa FASTA           fasta file
  -n NAME             sample name
  -reg REGION         genome region, format: <chr:start-end>
  -samtools SAMTOOLS  samtools bin
  -insP1 INSP1        left insertion site
  -insP2 INSP2        right insertion site
  -insLen INSLEN      ins len
  -vaf VAF            mut allele percentages
  -num NUM            how many reads
  -rlen RLEN          read len
  -flen FLEN          fragment len
  -outdir OUTDIR      outdir to store files
```



> at present, only support Insertion & Tandem Duplication type. SNV and Del modules have not been implemented.


## Testing
1. `cd /path/SimVar2Fastq/test`
2. `sh cmd.sh` 
3. `sh run.sh`

after `sh run.sh`, you will get the following files:

```
lffu@bogon:~/git_repo/SimVar2Fastq/test$ ll
total 6456
-rw-r--r--  1 lffu  staff   107B  9 16 16:32 cmd.sh
-rw-r--r--  1 lffu  staff   1.2K  9 16 16:32 run.sh
-rw-r--r--  1 lffu  staff   2.1K  9 16 16:32 sim.log
-rw-r--r--  1 lffu  staff   2.1K  9 16 16:32 tandemDup.log
-rw-r--r--  1 lffu  staff   1.5M  9 16 16:32 test.R1.fastq
-rw-r--r--  1 lffu  staff   1.5M  9 16 16:32 test.R2.fastq
-rw-r--r--  1 lffu  staff   108K  9 16 16:32 test.sort.TD.50.bam
-rw-r--r--  1 lffu  staff    14K  9 16 16:32 test.sort.TD.50.bam.bai
```

### File content specification
1. `cmd.sh` is the main shell
2. `run.sh` is the generated file after runing `sh cmd.sh`
3. `sim.log` is the log info
4. `test.R1.fastq` is the read1 file
5. `test.R2.fastq` is the read2 file
6. `tandemDup.log` is the log info
7. `test.sort.TD.50.bam` is the bwa aligned bam file
8. `test.sort.TD.50.bam.bai` is the bwa index file

### Notes
* the `-reg` format should same as your fasta file. if your fasta file's chr naming is `chr1`, then you should specify `-reg chr1:100-200`


