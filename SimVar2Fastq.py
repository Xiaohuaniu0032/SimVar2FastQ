import os
import sys
import argparse




def parse_args():

    AP = argparse.ArgumentParser("Simulate SNV/InDel/TandemDup into Fastq files")
    AP_subparsers = AP.add_subparsers()

    P_TandemDup = AP_subparsers.add_parser('TD')
    P_TandemDup.add_argument('-fa',help="fasta file",dest="fasta",default='/data1/database/b37/human_g1k_v37.fasta')
    P_TandemDup.add_argument('-n',help="sample name",dest="name")
    P_TandemDup.add_argument('-reg',help="genome region, format: <chr:start-end>",dest="region")
    P_TandemDup.add_argument('-samtools',help="samtools bin",dest="samtools",default='/home/fulongfei/miniconda3/bin/samtools')
    P_TandemDup.add_argument('-dupStart',help="dup start pos",dest="dupStart")
    P_TandemDup.add_argument('-dupLen',help="dup length",dest="dupLen")
    P_TandemDup.add_argument('-repTimes',help="how many times this DUP exists",default=1)
    P_TandemDup.add_argument('-vaf',help="mut allele percentages",default=0.2)
    P_TandemDup.add_argument('-num',help="how many reads",default=5000)
    P_TandemDup.add_argument('-rlen',help="read len",default=150)
    P_TandemDup.add_argument('-flen',help="fragment len",default=260)
    P_TandemDup.add_argument('-outdir',help="outdir to store files",dest="outdir")

    P_Ins = AP_subparsers.add_parser('Ins')
    P_Ins.add_argument('-fa',help="fasta file",dest="fasta",default='/data1/database/b37/human_g1k_v37.fasta')
    P_Ins.add_argument('-n',help="sample name",dest="name")
    P_Ins.add_argument('-reg',help="genome region, format: <chr:start-end>",dest="region")
    P_Ins.add_argument('-samtools',help="samtools bin",dest="samtools",default='/home/fulongfei/miniconda3/bin/samtools')
    P_Ins.add_argument('-insP1',help="left insertion site",dest='insP1')
    P_Ins.add_argument('-insP2',help="right insertion site",dest='insP2')
    P_Ins.add_argument('-insLen',help="ins len",dest='insLen')
    P_Ins.add_argument('-vaf',help="mut allele percentages",default=0.2)
    P_Ins.add_argument('-num',help="how many reads",default=5000)
    P_Ins.add_argument('-rlen',help="read len",default=150)
    P_Ins.add_argument('-flen',help="fragment len",default=260)
    P_Ins.add_argument('-outdir',help="outdir to store files",dest="outdir")



    return AP.parse_args()



if __name__ == '__main__':
    args = parse_args()
    runsh = args.outdir + '/run.sh'
    of = open(runsh,'w')
    abs_path = os.path.split(os.path.realpath(__file__))[0]
    
    if sys.argv[1] == "TD":
        cmd = 'perl %s/bin/makeTandemDupSeq.pl -fa %s -reg %s -samtools %s -dupStart %s -dupLen %s -repTimes %s -outdir %s' % (abs_path,args.fasta,args.region,args.samtools,args.dupStart,args.dupLen,args.repTimes,args.outdir)
        of.write(cmd+'\n')

        dup_log = "%s/tandemDup.log" % (args.outdir)
        cmd = 'perl %s/bin/SimTandemDup.pl -in %s -n %s -vaf %s -num %s -rlen %s -flen %s -outdir %s' % (abs_path,dup_log,args.name,args.vaf,args.num,args.rlen,args.flen,args.outdir)
        of.write(cmd+'\n')

    if sys.argv[1] == "Ins":
        cmd = 'perl %s/bin/makeInsSeq.pl -fa %s -reg %s -samtools %s -insP1 %s -insP2 %s -od %s' % (abs_path,args.fasta,args.region,args.samtools,args.insP1,args.insP2,args.outdir)
        of.write(cmd+'\n')

        ins_log = "%s/ins.log" % (args.outdir)
        cmd = 'perl %s/bin/SimIns.pl -in %s -n %s -vaf %s -num %s -rlen %s -flen %s -minlen %s -od %s' % (abs_path,ins_log,args.name,args.vaf,args.num,args.rlen,args.flen,args.minlen,args.outdir)
        of.write(cmd+'\n')


    # aln
    bam = "%s/%s.bam" % (args.outdir,args.name)
    r1 = "%s/%s.R1.fastq" % (args.outdir,args.name)
    r2 = "%s/%s.R2.fastq" % (args.outdir,args.name)
    cmd = 'bwa mem -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:illumina\" %s %s %s | %s view -b -o %s -' % (args.name,args.name,args.fasta,r1,r2,args.samtools,bam)
    of.write(cmd+'\n')

    bam_sort = "%s/%s.sort.bam" % (args.outdir,args.name)
    cmd = '%s sort %s -o %s' % (args.samtools,bam,bam_sort)
    of.write(cmd+'\n')

    cmd = '%s index %s' % (args.samtools,bam_sort)
    of.write(cmd+'\n')

    cmd = 'rm %s' % (bam)
    of.write(cmd+'\n')

  

    of.close()

    

