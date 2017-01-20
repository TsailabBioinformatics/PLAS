#!/bin/bash

platform="Zcluster"
mode="paired-end"

###Only run AFTER preRun.sh has been successfully run at least once for data set!

####################### this starts with proteome #######################
## this part only needs to be run once, it prepares and formats the data
## manually run
# WU-BLAST, construct similarity matrix for ortholog group clustering
#/usr/local/wublast/latest/xdformat -p -o 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa
#/usr/local/wublast/latest/blastp 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa -o 01.data/03.MCL/01.blast/wu.blast.all.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg

## manually run
# mcl, cluster ortholog groups
#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/wu.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt wu
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/wu.mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5

## manually run
# construct meta-group, combine ortholog groups into meta-group, each group contains 1000 genes
#module load R/3.1.2    # this is for Sapelo
#time /usr/local/apps/R/3.1.2/bin/Rscript 00.script/a3.geneSelection.R 01.data/03.MCL/02.mcl/mcl.out.txt 01.data/04.GeneOfInterest/GeneID.txt 1000 Potri

## manually run
# split gene, based on the meta-group, split gene sequences accordingly
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/01.Protein/run.0 1000
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/02.Transcript/run.0 1000

## manually run
# get meta-data for the meta-group, eg. gene/protein length, which group each gene belongs to
#time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/00.PriorData/gene.gff3 01.data/04.GeneOfInterest/GeneID.v1.txt 1000

###### preprocess data

## manually run
# construct DIAMOND and blast database for the meta-group

time perl 00.script/02.makeblastdb.folder.pl 01.data/05.SplitGenes/01.Protein/run.0 prot DIAMOND $platform
wait
time perl 00.script/02.makeblastdb.folder.pl 01.data/05.SplitGenes/02.Transcript/run.0 nucl FALSE $platform
wait
# separated paired reads and single reads
time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/01.Fastq " " $platform
wait
# convert fastq file to fasta
time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta $mode $platform
wait
# simplify read IDs
time perl 00.script/01.folder.IDConverter.pl 01.data/02.Fasta $mode $platform
wait
## automatic run
:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

###################################################l
#: '
a=0
evalue=1e-3
while [ $a -le 14 ]
do
    b=`expr $a + 1`
    echo "Run number: $b" >> job.monitor.txt

    if [ $b -eq 0 ];then

        time perl 00.script/03.diamond.folder.pl 01.data/02.Fasta 01.data/05.SplitGenes/01.Protein/run.$b 03.blast/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $evalue $mode $platform 
        time perl 00.script/040.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b 1000 genome $mode $platform 20
        echo "Already been done"
    else
        time perl 00.script/021.makebowtiedb.folder.pl 01.data/05.SplitGenes/02.Transcript/run.$b $platform 10
        time perl 00.script/03.bowtie.folder.pl 01.data/02.Fasta 01.data/05.SplitGenes/02.Transcript/run.$b 03.blast/03.bowtie.nucl/run.$b nucl local bowtie.log/bowtie.run.$b $mode $platform 10

        time perl 00.script/04.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $mode $platform 20
    fi

    time perl 00.script/06.assembly.trinity.folder.pl 04.retrieve.reads/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b genome $mode $platform 20
    time perl 00.script/06.truncate.header.folder.pl 06.assembly/03.bowtie.nucl/run.$b $platform 20
    time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.$b $platform 10
    time perl 00.script/07.blastn.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/02.Transcript/run.0 07.map.back/02.blastn/run.$b $platform 10    

    c=`expr $b + 1`

    time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b 07.map.back/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
    time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/01.Protein/run.0 01.data/05.SplitGenes/03.Full.Length/run.$c pct 0.1
    cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
    cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

    a=`expr $a + 1`
done


#### summarize all runs

grep ">" 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta > 01.data/05.SplitGenes/03.Full.Length/count1
sed -i "s/>//" 01.data/05.SplitGenes/03.Full.Length/count1
perl 00.script/b3.full.length.format.pl 01.data/04.GeneOfInterest/GeneID.v1.txt 01.data/05.SplitGenes/03.Full.Length/count1 01.data/05.SplitGenes/03.Full.Length/count2 01.data/05.SplitGenes/03.Full.Length/count3
#'
#### assemble unmapped reads
cp 01.data/05.SplitGenes/03.Full.Length/count3 08.full.length/
cp 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta 08.full.length/
time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.nucl.fasta 08.full.length/Final.v1.fasta

time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1
wait
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
wait
time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1
wait
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1
wait
time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02
wait

cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../

#module load bowtie2/2.2.4
export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final

wait
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length unmap $platform
wait
time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform
wait
time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta
wait

### combine full length contigs from local assembly and assembly of unmapped reads
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.fasta -dbtype nucl
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/Trinity.new.fasta -out 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 1
wait
time perl 00.script/c11.remove.redundancy.pl 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/temp.fasta query > 10.unmapped.reads.trinity/remove.redundancy.log
wait
cat 08.full.length/Final.fasta 10.unmapped.reads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.v1.fa
#rm 10.unmapped.reads.trinity/temp.fasta


