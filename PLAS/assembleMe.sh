#PBS -S /bin/bash
#PBS -q batch
#PBS -N PLAS
#PBS -l nodes=1:ppn=48:AMD
#PBS -l walltime=01:00:00
#PBS -l mem=24gb

#######################################
	time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$run 07.map.back/03.bowtie.nucl/run.$run 01.data/05.SplitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
#######################################
	time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$run 06.assembly/03.bowtie.nucl/run.$run 01.data/05.SplitGenes/01.Protein/run.0 01.data/05.SplitGenes/03.Full.Length/run.$c pct 0.1
#######################################

cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta