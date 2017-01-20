#!/bin/bash

#for f in 00.select.ncRNA.pl 00.folder.fastqc.pl 01.folder.CutAdapter.pl 00.download.sra.pl 02.fastq2fas.pl 02.folder.fas.pl 03.folder.blat.pl 04.clrna.pl 04.folder.ncRNA.pl 05.trimmomatic.pl 06.bowtie.folder.pl 06.folder.tophat.pl 07.folder.htseq.pl 08.deseq2.R 08.folder.deseq.pl 11.ooc a1.extract.read.number.pl a1.extract.read.number.R a1.folder.extract.read.pl a2.identify.quality.scheme.pl a3.read.number.summary.pl a6.folder.dustMasker.fasta.pl a7.folder.retrieve.dustMasker.fasta.pl a7.retrieve.dustMasker.fasta.pl 02.folder.fastq2fasta.pl
for f in 01.folder.fastaCombinePairedEnd.pl 01.fastaCombinePairedEnd.py 01.fastq2fasta.folder.pl 01.folder.IDConverter.pl 021.makebowtiedb.folder.pl 02.makeblastdb.folder.pl 02.makebowtiedb.folder.pl 03.diamond.folder.pl 03.bowtie.folder.pl 040.folder.retrievebowtie.reads.pl 040.retrievebowtie.reads.pl 04.folder.retrievebowtie.reads.pl 04.retrievebowtie.reads.pl 06.assembly.trinity.folder.pl 06.truncate.header.folder.pl 06.truncate.header.pl 07.blastn.back.pl 07.blastx.back.pl 10.transfer.saturate.seq.pl 100.transfer.saturate.seq.pl 101.transfer.saturate.seq.pl a1.mcl.prepare.graph.pl a5.releventInfo.pl c11.combine.full.length.pl c12.filter.by.length.pl d11.remove.redundant.contigs.pl d12.mcl.prepare.graph.pl d13.remove.redundant.mcl.pl d14.blast.batch.pl d15.folder.get.unplant.seq.pl d15.get.unplant.seq.pl d16.remove.nonplant.seq.pl d16.remove.hostplant.seq.pl d17.transform.ID.pl c13.compare.Trinity.local.pl d18.annotation.pl z1.transform.ID.pl d19.clrna.pl 13.estimate.abundance.pl a8.summarize.expression.R b2.extract.gene.of.interest.pl b3.full.length.format.pl c9.get.full.length.seq.pl c10.folder.bowtie.full.length.pl c10.unmapped.reads.trinity.pl c11.remove.redundancy.pl 10.folder.detect.full.length.seq.pl c13.compare.Trinity.Local.pl c14.retrieveseq.Trinity.Local.pl c16.find.hybridicity.pl 10.detect.full.length.seq.pl 101.transfer.saturate.seq.nucl.pl
#for f in 101.transfer.saturate.seq.nucl.pl
do
    #cp ../../../12.Parasite/04.StHe/03.Protein.Clean/00.script/$f ./
    #unlink $f
    ln -sf ../../01.Ptrichocarpa/00.script/$f 
done
