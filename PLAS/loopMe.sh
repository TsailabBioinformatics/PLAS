#PBS -S /bin/bash
#PBS -q batch
#PBS -N PLAS
#PBS -l nodes=1:ppn=48:AMD
#PBS -l walltime=05:00:00
#PBS -l mem=48gb

###################
###LOOP OF DEATH###
###################
run=${RUN}
evalue=1e-3

echo "Beginning run $run" >> 00.script/01.log/job.monitor_$sub.txt

	if [ $run -gt "0" ]; then
		tgtfolder="01.data/05.SplitGenes/02.Transcript/run.$run"
		sleeptime="10"
		let prerun="$run-1"

		mkdir $tgtfolder/$sub
        cd $tgtfolder/$sub
		bowtie2-build -f -q $sub.fasta $sub                        
		cd ../../../../../

	echo "Finished 021.makebowtiedb.folder.pl!" >> 00.script/01.log/job.monitor_$sub.txt
		
################
    qryfolder="01.data/02.Fasta"
    dbfolder="01.data/05.SplitGenes/02.Transcript/run.$run"
    tgtfolder="03.blast/03.bowtie.nucl/run.$run"
    reffolder="01.data/05.SplitGenes/01.Protein/run.0"
    seqtype="nucl"
    logfolder="bowtie.log/bowtie.run.$run"
    thread="24"
    R1=()
    R2=()
    R3=()
    str2=""
    str3=""
		
   	mkdir $tgtfolder 
		
	mkdir -p $tgtfolder/$sub

		if [ "$mode" == "paired-end" ]; then
		##Construct R1 input string
			str1=""
			R1=()
			for sam in $qryfolder/*; do cd $sam
				for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
				item=$(basename ${item}); R1+=("$sam/$item")
				done; cd ../../../
			done

			for elem in "${R1[@]}"; do str1="$str1,$elem"; done
			str1=${str1#?}
			
		##Construct R2 input string
                        str2=""
			R2=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"R2"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R2+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R2[@]}"; do str2="$str2,$elem"; done
                        str2=${str2#?}
		
		##Construct singles input string
                        str3=""
			R3=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"singles"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R3+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R3[@]}"; do str3="$str3,$elem"; done
                        str3=${str3#?}
							#-unmap
			##Make a run bowtie script
			time bowtie2 -f -x $dbfolder/$db/$db -p $thread --local -1 "$str1" -2 "$str2" -S $tgtfolder/$db/bowtie.out.$db.sam
       

	         elif [ "$mode" == "Single-end" ]; then
		##Construct R1 input string
                        str1=""
			R1=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R1+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R1[@]}"; do str1="$str1,$elem"; done
                        str1=${str1#?}
						
						##FIX COMMAND GET RUNME SCRIPT MADE
							#$unmap
                        time bowtie2 -f -x $dbfolder/$db/$db -p $thread --local -U "$str1" -S $tgtfolder/$db/bowtie.out.$db.sam
                else
                        echo "No read mode selected!" >> 00.script/01.log/job.monitor_$sub.txt
                        exit
        fi


    echo 'Finished 03.bowtie.folder.pl!' >> 00.script/01.log/job.monitor_$sub.txt
						
#############################

	srcfolder="03.blast/03.bowtie.nucl/run.$run"
	tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$run"
	seqtype="nucl"
	logfolder="bowtie.log/bowtie.run.$run"
	sleeptime="20"
	unmap="map"

            mkdir -p $tgtfolder/$sub
			if [ "$mode" == "paired-end" ]; then
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.R1.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta $tgtfolder/$sub/retrieved.$sub.R2.fasta $tgtfolder/$sub/unmap.$sub.R2.fasta
				error_check "Bowtie broke at $sub"
			elif [ "$mode" == "single-end" ]; then                        
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta
			else
				echo "Error, specify read mode" >> 00.script/01.log/job.monitor_$sub.txt
				exit
            fi
			touch "00.script/04.retrieve.script/run.$run/$sub.done.log"

echo 'Finished 04.folder.retrievebowtie.reads.pl!' >> 00.script/01.log/job.monitor_$sub.txt
chmod 777 -R 00.script/04.retrieve.script/run.$run

fi

#######################################
echo "running trinity.runMe.sh ..." >> 00.script/01.log/job.monitor_$sub.txt
thread="2"
memory="24"

tgtfolder=06.assembly/03.bowtie.nucl/run."$run"/"$sub".trinity
srcfolder=04.retrieve.reads/03.bowtie.nucl/run."$run"/"$sub"
mkdir $tgtfolder

time Trinity --seqType fa --CPU $thread --max_memory "$memory"G \
		     --left "$srcfolder"/retrieved."$sub".R1.fasta \
			 --right "$srcfolder"/retrieved."$sub".R2.fasta \
			 --output $tgtfolder --min_contig_length 100
		
echo "Trinity output to $tgtfolder/$sub.trinity"

#######################################
###06.truncate.header.folder.pl

srcfolder="06.assembly/03.bowtie.nucl/run.$run"

echo "Running 06.truncate.header.folder.pl ..." >> 00.script/01.log/job.monitor_$sub.txt

time perl 00.script/06.truncate.header.pl $srcfolder/$sub/Trinity.fasta $srcfolder/$sub/Trinity.new.fasta

echo "Finished running 06.truncate.header.folder.pl!" >> 00.script/01.log/job.monitor_$sub.txt
#######################################

echo "Running 07.blastx.back.pl ..." >> 00.script/01.log/job.monitor_$sub.txt
srcfolder="06.assembly/03.bowtie.nucl/run.$run"
dbfolder="01.data/05.SplitGenes/01.Protein/run.0"
tgtfolder="07.map.back/03.bowtie.nucl/run.$run"
reffolder="01.data/05.SplitGenes/01.Protein/run.0"
thread="24"

mkdir $tgtfolder

	mkdir -p $tgtfolder/$sub
	blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-5  -outfmt 6 -num_threads 1 -max_target_seqs 1
	blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-5  -outfmt 5 -num_threads 1 -max_target_seqs 1

echo 'Finished 07.blastx.back.pl!' >> 00.script/01.log/job.monitor_$sub.txt

#######################################


echo "Running 07.blastn.back.pl ..." >> 00.script/01.log/job.monitor_$sub.txt

srcfolder="06.assembly/03.bowtie.nucl/run.$run"
dbfolder="01.data/05.SplitGenes/02.Transcript/run.0"
tgtfolder="07.map.back/02.blastn/run.$run"


mkdir $tgtfolder

	mkdir -p $tgtfolder/$sub
	blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-10  -outfmt 6 -num_threads 1 -max_target_seqs 1
	blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-10  -outfmt 5 -num_threads 1 -max_target_seqs 1

echo "Finished 07.blastn.back.pl!" >> 00.script/01.log/job.monitor_$sub.txt

#######################################
