#PBS -S /bin/bash
#PBS -q batch
#PBS -N Master_PLAS
#PBS -l nodes=1:ppn=48:AMD
#PBS -l walltime=24:00:00
#PBS -l mem=48gb
cd $PBS_O_WORKDIR
repeats=5
counter=0
./setUp.sh
./preRun.sh

for sub in 01.data/05.SplitGenes/01.Protein/run.0/*; do 
qsub -N PLAS_$sub -v SUB=$sub, runMe.sh
done

while [ $repeats -ge $counter ]; do
	for sub in 01.data/05.Splitgenes/01.Protein/run.0/*; do
		qsub -hold_jid PLAS_* -N PLASL_$sub -v RUN=$counter, loopMe.sh
	done
	
	qsub -hold_jid PLASL* -N PLAS_$counter -v COUNTER=$counter, assembleMe.sh
	let counter=$counter+1
done

qsub -hold_jid PLAS_* -v COUNTER=$counter, assembly.sh
echo "All done!" >> 00.script/01.log/job_monitor_master.sh