#!/bin/bash
#PBS -A ONRDC17423173
#PBS -l select=2:ncpus=36:mpiprocs=36
#PBS -l walltime=01:00:00
#PBS -q standard
#PBS -j oe
#PBS -V
#PBS -N mcnpu0t4

module swap mpi/sgimpt/2.15 mpi/intelmpi/17.0.1

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE > nodes
num_per_job=72
pre="A_"

for (( i=1 ; i<2 ; i++ ))
do
	n1=$( expr $num_per_job \* \( $i + 1 \) )  #specify lines
	if (( $i <= 9 )) ; then
		j=00$i
	elif (( $i <= 99 )) ; then
		j=0$i
	else
		j=$i
	fi
	folder=$pre$j
	mkdir -p $folder
	cp -r bin  bv.py el_list.txt  io.py  main.py  mc.py  PSPs  structure.xsf  templates $folder
	cd $folder

	cat $PBS_O_WORKDIR/nodes | head -n $n1 | tail -n $num_per_job > Nodes
	( export PBS_NODEFILE=Nodes; mkdir -p temp; cp Nodes temp; python main.py structure.xsf el_list.txt >& err_info.log ) &
	sleep 2s

	cd $PBS_O_WORKDIR
done
wait
