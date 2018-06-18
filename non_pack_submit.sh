#!/bin/bash

pre="B_"
for (( i=$1 ; i<=$2 ; i++ ))
do
	wait_time=$( expr $i \* 2 )
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
	cat > runscript <<EOF
#!/bin/bash
#PBS -A ONRDC17423173
#PBS -l select=8:ncpus=40:mpiprocs=36
#PBS -l walltime=24:00:00
#PBS -q standard
#PBS -j oe
#PBS -V
#PBS -N gcmc-$i

cd \$PBS_O_WORKDIR
sleep ${wait_time}s
python main.py structure.xsf el_list.txt >& err_info.log
EOF
	qsub runscript
	cd ../
done
