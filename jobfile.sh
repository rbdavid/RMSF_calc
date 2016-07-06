#!/bin/bash
#SBATCH --job-name=RMSF.a_apo
#SBATCH --output=RMSF.a_apo.output
#SBATCH --time=96:00:00 
#SBATCH --nodes=1
#SBATCH --exclusive

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/software/usr/gcc-4.9.2/lib64"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/software/usr/hpcx-v1.2.0-292-gcc-MLNX_OFED_LINUX-2.4-1.0.0-redhat6.6/ompi-mellanox-v1.8/lib"
export PYTHON_EGG_CACHE="./"

AVG_LOC=/mnt/lustre_fs/users/rbdavid/molec_machines/dns3h/Analysis/Avg_structure/AMBER_apo/newAvg/
PDB_LOC=/mnt/lustre_fs/users/rbdavid/molec_machines/dns3h/AMBER_apo/truncated.pdb
TRAJ_LOC=/mnt/lustre_fs/users/rbdavid/molec_machines/dns3h/AMBER_apo/Truncated/
NPRODS=150
NCPUS=20

prod=1
for ((i=1;i<=2;i++))
do
	j=1
	while ((j <= $NCPUS)) && ((prod <= $NPRODS))
	do
		echo $j $i $prod
		((a=$prod+4))
		printf -v x "%03d" $prod
		printf -v y "%03d" $a
		mkdir $x.$y.RMSF
		cd $x.$y.RMSF
		time ../rmsf_calc.py $AVG_LOC $PDB_LOC $TRAJ_LOC $prod $a AMBER_apo > RMSF.output &
		cd ../
		((j=$j+1))
		((prod=$prod+5))
	done
	wait
done

mkdir 021.150.RMSF
cd 021.150.RMSF
time ../rmsf_calc.py ../../Avg_structure/ ../../truncated.pdb ../../Trajectories/ 21 150 AMBER_apo > RMSF.output 

