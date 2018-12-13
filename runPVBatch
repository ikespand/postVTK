#!/bin/bash
#PBS -l nodes=01:ppn=24
#PBS -N FindMinMax
#PBS -l walltime=00:05:00
#PBS -j oe

cd $PBS_O_WORKDIR

#source ${MODULESHOME}/init/bash
myenv=$(echo $PE_ENV | tr '[:upper:]' '[:lower:]')
module sw PrgEnv-${myenv} PrgEnv-gnu
module sw cray-mpich cray-mpich/7.7.0
module unload craype-hugepages16M
module load tools/paraview/5.5-git-master-parallel-Mesa

# process the TimeAvgEuu.py script in single process mode.
#aprun -n 1 pvbatch FindMinMax.py > FindMinMax_$PBS_JOBID.out 2>&1
#aprun -n 24 pvbatch RMSMeanRSS.py > RMSMeanRSS_$PBS_JOBID.out 2>&1

#for i in {41..128}
#do aprun -n 1 pvbatch 2TimeAvgEuu.py Channel_al_$i.vtu > 2TimeAvgEuu_$i.out 2>&1 &
#done
#wait
#exit 


aprun -n 1 pvbatch AveragingSpectrum.py > AveragingSpectrum_$PBS_JOBID.out 2>&1