#PBS -l nodes=1:ppn=16
#PBS -l walltime=600:00:00

#! change the working directory (default is home directory)
cd $PBS_O_WORKDIR
#! Record some useful job details in the output file
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following nodes:
echo `cat $PBS_NODEFILE | uniq`
#! add the MATLAB module (as per BCp3)
module add apps/matlab-r2016b
module add cuda/toolkit/7.5.18
options="-nodesktop -nosplash -noFigureWindows"
#! Run MATLAB in batch mode
matlab $options -r runOptimisationIgnotus
