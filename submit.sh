#!/bin/bash
#PBS -N Corr_MKL
#PBS -P ee
#PBS -m bea
#PBS -M ee3140506@iitd.ac.in
#PBS -l select=1:ncpus=12
#PBS -o stdout_file
#PBS -e stderr_file
#PBS -l walltime=04:00:00
#PBS -V 


echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
#cd $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
module load compiler/gcc/4.9.3/compilervars 
./fconn.o -i ../scratch/Analysis/funcAnalysis/Preprocessing_group0/ProcessedFile_sub0.nii.gz -o Corr_sub0 -r ../scratch/Analysis/atlases/fullbrainatlas_3mm/fullbrain_atlas_3mm.nii 274
#./fconn.o -i  -o Corr_sub2 -r fullbrain_atlas_thr0-2mm.nii.gz 274 -f -m MNI152_T1_2mm_brain_mask.nii.gz
#python /home/ee/btech/ee3140506/funcConn/funcConnGUI/Scripts/read_json.py funcAnalysisDesign.json
#bash prot.sh
#module load apps/Caffe
#time -p mpirun -n {1*2} bash /home/ee/btech/ee3140506/caffe/examples/mnist/train_lenet.sh
