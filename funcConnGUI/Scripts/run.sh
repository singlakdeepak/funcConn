START=3
END=12

for i in `seq $START  $END`;
do

echo $i 
#bash XGBoostSubmission.sh  $i  &
fslmaths file$i -mul MNI152_T1_2mm_brain_mask.nii.gz file$i
done

