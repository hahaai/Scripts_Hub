#/bin/bash
# This will correct the auto defaced T1 images. 
# input: original T1 image, defaced image
# usuage: correct_AutoDeface.sh T1.nii.gz T1_defaced.nii.gz
# Aurthor: Lei Ai 06/2018


#goodhead='/data2/HBNcore/CMI_HBN_Data/MRI/RU/Release/R4_20171201_20180531/data/sub-5000820/anat/sub-5000820_acq-HCP_T1w.nii.gz'

#gooddeface='/data2/HBNcore/CMI_HBN_Data/MRI/RU/Release/R4_20171201_20180531/data/sub-5000820/anat/sub-5000820_acq-HCP_T1w_defaced.nii.gz'

goodhead=$1
gooddeface=$2

if [[ ! -f ${gooddeface/.nii.gz/_blocked.nii.gz} ]];then

##### 3dskullastirp to generate mask
mask=$(dirname $goodhead)'/brain_mask_'$RANDOM'.nii.gz'
3dSkullStrip -input $goodhead -mask_vol -prefix $mask

fslmaths $mask -kernel sphere 3 -dilM -bin $mask

# call python script to remove any remainign nose!!! 
# and also put the brain back. 

cd  /data2/HBNcore/CMI_HBN_Data/Scripts
ain=$goodhead'^^'$gooddeface'^^'$mask python - <<END
from deface_manually import correct_autodeface
import os
tmp=os.environ["ain"]
head=tmp.split('^^')[0]
anat=tmp.split('^^')[1]
mask=tmp.split('^^')[2]
correct_autodeface(head,anat,mask)
END

rm $mask

else
echo 'Already corrected'

fi




