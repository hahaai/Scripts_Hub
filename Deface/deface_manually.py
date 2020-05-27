# this removes any face (particular nose) remaining from the auto deface.
def correct_autodeface(f_head,f_anat,f_mask):
    import nibabel as nb
    import numpy as np

    #f_mask='/Users/lei.ai/Downloads/deface/anat_new/brain_mask.nii.gz'
    #f_anat='/Users/lei.ai/Downloads/deface/anat_new/sub-5109687_acq-HCP_T1w_defaced.nii.gz'

    img_mask=nb.nifti1.load(f_mask)
    img=nb.nifti1.load(f_anat)
    img_head=nb.nifti1.load(f_head)
    
    # for anatomical image. need it later
    img_header=img.get_header()
    img_aff=img.get_affine()

    img_shape=img.shape
    img_voxels=img.get_header().get_zooms()

    img_data_type=img.get_data_dtype()
    img_data=img.get_data()
    img_all=np.array(img_data)

    # mask:
    img_mask_data=img_mask.get_data()
    img_mask_all=np.array(img_mask_data)


    # anat head,then brain:
    img_head_data=img_head.get_data()
    img_head_all=np.array(img_head_data)
    img_brain_all=img_head_all*img_mask_all


    # add any brain cut back into the defaced one.
    temp1=img_all==0
    temp2=img_mask_all==1
    temp=temp1 & temp2

    img_all[temp]=img_brain_all[temp]

    # loop through each sagital slice
    for i in range(img_mask_all.shape[0]):
        img_2d=img_mask_all[i,:,:]

        # img_2d=img_mask_all.mean(axis=0) # average axially, to get a 2d image. Not using mean, actuall do this for each sagital slice

    # find the most right-index that are non-zeor
        xx=np.where(img_2d>0)
        if len(xx[0])>0:
            area=len(xx[0])

            tmp_idx=np.where(xx[0]==max(xx[0]))
            AP_idx=xx[0][tmp_idx[0][np.argmin(xx[1][tmp_idx])]]
            IS_idx=xx[1][tmp_idx[0][np.argmin(xx[1][tmp_idx])]]


            IS_idx_final=IS_idx-20*area/27132 # move the postion down (S --> I) 27132 is a arbitory number!!!!!

            AP_idx_final=AP_idx-7*area/27132 # move the positon left (A --> P), but make sure not cutting in teh brain.
            while img_2d[AP_idx_final,IS_idx_final]>0:
                AP_idx_final += 1
            # Update the anat image based the above coordinates.
            img_all[i,AP_idx_final:,:IS_idx_final]=0



    op_nifti=nb.nifti1.Nifti1Image(img_all,img_aff,header=img_header)    

    nb.nifti1.save(op_nifti, f_anat.replace('.nii.gz','_blocked.nii.gz'))


    #coords=':,' + str(AP_idx_final) + ':,:' + str(IS_idx_final)

    #newimg=block_deface_new(f_anat,coords)

