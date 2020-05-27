# How to run CPAC using Singularity

# path to singularity file
simgfile='/data3/cnl/fmriprep/Lei_working/FCP-INDI-C-PAC-master-latest.simg'

## cpac output folder
output='/data3/cdb/Lei/working/CPAC1.6.2_Singularity/cpac_out/rest_test_CPACv162'
mkdir -p $output

# datafolder where the data in BIDS format is
data='/data3/cdb/Lei/working/CPAC1.6.2_Singularity/data'

# config folder where the pipeline config file is
config='/data3/cdb/Lei/working/CPAC1.6.2_Singularity/config/rest_test_CPACv1.6.2_nuis1-3a_test.yml'

## creat binding folder for inside the container, this should be in your home directory where you have write permission.
Binddir='/home/lai/cpac_test_output'
mkdir -p $Binddir'/data'
mkdir -p $Binddir'/output'
mkdir -p $Binddir'/config'

singularity run \
        -B $data:$Binddir'/data' \
        -B $output:$Binddir'/output' \
        -B $(dirname $config):$Binddir'/config' \
        $simgfile \
        $Binddir'/data' \
        $Binddir'/output' \
        participant \
        --pipeline_file $Binddir'/config/'$(basename $config) \
        --n_cpus 3 --pipeline_override "maximumMemoryPerParticipant: 15" \
        --pipeline_override "numParticipantsAtOnce: 5"




# go in to singularity
#singularity shell /data3/cnl/fmriprep/Lei_working/FCP-INDI-C-PAC-master-latest.simg 




### working one, simple:

#singularity run \
#        -B /data3/cdb/Lei/working/CPAC1.6.2_Singularity:/home/lai/cpac_test_output \
#        /data3/cnl/fmriprep/Lei_working/FCP-INDI-C-PAC-master-latest.simg \
#        /home/lai/cpac_test_output/data \
#        /home/lai/cpac_test_output/rest_test_CPACv162 \
#        participant \
#        --pipeline_file /home/lai/cpac_test_output/rest_test_CPACv1.6.2_nuis1-3a_test.yml \
#        --n_cpus 3 --pipeline_override "maximumMemoryPerParticipant: 15" \
#        --pipeline_override "numParticipantsAtOnce: 5"



