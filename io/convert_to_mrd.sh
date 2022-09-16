#!/bin/bash 

idsCartesian=(82 83)
idsSpiral=(68 70 72 74 76 78 99 101 103 105 107 109)

pathData="/mnt/e/SPIDI/data/SPIDI_0007/Phantom/rawdata"
# standard conversion of .dat files
for id in ${idsCartesian[@]}; do
    printf -v filepattern "meas_MID%05d" "${id}"
    echo $filepattern 
    foundFiles=$(ls ${pathData}/${filepattern}*.dat)
    echo "Converting $foundFiles to mrd using IsmrmrdParameterMap_Siemens_NX.xsl"
    siemens_to_ismrmrd -f $foundFiles -z 2 -x IsmrmrdParameterMap_Siemens_NX.xsl
done

for id in ${idsSpiral[@]}; do
    printf -v filepattern "meas_MID%05d" "${id}"
    echo $filepattern 
    foundFiles=$(ls ${pathData}/${filepattern}*.dat)
    echo "Converting $foundFiles to mrd using IsmrmrdParameterMap_Siemens_NX_ReadinSpiral.xsl"
    siemens_to_ismrmrd -f $foundFiles -z 2 -x IsmrmrdParameterMap_Siemens_NX_ReadinSpiral.xsl
done