#!/bin/bash 

function print_usage {
    echo "Usage: $0 -p DataPath -c CartesianMIDs -s SpiralMIDs"
    echo "  -p DataPath             The directory of Cartesian and spiral measurement data, in double quotes. e.g. -p \"OneDrive - UHN/SPIDI003/dat\" "
    echo "  -c CartesianMIDs        The list of measurement IDs of Cartesian data, in double quotes and separated by space. e.g. -c \"82 85 86\" "
    echo "  -s SpiralMIDs           The list of measurement IDs of spiral data, in double quotes and separated by space. e.g. -s \"82 85 86\" "
    echo "All arguments should be in double quotes."
    exit 1
}

if [[ $# == 0 ]]; then
    print_usage
    exit 1
fi

while getopts p:c:s:h flag
do
    case "${flag}" in
        p)
            pathData="$OPTARG"
            ;;
        c)
            set -f # disable globbing temporarily
            IFS=' ' # split on space characters
            idsCartesian=($OPTARG) # use the split+glob operator
            ;;
        s)
            set -f # disable globbing temporarily
            IFS=' ' # split on space characters
            idsSpiral=($OPTARG)
            ;;
        h)
            print_usage
            ;;
        *)
            print_usage
            ;;
    esac
done

set +f # Turn globbing on

## standard conversion of .dat files
for id in ${idsCartesian[@]}; do
    printf -v filepattern "meas_MID%05d" "${id}"
    echo $filepattern 
    foundFiles=$(ls "${pathData}"/${filepattern}*.dat)
    echo "Converting $foundFiles to mrd using IsmrmrdParameterMap_Siemens_NX.xsl"
    siemens_to_ismrmrd -f "$foundFiles" -z 2 -x IsmrmrdParameterMap_Siemens_NX.xsl --skipSyncData
done

for id in ${idsSpiral[@]}; do
    printf -v filepattern "meas_MID%05d" "${id}"
    echo $filepattern 
    foundFiles=$(ls "${pathData}"/${filepattern}*.dat)
    echo "Converting $foundFiles to mrd using IsmrmrdParameterMap_Siemens_NX_ReadinSpiral.xsl"
    siemens_to_ismrmrd -f "$foundFiles" -z 2 -x IsmrmrdParameterMap_Siemens_NX_ReadinSpiral.xsl --skipSyncData
done