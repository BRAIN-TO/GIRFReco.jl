using Dates

paramsGeneral = Dict{Symbol,Any}()

# UHN work
# paramsGeneral[:pathProject] = "C:\\Users\\Lars Kasper\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Laptop home, one drive sync
# paramsGeneral[:pathProject] = "C:\\Users\\kasperla\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Gadgetron Server
# laptop home, external drive
# paramsGeneral[:pathData] = "e:\\SPIDI\\data\\SPIDI_0007\\Phantom\\rawdata"

paramsGeneral[:pathProject] = "/home/kasperl/SPIDI"
paramsGeneral[:pathData] = joinpath(paramsGeneral[:pathProject], "data", "SPIDI_0007", "Human", "dat")
paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "SPIDI_0007", "Human")
paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon")


# paramsGeneral[:fileNameMultiEcho] = "meas_MID00083_FID06181_GRE_FieldMap_DualEcho_2mm.mrd"
paramsGeneral[:fileNameMultiEcho] = "field_map_132_2.h5"
paramsGeneral[:fullPathMultiEcho] = joinpath(paramsGeneral[:pathData], paramsGeneral[:fileNameMultiEcho])

paramsGeneral[:fileNameProcessedCartesian] = "processedCartesianData.h5"
paramsGeneral[:fullPathProcessedCartesian] = joinpath(paramsGeneral[:pathData], paramsGeneral[:fileNameProcessedCartesian])

paramsGeneral[:timeStamp] = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS")
paramsGeneral[:fileNameSaveRecon] = splitext(paramsGeneral[:fileNameMultiEcho])[1] * "_recon_" * paramsGeneral[:timeStamp] * ".h5"
paramsGeneral[:fileNameSaveNifti] = splitext(paramsGeneral[:fileNameMultiEcho])[1] * "_recon_" * paramsGeneral[:timeStamp] * ".nii"
paramsGeneral[:fileNameSaveSense] = splitext(paramsGeneral[:fileNameMultiEcho])[1] * "_sensemap_" * paramsGeneral[:timeStamp] * ".nii"
paramsGeneral[:fullPathSaveRecon] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveNifti] )
paramsGeneral[:fullPathSaveSense] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveSense] )
paramsGeneral[:doSaveRecon] = true