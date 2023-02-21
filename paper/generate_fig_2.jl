#-----------------------------------------------------------------------------------
# # [Figure 3 Generation](@id generate_fig_3)
#-----------------------------------------------------------------------------------

#=
This page provides the ability to generate figure 3 as a result of the reconstructions with successive corrections
This page was generated from the following Julia file: [`generate_fig_3.jl`](@__REPO_ROOT_URL__/docs/lit/examples/generate_fig_3.jl)
=#

#=
## 1. Setup
The necessary Julia packages for reading NIfTI data
=#

#Base packages for file loading and plotting
using Measures
using NIfTI

#=
## 2. Where is my data?

We need to provide a location for the NIfTI files with each successive correction. Start by making sure we have the same file structure as in recon

=#
# rootProjPath = "/home/kasperl/SPIDI" # Root path of the project needs to be defined
rootProjPath = "/home/wuz/spiralDiffusion/data/demo_data"
# rootProjPath = "/Users/ajaffray/Documents/PhD/Data/SPIDI/"
include("ReconConfig_joss_demo.jl")

#=

## 3. Need to use a set of tags to locate the files corresponding to the reconstructions with different corrections

=#

tags = ["vNoCorr","vB0Corr","vB0_GIRF","vB0_GIRF_k0"]
paths = joinpath.(paramsGeneral[:pathResults], "recon",tags)
filename_magn = splitext(paramsGeneral[:fileNameSaveRecon])[1] * "_magn.nii"
files = joinpath.(paths,filename_magn)

#=

## 4. Load the NIfTI files

=#

noCorr = niread(files[1]).raw[:,:,5,1,1,1]
b0Corr = niread(files[2]).raw[:,:,5,1,1,1]
b0GirfCorr = niread(files[3]).raw[:,:,5,1,1,1]
b0GirfK0Corr = niread(files[4]).raw[:,:,5,1,1,1]
#=

## 5. Plot!

=#

gr()

p1 = heatmap(noCorr,color=:grays,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title="No Correction",titlefontsize=10,top_margin=0mm)
p2 = heatmap(b0Corr,color=:grays,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title="B₀ Correction",titlefontsize=10,top_margin=0mm)
p3 = heatmap(b0GirfCorr,color=:grays,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title="B₀ + GIRF",titlefontsize=10,top_margin=0mm)
p4 = heatmap(b0GirfK0Corr,color=:grays,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title = "B₀ + GIRF + k₀",titlefontsize=10,top_margin=0mm)

pc1 = heatmap(abs.((b0Corr - noCorr)./maximum(b0Corr)),color=:viridis,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title = "Δ: +B₀",titlefontsize=10,top_margin=0mm)
pc2 = heatmap(abs.((b0GirfCorr - b0Corr)./maximum(b0GirfCorr)),color=:viridis,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title = "Δ: +GIRF",titlefontsize=10,top_margin=0mm)
pc3 = heatmap(abs.((b0GirfK0Corr - b0GirfCorr)./maximum(b0GirfK0Corr)),color=:viridis,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,title = "Δ: +k₀",titlefontsize=10,top_margin=0mm)

p_tot = plot(p1,p2,p3,p4,heatmap(ones(200,200);xgrid=false,ygrid=false,aspectratio = 1,xlims=(0,200),ylims=(0,200),xshowaxis=false,yshowaxis=false,colorbar=:none,color=:grays),pc1,pc2,pc3;layout=(2,4),aspectratio=1,top_margin=0mm,bottom_margin=0mm)

savefig(p_tot,"./paper/figure/figure2.pdf")

# TODO keep working on the plots to generate something useable for the figure 3