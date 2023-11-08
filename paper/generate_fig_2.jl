#-----------------------------------------------------------------------------------
# # [Figure 2 Generation](@id generate_fig_2)
#-----------------------------------------------------------------------------------

#=
This page provides the ability to generate figure 2 as a result of the reconstructions with successive corrections
This page was generated from the following Julia file: [`generate_fig_2.jl`](@__REPO_ROOT_URL__/docs/lit/examples/generate_fig_2.jl)
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
root_project_path = "/Your/Project/Path/Here"
include("Your/Project/Path/docs/lit/examples/recon_config_joss_demo.jl")

#=

## 3. Need to use a set of tags to locate the files corresponding to the reconstructions with different corrections

=#

tags = ["vNoCorr_Sing", "vB0Corr_Sing", "vB0_GIRF_Sing", "vB0_GIRF_k0_Sing"]
paths = joinpath.(params_general[:results_path], "recon", tags)
filename_magnitude = splitext(params_general[:recon_save_filename])[1] * "_magn.nii"
files = joinpath.(paths, filename_magnitude)

#=

## 4. Load the NIfTI files

=#

uncorrected = niread(files[1]).raw[:, :, 8, 1, 1, 1]
b0_corrected = niread(files[2]).raw[:, :, 8, 1, 1, 1]
b0_girf_k1_corrected = niread(files[3]).raw[:, :, 8, 1, 1, 1]
b0_girf_k0_corrected = niread(files[4]).raw[:, :, 8, 1, 1, 1]
#=

## 5. Plot!

=#

gr()

p1 = heatmap(
    uncorrected,
    color = :grays,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "No Correction",
    titlefontsize = 10,
    top_margin = 0mm,
)
p2 = heatmap(
    b0_corrected,
    color = :grays,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "B₀ Correction",
    titlefontsize = 10,
    top_margin = 0mm,
)
p3 = heatmap(
    b0_girf_k1_corrected,
    color = :grays,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "B₀ + GIRF",
    titlefontsize = 10,
    top_margin = 0mm,
)
p4 = heatmap(
    b0_girf_k0_corrected,
    color = :grays,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "B₀ + GIRF + k₀",
    titlefontsize = 10,
    top_margin = 0mm,
)

pc1 = heatmap(
    abs.((b0_corrected - uncorrected) ./ maximum(b0_corrected)),
    color = :viridis,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "Δ: +B₀",
    titlefontsize = 10,
    top_margin = 0mm,
)
pc2 = heatmap(
    abs.((b0_girf_k1_corrected - b0_corrected) ./ maximum(b0_girf_corrected)),
    color = :viridis,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "Δ: +GIRF",
    titlefontsize = 10,
    top_margin = 0mm,
)
pc3 = heatmap(
    abs.((b0GirfK0k1_Corr - b0_girf_corrected) ./ maximum(b0_girf_k0_corrected)),
    color = :viridis,
    aspectratio = 1,
    xlims = (0, 200),
    ylims = (0, 200),
    xshowaxis = false,
    yshowaxis = false,
    colorbar = :none,
    title = "Δ: +k₀",
    titlefontsize = 10,
    top_margin = 0mm,
)

p_tot = plot(
    p1,
    p2,
    p3,
    p4,
    heatmap(
        ones(200, 200);
        xgrid = false,
        ygrid = false,
        aspectratio = 1,
        xlims = (0, 200),
        ylims = (0, 200),
        xshowaxis = false,
        yshowaxis = false,
        colorbar = :none,
        color = :grays,
    ),
    pc1,
    pc2,
    pc3;
    layout = (2, 4),
    aspectratio = 1,
    top_margin = 0mm,
    bottom_margin = 0mm,
)

savefig(p_tot, "./paper/fig2.pdf")
