push!(LOAD_PATH,"../src")
push!(LOAD_PATH,"../recon/")
push!(LOAD_PATH,"../utils/")

using Documenter, Literate, GIRFReco

# Generating Example Markdown files using Literate
lit = joinpath(@__DIR__, "lit", "examples")
gen = joinpath(@__DIR__, "src", "generated")

ipath = joinpath(lit, "joss_demo.jl")
opath = gen

Literate.markdown(ipath, opath; documenter = true)

makedocs(sitename="GIRFReco Documentation",
    modules = [GIRFReco],
    pages = [
        "Home" => "index.md",
        "Utilities" => "Utilities.md",
        "Examples" => joinpath("generated/", "joss_demo.md")
    ],
)

deploydocs(
    repo="github.com/BRAIN-TO/GIRFReco.git",
    push_preview = true,
    # deploy_config = Documenter.GitHubActions(),
    devbranch = "PaperPreparation",
    devurl = "dev",
    versions = ["stable" => "v^", "dev" => "dev"],
    forcepush = true,
)
