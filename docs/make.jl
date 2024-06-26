push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "../src/io/")
push!(LOAD_PATH, "../src/utils/")

using Documenter, Literate, GIRFReco

# Generating Example Markdown files using Literate
lit = joinpath(@__DIR__, "lit", "examples")
gen = joinpath(@__DIR__, "src", "generated")

ipath_script = joinpath(lit, "joss_demo.jl")
ipath_config = joinpath(lit, "recon_config_joss_demo.jl")
opath = gen
Literate.markdown(ipath_script, opath; documenter = true)
Literate.markdown(ipath_config, opath; documenter = true)

makedocs(
    sitename = "GIRFReco.jl Documentation",
    modules = [GIRFReco],
    pages = [
        "Home" => "index.md",
        "API" => "Utilities.md",
        "Examples" => [joinpath("generated/", "joss_demo.md"), joinpath("generated/", "recon_config_joss_demo.md")],
    ],
)

deploydocs(
    repo = "github.com/BRAIN-TO/GIRFReco.jl.git",
    push_preview = true,
    # deploy_config = Documenter.GitHubActions(),
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "dev" => "dev"],
    forcepush = true,
)
