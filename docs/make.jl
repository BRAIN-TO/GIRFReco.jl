push!(LOAD_PATH,"../src")
push!(LOAD_PATH,"../recon/")
push!(LOAD_PATH,"../utils/")

using Documenter, GIRFReco

makedocs(sitename="GIRFReco Documentation",
    modules = [GIRFReco],
    pages = [
        "Home" => "index.md",
        "Utilities" => "Utilities.md"
    ]
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
