push!(LOAD_PATH,"../")

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
    deploy_config = Documenter.GitHubActions(),
)
