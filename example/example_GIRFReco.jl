using Pkg

Pkg.activate(".")
Pkg.add("ZipFile")
Pkg.add("Downloads")
Pkg.add("ZipFile")
Pkg.add("Flux")
Pkg.add("HTTP")

using HTTP, Downloads, ZipFile, Flux

example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"
example_config_url = "https://github.com/BRAIN-TO/GIRFReco.jl/blob/example-workflow-updates-readme/example/recon_config_joss_demo.jl"
example_demo_script_url = "https://github.com/BRAIN-TO/GIRFReco.jl/blob/example-workflow-updates-readme/example/joss_demo.jl" 
example_cartesian_reco_script_url = "https://github.com/BRAIN-TO/GIRFReco.jl/blob/example-workflow-updates-readme/example/cartesian_recon.jl"
example_fieldmap_estimator_url = "https://github.com/BRAIN-TO/GIRFReco.jl/blob/example-workflow-updates-readme/example/fieldmap_estimator.jl"
example_env_url = "https://github.com/BRAIN-TO/GIRFReco.jl/blob/example-workflow-updates-readme/example/Project.toml"

HTTP.download(example_zenodo_data_url, "./data.zip") #Then uncompress

unzip("./data.zip", exdir="./data")

# Run the Demo
include("joss_demo.jl")
