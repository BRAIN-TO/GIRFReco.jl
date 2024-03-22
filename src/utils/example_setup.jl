# Getting necessary files
using Downloads

function download_example()

    example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"
    example_config_url = ""
    example_demo_script_url = "" 
    example_cartesian_reco_script_url = ""
    example_fieldmap_estimator_url = ""
    example_env_url = ""

    Downloads.download(example_zenodo_data_url, "./data.zip") #Then uncompress
    Downloads.download(example_config_url, "./recon_config_joss_demo.jl")
    Downloads.download(example_demo_script_url,"./joss_demo.jl")
    Downloads.download(example_cartesian_reco_script_url,"./cartesian_recon.jl")
    Downloads.download(example_fieldmap_estimator_url,"./fieldmap_estimator.jl")
    Downloads.download(example_env_url, "./Project.toml")

end