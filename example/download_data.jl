# using Pkg

# Pkg.activate(".")
# Pkg.add("HTTP")
# Pkg.instantiate()
# using HTTP

example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"

# Zip data file download
if ~isfile("./data.zip") && ~isfile("./joss_data_zenodo/GIRF/GIRF_ISMRMR2022/2021Nov_PosNeg_Gx") && ~isfile("./joss_data_zenodo/Gradients/gradients508.txt")
    @info "Downloading demo data from Zenodo."
    HTTP.download(example_zenodo_data_url, "./data.zip")
end

# Unzipping data files
if ~isdir("./joss_data_zenodo")
    @info "Unzipping demo data."
    run(`unzip data.zip -d ./`)
end