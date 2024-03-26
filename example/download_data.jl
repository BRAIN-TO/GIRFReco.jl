using Pkg
Pkg.add("HTTP")
using HTTP

example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"

# Zip data file download
if ~isfile("./data.zip")
    HTTP.download(example_zenodo_data_url, "./data.zip")
end

# Unzipping data files
if ~isdir("./joss_data_zenodo")
    mkpath("./joss_data_zenodo")
    run(`unzip data.zip -d ./joss_data_zenodo/`)
end