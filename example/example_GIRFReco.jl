using Pkg

Pkg.activate(".")
Pkg.add("Flux")
Pkg.add("HTTP")

using HTTP, Downloads, ZipFile, Flux

example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"

if ~isfile("./data.zip")
    HTTP.download(example_zenodo_data_url, "./data.zip") #Then uncompress
end

if ~isdir("./joss_data_zenodo") # need to download
    Base.prompt("Please navigate to $(pwd()) and unzip data.zip using your favourite unzip tool. Once complete, press any key to proceed")
end

if isdir("./joss_data_zenodo") # need to download # Run the Demo
    include("joss_demo.jl")
end