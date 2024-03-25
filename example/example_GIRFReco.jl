using Pkg

Pkg.activate(".")
Pkg.add("ZipFile")
Pkg.add("Downloads")
Pkg.add("ZipFile")
Pkg.add("Flux")
Pkg.add("HTTP")

using HTTP, Downloads, ZipFile, Flux

example_zenodo_data_url = "https://zenodo.org/records/7779045/files/data_girfreco_03_28_2023.zip?download=1"

if ~isfile("./data.zip")
    HTTP.download(example_zenodo_data_url, "./data.zip") #Then uncompress
end

if ~isdir("./data") # need to download
    Base.prompt("Please navigate to $(pwd()) and unzip data.zip")
else # Run the Demo
    include("joss_demo.jl")
end