using Ecophys
using Documenter

makedocs(;
    doctest = false,
    modules = [Ecophys],
    authors = "Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo = "https://github.com/VirtualPlantLab/Ecophys.jl/blob/{commit}{path}#{line}",
    sitename = "Ecophys.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "master",
        assets = String[]),
    pages = [
        "Photosynthesis" => "Photosynthesis.md",
        "Growth" => "Growth.md"
    ])

deploydocs(;
    repo = "github.com/VirtualPlantLab/Ecophys.jl.git",
    devbranch = "master")
