using AcousticsToolbox
using Documenter

DocMeta.setdocmeta!(AcousticsToolbox, :DocTestSetup, :(using AcousticsToolbox); recursive=true)

makedocs(;
    modules=[AcousticsToolbox],
    authors="Mandar Chitre <mandar@nus.edu.sg> and contributors",
    repo="https://github.com/mchitre/AcousticsToolbox.jl/blob/{commit}{path}#{line}",
    sitename="AcousticsToolbox.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mchitre.github.io/AcousticsToolbox.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mchitre/AcousticsToolbox.jl",
    devbranch="main",
)
