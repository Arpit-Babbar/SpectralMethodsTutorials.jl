using SpectralMethodsTutorials
using Documenter

DocMeta.setdocmeta!(SpectralMethodsTutorials, :DocTestSetup, :(using SpectralMethodsTutorials); recursive=true)

makedocs(;
    modules=[SpectralMethodsTutorials],
    authors="Arpit Babbar",
    repo="https://github.com/arpit-babbar/SpectralMethodsTutorials.jl/blob/{commit}{path}#{line}",
    sitename="SpectralMethodsTutorials.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://arpit-babbar.github.io/SpectralMethodsTutorials.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/arpit-babbar/SpectralMethodsTutorials.jl",
    devbranch="main",
)
