using bigleaf
using Documenter

DocMeta.setdocmeta!(bigleaf, :DocTestSetup, :(using bigleaf); recursive=true)

makedocs(;
    modules=[bigleaf],
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de> and contributors",
    repo="https://github.com/bgctw/bigleaf.jl/blob/{commit}{path}#{line}",
    sitename="bigleaf.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/bigleaf.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bgctw/bigleaf.jl",
    devbranch="main",
)
