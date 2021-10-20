using bigleaf
using Documenter, Latexify

# allow plot to work without display
# https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988/2
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(bigleaf, :DocTestSetup, :(using bigleaf, Latexify); recursive=true, warn=false)
doctest(bigleaf, manual = false)

makedocs(;
    modules=[bigleaf],
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de>. Jürgen Knauer <Juergen.Knauer@csiro.au> and contributors",
    repo="https://github.com/bgctw/bigleaf.jl/blob/{commit}{path}#{line}",
    sitename="bigleaf.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/bigleaf.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        #"Meteorological variables" => "metorological_variables.md",
        #"Unit conversions" => "unit_conversions.md",
        "Walkthrough" => "walkthrough.md",
        hide("metorological_variables.md"),
        hide("unit_conversions.md"),
        hide("bigleaf_constants.md"),
        "Index" => "autodocs.md",
        ],
)

deploydocs(;
    repo="github.com/bgctw/bigleaf.jl",
    devbranch="main",
)
