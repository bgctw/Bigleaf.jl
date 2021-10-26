using Bigleaf
using Documenter, Latexify

# allow plot to work without display
# https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988/2
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(Bigleaf, :DocTestSetup, :(using Bigleaf, Latexify); recursive=true, warn=false)
doctest(Bigleaf, manual = false)

makedocs(;
    modules=[Bigleaf],
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de>. Jürgen Knauer <Juergen.Knauer@csiro.au> and contributors",
    repo="https://github.com/bgctw/Bigleaf.jl/blob/{commit}{path}#{line}",
    sitename="Bigleaf.jl",
    doctestfilters=[r".*Info.*"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/Bigleaf.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        #"Meteorological variables" => "metorological_variables.md",
        #"Unit conversions" => "unit_conversions.md",
        "Walkthrough" => "walkthrough.md",
        hide("metorological_variables.md"),
        hide("evapotranspiration.md"),
        hide("global_radiation.md"),
        hide("unit_conversions.md"),
        hide("bigleaf_constants.md"),
        "Index" => "autodocs.md",
        ],
)

deploydocs(;
    repo="github.com/bgctw/Bigleaf.jl",
    devbranch="main",
)
