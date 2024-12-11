using TimeVaryingPeriodograms
using Documenter

DocMeta.setdocmeta!(TimeVaryingPeriodograms, :DocTestSetup, :(using TimeVaryingPeriodograms); recursive=true)

makedocs(;
    modules=[TimeVaryingPeriodograms],
    authors="OskarGustafsson",
    sitename="TimeVaryingPeriodograms.jl",
    format=Documenter.HTML(;
        canonical="https://OskarGU.github.io/TimeVaryingPeriodograms.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/OskarGU/TimeVaryingPeriodograms.jl",
    devbranch="master",
)
 