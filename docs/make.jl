using EOptInterface
using Documenter

DocMeta.setdocmeta!(EOptInterface, :DocTestSetup, :(using EOptInterface); recursive=true)

makedocs(;
    modules=[EOptInterface],
    authors="Joseph Choi <jsphchoi@mit.edu>",
    sitename="EOptInterface.jl",
    format=Documenter.HTML(;
        canonical="https://PSORLab.github.io/EOptInterface.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PSORLab/EOptInterface.jl",
    devbranch="master",
)
