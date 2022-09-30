using StirredReactor
using Documenter

DocMeta.setdocmeta!(StirredReactor, :DocTestSetup, :(using StirredReactor); recursive=true)

makedocs(;
    modules=[StirredReactor],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/StirredReactor.jl/blob/{commit}{path}#{line}",
    sitename="StirredReactor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/StirredReactor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/StirredReactor.jl",
    devbranch="main",
)
