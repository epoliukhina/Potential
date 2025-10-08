using Documenter, Potential

makedocs(;
sitename="Potential",
authors="Ekaterina Poliukhina",
clean=true,
modules=[Potential],
checkdocs=:exports,
pages = [
    "Home" => "index.md",
    "Theory" => "theory.md"
    ]
)

deploydocs(;
    repo="github.com/epoliukhina/Potential.jl",
    devbranch="main",
)
