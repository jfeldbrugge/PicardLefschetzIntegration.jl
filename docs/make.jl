using PicardLefschetzIntegration
using Documenter

DocMeta.setdocmeta!(PicardLefschetzIntegration, :DocTestSetup, :(using PicardLefschetzIntegration); recursive=true)

makedocs(;
    modules=[PicardLefschetzIntegration],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="PicardLefschetzIntegration.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/PicardLefschetzIntegration.jl",
    devbranch="main",
)
