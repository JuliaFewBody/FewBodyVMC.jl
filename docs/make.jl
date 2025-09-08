using FewBodyVMC
using Documenter

DocMeta.setdocmeta!(FewBodyVMC, :DocTestSetup, :(using FewBodyVMC); recursive=true)

makedocs(;
    modules=[FewBodyVMC],
    authors="Shuhei Ohno",
    sitename="FewBodyVMC.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaFewBody.github.io/FewBodyVMC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaFewBody/FewBodyVMC.jl",
    devbranch="main",
)
