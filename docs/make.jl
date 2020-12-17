using RayTracing
using Documenter

makedocs(;
    modules=[RayTracing],
    authors="Ramiro Vignolo <ramirovignolo@gmail.com>",
    repo="https://github.com/rvignolo/RayTracing.jl/blob/{commit}{path}#L{line}",
    sitename="RayTracing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
