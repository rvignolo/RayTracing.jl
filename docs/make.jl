pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using RayTracing
using Documenter

makedocs(;
    modules=[RayTracing],
    checkdocs=:exports,
    authors="Ramiro Vignolo <ramirovignolo@gmail.com>",
    repo="https://github.com/rvignolo/RayTracing.jl/blob/{commit}{path}#L{line}",
    sitename="RayTracing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        repolink="https://github.com/rvignolo/RayTracing.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
    ],
)
