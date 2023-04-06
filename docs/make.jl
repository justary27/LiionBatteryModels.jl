using Documenter
using CaseStudy

makedocs(
    sitename = "CaseStudy.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
    ),
    modules = [CaseStudy],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "manual/guide.md",
            "Examples" => "manual/examples.md",
        ],
        "Library" => [
            "Core" => "library/core.md",
            "Utils" => "library/utils.md",
        ],
        "License" => "license.md"
    ]
)

deploydocs(
    repo = "github.com/just-ary27/CaseStudy.git",
    versions = nothing,
)