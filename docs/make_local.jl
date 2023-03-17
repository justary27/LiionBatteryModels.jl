using Documenter
using CaseStudy, CaseStudy.SenthilModel, CaseStudy.AshwiniModel

makedocs(
    sitename = "CaseStudy.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
    ),
    modules = [
        CaseStudy,
        CaseStudy.SenthilModel, 
        CaseStudy.AshwiniModel,
    ],
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