using Documenter
using LiionBatteryModels

makedocs(
    sitename="LiionBatteryModels.jl",
    format=Documenter.HTML(
        assets=["assets/favicon.ico"]
    ),
    modules=[LiionBatteryModels],
    pages=[
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
    repo="github.com/just-ary27/LiionBatteryModels.jl.git",
    versions=nothing,
)
