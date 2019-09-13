using Documenter, SDPAFamily

makedocs(;
    modules = [SDPAFamily],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Examples" => "examples.md",
        "Benchmarks" => "benchmarks.md",
        "Problematic problems & troubleshooting" => "problematicproblems.md",
        "Developer reference" => "reference.md"
    ],
    repo = "https://github.com/ericphanson/SDPAFamily.jl/blob/{commit}{path}#L{line}",
    sitename = "SDPAFamily.jl")


deploydocs(repo = "github.com/ericphanson/SDPAFamily.jl.git")
