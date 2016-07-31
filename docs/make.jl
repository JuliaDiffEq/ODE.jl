using Documenter, ODE

makedocs()

deploydocs(
           repo = "github.com/JuliaODE/ODE.jl.git",
           deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math")
           )
