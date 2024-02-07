using Documenter, SynthControl

makedocs(sitename="SynthControl",
    pages = ["Home" => "index.md", 
             "Methods" => [
                "SimpleSCM" => "methods/SimpleSCM.md"
                "ADH2010" => "methods/ADH2010.md"
                ]
             ])

deploydocs(
    repo = "github.com/nilshg/SynthControl.jl.git",
)