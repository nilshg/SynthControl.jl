using Documenter, SynthControl

makedocs(sitename="SynthControl",
    pages = ["Home" => "index.md", 
             "Methods" => [
                "SimpleSCM" => "methods/SimpleSCM.md"
             ]])

deploydocs(
    repo = "github.com/nilshg/SynthControl.jl.git",
)