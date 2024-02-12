using Documenter, SynthControl

makedocs(sitename="SynthControl",
    pages = ["Home" => "index.md", 
             "Introduction to Synthetic Controls" => "sc_intro.md",
             "Methods" => [
                "SimpleSCM" => "methods/SimpleSCM.md",
                "ADH2010" => "methods/ADH2010.md",
                "AugmentedSCM" => "methods/AugmentedSCM.md",
                "PenalizedSCM" => "methods/PenalizedSCM.md",
                "SyntheticDiD" => "methods/SyntheticDiD.md",
                "MC-NNM" => "methods/MC-NNM.md"
                ]
             ])

deploydocs(
    repo = "github.com/nilshg/SynthControl.jl.git",
)