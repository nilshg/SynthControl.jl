
# Define plotting recipe
@recipe function f(s::SynthControlModel; kind = "overall")
    legend --> :topleft

    if kind == "weights"

        wp = sort(DataFrame(comp = s.treatment_panel.comp_labels,
                  weight = s.w), :weight, rev = true)

        @series begin
            seriestype := :bar
            label --> ""
            xguide --> "Donor"
            yguide --> "Weight"
            linecolor := "white"
            wp[wp.weight .> 0.05, :comp], wp[wp.weight .> 0.05, :weight]
        end

    elseif kind == "overall"

        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁]
        end

        @series begin
            label --> "Control"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), s.ŷ₁
        end

        @series begin
            label --> "Impact"
            seriestype --> :bar
            seriesalpha --> 0.5
            linecolor --> "white"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))[findfirst(x -> x >= s.treatment_panel.treatment[2],
                                    sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))):end], s.τ̂
        end

        @series begin
            label --> "Intervention"
            seriestype := :vline
            linestyle := :dash
            seriescolor := "black"
            xguide := "Time"
            yguide := "Outcome"
            [s.treatment_panel.treatment[2]]
        end

    elseif kind == "diffplot"
        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁] .- s.ŷ₁
        end

        @series begin
            label --> ""
            seriestype := :scatter
            seriescolor := 1
            markerstrokecolor := "white"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁] .- s.ŷ₁
        end

        @series begin
            label --> "Intervention"
            seriestype := :vline
            linestyle := :dash
            seriescolor := "black"
            xguide := "Time"
            yguide := "Outcome"
            [s.treatment_panel.treatment[2]]
        end
    end

     nothing
end
