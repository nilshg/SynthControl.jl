
# Define plotting recipe
@recipe function f(s::SimpleSCM; kind = "overall")
    
    tp = s.treatment_panel
    
    legend --> :topleft

    if kind == "weights"

        wp = sort(DataFrame(comp = tp.is[Not(treated_ids(tp))], 
                            weight = s.w), 
                  :weight, rev = true)

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
            label --> treated_labels(tp)
            tp.ts, vec(tp.Y[treated_ids(tp), :])
        end

        @series begin
            label --> "Control"
            tp.ts, s.ŷ₁
        end

        @series begin
            label --> "Impact"
            seriestype --> :bar
            seriesalpha --> 0.5
            linecolor --> "white"
            tp.ts[only(length_T₀(tp)):end], s.τ̂
        end

        @series begin
            label --> "Intervention"
            seriestype := :vline
            linestyle := :dash
            seriescolor := "black"
            xguide := "Time"
            yguide := "Outcome"
            [tp.ts[only(length_T₀(tp))]]
        end

    elseif kind == "diffplot"
        @series begin
            label --> tp.treatment[1]
            sort!(unique(tp.data[!, tp.t_var])), [tp.y₁₀; tp.y₁₁] .- s.ŷ₁
        end

        @series begin
            label --> ""
            seriestype := :scatter
            seriescolor := 1
            markerstrokecolor := "white"
            sort!(unique(tp.data[!, tp.t_var])), [tp.y₁₀; tp.y₁₁] .- s.ŷ₁
        end

        @series begin
            label --> "Intervention"
            seriestype := :vline
            linestyle := :dash
            seriescolor := "black"
            xguide := "Time"
            yguide := "Outcome"
            [tp.treatment[2]]
        end
    end

     nothing
end
