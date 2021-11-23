
# Data loading for example
"""
    load_brexit()

Returns a `DataFrame` with quarterly GDP data on 30 OECD countries, borrowed from the
analysis of the effect of Brexit on UK GDP undertaken by the [Centre for European Reform](
https://www.cer.eu/insights/cost-brexit-june-2019)

The `TreatmentPanel` for Brexit should be specified as
```
TreatmentPanel("United Kingdom" => Date(2016, 7, 1), df;
                       outcome = :realgdp, id_var = :country, t_var = :quarter)
```
"""
function load_brexit()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"), DataFrame; stringtype = String)
end

function load_brexit_panel()
    df_brexit = load_brexit()
    BalancedPanel(df_brexit, "United Kingdom" => Date(2016, 7, 1); 
        id_var = :country, t_var = :quarter, outcome_var = :realgdp)
end

function load_germany()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","germany.csv"), DataFrame; stringtype = String)
end

function load_germany_panel()
    df_germany = load_germany()
    BalancedPanel(df_germany, "West Germany" => 1990; 
        id_var = "country", t_var = "year", outcome_var = "gdp")
end

function load_basque()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","basque.csv"), DataFrame; 
    stringtype = String, missingstring = "NA")
end

function load_basque_panel()
    df_basque = load_basque()
    BalancedPanel(df_basque, "Basque Country (Pais Vasco)" => 1970;
        id_var = "regionname", t_var = "year", outcome_var = "gdpcap")
end

function load_smoking()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","smoking.csv"), DataFrame; 
    stringtype = String, missingstring = "NA")
end
