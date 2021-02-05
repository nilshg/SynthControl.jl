
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
    CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"), DataFrame)
end

function load_germany()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","germany.csv"), DataFrame)
end

function load_basque()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","basque.csv"), DataFrame)
end


function load_smoking()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","smoking.csv"), DataFrame)
end
