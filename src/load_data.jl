
# Data loading for example
"""
    load_brexit()

Returns a `DataFrame` with quarterly GDP data on 30 OECD countries, borrowed from the analysis of
the effect of Brexit on UK GDP undertaken by the [Centre for European Reform](
https://www.cer.eu/insights/cost-brexit-june-2019)

With the Brexit referendum having taken place on 23 June 2016, the `TreatmentPanel` for Brexit
should be specified as
    ```
    TreatmentPanel(df, "United Kingdom" => Date(2016, 7, 1);
        outcome = :realgdp, id_var = :country, t_var = :quarter)
    ```
"""
function load_brexit()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"), DataFrame; stringtype = String)
end

"""
    load_brexit_panel()

Returns a `BalancedPanel` object based on the data set available through load_brexit(), with
treatment assignment specified as `"United Kingdom" => Date(2016, 7, 1)`

"""
function load_brexit_panel()
    df_brexit = load_brexit()
    BalancedPanel(df_brexit, "United Kingdom" => Date(2016, 7, 1); 
        id_var = :country, t_var = :quarter, outcome_var = :realgdp)
end


"""
    load_germany()

    Returns a `DataFrame` with annual observations of per-capita GDP and covariates between 1960 and
    2003. The data is used in Abadie, Diamond and Hainmueller (2015) to study the effect of German
    unification on per-capita GDP. 

    As German reunification occurred in 1990, the `TreatmentPanel` for use in synthetic control
    estimation should be specified as:

    ```
    TreatmentPanel(df, "West Germany" => 1990,
        id_var = "country", t_var = "year", outcome_var = "gdp")
    ```
"""
function load_germany()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","germany.csv"), DataFrame; stringtype = String)
end

"""
    load_germany_panel()

    Returns a `BalancedPanel` object based on the data available from `load_germany()`, with
    treatment assignment specified as `"West Germany" => 1990`
"""
function load_germany_panel()
    df_germany = load_germany()
    BalancedPanel(df_germany, "West Germany" => 1990; 
        id_var = "country", t_var = "year", outcome_var = "gdp")
end

"""
    load_basque()

    Returns a `DataFrame` with annual observations of per-capita GDP and covariates between 1955 and
    1997. The data is used in Abadie and Gardeazabal (2003) to study the effect of ETA terrorism  
          on per-capita GDP in the Basque country. 

    To undertake analysis simliar to that in Abadie and Gardeazabal (2003), the `TreatmentPanel` for
    use in synthetic control estimation should be specified as:

    ```
    TreatmentPanel(df, "Basque Country (Pais Vasco)" => 1970,
        id_var = "regionname", t_var = "year", outcome_var = "gdpcap")
    ```
"""
function load_basque()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","basque.csv"), DataFrame; 
    stringtype = String, missingstring = "NA")
end

"""
    load_basque_panel()

    Returns a `BalancedPanel` object based on the data available from `load_basque()`, with
    treatment assignment specified as `"Basque Country (Pais Vasco)" => 1970`
"""
function load_basque_panel()
    df_basque = load_basque()
    BalancedPanel(df_basque, "Basque Country (Pais Vasco)" => 1970;
        id_var = "regionname", t_var = "year", outcome_var = "gdpcap")
end

"""
    load_smoking()

    Returns a `DataFrame` with annual observations of cigarette sales and covariates for 39 US
    states between 1970 and 2000. The data is used in Abadie, Diamond and Hainmueller (2010) to
    study the effect of Proposition 99 on cigarette sales in California

    To undertake analysis simliar to that in Abadie, Diamond and Hainmueller (2010), the
    `TreatmentPanel` for use in synthetic control estimation should be specified as:

    ```
    TreatmentPanel(df, 3 => 1989,
        id_var = "state", t_var = "year", outcome_var = "cigsale")
    ```
"""
function load_smoking()
    df_smoking = CSV.read(joinpath(dirname(@__FILE__),"..","data","smoking.csv"), DataFrame; 
    stringtype = String, missingstring = "NA")
    select!(df_smoking, Not(:Column1))
end

"""
    load_smoking_panel()

    Returns a `BalancedPanel` object based on the data available from `load_smoking()`, with
    treatment assignment specified as `3 => 1989`
"""
function load_smoking_panel()
    df_smoking = load_smoking()
    BalancedPanel(df_smoking, 3 => 1989;
        id_var = :state, t_var = :year, outcome_var = :cigsale)
end