include("./src/gravity-tools.jl")
include("./src/trade-environment.jl")
include("./src/simmulate-eaton-kortum.jl")
using CSV
using DataFrames
using Plots
using MINPACK

################################################################
# builds the EK dataset

dftrade = DataFrame(CSV.File("./ek-data/ek-data.csv"))

dflang = DataFrame(CSV.File("./ek-data/ek-language.csv"))

dflabor = DataFrame(CSV.File("./ek-data/ek-labor.csv"))

filter!(row -> ~(row.trade ≈ 1.0), dftrade);

filter!(row -> ~(row.trade ≈ 0.0), dftrade);

dftrade = hcat(dftrade, dflang);

#dfcntryfix = select(dftrade,Not("trade"))
dfcntryfix = DataFrame(CSV.File("./ek-data/ek-cntryfix.csv"));

################################################################

grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix)

grvdata = gravity(dftrade, display = true);

# ################################################################
# # Recover the trade costs and technology parameters

Ncntry = 19

d = zeros(Ncntry,Ncntry)
T = zeros(Ncntry)
W = ones(Ncntry)

make_trade_costs!(grvdata, d, grv_params)

make_technology!(grvdata, T, W, grv_params)

# ################################################################
# # Recover the trade costs and technology parameters


@time πshares, foo = sim_trade_pattern_ek(exp.(grvdata.S), d, grv_params.θ, 1.5);

@time πshares_BEKK, foo = sim_trade_pattern_BEJK(exp.(grvdata.S), d, grv_params.θ, 1.5);

# πshares = average_trade_pattern(exp.(grvdata.S), d, grv_params.θ, 1.5, Nruns = 30);

# ################################################################
# # Recover the trade costs and technology parameters

dfmodel = plot_trade(πshares, Ncntry);

dfmodel_BEJK = plot_trade(πshares_BEKK, Ncntry);

plot(dfmodel.trade, dfmodel_BEJK.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)