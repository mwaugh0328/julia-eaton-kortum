include("./src/gravity-tools.jl")
include("./src/trade-environment.jl")
include("./src/simmulate-eaton-kortum.jl")
using CSV
using DataFrames
using Plots
using MINPACK

################################################################
# builds the EK dataset

dftrade, dfcntryfix, dflabor = make_ek_dataset()
# this one has the country numbers which allows for the construction of the 
# trade costs given the estimated fixed effects from the gravity regression

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

τ = zeros(Ncntry,Ncntry)

trd_prm = trade_params(θ = grv_params.θ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# @time πshares, foo = sim_trade_pattern_ek(exp.(grvdata.S), d, τ, grv_params.θ, 1.5);

πshares, foo = sim_trade_pattern_ek(trd_prm);

@time πshares_BEKK, foo = sim_trade_pattern_BEJK(exp.(grvdata.S), d, grv_params.θ, 1.5);

# πshares = average_trade_pattern(exp.(grvdata.S), d, grv_params.θ, 1.5, Nruns = 30);

# ################################################################
# # Recover the trade costs and technology parameters

dfmodel = plot_trade(πshares, Ncntry);

dfmodel_BEJK = plot_trade(πshares_BEKK, Ncntry);

plot(dfmodel.trade, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)