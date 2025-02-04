include("sim_trade_pattern_ek.jl")

Ncntry = 2
θ = 4
code = 1
σ = 1.5
τ = ones(Ncntry,Ncntry)
λ = [3.0; 3.0]

m, rec_low_price = sim_trade_pattern_ek(λ, τ, θ, σ, code)
m

τ2 = [1.0 1.5; 1.5 1.0]
m2, rec_low_price = sim_trade_pattern_ek(λ, τ2, θ, σ, code)
m2