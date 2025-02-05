include("./src/sim_trade_pattern_ek.jl")

Ncntry = 2
θ = 4
code = 1
σ = 1.5
τ = ones(Ncntry,Ncntry)
λ = [3.0; 3.0]

@time m, _ = sim_trade_pattern_ek(λ, τ, θ, σ, code);
@time m_fast, _ = sim_trade_pattern_ek_fast(λ, τ, θ, σ, code);
m
m_fast

τ2 = [1.0 1.5; 1.5 1.0]
@time m2, _ = sim_trade_pattern_ek(λ, τ2, θ, σ, code);
@time m2_fast, _ = sim_trade_pattern_ek_fast(λ, τ2, θ, σ, code);
m2
m2_fast

