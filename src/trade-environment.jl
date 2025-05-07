
using BenchmarkTools, SpecialFunctions
using Statistics
using Parameters

##########################################################################
##########################################################################
# some structures that I use

@with_kw struct trade_params
    θ::Float64 = 4.0
    σ::Float64 = 1.5
    d::Array{Float64} = [1.0  1.95; 1.95 1.0]
    N::Array{Float64} = [1.0, 1.0]
    S::Array{Float64} = [1.0, 1.0]
    T::Array{Float64} = ones(length(S))
    Wguess::Array{Float64} = ones(length(S))
    Ncntry::Int64 = length(S)
    τ::Array{Float64} = zeros(length(T), length(T))
end

##########################################################################
##########################################################################

function log_wages(xxx, Ncntry)
    # function to compute log wages
    # where xxx is the wage and tariff revene vector

    return vcat(log.(xxx[1:Ncntry - 1]), xxx[Ncntry:end])

end

function exp_wages(xxx, Ncntry)
    # function to compute exp wages
    # where xxx is the wage and tariff revene vector

    return vcat(exp.(xxx[1:Ncntry - 1]), xxx[Ncntry:end])

end

##########################################################################
##########################################################################

function eaton_kortum(W, d, T, θ)
    # constructs pattern of trade for eaton and kortum model

    Ncntry = size(d)[1]
    
    πshares = Array{eltype(d)}(undef, size(d))
    
    Φ = similar(T)

    for importer = 1:Ncntry

        for exporter = 1:Ncntry

            πshares[importer, exporter] = T[exporter] * (W[exporter] * d[importer, exporter])^(-θ)
            # equation (8)-like 

        end

        Φ[importer] = sum( πshares[importer, :])
        # equation (7)

        πshares[importer, : ] .= πshares[importer, : ] / Φ[importer]
        # complete equation (8)

    end

    return πshares, Φ

end

# ################################################################
# ################################################################

# function trade_balance(W, L, πshares; method = "solver")
#     # function to compute trade (im)balance
#     trade_balance = similar(L)
#     Ncntry = length(L)

#     for importer = 1:Ncntry

#         trade_balance[importer] = W[importer]*L[importer] - sum(πshares[:, importer] .* W .* L)
#         # equation (20) in Eaton and Kortum (no exogenous income, β = 1)
#         # need better explanation of this

#     end

#     if method == "solver"

#         return trade_balance

#     else

#         return trade_balance

#     end


# end

##########################################################################
##########################################################################

function trade_equilibrium(w, τrev, trade_parameters; display = false)
    # constructs zero function, takes wages, demand, tariff revenue
    # returns diffrence between procution and demand and guessed tariff transfer
    # and relized tariff transfer
    #
    # If output = "all" then returns all trade statistics

    @assert length(w) == length(τrev)

    @unpack d, τ, θ, σ, T, N = trade_parameters

    τ_zero = similar(w)


    # Step 1. Compute aggregate demand = labor income + tariff revenue
    AD = w.*N .+ τrev

    @assert length(w) ≈ length(AD)

    S = make_S(w, trade_parameters)

    #@assert exp.(S) ≈ trade_parameters.S

    πshares, τ_revenue, Φ = average_trade_pattern(exp.(S), d, τ, θ, σ)
        
    # Note too, that this will use the S, not the T....

    τ_revenue = sum(τ_revenue, dims = 2)[:].*AD
    # this is a place to  be carefull about...tariff revenue as poped out is 
    # on a per-unit of AD basis, so we need to multiply it by AD to get the
    # actual tariff revenue


    ##########################################################################
    ##########################################################################

    # Step 4. Compute trade balance
    trd_blnce = trade_balance(AD, πshares)

    # Step 5. Compute zero function for tarff revene, i.e. how does the 
    # guess of tariff revenue compare to the realized tariff revenue
    
    τ_zero .= τrev .- τ_revenue

    if display == false

        return vcat(τ_zero, trd_blnce)[1:end-1]
        # take the two tariff revenue conditions, 
        # then one of the market clearing conditions

    else

        return πshares, τ_revenue, Φ

    end

end

##########################################################################
##########################################################################

function trade_equilibrium_gravity(xxx, gravity_results, trade_parameters; display = false)
    # multiple dispacth function to compute trade equilibrium
    # takes a vector of wages and tariff revenue and trade parameters
    # returns the difference between production and demand and the guessed tariff transfer

    @unpack Ncntry = trade_parameters

    w = xxx[1:Ncntry - 1]

    push!(w, 1.0) # add the numeraire

    w = w ./ ( sum(w) / Ncntry)

    # println(w)

    T = make_technology(gravity_results, w, trade_parameters)

    # @assert T ≈ trade_parameters.T

    τrev = xxx[Ncntry:end]

    @assert length(w) == length(τrev)

    foo_trade_params = trade_params(trade_parameters, T = T)

    if display == false

        return trade_equilibrium(w, τrev, foo_trade_params, display = display)

    else

        return trade_equilibrium(w, τrev, foo_trade_params, display = display), T

    end

end

##########################################################################
##########################################################################

function trade_balance(AD, trade_share)
    # function to compute trade (im)balance
    trade_balance = similar(AD)
    Ncntry = length(AD)

    for importer = 1:Ncntry

        trade_balance[importer] = AD[importer]- sum(trade_share[:, importer] .* AD)
        # equation (20) in Eaton and Kortum (no exogenous income, β = 1)
        # need better explanation of this

    end

    return trade_balance

end

##########################################################################
##########################################################################

function trade_equilibrium(xxx, trade_parameters)
    # multiple dispacth function to compute trade equilibrium
    # takes a vector of wages and tariff revenue and trade parameters
    # returns the difference between production and demand and the guessed tariff transfer

    @unpack Ncntry = trade_parameters

    w = xxx[1:Ncntry - 1]

    push!(w, 1.0) # add the numeraire

    w = w ./ ( sum(w) / Ncntry)

    τrev = xxx[Ncntry:end]

    @assert length(w) == length(τrev)

    return trade_equilibrium(w, τrev, trade_parameters)

end

##########################################################################
##########################################################################

function trade_equilibrium(trade_parameters; display = true)
    # multiple dispatch function to compute trade equilibrium
    # just needs trade parameters and returns
    # the equilibrium wages, tariff revenue, and trade statistics

    f(x) = trade_equilibrium(exp_wages(x, trade_parameters.Ncntry), trade_parameters)

    function f!(fvec, x)

        fvec .= f(x)

    end

    wguess = trade_parameters.Wguess[1:trade_parameters.Ncntry-1]

    xguess = vcat(wguess, zeros(trade_parameters.Ncntry))

    xguess = log_wages(xguess, trade_parameters.Ncntry)

    n = length(xguess)
    diag_adjust = n - 1

    sol = fsolve(f!, xguess, show_trace = display, method = :hybr;
        ml=diag_adjust, mu=diag_adjust,
        diag=ones(n),
        mode= 1,
        tol=1e-3,)

    if sol.converged != true

        println("Convergence failed")

    end


    w = [exp.(sol.x[1:(trade_parameters.Ncntry - 1)]) ; 1.0]

    w = w ./ ( sum(w) / trade_parameters.Ncntry)

    τrev = sol.x[trade_parameters.Ncntry:end]

    # println(w)
    # println(τrev)

    πshares, τ_rev, Φ = trade_equilibrium(w, τrev, trade_parameters; display = true)

    # println(out.Qindex[19])

    return πshares, τ_rev, Φ, w
end