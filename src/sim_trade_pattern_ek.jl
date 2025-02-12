using Random

function sim_trade_pattern_ek(λ, τ, θ, σ, code)

    # Parameters for goods and countries and the sample size for prices
    Ngoods = 1000000 # Adjust if too slow
    Ncntry = length(λ)

    # Parameters for technologies
    η = σ
    inv_η = 1 / (1 - η)
    inv_Ngoods = 1 / Ngoods
    low_price = 1.0e7

    # Draw productivities and and compute unit costs to produce each good in
    # each country
    Random.seed!(03281978 + code)
    pconst = Array{eltype(λ)}(undef, Ngoods, Ncntry)
 

    Threads.@threads for j in 1:Ncntry

        u = rand(MersenneTwister(03281978 + j), Ngoods)

        pconst[:, j] .= (log.(u) ./ (-λ[j])) .^ (1/θ) 

    end

    # Loop to calculate the low price and country suppliers
    m = Array{eltype(λ)}(undef,Ncntry, Ncntry)
    sum_price = Array{eltype(λ)}(undef,Ncntry)
    rec_low_price = Array{eltype(λ)}(undef,Ngoods, Ncntry)

    for gd in 1:Ngoods # This is the good
        
        for im in 1:Ncntry # Importing country

            cif_price = Array{eltype(λ)}(undef, Ncntry) 

            for ex in 1:Ncntry # This is the country (potentially) exporting the good

                cif_price[ex] = τ[ex, im] * pconst[gd, ex]

                if ex != 1
                    low_price = min(cif_price[ex], low_price)
                else
                    low_price = cif_price[ex]
                end

            end

            if low_price == cif_price[im]

                m[im, im] += low_price^(1 - η)

            else

                for ex in 1:Ncntry

                    if low_price == cif_price[ex]

                        m[ex, im] += low_price^(1 - η)
                        break

                    end 

                end

            end

            # Now record the low price
            sum_price[im] += low_price^(1 - η)
            rec_low_price[gd, im] = low_price 
        end

    end

    # Loop to calculate aggregate price index and the trade shares.
    g = zeros(Ncntry)

    for im in 1:Ncntry

        g[im] = (sum_price[im] * inv_Ngoods)^inv_η

        for ex in 1:Ncntry

            m[ex, im] = (inv_Ngoods*m[ex, im]) / g[im]^(1 - η)

        end

    end

    return m, rec_low_price

end

###############################################################
###############################################################

function average_trade_pattern(λ, τ, θ, σ; Ngoods = 100000, Nruns = 30)

    πshares = zeros(size(τ))    

    Threads.@threads for xxx = 1:Nruns

        πshares = πshares .+ (1.0 / Nruns)*sim_trade_pattern_ek_fast(λ, τ, θ, σ, Ngoods = Ngoods, code = xxx)[1]

    end

    return πshares

end

###############################################################
###############################################################


function sim_trade_pattern_ek_fast(λ, τ, θ, σ; Ngoods = 100000, code = 1)
    # Parameters for goods and countries and the sample size for prices
    Ngoods = Ngoods # Adjust if too slow
    Ncntry = length(λ)

    # Parameters for technologies
    
    inv_Ngoods = 1.0 / Ngoods

    ###############################################################
    # Draw productivities and and compute unit costs to produce each good in
    # each country
    
    p = Array{Float64}(undef, Ncntry, Ngoods)

    u = Array{Float64}(undef, Ncntry, Ngoods)

    rand!(MersenneTwister(03281978 + code ), u)

    println(code)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p[j, :] .= marginal_cost.(u[j,:], λ[j], θ) 
        #not sure that this helped... may need to return back to the original code

    end

    ###############################################################

    # Loop to calculate the low price and country suppliers
    m = zeros(Ncntry, Ncntry) # need to be zero as I'm summing over stuff
    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods  # Loop over goods # threading here messes stuff up  

        @inbounds for im in 1:Ncntry  # Loop over importing countries

            low_price = p[im, gd]
            min_ex = im

            @inbounds for ex in 1:Ncntry

                cif_price = τ[im, ex] * p[ex, gd] # price of exporter

                if cif_price < low_price # if the price is lower than the current low price

                    low_price = cif_price # it is the low price
                    
                    min_ex = ex # and the exporter is the one with the lowest price
                end

            end

            # Update trade matrix `m`
            m[im, min_ex] += low_price^(1.0 - σ)  

            # Update sum price and record lowest price
            sum_price[im] += low_price^(1.0 - σ) 

            rec_low_price[im, gd] = low_price
        end

    end

    # Loop to calculate aggregate price index and the trade shares.

    for im in 1:Ncntry

        g_val = (sum_price[im] * inv_Ngoods)

        for ex in 1:Ncntry

            m[im, ex] = inv_Ngoods*( m[im, ex] ) / g_val

        end

    end

    return m, rec_low_price

end

###############################################################
###############################################################

function marginal_cost(u, λ, θ)
    # takes random number u, productivity λ and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-λ) )^ ( one(θ) / θ )

    # (log.(u) ./ (-λ[j])) .^ (1/θ) 

end

###############################################################

function plot_trade(πshares, Ncntry)

    trademodel = log.(vec(normalize_by_home_trade(πshares, Ncntry)'))

    dfmodel = DataFrame(trade = trademodel)

    filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

    filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

    return dfmodel

end

