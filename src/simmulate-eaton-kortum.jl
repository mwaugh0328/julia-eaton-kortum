using Random

###############################################################
###############################################################

function average_trade_pattern(S, d, θ, σ; Ngoods = 100000, Nruns = 30)
    # Computes the trade shares when avergaged over diffrent simmulaitons
    # of the trade pattern. The function returns the average trade shares.

    πshares = zeros(size(d))    

    Threads.@threads for xxx = 1:Nruns

        πshares = πshares .+ (1.0 / Nruns)*sim_trade_pattern_ek(S, d, θ, σ, Ngoods = Ngoods, code = xxx)[1]

    end

    return πshares

end

###############################################################
###############################################################

function sim_trade_pattern_ek(S, d,  θ, σ; Ngoods = 100000, code = 1)
    # Constructs pattern of trade for the perfectly competitive model with Frechet
    # distributed productivity shocks. The function returns the trade shares and the
    # lowest price of each good in each country.
    #
    # S are parameters from gravity equation and are sufficient to simmulate marginal costs
    # d is the trade costs matrix with rows being imports, columns being exports
    # θ is the Frechet shape parameter
    # σ is the elasticity of substitution
    # options include number of goods and a code for the random number generator
    
    Ncntry = length(S)

    inv_Ngoods = 1.0 / Ngoods

    ###############################################################
    # Draw productivities and and compute unit costs to produce each good in
    # each country
    
    p = Array{Float64}(undef, Ncntry, Ngoods)

    u = Array{Float64}(undef, Ncntry, Ngoods)

    rand!(MersenneTwister(03281978 + code ), u)

    println(code)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
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

                cif_price = d[im, ex] * p[ex, gd] # price of exporter

                # low_price, min_ex = ifelse(cif_price < low_price, (cif_price, ex), (low_price, min_ex)) 

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

function marginal_cost(u, S, θ)
    # takes random number u, productivity S and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-S) )^ ( one(θ) / θ )

    # (log.(u) ./ (-λ[j])) .^ (1/θ) 

end

###############################################################

function plot_trade(πshares, Ncntry)
    # helps for plotting the trade shares

    trademodel = log.(vec(normalize_by_home_trade(πshares, Ncntry)'))

    dfmodel = DataFrame(trade = trademodel)

    filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

    filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

    return dfmodel

end

