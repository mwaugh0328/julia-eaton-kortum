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
    pconst = zeros(Ngoods, Ncntry)

    for j in 1:Ncntry

        u = rand(MersenneTwister(03281978 + j), Ngoods)

        # should be p = (log(u) / (-T_i))^(-1/θ) ?
        pconst[:, j] .= (log.(u) ./ (-λ[j])) .^ (1/θ) 

    end

    # Loop to calculate the low price and country suppliers
    m = zeros(Ncntry, Ncntry)
    sum_price = zeros(Ncntry, 1)
    rec_low_price = zeros(Ngoods, Ncntry)

    for gd in 1:Ngoods # This is the good
        
        for im in 1:Ncntry # Importing country

            cif_price = zeros(Ncntry) 

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