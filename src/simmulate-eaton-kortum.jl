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

function average_trade_pattern(S, d, τ, θ, σ; Ngoods = 50000, Nruns = 5)
    # Computes the trade shares when avergaged over diffrent simmulaitons
    # of the trade pattern. The function returns the average trade shares.

    πshares = Array{Float64}(undef, length(S), length(S), Nruns)
    τ_revenue = Array{Float64}(undef, length(S), length(S), Nruns)
    Φ = Array{Float64}(undef, length(S), Nruns)

    avg_πshares = zeros(size(d))
    avg_τ_revenue = zeros(size(S))
    avg_Φ = zeros(size(S))

    Threads.@threads for xxx = 1:Nruns

       πshares[:,:, xxx], Φ[:,xxx], τ_revenue[:, :, xxx] = sim_trade_pattern_ek(S, d, τ, θ, σ, Ngoods = Ngoods, code = xxx)

    end

    for xxx = 1:Nruns

        avg_πshares = avg_πshares .+ (1.0 / Nruns)*πshares[:,:, xxx]

        avg_τ_revenue = avg_τ_revenue .+ (1.0 / Nruns).*τ_revenue[:,:,xxx]

        avg_Φ = avg_Φ .+ (1.0 / Nruns).*Φ[:,xxx]

    end

    return avg_πshares, avg_τ_revenue, avg_Φ

end

###############################################################
###############################################################

function sim_trade_pattern_ek(trade_parameters; Ngoods = 100000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_ek function
    # this allows me to pass the trade_parameters structure and it will work

    return sim_trade_pattern_ek(trade_parameters.S, trade_parameters.d, trade_parameters.τ, 
    trade_parameters.θ, trade_parameters.σ; Ngoods = Ngoods, code = code)

end


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

    #println(code)

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

            ###############################################################
            # This is an alternative way to find the low cost exporter

            # cif_price = d[im, :] .* p[:, gd]

            # sorted_price = sort(cif_price)

            # low_price = sorted_price[1]

            # min_ex = findfirst(==(low_price), cif_price)

            # # This ==(low_price) creates an anonymous function that checks if an element is equal to low_price.

            # ###############################################################

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

function marginal_cost_second(u, first_price, S, θ)
    # takes random number u, productivity S and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-S) + (first_price)^ (θ) )^ (one(θ) / θ)

    # (log(u2)./(-S(j)) + (p1const(:,j).^-1).^(-theta)).^(-1./theta);

    # Second draw, second best productivity, this comes from
    # 1-exp(-S*z_two^(-θ) + S*z_one^(-θ)) 

end

###############################################################
###############################################################

function sim_trade_pattern_BEJK(S, d, θ, σ; Ngoods = 100000, code = 1)
    # A function to simmulate a pattern of trade and then generate a trade
    # share matrix and a random sample of final goods prices given bertrand
    # pricing from BEJK (2003). 

    Ncntry = length(S)

    inv_Ngoods = 1 / Ngoods

    markup = σ / (σ - one(σ))

    ###########################################################
    # Draw productivities and compute unit costs to produce each good in 
    # each country

    p1const = Array{Float64}(undef, Ncntry, Ngoods)

    p2const = Array{Float64}(undef, Ncntry, Ngoods)

    u = Array{Float64}(undef, Ncntry, Ngoods)

    rand!(MersenneTwister(03281978 + code ), u)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p1const[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
        # Invert to convert to unit costs. Here assume w[j] = 1 ?

    end

    rand!(MersenneTwister(02071983 + code ), u)
    
    @inbounds @views Threads.@threads for j in 1:Ncntry

        p2const[j, :] .= marginal_cost_second.(u[j,:], p1const[j, :], S[j], θ)

    end

    ###########################################################
    # Loop to calculate the low price and country suppliers

    m = zeros(Ncntry, Ncntry)

    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    # rec_cost = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods # This is the good

        @inbounds for im in 1:Ncntry # This is the country importing the good

            # In BEJK there are two seperate issues:
            # 1. Needs to find who is the low cost producer AND second lowest cost producer.
            # 2. Need to find the price charged by the low cost producer which
            # is either the 2nd domestic low cost producer, the 2nd foreign low cost producer or the monopolist price.

            # this is how the old matlab code worked, the first julia BEJK version put together earlier was not correct (unsure why)
            # The first step is to find the two lowest international cost producers
        
            cif_price1 = d[im, 1] * p1const[1,gd]
            
            cif_price2 = d[im, 2] * p1const[2,gd]
           
            if cif_price1 < cif_price2 # if the first country is lower, its low cost

               low_cost = cif_price1

               low2_cost = cif_price2 # second country is second lowest
                
               min_ex = 1
            else # if the second country is lower, its low cost 1 is second lowest
               
                low_cost = cif_price2
               
                low2_cost = cif_price1 # first country is second lowest
               
                min_ex = 2
            end


            for ex in 3:Ncntry
                # now walk through remaining potential exporters

                cif_price = d[im, ex] * p1const[ex, gd]
                # this is the price of the exporter

                low2_cost = min(max(cif_price, low_cost), low2_cost)

                if cif_price < low_cost # if the exporter price is lower than the current low price

                    low_cost = cif_price # the low costs is that exporters price
                    
                    min_ex = ex # and the exporter is the one with the lowest price

                end

            end

            price_charged = min(min(d[im, min_ex] * p2const[min_ex, gd], low2_cost), markup * low_cost)

            m[im, min_ex] += (price_charged )^(one(σ) - σ)

            sum_price[im] += (price_charged )^(one(σ) - σ)

            rec_low_price[im, gd] = price_charged

            # rec_cost[im, gd] = low_cost

        end

    end

    ###########################################################
    # Loop to calculate aggregate price index and the trade shares

    for im in 1:Ncntry

        g_val = (sum_price[im] * inv_Ngoods)

        for ex in 1:Ncntry

            m[im, ex] = inv_Ngoods*( m[im, ex] ) / g_val

        end

    end

    return m, rec_low_price

end

###############################################################
# this is the version with tariffs using muliple dispatch so when I call the function
# I can just pass the tariff matrix and it will work

function sim_trade_pattern_ek(S, d, τ, θ, σ; Ngoods = 100000, code = 1)
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

    #println(code)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
        #not sure that this helped... may need to return back to the original code

    end

    ###############################################################

    # Loop to calculate the low price and country suppliers
    m = zeros(Ncntry, Ncntry) # need to be zero as I'm summing over stuff

    τ_rev = zeros(Ncntry, Ncntry)

    Φ = zeros(Ncntry)

    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods  # Loop over goods # threading here messes stuff up  

        @inbounds for im in 1:Ncntry  # Loop over importing countries

            low_price = p[im, gd]
            min_ex = im

            @inbounds for ex in 1:Ncntry

                cif_price = (1.0 + τ[im, ex] ) * d[im, ex] * p[ex, gd] # price of exporter

                # low_price, min_ex = ifelse(cif_price < low_price, (cif_price, ex), (low_price, min_ex)) 

                if cif_price < low_price # if the price is lower than the current low price

                    low_price = cif_price # it is the low price
                    
                    min_ex = ex # and the exporter is the one with the lowest price
                end

            end

            ###############################################################
            # This is an alternative way to find the low cost exporter

            # cif_price = d[im, :] .* p[:, gd]

            # sorted_price = sort(cif_price)

            # low_price = sorted_price[1]

            # min_ex = findfirst(==(low_price), cif_price)

            # # This ==(low_price) creates an anonymous function that checks if an element is equal to low_price.

            # ###############################################################

            # Update trade matrix `m`
            m[im, min_ex] += low_price^(1.0 - σ)
            
            τ_rev[im, min_ex] += τ[im, min_ex] * ( low_price^(1.0 - σ) / (1.0 + τ[im , min_ex]) )
            # this then computes tariff revenue with two caveats...
            # need to divide through by the price index in the importing country, just like the m matrix
            # this is on a per-unit basis, so need to multiply by total demand to get total tariff revenue

            # Update sum price and record lowest price
            sum_price[im] += low_price^(1.0 - σ) 

            rec_low_price[im, gd] = low_price
        end

    end

    # Loop to calculate aggregate price index and the trade shares.

    for im in 1:Ncntry

        Φ[im] = (sum_price[im] * inv_Ngoods)

        for ex in 1:Ncntry

            m[im, ex] = inv_Ngoods*( m[im, ex] ) / Φ[im]

            τ_rev[im, ex] = inv_Ngoods*( τ_rev[im, ex] ) / Φ[im]

        end

    end

    P = Φ.^(1.0 / (1.0 - σ)) # this is the CES price index

    return m, P, τ_rev

end



