using FixedEffectModels
using Parameters
using DataFrames

struct gravity_results{T}
    dist_coef::Array{T} # distance bins
    lang_coef::Array{T} # border, language, ec, efta
    S::Array{T} # source technology part
    θm::Array{T} # asymetric part
end

##########################################################################
##########################################################################

@with_kw struct gravity_params
    Ncntry::Int = 19
    θ::Float64 = 4.0
    L::Array{Float64} = ones(Ncntry)
    dfcntryfix::DataFrame = DataFrame(foo = ones(Ncntry))
end

##########################################################################
##########################################################################

struct trade_costs{T}
    dist_coef::Array{T} # distance bins
    lang_coef::Array{T} # border, language, ec, efta
    θm::Array{T} # asymetric part
end

##########################################################################
##########################################################################

# function gravity!(tradedata, d, T, W, gravity_params)
#     #mulitple dispatch version, runs gravity regression
#     # then fills matrix d with trade costs

#     grv  = gravity(tradedata)

#     make_trade_costs!(tradedata, grv, d, gravity_params)

#     make_technology!(tradedata, grv, d, gravity_params)

# end

function replace_tariff!(tariff_matrix, row_index, new_row)
    tariff_matrix[row_index, :] = new_row

    if tariff_matrix[row_index, row_index] != 0.0
        tariff_matrix[row_index, row_index] = 0.0;
    end

end



function make_ek_dataset(data_path = ".\\ek-data\\")

    # Construct file paths using joinpath
    ek_data_path = joinpath(data_path, "ek-data.csv")
    ek_language_path = joinpath(data_path, "ek-language.csv")
    ek_labor_path = joinpath(data_path, "ek-labor.csv")
    ek_cntryfix_path = joinpath(data_path, "ek-cntryfix.csv")

    # Read CSV files into DataFrames
    dftrade = DataFrame(CSV.File(ek_data_path))
    dflang = DataFrame(CSV.File(ek_language_path))
    dflabor = DataFrame(CSV.File(ek_labor_path))

    # Filter out rows where trade is exactly 1.0 or 0.0
    filter!(row -> ~(row.trade ≈ 1.0), dftrade)
    filter!(row -> ~(row.trade ≈ 0.0), dftrade)

    # Combine dftrade and dflang DataFrames horizontally
    dftrade = hcat(dftrade, dflang)

    # Read the country fix DataFrame
    dfcntryfix = DataFrame(CSV.File(ek_cntryfix_path))

    return dftrade, dfcntryfix, dflabor

end

##########################################################################
##########################################################################

##########################################################################
##########################################################################

function make_technology!(gravity_results, T, W, gravity_params)
    
    @unpack θ, Ncntry = gravity_params
    @unpack S = gravity_results

    for importer = 1:Ncntry

        T[importer] = exp( (S[importer] + θ*log(W[importer])) )
        #equation (27) from EK, but with β = 1 (no round about)

    end

end


function make_technology(gravity_results, W, gravity_params)
    
    @unpack θ, Ncntry = gravity_params
    @unpack S = gravity_results

    T = similar(W)

    for importer = 1:Ncntry

        T[importer] = exp( (S[importer] + θ*log(W[importer])) )
        #equation (27) from EK, but with β = 1 (no round about)

    end

    return T

end

function make_S(W, trade_parameters)
    
    @unpack θ, Ncntry, T = trade_parameters

    S = similar(W)

    for importer = 1:Ncntry

        S[importer] =  log(T[importer]) - θ*log(W[importer])
        #equation (27) from EK, but with β = 1 (no round about)

    end

    return S

end


##########################################################################
##########################################################################

function make_trade_costs!(gravity_results, d, gravity_params; trade_cost_type = "ek")
    # makes the trade costs given fixed country characteristics
    # this it he dffix
    @unpack θ, Ncntry, dfcntryfix = gravity_params
    @unpack dist_coef, lang_coef, θm  = gravity_results

    inv_θ = (1.0 / θ)

    for importer = 1:Ncntry

        foo = dfcntryfix[dfcntryfix.importer .== importer, :]

        for exporter = 1:Ncntry

            if exporter != importer

                get_exporter = foo.exporter .== exporter

                distance_effect = exp(-inv_θ * dist_coef[Int(foo[get_exporter, :].distbin[1])]) 

                border_effect = exp(-inv_θ * lang_coef[1] * (foo[get_exporter, :].border[1])) 

                language_effect = exp(-inv_θ * lang_coef[2] * foo[get_exporter, :].sharedlanguage[1]) 

                europeancom_effect = exp(-inv_θ * lang_coef[3] * foo[get_exporter, :].europeancom[1]) 

                efta_effect = exp(-inv_θ * lang_coef[4] * foo[get_exporter, :].efta[1])

                if trade_cost_type == "ek"

                    asym_effect = exp( -inv_θ * θm[importer] )   
                    
                elseif trade_cost_type == "waugh"

                    asym_effect = exp( -inv_θ * θm[exporter] )

                end  

                d[importer, exporter] =(distance_effect * border_effect  * language_effect
                                         * europeancom_effect * efta_effect * asym_effect)
                                         # equation (29) exponentiated

                d[importer, exporter] = max(d[importer, exporter], 1.0)

            elseif exporter == importer

                d[importer, exporter] = 1.0

            end

        end

    end

end



##########################################################################
##########################################################################

function gravity(tradedata; trade_cost_type = "ek", display = false)
    #function to perform basic gravity regression
    #assumes tradedata is a DataFrame, takes on the structure of EK dataset

    outreg = reg(tradedata, @formula(trade ~ fe(importer) + fe(exporter) +
         bin375 + bin750 + bin1500 + bin3000 + bin6000 + binmax  + border + sharedlanguage +
                europeancom + efta), save = true, tol = 1e-10)

    lang_coef = outreg.coef[7:end]

    if trade_cost_type == "ek"

        S, θm, dist_bins = eaton_kortum_trade_costs(outreg)

        if display == true

            println(outreg)
            println(" ")

            println("Compare to Table III (1762)")
            println(" ")
            println("Distance Effects")
            dffoo = DataFrame(distance_effects = dist_bins);
            println(dffoo)
            println(" ")
            println("Border, language, Eupope, etc. Effects")
            dffoo = DataFrame(boder_lang_effects = lang_coef);
            println(dffoo)
            println(" ")
            println("Source and Destination Effects (The S's and θm's)")
            dffoo = DataFrame(source_effects = S, destination_effects = θm);
            println(dffoo)

        end

    elseif trade_cost_type == "waugh"

        S, θm, dist_bins = waugh_trade_costs(outreg)

        if display == true

            println(outreg)
            println(" ")

            println("Waugh (2010) Formulation")
            println(" ")
            println("Distance Effects")
            dffoo = DataFrame(distance_effects = dist_bins);
            println(dffoo)
            println(" ")
            println("Border, language, Eupope, etc. Effects")
            dffoo = DataFrame(boder_lang_effects = lang_coef);
            println(dffoo)
            println(" ")
            println("Source and Exporter Effects (The S's and θex's)")
            dffoo = DataFrame(source_effects = S, exporter_effects = θm);
            println(dffoo)

        end
    end

    return gravity_results(dist_bins, lang_coef, S, θm)

end

##########################################################################
##########################################################################

function eaton_kortum_trade_costs(outreg)

    grp = groupby(outreg.fe, "exporter")

    S = get_first(grp, "fe_exporter")

    grp = groupby(outreg.fe, "importer")

    θm = get_first(grp, "fe_importer")

    θm = θm .+ S 

    norm_fe = sum(θm) / length(θm)

    θm = θm .- norm_fe

    dist_bins = outreg.coef[1:6] .+ norm_fe

    return S, θm, dist_bins

end

function waugh_trade_costs(outreg)

    grp = groupby(outreg.fe, "importer")

    S = get_first(grp, "fe_importer")

    norm_fe = sum(S) / length(S)

    S = -(S .- norm_fe)

    grp = groupby(outreg.fe, "exporter")

    θm = get_first(grp, "fe_exporter")

    θm = θm .- S
    
    dist_bins = outreg.coef[1:6] .+ norm_fe

    return S, θm, dist_bins

end

##########################################################################
##########################################################################

function gravity_as_guide(xxx, grv_data, gravity_params; trade_cost_type = "ek")
    # mulitple dispatch version to use in solver

    @unpack Ncntry, θ, dfcntryfix = gravity_params 

    T = [exp.(xxx[1:(Ncntry - 1)]); 1] # S's are normalized -> only have 18 degrees of freedom on Ts
    
    θm = [xxx[Ncntry:((Ncntry - 1)*2)]; -sum(xxx[Ncntry:((Ncntry - 1)*2)])] # same with this, they sum zero

    dist_coef = xxx[((Ncntry - 1)*2 + 1):((Ncntry - 1)*2 + 6)] # six distance bins
    
    lang_coef = xxx[((Ncntry - 1)*2 + 7):end] # the language stuff

    # build the trade cost structure 
    trc = trade_costs(dist_coef, lang_coef, θm)

    grv = gravity_as_guide(trc, T, dfcntryfix, gravity_params, trade_cost_type = trade_cost_type)
    # multiple dispacth calls the base file

    outvec = [grv.S[1:end-1] .- grv_data.S[1:end-1] ; 
                grv.θm[1:end-1] .- grv_data.θm[1:end-1] ;
                grv.dist_coef .- grv_data.dist_coef;
                grv.lang_coef .- grv_data.lang_coef]

    return outvec

end

##########################################################################
##########################################################################

function gravity_as_guide(trade_costs, T, dfcntryfix, gravity_params; 
                            solver = true, trade_cost_type = "ek")

    @unpack Ncntry, θ = gravity_params 

    # construct trade costs
    d = zeros(Ncntry,Ncntry)

    make_trade_costs!(trade_costs, d, gravity_params, trade_cost_type = trade_cost_type)

    # then given T's, d's, θ...we can make trade flows
    # and solver for balanced trade

    W = solve_trade_balance(d, T, gravity_params)

    # recover the trade flows to run gravity eq on.
    πshares = eaton_kortum(W, d, T, θ)[1]

    # Organize "model" dataset

    trademodel = log.(normalize_by_home_trade(πshares, Ncntry)')

    dfmodel = DataFrame(trade = vec(drop_diagonal(trademodel, Ncntry)))

    dfmodel = hcat(dfmodel, dfcntryfix)

    # run the gavity regression and return 
    # the gravity structure

    if solver == true

        return gravity(dfmodel, trade_cost_type = trade_cost_type)

    else

        return gravity(dfmodel, trade_cost_type = trade_cost_type), W, πshares, dfmodel

    end

end

##########################################################################
##########################################################################

function solve_trade_balance(xxx, L, d, T, θ)

    W = [ exp.(xxx) ; 1.0]
    # build the wage vector

    πshares = eaton_kortum(W, d, T, θ)[1]
    # make the trade flows

    return trade_balance(W, L, πshares)
    # return trade balance

end

function solve_trade_balance(d, T, gravity_params)
    # muliple dispatch version to find zero

    @unpack L, θ, Ncntry = gravity_params

    f(x) = solve_trade_balance(x, L, d, T, θ);

    function f!(fvec, x)

        fvec .= f(x)

    end

    initial_x = zeros(Ncntry - 1)
    n = length(initial_x)

    sol = fsolve(f!, initial_x, method = :hybr;
        ml= (n - 1), mu= (n - 1),
        diag=ones(n),
        mode= 1,
        tol=1e-15,
        )

    return [exp.(sol.x); 1.0]

end

##########################################################################
##########################################################################

function normalize_by_home_trade(πshares, Ncntry)

    norm_πshares = similar(πshares)

    for importer = 1:Ncntry

        norm_πshares[importer, :] .= πshares[importer, : ] / πshares[importer, importer]

    end
    
    return norm_πshares
    
end 

function drop_diagonal(πshares, Ncntry)

    nodiag = Array{Float64}(undef, Ncntry - 1, Ncntry)

    for exporter = 1:Ncntry

            nodiag[:, exporter] .= deleteat!(πshares[:, exporter], exporter)

    end

    return nodiag

end

##########################################################################
##########################################################################

function make_trade_shares(tradedata, Ncntry)
    # function to go from log, normalized trade data
    # back into the tradeshare matrix

    πdata = Array{Float64}(undef, Ncntry, Ncntry)
    fill!(πdata, 0.0)

    for importer = 1:Ncntry

        foo = tradedata[tradedata.importer .== importer, :]

        for exporter = 1:Ncntry

            if exporter != importer

            get_exporter = foo.exporter .== exporter

            πdata[importer, exporter] = exp(foo[get_exporter, : ].trade[1])

            end

        end

        hometrade = (1.0 + sum(πdata[importer, :]))^(-1.0)

        πdata[importer, :] = πdata[importer, :]*hometrade

        πdata[importer, importer] = hometrade

    end

    return πdata

end

##########################################################################
##########################################################################

function get_first(grp, variable)
       
    var = Array{Float64}(undef, length(grp))

    for i ∈ eachindex(grp)

        var[i[1]] = grp[i][1, variable]

    end
    
    return var

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