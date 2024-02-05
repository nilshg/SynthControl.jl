using GLM: lm, @formula, DummyCoding

function inter_fe_mc(Y, Y0, X, I, W, β₀, r, λ, force, tol=1e-5, max_iter=1_000)

    T, N = size(Y) #!# Double check - this is inconsistent with usual format
    p = size(X, 3)
    obs = sum(I)
    n_iter = 0
    validF = 1
    μ = 0
    μ_Y = 0
    
    use_weight = 0
    
    if size(Y) == size(W)
        use_weight = 1
        W ./= max(W)
        WI = W .* I
    end

    α = zeros(N, 1)
    xᵢ = zeros(T, 1)
    fit = zeros(T, N)

    # Duplicate data
    YY = copy(Y)
    XX = copy(X)

    # Check if XX has enough variation
    p1 = p
    X_invar = sum.(eachslice(abs.(XX), dims=3)) .< 1e-5
    
    if any(X_invar)
        XX = X[:, :, filter(!in(X_invar))]
        p1 -= sum(X_invar)
    end

    validX = 1
    if p1 == 0 
        validX = 0
        if force == 0 && r == 0 # no covariate and force == 0 and r == 0
            if use_weight
                μ_Y = sum(YY .* WI)/sum(WI)
            else
                μ_Y = sum(YY) / obs
            end
            
            μ = μ_Y
            YY = FE_adj(YY .- μ_Y, I)
        end
    end

    # Main algorithm
    if p1 == 0
        if r > 0
            # Add FE; inter FE; iteration

            # Soft impute as starting value
            fe_ad_inter = fe_ad_inter_iter(YY, Y0, I, W, force, 1, 1, false, λ, tol, max_iter)

            # Hard impute
            (; fit, μ, niter, validF) = fe_ad_inter
            U = fe_ad_inter.e
            
            if (force == 1) || (force == 3)
                α = fe_ad_inter.α
            elseif (force == 2) || (force == 3)
                xᵢ = fe_ad_inter.xᵢ
            end
        
        else

            if force == 0
                U = YY
                fit .= μ
                validF = 0
            else
                # Add FE, iteration
                fe_ad = fe_ad_iter(YY, Y0, I, W, force, tol, max_iter)
                (; μ, fit, niter) = fe_ad
                U = fe_ad.e
                
                if (force == 1) || (force = 3)
                    α = fe_ad.α
                elseif (force == 2) || (force == 3)
                    xᵢ = fe_ad.xᵢ
                end
            end
        end
    else
        # Starting value: the OLS estimator
        invXX = use_weight ? wXXinv(XX, W) : XXinv(XX) #!# is use_weight an Int or  a Bool?
            
        if r == 0 #!# r should be an Int
            # Add FE, covar, iteration
            fe_ad = fe_ad_covar_iter(XX, invXX, YY, Y0, I, β₀, W, force, tol, max_iter)
            (; μ, β, U, fit, niter) = fe_ad
            if (force == 1) || (force == 3)
                α = fe_ad.α
            end
            if (force == 2) || (force == 3)
                xᵢ = fe_ad.xᵢ
            end
            validF = 0
        elseif r > 0
            # Add covar, interactive, iteration
            # Soft impute as starting value
            fe_ad_inter_covar = fe_ad_inter_covar_iter(XX, invXX, YY, Y0, I, W, β₀, force, 1, 1, 0, λ, tol, max_iter)
            (; β, fit, μ, niter, validF) = fe_ad_inter_covar
            U = fe_ad_inter_covar.e
            if (force == 1) || (force == 3)
                α = fe_ad_inter_covar.α
            end
            if (force == 2) || (force == 3)
                xᵢ = fe_ad_inter_covar.xᵢ
            end
        end
    end

    #############
    ## Storage ##
    #############

    β̂ = zeros(p)

    if p > 0
        if p > p1 # Some covariates had insufficient variation, set their coefficients to NaN
            j4 = 0
            for i ∈ 1:p
                if X_invar[i]
                    β̂[i] = NaN
                else
                    β̂[i] = β[j4]
                    j4 += 1
                end
            end
        else # All covariates estimated
            β̂ = β
        end
    end
    
    return (; β̂ = β̂, μ = μ, fit = fit, validF = validF, niter = niter, α = α, xᵢ = xᵢ, residuals = U, validX = validX)
end

function XXinv(X)
    # Three dimensional matrix inverse
    p = size(X, 3)
    xx = zeros(p, p)
    for k ∈ 1:p
        for m ∈ 1:p
            if k > m
                xx[k, m] = xx[m, k]
            else
                xx[k,m] = @views tr(X[:, :, k]' * X[:, :, m])
            end
        end
    end
    
    return inv(xx)
end

function wXXinv(X, w)
    p = size(X, 3)
    w_sr = (√).(w)
    xx = zeros(p, p)
    for k ∈ 1:p
        for m ∈ 1:p
            if k > m
                xx[k, m] = xx[m, k]
            else
                xx[k, m] = tr((w_sr .* X[:, :, k])' * (w_sr .* X[:, :, m])) 
            end
        end
    end

    return inv(cholesky(xx)) #!# cholesky as we assume matrix is symmetric PD here - inv_sympd in C++
end

function FE_adj(FE, I)
    # Reset FEᵢⱼ = 0 if Iᵢⱼ == 0
    x = copy(FE)
    x[I .== 0] .= 0.0
    return x
end

function fe_ad_iter(Y, Y0, I, W, force, tolerate, max_iter = 500)
    T, N = size(Y)
    μ = 0.0
    use_weight = size(Y) == size(W) #!# this is an Int in the C++ code

    fit = Y0 # Initial value
    fit_old = Y0
    α = zeros(N, 1)
    xᵢ_Y = zeros(T, 1)
    ε = zeros(T, N) # residual
    YY = copy(Y)

    dif = 1.0; niter = 0
    while (dif > tolerate && niter <= max_iter)
        if use_weight 
            YY = wE_adj(Y, fit, W, I) # e step: expectation
        else
            YY = E_adj(Y, fit, I)
        end

        Y_ad = Y_demean(YY, force)
        μ_Y = Y_ad.μ_Y
        if in((1, 3))(force)
            α_Y = Y_ad.α_Y
        elseif in((2, 3))(force)
            xᵢ_Y = Y_ad.xᵢ_Y
        end
        Y_fe_ad = fe_add(α_Y, xᵢ_Y, μ_Y, T, N, force) # m step: estimate FE

        fit = Y_fe_ad.FE_ad

        if use_weight
            dif = norm(W .* (fit .- fit_old))/norm(W .* fit_old)
        else
            dif = norm(fit .- fit_old)/norm(fit_old)
        end

        fit_old = fit

        niter += 1
    end

    e = FE_adj(YY .- fit, I)

    return (; μ = Y_fe_ad.μ, fit, niter, e, α, xᵢ)
end

function Y_demean(Y, force) ##test## Checked against R
    T, N = size(Y)
    μ_Y = 0.0
    α_Y = zeros(N)
    xᵢ_Y = zeros(T)
    YY = copy(Y)

    μ_Y = sum(YY)/(N*T) #!# that should just be mean(Y)?
    if force == 0
        YY = YY .- μ_Y
    end

    # Unit fixed effects
    if force == 1
        α_Y .= vec(mean(YY, dims = 1)) # colMeans, N×1 matrix
        YY .-= α_Y'
    end

    # Time fixed effects
    if force == 2
        xᵢ_Y .= vec(mean(YY, dims = 2)) # rowMeans, N×1 matrix (#!# should be T×1?)
        YY .-= xᵢ_Y
    end

    # Unit and time fixed effects
    if force == 3
        α_Y .= vec(mean(YY, dims = 1))
        xᵢ_Y .= vec(mean(YY, dims = 2))
        YY .=  (YY .- α_Y' .- xᵢ_Y) .+ μ_Y 
    end

    return (; YY, μ_Y, α_Y, xᵢ_Y)
end

function fe_add(α_Y, xᵢ_Y, μ_Y, T, N, force)##test note## tested against R
    FE_ad = zeros(T, N)
    μ = 0.0
    α = zeros(N, 1)
    xᵢ = zeros(T, 1)
    μ = μ_Y

    if in((1, 3))(force)
        α = α_Y .- μ_Y
    end

    if in((2, 3))(force)
        xᵢ = xᵢ_Y .- μ_Y
    end

    FE_ad .+= μ

    if in((1, 3))(force)
        FE_ad .+= α'
    end

    if in((2, 3))(force)
        FE_ad .+= xᵢ
    end

    return (; μ, FE_ad, α, xᵢ)
end

function E_adj(E, FE, I) 
    # Expectation :E if Iᵢⱼ == 0 then Eᵢⱼ = FEᵢⱼ
    @assert size(E) == size(FE) == size(I)
    EE = copy(E)
    
    for i ∈ eachindex(E)
        if I[i] == 0
            EE[i] = FE[i]
        end
    end

    return EE
end

function wE_adj(E, FE, W, I)
    @assert size(E) == size(FE) == size(W) == size(I)
    EE = copy(E)

    for i ∈ eachindex(E)
        if I[i] == 0
            EE[i] = FE[i]
        else
            EE[i] = (1-W[i])*FE[i] + W[i]*E[i]
        end
    end

    return EE
end

# Obtain additive FE for ub data; assume r>0 but p=0
function fe_ad_inter_iter(Y, Y0, I, W, force, mc, r, hard, λ, tolerate, max_iter)
    T, N = size(Y)
    μ = 0.0
    dif = 1.0
    niter = 0
    validF = 1 # whether has a factor structure
    use_weight = size(Y) == size(W)
    r_burnin = 0
    d = T <= N ? T : N

    fit = Y0
    fit_old = copy(fit)
    stop_burnin = 0
    YY = use_weight ? wE_adj(Y, fit, W, I) : E_adj(Y, fit, I)
    ife_inner = ife(YY, force, 1, 1, hard, λ)

    while dif > tolerate && niter <= max_iter
        YY = use_weight ? wE_adj(Y, fit, W, I) : E_adj(Y, fit, I)

        if mc == 0
            if use_weight && stop_burnin == 0
                r_burnin = d - niter
                if r_burnin <= r
                    r_burnin = r
                end
                ife_inner = ife(YY, force, 0, r_burnin, 0, 0)
            else
                ife_inner = ife(YY, force, 0, r, 0, 0)
            end
        else
            ife_inner = ife(YY, force, 1, 1, hard, λ)
        end
        
        fit = ife_inner.FE
        
        if use_weight
            dif = norm(W .* (fit .- fit_old))/norm(W .* fit_old)
        else
            dif = norm(fit .- fit_old)/norm(fit_old)
        end

        fit_old = fit
        niter += 1

        if (dif <= tolerate) && (niter <= d) && use_weight && (stop_burnin == 0) && (mc == 0)
            stop_burnin = 1
            dif = 1.0
            niter = 0
            fit = Y0
            fit_old = fit
        end
    end
    
    e = FE_adj(YY .- fit, I)
    FE_inter_use = ife_inner.FE_inter_use

    if sum(abs, FE_inter_use) < 1e-10
        validF = 0
    end

    α = in((1, 3))(force) ? ife_inner.α : zeros(N, 1)
    xᵢ = in((2, 3))(force) ? ife_inner.xᵢ : zeros(T, 1)

    return (; μ = ife_inner.μ, niter, fit, e, validF, α, xᵢ)
end

# Factor analysis: μ add ife
function ife(E, force, mc, r, hard, λ) # Third parameter to set pac or mc method
    T, N = size(E)

    E_ad = Y_demean(E, force)
    EE = E_ad.YY
    μ_E = E_ad.μ_Y

    if in((1, 3))(force)
        α_E = E_ad.α_Y
    end
    if in((2, 3))(force)
        xᵢ_E = E_ad.xᵢ_Y
    end

    E_fe_ad = fe_add(α_E, xᵢ_E, μ_E, T, N, force)

    FE_add_use = E_fe_ad.FE_ad # additive FE

    FE_inter_use = zeros(T, N)

    if r > 0
        if mc == 0
            pf = panel_factor(EE, r)
            F = pf.factor
            L = pf.λ
            VNT = pf.VNT
            FE_inter_use .= F * L' # interactive FE
        else
            FE_inter_use .= panel_FE(EE, λ, hard)
        end
    end

    FE = FE_add_use .+ FE_inter_use

    α = in((1, 3))(force) ? E_fe_ad.α : zeros(N, 1)
    xᵢ = in((2, 3))(force) ? E_fe_ad.xᵢ : zeros(T, 1)

    return (; μ = E_fe_ad.μ, FE, α, xᵢ, FE_inter_use)

end

# Obtain factors and loading given error
function panel_factor(E, r)
    T, N = size(E)

    if T < N
        EE = (E * E') / (N*T)
        U, s = svd(EE) # Don't need V
        factor = U[:, 1:r] .* √(T)
        λ = E' * factor ./ T
    else
        EE = (E' * E) / (N*T)
        U, s = svd(EE)
        λ = U[:, 1:r] .* √(N)
        factor = E * λ ./ N
    end

    VNT = Diagonal(s[1:r])

    FE = factor * λ'

    return (; λ, factor, VNT, FE)
end

# Obtain interactive fe directly
function panel_FE(E, λ, hard)
    T, N = size(E)
    r = (T >= N) ? N : T #!# I think this is basically r = min(T, N)
    n_obs = T*N

    D = zeros(r, r) #!# Should probably be Diagonal(zeros(r))
    U, s, V = svd(E ./ n_obs)

    for i ∈ 1:r
        if s[i] > λ
            if hard 
                D[i, i] = s[i]
            else
                D[i, i] = s[i] - λ
            end
        else
            D[i, i] = 0.0
        end
    end

    if T >= N
        UU = U[:, 1:r]
        FE = UU * D * V' * n_obs
    else
        VV = V[:, 1:r] #!# check whether this is correct given transposed V in Julia
        FE = U * D * VV' * n_obs
    end

    return FE
end

function initialFit(data, force, oci; weight = nothing)
    N = length(unique(data[:, 2]))
    T = length(unique(data[:, 3]))
    p = size(data, 2) - 3

    x = nothing; x_sub = nothing
    if p > 0
        x = data[:, 4:end]
        x_sub = x[oci, :]
    end

    y = data[:, 1]
    y_sub = y[oci]
    β₀ = zeros(1, 1)

    use_weight = !isnothing(weight)
    if use_weight
        w = weight
        w_sub = w[oci]
    end

    ind = if force == 0
        nothing
    elseif force == 1
        data[:, 2]
    elseif force == 2
        data[:, 3]
    elseif force == 3
        data[:, [2, 3]]
    end

    Y₀ = zeros(T, N)
    if force == 0
        if p == 0 # No FEs, no covariates - prediction is simply (weighted) average Y
            if use_weight
                μ = sum(y_sub .* w_sub)/sum(w_sub)
                Y₀ .= μ
            else
                μ = mean(y_sub)
                Y₀ .= μ
            end
        else # No FEs, but covariates - need regression on covariates without FEs
            ols = if use_weight
                reg([DataFrame(y = y_sub) DataFrame(x_sub, :auto) DataFrame(w = w_sub)],
                    term("y") ~ sum(term.("x$i" for i ∈ 1:p)), weights = :w)
            else
                reg([DataFrame(y = y_sub) DataFrame(x_sub, :auto)],
                    term("y") ~ sum(term.("x$i" for i ∈ 1:p)))
            end

            cf = coef(ols)
            μ = first(cf)
            β₀ = cf[2:end]

            y0 = μ .+ x * β₀
            #!# to complete - this doesn't really make sense 
            #!# this whole section should really just be
            #!# formula = term("y") ~ sum(term.("x$i" for i ∈ 1:p)) + fe(:t) + fe(:id)
            #!# fe_model = reg(df[oci, :], formula, weights = :w)
            #!# Y₀ = predict(fe_model, df)
            #!# β₀ = coef(fe_model)
        end
    else
        df = DataFrame(y = y, id = data[:, 2], t = data[:, 3])
        if p > 0
            df = [df DataFrame(x, :auto)]
        end

        #!# This should be FixedEffectModels but that is buggy for only FE models

        f = @formula(y ~ id + t)
        cs = Dict(:id => DummyCoding(), :t => DummyCoding())

        if p > 0
            f += sum(term.("x$i" for i ∈ 1:p))
        end

        fe_model = if use_weight
            lm(f, df[oci, :], wts = w_sub; contrasts = cs) 
        else 
            lm(f, df[oci, :]; contrasts = cs)
        end

        Y₀ .= reshape(predict(fe_model, df), T, N)
    end

    return (; Y₀, β₀)
end

function cv_sample(I, D, my_count; cv_count = 3, cv_treat = false, cv_donut = 1)
    ## Cross validation sampling

    TT, N = size(I)
    tr_pos = findall(any.(eachcol(D)))
    D_fake = zeros(size(D)); D_fake[:, tr_pos] .= 1
    prop = TT*N/sum(I)

    if !cv_treat
        oci = findall(vec(I) .== 1)
    else
        oci = findall((vec(I) .== 1) .&& (vec(D) .== 1))
        lenght(oci) > my_count || error("Too few observations are valid for cross validation. Try setting cv_treat to false or a smaller cv_prop.")
    end

    if cv_count <= 2 # randomly missing
        cv_id = sample(oci, my_count; replace = false)
        rm_id_use = cv_id
    else # Not missing at random - take consecutive observations (this is cv_count = 3)
        res = TT % cv_count
        int = Int(floor(TT / cv_count))

        rm_pos = []

        ## randomly select 1/cv.count
        subcount = Int(floor(my_count/cv_count * prop))
        subcount = subcount == 0 ? 1 : subcount

        seq = range(start = 1, step = cv_count, length = int)
        
        if !cv_treat
            rm_pos = reduce(vcat, TT*(i-1) + res .+ seq for i ∈ 1:N)
        else
            error("cv_treat = true not implemented")
        end
        
        rm_id = sample(rm_pos, subcount; replace = false)
        rm_id_all = rm_id
        rm_id_use = Int[]

        for i ∈ 0:(cv_count - 1)
            rm_id_all = vcat(rm_id_all, rm_id .+ i)
            if i >= cv_donut && i <= cv_count - cv_donut - 1
                rm_id_use = vcat(rm_id_use, rm_id .+ i)
            end
        end
        
        intersect!(rm_id_all, oci) ## remove missing values #!# this also removes duplicate indices
        sort!(rm_id_all)

        intersect!(rm_id_use, oci)
        sort!(rm_id_use)

        if length(rm_id_all) >= my_count
            cv_id = rm_id_all[1:my_count]
            pos_cv = findall(rm_id_use .>= minimum(cv_id) .&& rm_id_use .<= maximum(cv_id))
            rm_id_use = rm_id_use[pos_cv]
        else
            cv_id2 = sample(setdiff(oci, rm_id_all), my_count - length(rm_id_all); replace = false)
            cv_id = sort(vcat(rm_id_all, cv_id2))
        end
    end

    return (; cv_id, est_id = rm_id_use)
end

function get_term(d, ii; type = "on")
    if any(d)
        treat_pos = findfirst(d)
        return collect((1:length(d)) .- treat_pos .+ 1)
    else
        return fill(NaN, length(d))
    end
end

function fect_cv(Y, X, D, W, I, II, T_on, force, tol; 
    T_off = nothing, T_on_carry = nothing, T_on_balance = nothing, balance_period = nothing,
    method = "mc", criterion = "mspe", k = 5, cv_prop = 0.1, cv_treat = false, cv_nobs = 3,
    cv_donut = 1, min_T0 = 5, r = 0, r_end, proportion = 0, n_λ = 10, λ = nothing,
    max_iteration = 1_000, norm_para = nothing, group_level = nothing, group = nothing, verbose = true)

    placebo_pos = na_pos = nothing
    TT, N = size(Y)
    p = 0
    X = zeros(1, 1, 0)
    W_use = [0;;]
    use_weight = false
    YY = copy(Y)
    YY[II .== 0] .= 0.0

    T0_min = minimum(sum(II, dims = 1))
    trX = findall(any.(eachcol(D)))
    Ntr = sum(trX .> 0)
    co = findall(all.(eachcol(.!(D))))
    Nco = sum(co .> 0)

    T_on = reduce(hcat, get_term.(eachcol(D), ones(N))) #!# Double check how this works for different treatment patterns
    t_on = vec(T_on)

    # Initial fit using fastplm
    #!# R code uses long data here, i.e. a (N×T)-by-3 matrix 
    data_ini = [vec(Y) repeat(1:N, inner = TT) repeat(1:TT, outer = N)] #!# Double check this
    if p > 0 #!# Check that this works i.e. X is the correct shape
        data_ini = [data_ini reshape(X, N*TT, p)]
    end
    oci = findall(==(1), vec(II)) #!# Again check that this makes sense - checked against R

    Y₀, β₀ = initialFit(data_ini, force, oci)

    ## Restrictions on candidate hyperparameters
    obs_con = (sum(II) - r_end * (N + TT) + r_end^2 - p) <= 0

    if obs_con
        while (sum(II) - r_end * (N + TT) + r_end^2 - p) <= 0
            r_end -= 1
        end
    end

    if r_end >= T0_min
        r_end = T0_min - 1
    end

    r_max = min(TT, r_end)
    r_cv = 0

    ############################
    ####   Main algorithm   ####
    ############################

    validX = 1 ## no multi-colinearity 
    CV_out_ife = nothing; CV_out_mc = nothing

    ## Cross validation of r and λ
    r_old = r
    r_max = min(TT, r_end)
    r_cv = 0
    if verbose
        @info "Cross-validating..."
    end

    #### ----- Initialize ----- ####
    cv_pos = findall(t_on .<= 0)
    t_on_cv = t_on[cv_pos]
    count_on_cv = ones(length(t_on_cv))
    rm_count = floor(Int, sum(II) * cv_prop)
    cv_count = sum(II) - rm_count

    Y0CV = fill(NaN, TT, N, k) ## store initial Y0
    β₀_cv = if p > 0 
        fill(NaN, p, 1, k) 
    else
        fill(0, 1, 0, k) ## store initial beta0
    end

    flag = false

    rmCV = Vector{Vector{Int}}(undef, k) ## removed indicator
    ociCV = Vector{Vector{Int}}(undef, k) ## store indicator
    estCV = Vector{Vector{Int}}(undef, k)
    #!# The following is not in the R code - scoping issue?
    get_cv = cv_sample(II, D, rm_count; cv_count = cv_nobs, cv_treat = cv_treat, cv_donut = cv_donut)
    cv_id = get_cv.cv_id

    
    for i ∈ 1:k
        cv_n = 0
        while true
            cv_n += 1

            get_cv = cv_sample(II, D, rm_count; cv_count = cv_nobs, 
                cv_treat = cv_treat, cv_donut = cv_donut)
            
            cv_id = get_cv.cv_id
            
            II_cv_valid = copy(II); II_cv = copy(II) #!# These should be pre-allocated and re-used
            II_cv[cv_id] .= 0
            #II_cv_valid[cv_id] .= -1 #!# What does this do? Fails if II_cv_valid is a BitMatrix

            con1 = sum(>(0), sum.(eachrow(II_cv))) == TT
            con2 = sum(>(0), sum.(eachcol(II_cv))) == N

            if con1 && con2 
                break
            end

            if cv_n >= 200
                flag = true
                keep1 = findall(<(1), sum.(eachrow(II_cv)))
                keep2 = findall(<(min_T0), sum.(eachcol(II_cv)))
                II_cv[keep1, :] .= II[keep1, :]
                II_cv[:, keep2] .= II[:, keep2]
                cv_id = findall(II_cv_valid .!= II)
                break
            end
        end

        length(cv_id) == 0 && error("Some units have too few pre-treatment observations. Set a larger `cv_prop` or set `cv_treat` to false")

        rmCV[i] = cv_id
        ocicv = setdiff(oci, cv_id)
        ociCV[i] = ocicv

        if cv_n < 200
            estCV[i] = get_cv.est_id
        else
            cv_id_old = get_cv.cv_id
            cv_diff = setdiff(cv_id_old, cv_id)
            estCV[i] = setdiff(get_cv.est_id, cv_diff)
        end

        if !use_weight
            initialOutCv = initialFit(data_ini, force, ocicv)
        else
            initialOutCv = initialFit(data_ini, force, oci; weight = W)
        end

        Y0CV[:, :, i] .= initialOutCv.Y₀

        if p > 0
            β₀_now = initialOutCv.β₀
            replace!(β₀_now, NaN => 0.0)
            β₀_cv[:, :, i] .= β₀_now
        end
    end

    if flag
        @info "Some units have too few pre-treatment observations. Remove them automatically in Cross-Validation.\n"
    end

    # Get count matrix
    if !use_weight
        count_T_cv = T_on
    end
    #!# To complete

    ### Cross validation for MC-NNM ###

    # Create the hyperparameter sequence
    
    # Biggest candidate λ is largest eigenvalue of Y after removing FEs
    Y_λ = YY .- Y₀
    Y_λ[II .== 0] .= 0.0
    eigen_all = svd(Y_λ ./ (TT*N)).S

    if isnothing(λ) || length(λ) == 1
        λ_max = log(10, maximum(eigen_all))
        λ_step = 3/(n_λ-2)
        λ = [[10^(λ_max - i*λ_step) for i ∈ 0:n_λ-2]; 0.0]
    end

    # Store all MSPE
    CV_out_mc = (; λ_norm = λ ./ maximum(eigen_all), MSPE = fill(1e20, length(λ)))
    break_count = 0
    break_check = 0
    MSPE_best = 0
    λ_cv = maximum(λ) #!# This isn't being defined in the R code - might be a scoping difference

    for i ∈ 1:length(λ) # Candidate λs
        SSE = 0.0

        for ii ∈ 1:k # CV folds
            II_cv = copy(II)
            II_cv[rmCV[ii]] .= 0
            YY_cv = copy(YY)
            YY_cv[rmCV[ii]] .= 0
            if use_weight
                #!# not implemented
            else
                W_use2 = zeros(1, 1)
            end
            est_cv_fit = inter_fe_mc(YY_cv, Y0CV[:, :, ii], X, II_cv, W_use2, β₀_cv[:, :, ii], 1, λ[i], force, 
                tol, max_iteration).fit

            if use_weight
                #!# not implemented
            else
                SSE += sum(abs2, YY[estCV[ii]] .- est_cv_fit[estCV[ii]])
            end
        end

        if use_weight 
            #!# not implemented
        else
            MSPE = SSE/sum(length, estCV) #!# This is length(unlist(estCV))
        end

        est_cv = inter_fe_mc(YY, Y₀, X, II, W_use, β₀, 1, λ[i], force, tol, max_iteration)

        eff_v_cv = vec(Y .- est_cv.fit)[cv_pos] #!# cv_pos should just get the right positions in the Y matrix instead
        MSE = sum(abs2, eff_v_cv)/length(eff_v_cv)

        #!# if !isnothing(norm_para) -> Check for normalization 

        if criterion == "mspe"
            if (minimum(CV_out_mc.MSPE) - MSPE) > 0.01*minimum(CV_out_mc.MSPE) # at least 1% MSPE improvement
                MSPE_best = MSPE
                est_best = est_cv
                λ_cv = λ[i]
                break_count = 0
                break_check = 0
            else
                if i > 1
                    if λ_cv == λ[i - 1]
                        break_check = 0
                        break_count = 1
                    end
                end
            end
        else
            #!# not implemented
        end

        if break_check == 1
            break_count += 1
        end

        CV_out_mc.MSPE[i] = MSPE
        if verbose
            println("\tλ_norm = $(round(λ[i]/maximum(eigen_all), digits=5)); MSPE = $(round(MSPE, digits = 2))")
        end
            
        break_count == 3 && break
    end

    est_best_mc = est_best
    MSPE_best_mc = MSPE_best
    if verbose
        println("\n\tλ_norm* = $(round(λ_cv/maximum(eigen_all), digits = 5))")
    end
    
    ## End of cross validation

    est_best = est_best_mc
    validF = est_best.validF
    validX = est_best.validX

    eff = Y .- est_best.fit
    τ̂ᵃᵗᵗ = mean(eff[D])

    return (; τ̂ᵃᵗᵗ, eff)
end

function fect_default(cp; 
    force = "two-way", # Fixed effects
    λ = nothing, n_λ = 10, # λ
    r = (0, 5), # Number of factors - if CV is true then CV procedure will select optimal r ∈ [r, 5]
    CV = true, # whether to CV
    k = 5, # CV folds)
    cv_prop = 0.1, # proportion of CV counts
    cv_treat = false, # whether CV targets treated units
    cv_nobs = 3, # CV taking consecutive units
    cv_donut = 0, # CV MSPE
    criterion = "mspe", 
    method = :mc,
    se = false, # whether to construct standard error estimates,
    vartype = :bootstrap, # or :jackknife
    quantile_CI = false,
    nboots = 200, # number of bootstraps
    α = 0.05, # significance level
    parallel = false, # whether to use parallel computing
    cores = nothing, # number of cores
    tol = 0.001, # tolerance levels
    max_iteration = 1_000,
    min_T0 = nothing, # minimum untreated periods
    max_missing = nothing, # maximum number of missing observations
    verbose = true)    

    Y = collect(cp.Y') #!# Whole routine is written with a T×N matrix
    D = collect(cp.W')    
    TT_old, N_old = size(Y)
    TT, N = TT_old, N_old
    X = zeros(TT, N, 0) # Covariates currently not supported
    

    force = Dict("none" => 0, "unit" => 1, "time" => 2, "two-way" => 3)[force]

    if isnothing(min_T0)
        if method == :mc
            min_T0 = 5
        end
    end

    r_end = length(r) == 1 ? r : maximum(r)
    r = minimum(r)

    TT_old, N_old = size(Y)
    TT, N = TT_old, N_old
    p = 0
    max_missing = isnothing(max_missing) ? TT : max_missing
    prod(size(Y)) == TT * N

    I_D = ones(TT, N)
    I = ones(TT, N)
    Y_ind = Y
    D_ind = D
    D_origin = copy(D)
    II = copy(I)
    II[D] .= 0 # Regard treated values as missing
    
    T_on = fill(NaN, TT, N)
    for i ∈ 1:N
        T_on[:, i] = get_term(D[:, i], I_D[:, i]; type = "on")
    end

    W = zeros(1, 1)

    return fect_cv(Y, X, D, W, I, II, T_on, force, tol; 
        method = method, criterion = criterion, r=r, r_end = r_end, cv_treat = cv_treat, 
        cv_nobs = cv_nobs, cv_donut = cv_donut, min_T0 = min_T0, n_λ = n_λ, λ = λ, 
        max_iteration = max_iteration, verbose = verbose)
end
