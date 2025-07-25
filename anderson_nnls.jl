using Octavian

using LinearAlgebra
using Printf, HDF5

# datafit_and_penalty.jl

# Abstract types to define the interface
abstract type AbstractDatafit end
abstract type AbstractPenalty end

#=
--------------------------------------------------------------------------------
                                DATAFIT
--------------------------------------------------------------------------------
=#

"""
Quadratic datafit (i.e., Ordinary Least Squares).
The datafit is: 0.5 * ||y - Xw||² / n_samples
"""
mutable struct QuadraticDatafit{T <: AbstractFloat} <: AbstractDatafit
    Xty::Vector{T}

    # Constructor for initialization
    QuadraticDatafit{T}() where T = new{T}(T[])
end

# Alias for the example file
const NNLSDatafit = QuadraticDatafit

"""
Initialize the datafit by pre-computing X' * y.
"""
function initialize!(
    df::QuadraticDatafit{T},
    X::AbstractMatrix{T},
    y::AbstractVector{T}
) where T <: AbstractFloat
    df.Xty = X' * y
end

"""
Compute the value of the datafit.
"""
function value(
    ::QuadraticDatafit{T},
    y::AbstractVector{T},
    Xw::AbstractVector{T}
) where T <: AbstractFloat
    n_samples = length(y)
    return sum((y .- Xw).^2) / (2 * n_samples)
end

"""
Compute the gradient for a single coordinate j.
"""
@inline function gradient_scalar(
    df::QuadraticDatafit{T},
    X::AbstractMatrix{T},
    Xw::AbstractVector{T},
    j::Int
) where T <: AbstractFloat
    n_samples = size(X, 1)
    # Use view to avoid allocating a new vector for X[:, j]
    return (dot(view(X, :, j), Xw) .- df.Xty[j]) ./ n_samples
end


"""
Compute coordinate-wise Lipschitz constants.
"""
function get_lipschitz(
    ::QuadraticDatafit{T},
    X::AbstractMatrix{T}
) where T <: AbstractFloat
    n_samples, n_features = size(X)
    lipschitz = zeros(T, n_features)
    for j in 1:n_features
        lipschitz[j] = dot(view(X, :, j), view(X, :, j)) / n_samples
    end
    return lipschitz
end

"""
Compute the update for the intercept term.
"""
function intercept_update_step(
    ::QuadraticDatafit{T},
    y::AbstractVector{T},
    Xw::AbstractVector{T}
) where T <: AbstractFloat
    return mean(Xw .- y)
end


#=
--------------------------------------------------------------------------------
                                PENALTY
--------------------------------------------------------------------------------
=#

"""
Positivity constraint penalty (w >= 0).
Value is 0 if w >= 0, and ∞ otherwise.
"""
struct NNLSPenalty <: AbstractPenalty end

"""
Compute the value of the penalty.
"""
function value(::NNLSPenalty, w::AbstractVector{T}) where T <: AbstractFloat
    return any(w .< 0) ? T(Inf) : T(0)
end

"""
Compute the proximal operator for a single coordinate.
"""
function prox_1d(::NNLSPenalty, val::T, ::T, ::Int) where T <: AbstractFloat
    return max(zero(T), val)
end

"""
Compute the distance of the negative gradient to the subdifferential at w.
This is used to measure the violation of the optimality conditions.
"""
function subdiff_distance(
    ::NNLSPenalty,
    w::AbstractVector{T},
    grad::AbstractVector{T},
    ws::AbstractVector{Int}
) where T <: AbstractFloat
    dist = zeros(T, length(ws))
    for (idx, j) in enumerate(ws)
        if w[j] ≈ 0
            # subdifferential is (-∞, 0], distance of -grad[j] to it
            dist[idx] = max(zero(T), -grad[idx])
        else # w[j] > 0
            # subdifferential is {0}, distance of -grad[j] to it
            dist[idx] = abs(grad[idx])
        end
    end
    return dist
end

"""
Return a boolean vector indicating which features are penalized.
For NNLS, all features are constrained.
"""
is_penalized(::NNLSPenalty, n_features::Int) = trues(n_features)

"""
Return a boolean mask of the non-zero coefficients (the support).
"""
generalized_support(::NNLSPenalty, w::AbstractVector{T}) where T = w .!= 0

"""
Options for the AndersonCD solver.
"""
struct AndersonCDOptions{T <: AbstractFloat}
    max_iter::Int
    max_epochs::Int
    p0::Int
    tol::T
    ws_strategy::Symbol
    fit_intercept::Bool
    verbose::Bool
end

# Default constructor
function AndersonCDOptions{T}(;
    max_iter=50, max_epochs=50000, p0=10, tol=1e-4,
    ws_strategy=:subdiff, fit_intercept=true, verbose=false
) where T <: AbstractFloat
    return AndersonCDOptions{T}(
        max_iter, max_epochs, p0, T(tol), ws_strategy, fit_intercept, verbose
    )
end


"""
Helper for Anderson Acceleration.
"""
mutable struct AndersonAcceleration{T <: AbstractFloat}
    K::Int
    current_iter::Int
    arr_w::Matrix{T}
    arr_Xw::Matrix{T}
    # Initialize with dummy arrays, will be resized on first call
    AndersonAcceleration{T}(K::Int) where T = new(K, 0, zeros(T, 0, K+1), zeros(T, 0, K+1))
end


"""
Extrapolate a new solution based on previous iterates.
"""
function extrapolate(
    accel::AndersonAcceleration{T},
    w::AbstractVector{T},
    Xw::AbstractVector{T}
) where T <: AbstractFloat
    # Initialize arrays on first call
    if isempty(accel.arr_w)
        accel.arr_w = zeros(T, length(w), accel.K + 1)
        accel.arr_Xw = zeros(T, length(Xw), accel.K + 1)
    end

    if accel.current_iter <= accel.K
        accel.arr_w[:, accel.current_iter + 1] = w
        accel.arr_Xw[:, accel.current_iter + 1] = Xw
        accel.current_iter += 1
        return w, Xw, false # not extrapolated
    end

    # Form matrix of differences (residuals)
    U = diff(accel.arr_w, dims=2)

    # Solve linear system to find extrapolation coefficients
    try
        # U'U \ ones(K) is equivalent to solving (U'U)c = ones(K)
        inv_UTU_ones = matmul(U', U) \ ones(T, accel.K)
        C = inv_UTU_ones / sum(inv_UTU_ones)

        # Extrapolate and reset for next time
        w_acc = view(accel.arr_w, :, 2:(accel.K+1)) * C
        Xw_acc = view(accel.arr_Xw, :, 2:(accel.K+1)) * C
        accel.current_iter = 0
        return w_acc, Xw_acc, true # extrapolated
    catch e
        if isa(e, SingularException)
            accel.current_iter = 0
            return w, Xw, false # failed to extrapolate
        else
            rethrow(e)
        end
    end
end

"""
Run one epoch of coordinate descent on the features in the working set.
"""
function _cd_epoch!(
    X::AbstractMatrix{T},
    y::AbstractVector{T},
    w::AbstractVector{T},
    Xw::AbstractVector{T},
    lipschitz::AbstractVector{T},
    datafit::AbstractDatafit,
    penalty::AbstractPenalty,
    ws::AbstractVector{Int}
) where T <: AbstractFloat
    for j in ws
        stepsize = lipschitz[j] != 0 ? one(T) / lipschitz[j] : T(1000.0)
        old_w_j = w[j]
        grad_j = gradient_scalar(datafit, X, Xw, j)
        w[j] = prox_1d(penalty, old_w_j - stepsize * grad_j, stepsize, j)

        diff = w[j] - old_w_j
        if diff != 0
            # Use view to avoid allocation and update Xw in-place
            Xw .+= diff .* view(X, :, j)
        end
    end
end


"""
Compute the gradient of the datafit restricted to the working set `ws`.
"""
function construct_grad(
    X::AbstractMatrix{T},
    Xw::AbstractVector{T},
    datafit::AbstractDatafit,
    ws::AbstractVector{Int}
) where T <: AbstractFloat
    n_samples = size(X, 1)
    return (matmul(view(X, :, ws)', Xw) .- datafit.Xty[ws]) ./ n_samples
end

"""
Compute the fixed point distance for optimality check.
"""
function dist_fix_point_cd(
    w::AbstractVector{T},
    grad_ws::AbstractVector{T},
    lipschitz_ws::AbstractVector{T},
    penalty::AbstractPenalty,
    ws::AbstractVector{Int}
) where T <: AbstractFloat
    dist = zeros(T, length(ws))
    for (idx, j) in enumerate(ws)
        if lipschitz_ws[idx] == 0
            continue
        end
        step_j = 1 / lipschitz_ws[idx]
        dist[idx] = abs(
            w[j] - prox_1d(penalty, w[j] - step_j * grad_ws[idx], step_j, j)
        )
    end
    return dist
end


"""
Main solver function using Coordinate Descent with Anderson Acceleration.
"""
function solve(
    X::AbstractMatrix{T},
    y::AbstractVector{T},
    datafit::AbstractDatafit,
    penalty::AbstractPenalty,
    options::AndersonCDOptions{T},
    w_init::Union{Nothing, AbstractVector{T}} = nothing,
    Xw_init::Union{Nothing, AbstractVector{T}} = nothing,
) where T <: AbstractFloat
    n_samples, n_features = size(X)

    # Initialize coefficients and model fit
    w = isnothing(w_init) ? zeros(T, n_features + options.fit_intercept) : copy(w_init)
    Xw = isnothing(Xw_init) ? zeros(T, n_samples) : copy(Xw_init)
    
    # Pre-computations
    initialize!(datafit, X, y)
    lipschitz = get_lipschitz(datafit, X)
    
    pen = is_penalized(penalty, n_features)
    unpen = .!pen
    n_unpen = sum(unpen)
    all_feats = 1:n_features

    obj_out = T[]
    stop_crit = T(Inf)

    # Anderson acceleration buffers
    w_acc = zeros(T, n_features + options.fit_intercept)
    Xw_acc = zeros(T, n_samples)
    
    # Main loop
    for t in 1:options.max_iter
        grad = construct_grad(X, Xw, datafit, all_feats)

        # Optimality condition check
        if options.ws_strategy == :subdiff
            opt = subdiff_distance(penalty, view(w, 1:n_features), grad, all_feats)
        elseif options.ws_strategy == :fixpoint
            opt = dist_fix_point_cd(view(w, 1:n_features), grad, lipschitz, penalty, all_feats)
        else
            error("Unsupported ws_strategy: $(options.ws_strategy)")
        end

        intercept_opt = options.fit_intercept ? abs(intercept_update_step(datafit, y, Xw)) : zero(T)
        stop_crit = max(maximum(opt), intercept_opt)

        if options.verbose
            @printf("Iter %d, Stop crit: %.2e\n", t, stop_crit)
        end

        if stop_crit <= options.tol
            break
        end

        # Build working set (ws)
        g_support_size = sum(generalized_support(penalty, view(w, 1:n_features)))
        ws_size = max(
            min(options.p0 + n_unpen, n_features),
            min(2 * g_support_size - n_unpen, n_features)
        )

        # Force inclusion of unpenalized features and current support
        opt[unpen] .= T(Inf)
        opt[generalized_support(penalty, view(w, 1:n_features))] .= T(Inf)
        
        # Get indices of the `ws_size` largest values in `opt`
        ws = partialsortperm(opt, 1:ws_size, rev=true)

        if options.verbose
            println("Solving subproblem with $(length(ws)) features.")
        end
        
        accelerator = AndersonAcceleration{T}(5)

        # Inner loop (epochs on subproblem)
        for epoch in 1:options.max_epochs
            _cd_epoch!(X, y, view(w, 1:n_features), Xw, lipschitz, datafit, penalty, ws)

            if options.fit_intercept
                intercept_old = w[end]
                w[end] -= intercept_update_step(datafit, y, Xw)
                Xw .+= (w[end] - intercept_old)
            end

            # Anderson Acceleration step
            ws_intercept = options.fit_intercept ? vcat(ws, -1) : ws # -1 is a dummy for the intercept
            w_subset = options.fit_intercept ? vcat(view(w, ws), w[end]) : view(w, ws)
            
            # This part is tricky. A simpler approach is to accelerate the full vector.
            w_accel, Xw_accel, is_extrap = extrapolate(accelerator, w, Xw)

            if is_extrap
                p_obj = value(datafit, y, Xw) + value(penalty, view(w, 1:n_features))
                p_obj_acc = value(datafit, y, Xw_accel) + value(penalty, view(w_accel, 1:n_features))
                
                if p_obj_acc < p_obj
                    w, Xw = w_accel, Xw_accel
                end
            end

            # Inner loop stopping criterion
            if epoch % 10 == 0
                grad_ws = construct_grad(X, Xw, datafit, ws)
                opt_ws = subdiff_distance(penalty, view(w, 1:n_features), grad_ws, ws)
                stop_crit_in = maximum(opt_ws)

                if stop_crit_in < 0.3 * stop_crit
                    if options.verbose
                        println("   Epoch $epoch, inner stop crit: $(stop_crit_in), early exit.")
                    end
                    break
                end
            end
        end # end epoch loop
        
        current_obj = value(datafit, y, Xw) + value(penalty, view(w, 1:n_features))
        push!(obj_out, current_obj)

    end # end iter loop

    # Return solution, objectives, and final stopping criterion
    return w, obj_out, stop_crit
end

function nnls_solve(X::AbstractMatrix{T}, y::AbstractVector{T}, thread::Int=1, verbose::Bool=false) where T <: AbstractFloat
    BLAS.set_num_threads(thread)
    datafit = NNLSDatafit{T}()
    penalty = NNLSPenalty()

    options = AndersonCDOptions{T}(
        fit_intercept=false,
        verbose=verbose
    )

    w_sol, _, stop_crit = solve(X, y, datafit, penalty, options)

    return w_sol
end

function load_nnls_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A"]), read(file["B"]), read(file["B_depth_start"])
    end
end

function load_fail_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A_fail"])
    end
end

function nnls_solve_partial(X::AbstractMatrix{T}, y::AbstractVector{T},
                            idx::Vector{Int}, verbose::Bool=false) where T <: AbstractFloat
    BLAS.set_num_threads(1)
    datafit = NNLSDatafit{T}()
    penalty = NNLSPenalty()

    options = AndersonCDOptions{T}(
        fit_intercept=false,
        verbose=verbose
    )

    w_sol, _, stop_crit = solve(view(X, :, idx .+ 1), y, datafit, penalty, options)
    predict_suc_B = view(X, :, idx .+ 1) * w_sol
    
    return w_sol, predict_suc_B
end


using Base.Threads
function run_nnls_map(X::AbstractMatrix{T}, y::AbstractVector{T},
                      thread_data::Vector{Vector{Int}}, verbose::Bool=false) where T <: AbstractFloat
    tasks = [@spawn nnls_solve_partial(X, y, idx, verbose) for idx in thread_data]
    return [fetch(t) for t in tasks]
end
