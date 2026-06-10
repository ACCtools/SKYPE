using LinearAlgebra
using Statistics
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
    # One full pass over X; parallel when Julia is started with -t N
    Threads.@threads for j in 1:n_features
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
    return any(<(zero(T)), w) ? T(Inf) : T(0)
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
    gram_mem_gb::Float64   # memory budget (GB = 1e9 bytes) for Gram-path buffers
end

# Default constructor
function AndersonCDOptions{T}(;
    max_iter=50, max_epochs=50000, p0=10, tol=1e-4,
    ws_strategy=:subdiff, fit_intercept=true, verbose=false,
    gram_mem_gb=1.0
) where T <: AbstractFloat
    return AndersonCDOptions{T}(
        max_iter, max_epochs, p0, T(tol), ws_strategy, fit_intercept, verbose,
        Float64(gram_mem_gb)
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
        inv_UTU_ones = (U' * U) \ ones(T, accel.K)
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
    return (view(X, :, ws)' * Xw .- datafit.Xty[ws]) ./ n_samples
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


# Compute/cache cap for any Gram dimension. Beyond ~4096 the Gram matrix
# (dim²×8B) falls out of L3 cache and the per-update cost O(dim) approaches
# O(n_samples), so the Gram paths stop winning regardless of available RAM —
# measured on matrix2 (n_features=10091): full-Gram 376s vs ws-Gram 21s.
const GRAM_DIM_CAP = 4096

# Memory the Gram paths materialize: a Float64 copy of the involved columns
# plus the Gram matrix itself, i.e. 8·dim·(n_samples + dim) bytes where dim is
# n_features (`_solve_gram`) or the working-set size (`_solve_subproblem_gram!`).
# Both paths additionally require this to fit the user-facing `gram_mem_gb`
# budget in AndersonCDOptions; the default 1.0 GB saturates the dim cap at
# n_samples ≈ 27k. Raise it when n_samples is much larger (else the budget,
# not GRAM_DIM_CAP, becomes the binding limit); set ≤ 0 to disable Gram paths.
_gram_bytes(n_samples::Int, dim::Int) = 8.0 * dim * (n_samples + dim)

"""
Main solver entry point.

Dispatches to `_solve_gram` (coordinate descent on the precomputed Gram matrix
X'X, like skglm's GramCD) whenever problem shape and the `gram_mem_gb` budget
allow it: each coordinate update then costs O(n_features) instead of
O(n_samples), and the only n_samples-sized work left is the one-time
BLAS-threaded X'X / X'y products.

Otherwise (very wide X or tight budget), `_solve_xw` is used: X is never
densified into a full Gram matrix; instead each working-set subproblem gets its
own |ws|-sized Gram via `_solve_subproblem_gram!` while that fits the budget,
so memory stays bounded for any n_features.
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
    if !options.fit_intercept && n_features <= n_samples &&
       n_features <= GRAM_DIM_CAP &&
       _gram_bytes(n_samples, n_features) <= options.gram_mem_gb * 1e9
        return _solve_gram(X, y, datafit, penalty, options, w_init)
    end
    return _solve_xw(X, y, datafit, penalty, options, w_init, Xw_init)
end

"""
Run one epoch of coordinate descent using the scaled Gram matrix.
`grad` holds (G*w - X'y) / n_samples for ALL features and is kept up to date
in place: changing w[j] only requires adding diff * G[:, j] / n_samples.
"""
function _cd_epoch_gram!(
    G::Matrix{Float64},
    w::Vector{Float64},
    grad::Vector{Float64},
    lipschitz::Vector{Float64},
    penalty::AbstractPenalty,
    ws::AbstractVector{Int},
    inv_n::Float64,
)
    n_features = size(G, 1)
    @inbounds for j in ws
        lj = lipschitz[j]
        stepsize = lj != 0 ? 1.0 / lj : 1000.0
        old_w_j = w[j]
        w[j] = prox_1d(penalty, old_w_j - stepsize * grad[j], stepsize, j)
        diff = w[j] - old_w_j
        if diff != 0
            c = diff * inv_n
            @simd for i in 1:n_features
                grad[i] = muladd(c, G[i, j], grad[i])
            end
        end
    end
end

"""
Datafit value 0.5‖y - Xw‖²/n expressed through the maintained Gram gradient:
0.5‖y - Xw‖²/n = 0.5*(w'grad + (y'y - w'X'y)/n).
"""
@inline function _gram_objective(
    w::AbstractVector{Float64}, grad::AbstractVector{Float64},
    Xty::Vector{Float64}, y_sq::Float64, inv_n::Float64,
)
    return 0.5 * (dot(w, grad) + inv_n * (y_sq - dot(w, Xty)))
end

"""
Solve one working-set subproblem with a |ws|-sized Gram matrix, updating
`w[ws]` and `Xw` in place. Memory is bounded by |ws|² regardless of
n_features, so this also serves very wide problems where the full X'X is
impossible. Assumes no intercept and an index-uniform penalty (true for
NNLSPenalty: prox/subdiff do not depend on the coordinate index).
"""
function _solve_subproblem_gram!(
    w::AbstractVector{T},
    Xw::AbstractVector{T},
    X::AbstractMatrix{T},
    y64::Vector{Float64},
    penalty::AbstractPenalty,
    options::AndersonCDOptions{T},
    ws::AbstractVector{Int},
    stop_crit::Float64,
) where T <: AbstractFloat
    n_samples = size(X, 1)
    n_ws = length(ws)
    inv_n = 1.0 / n_samples
    local_feats = 1:n_ws

    # ws-sized one-time products (BLAS syrk/gemv). Manual gather+convert is
    # several times faster than broadcasting over a vector-indexed view.
    X_ws = Matrix{Float64}(undef, n_samples, n_ws)
    @inbounds for (c, j) in enumerate(ws)
        @simd for i in 1:n_samples
            X_ws[i, c] = Float64(X[i, j])
        end
    end
    G = X_ws' * X_ws
    lipschitz = diag(G) .* inv_n

    w0 = Float64.(view(w, ws))
    grad0 = X_ws' * (Float64.(Xw) .- y64) .* inv_n
    w_loc = copy(w0)
    grad = copy(grad0)

    # Subproblem objective up to an additive constant:
    # 0.5‖y - Xw‖²/n = const + 0.5*(w_loc - w0)'(grad(w_loc) + grad0),
    # enough for the Anderson acceptance comparison.
    local_obj = (wl, gl) -> begin
        s = 0.0
        @inbounds @simd for i in 1:n_ws
            s += (wl[i] - w0[i]) * (gl[i] + grad0[i])
        end
        return 0.5 * s + value(penalty, wl)
    end

    accelerator = AndersonAcceleration{Float64}(5)

    for epoch in 1:options.max_epochs
        _cd_epoch_gram!(G, w_loc, grad, lipschitz, penalty, local_feats, inv_n)

        w_accel, grad_accel, is_extrap = extrapolate(accelerator, w_loc, grad)
        if is_extrap && local_obj(w_accel, grad_accel) < local_obj(w_loc, grad)
            copyto!(w_loc, w_accel)
            copyto!(grad, grad_accel)
        end

        if epoch % 10 == 0
            stop_crit_in = maximum(subdiff_distance(penalty, w_loc, grad, local_feats))

            if stop_crit_in < 0.3 * stop_crit || stop_crit_in <= Float64(options.tol)
                if options.verbose
                    println("   Epoch $epoch, inner stop crit: $(stop_crit_in), early exit.")
                end
                break
            end
        end
    end

    # Push the subproblem solution back into the global iterate.
    delta = w_loc .- w0
    if any(!=(0.0), delta)
        Xw .+= X_ws * delta
    end
    w[ws] .= w_loc
    return nothing
end

"""
Gram-matrix variant of the solver (no intercept). All inner-loop state is
n_features-sized and kept in Float64 so the maintained gradient does not drift
below the stopping tolerance even for Float32 inputs.
"""
function _solve_gram(
    X::AbstractMatrix{T},
    y::AbstractVector{T},
    datafit::AbstractDatafit,
    penalty::AbstractPenalty,
    options::AndersonCDOptions{T},
    w_init::Union{Nothing, AbstractVector{T}} = nothing,
) where T <: AbstractFloat
    n_samples, n_features = size(X)
    inv_n = 1.0 / n_samples

    # One-time n_samples-sized work; X'X dispatches to BLAS syrk and is the
    # part that actually uses BLAS.set_num_threads.
    X64 = convert(Matrix{Float64}, X)
    y64 = convert(Vector{Float64}, y)
    G = X64' * X64
    Xty = X64' * y64
    y_sq = dot(y64, y64)

    w = isnothing(w_init) ? zeros(Float64, n_features) : Float64.(w_init)
    grad = (G * w .- Xty) .* inv_n
    lipschitz = diag(G) .* inv_n

    pen = is_penalized(penalty, n_features)
    unpen = .!pen
    n_unpen = sum(unpen)
    all_feats = 1:n_features

    obj_out = T[]
    stop_crit = Inf
    tol = Float64(options.tol)

    for t in 1:options.max_iter
        # Optimality condition check (grad is maintained, no matvec needed)
        if options.ws_strategy == :subdiff
            opt = subdiff_distance(penalty, w, grad, all_feats)
        elseif options.ws_strategy == :fixpoint
            opt = dist_fix_point_cd(w, grad, lipschitz, penalty, all_feats)
        else
            error("Unsupported ws_strategy: $(options.ws_strategy)")
        end
        stop_crit = maximum(opt)

        if options.verbose
            @printf("Iter %d, Stop crit: %.2e\n", t, stop_crit)
        end

        if stop_crit <= tol
            break
        end

        # Build working set (ws)
        g_support = generalized_support(penalty, w)
        ws_size = max(
            min(options.p0 + n_unpen, n_features),
            min(2 * sum(g_support) - n_unpen, n_features)
        )

        # Force inclusion of unpenalized features and current support
        opt[unpen] .= Inf
        opt[g_support] .= Inf

        # Get indices of the `ws_size` largest values in `opt`
        ws = partialsortperm(opt, 1:ws_size, rev=true)

        if options.verbose
            println("Solving subproblem with $(length(ws)) features.")
        end

        accelerator = AndersonAcceleration{Float64}(5)

        # Inner loop (epochs on subproblem)
        for epoch in 1:options.max_epochs
            _cd_epoch_gram!(G, w, grad, lipschitz, penalty, ws, inv_n)

            # Anderson Acceleration on (w, grad); grad is affine in w, so the
            # extrapolated pair stays consistent.
            w_accel, grad_accel, is_extrap = extrapolate(accelerator, w, grad)

            if is_extrap
                p_obj = _gram_objective(w, grad, Xty, y_sq, inv_n) + value(penalty, w)
                p_obj_acc = _gram_objective(w_accel, grad_accel, Xty, y_sq, inv_n) +
                            value(penalty, w_accel)

                if p_obj_acc < p_obj
                    copyto!(w, w_accel)
                    copyto!(grad, grad_accel)
                end
            end

            # Inner loop stopping criterion
            if epoch % 10 == 0
                opt_ws = subdiff_distance(penalty, w, view(grad, ws), ws)
                stop_crit_in = maximum(opt_ws)

                if stop_crit_in < 0.3 * stop_crit || stop_crit_in <= tol
                    if options.verbose
                        println("   Epoch $epoch, inner stop crit: $(stop_crit_in), early exit.")
                    end
                    break
                end
            end
        end # end epoch loop

        current_obj = _gram_objective(w, grad, Xty, y_sq, inv_n) + value(penalty, w)
        push!(obj_out, T(current_obj))

    end # end iter loop

    return T.(w), obj_out, T(stop_crit)
end

"""
Original X-based solver, kept for problems the Gram path does not cover
(fit_intercept=true or very wide X).
"""
function _solve_xw(
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
    y64 = Float64.(y)   # used by the ws-Gram subproblem solver

    # Pre-computations
    initialize!(datafit, X, y)
    lipschitz = get_lipschitz(datafit, X)
    
    pen = is_penalized(penalty, n_features)
    unpen = .!pen
    n_unpen = sum(unpen)
    all_feats = 1:n_features

    obj_out = T[]
    stop_crit = T(Inf)

    # All-or-nothing: enable ws-Gram subproblems only if the budget covers the
    # largest working-set Gram that can occur (the dim cap). Mixing Float64
    # Gram subproblems with Float32 X-based ones mid-solve can leave the
    # X-based inner target 0.3·stop_crit below Float32 noise and grind epochs.
    use_ws_gram = !options.fit_intercept &&
                  _gram_bytes(n_samples, min(n_features, GRAM_DIM_CAP)) <=
                  options.gram_mem_gb * 1e9

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
        g_support = generalized_support(penalty, view(w, 1:n_features))
        g_support_size = sum(g_support)
        ws_size = max(
            min(options.p0 + n_unpen, n_features),
            min(2 * g_support_size - n_unpen, n_features)
        )

        # Keep the subproblem on the Float64 Gram path while the forced-include
        # set (support + unpenalized) fits the dim cap: candidates beyond the
        # cap simply wait for a later iteration. Dropping the whole subproblem
        # to the Float32 X-based loop instead risks an unreachable inner target
        # (0.3·stop_crit below Float32 noise) once Gram iterations have made
        # stop_crit small.
        if use_ws_gram && g_support_size + n_unpen <= GRAM_DIM_CAP
            ws_size = min(ws_size, GRAM_DIM_CAP)
        end

        # Force inclusion of unpenalized features and current support
        opt[unpen] .= T(Inf)
        opt[g_support] .= T(Inf)

        # Get indices of the `ws_size` largest values in `opt`
        ws = partialsortperm(opt, 1:ws_size, rev=true)

        if options.verbose
            println("Solving subproblem with $(length(ws)) features.")
        end

        if use_ws_gram && length(ws) <= GRAM_DIM_CAP
            # Subproblem on a |ws|-sized Gram matrix: O(|ws|) per coordinate
            # instead of O(n_samples), memory bounded by the budget even for
            # huge n_features.
            _solve_subproblem_gram!(w, Xw, X, y64, penalty, options, ws, Float64(stop_crit))
        else

        accelerator = AndersonAcceleration{T}(5)

        # Progress tracking: Float32 noise can floor the inner criterion above
        # its target; bail back to the outer loop (which re-checks the global
        # criterion) instead of grinding to max_epochs.
        best_stop_crit_in = T(Inf)
        checks_no_improve = 0

        # Inner loop (epochs on subproblem)
        for epoch in 1:options.max_epochs
            _cd_epoch!(X, y, view(w, 1:n_features), Xw, lipschitz, datafit, penalty, ws)

            if options.fit_intercept
                intercept_old = w[end]
                w[end] -= intercept_update_step(datafit, y, Xw)
                Xw .+= (w[end] - intercept_old)
            end

            # Anderson Acceleration on the full (w, Xw) pair
            w_accel, Xw_accel, is_extrap = extrapolate(accelerator, w, Xw)

            if is_extrap
                p_obj = value(datafit, y, Xw) + value(penalty, view(w, 1:n_features))
                p_obj_acc = value(datafit, y, Xw_accel) + value(penalty, view(w_accel, 1:n_features))

                if p_obj_acc < p_obj
                    w, Xw = w_accel, Xw_accel
                end
            end

            # Inner loop stopping criterion (also exit at global tol: pushing a
            # subproblem below tol is wasted work and the relative target can
            # sit under the Float32 noise floor)
            if epoch % 10 == 0
                grad_ws = construct_grad(X, Xw, datafit, ws)
                opt_ws = subdiff_distance(penalty, view(w, 1:n_features), grad_ws, ws)
                stop_crit_in = maximum(opt_ws)

                if stop_crit_in < 0.99 * best_stop_crit_in
                    best_stop_crit_in = stop_crit_in
                    checks_no_improve = 0
                else
                    checks_no_improve += 1
                end

                if stop_crit_in < 0.3 * stop_crit || stop_crit_in <= options.tol ||
                   checks_no_improve >= 10
                    if options.verbose
                        println("   Epoch $epoch, inner stop crit: $(stop_crit_in), early exit.")
                    end
                    break
                end
            end
        end # end epoch loop

        end # end ws-gram / X-based branch
        
        current_obj = value(datafit, y, Xw) + value(penalty, view(w, 1:n_features))
        push!(obj_out, current_obj)

    end # end iter loop

    # Return solution, objectives, and final stopping criterion
    return w, obj_out, stop_crit
end

"""
Number of physical cores (Linux sysfs; falls back to logical count). BLAS is
never faster on SMT siblings, so `nnls_solve` clamps its thread argument to
this — e.g. asking for 16 threads on an 8-core/16-thread machine measured ~35%
slower than 8.
"""
function _physical_cores()
    try
        cores = Set{Tuple{String, String}}()
        base = "/sys/devices/system/cpu"
        for d in readdir(base)
            occursin(r"^cpu\d+$", d) || continue
            pkg = read(joinpath(base, d, "topology", "physical_package_id"), String)
            core = read(joinpath(base, d, "topology", "core_id"), String)
            push!(cores, (pkg, core))
        end
        return max(length(cores), 1)
    catch
        return Sys.CPU_THREADS
    end
end

"""
    nnls_solve(X, y, thread=1, verbose=false; gram_mem_gb=1.0, w_init=nothing)

Solve min ‖y - Xw‖² s.t. w ≥ 0. `thread` sets BLAS threads (clamped to the
physical core count — SMT oversubscription only slows BLAS down).
`gram_mem_gb` (GB = 1e9 bytes) caps the extra memory the solver may spend on
Gram-path buffers (Float64 column copy + Gram matrix, ≈ 8·dim·(n_samples+dim)
bytes). The default 1.0 already saturates the compute-optimal Gram size
(GRAM_DIM_CAP) at n_samples ≈ 27k — raise it only when n_samples is much
larger, lower it (or ≤ 0 to disable Gram paths) in memory-tight environments.
Speed is capped by GRAM_DIM_CAP, not by budget.

`w_init` warm-starts from a previous solution (length n_features, ideally
≥ 0): the closer it is to the new optimum, the fewer working-set iterations
and epochs are needed. The matching model fit Xw is derived internally, so any
consistent-length vector is safe.
"""
function nnls_solve(X::AbstractMatrix{T}, y::AbstractVector{T}, thread::Int=1, verbose::Bool=false;
                    gram_mem_gb::Real=1.0,
                    w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    BLAS.set_num_threads(min(thread, _physical_cores()))
    datafit = NNLSDatafit{T}()
    penalty = NNLSPenalty()

    options = AndersonCDOptions{T}(
        fit_intercept=false,
        verbose=verbose,
        gram_mem_gb=gram_mem_gb
    )

    if w_init === nothing
        w_sol, _, stop_crit = solve(X, y, datafit, penalty, options)
    else
        length(w_init) == size(X, 2) || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected n_features = $(size(X, 2))"))
        w0 = convert(Vector{T}, w_init)
        Xw0 = X * w0   # derive the matching fit; never trust a caller-supplied pair
        w_sol, _, stop_crit = solve(X, y, datafit, penalty, options, w0, Xw0)
    end

    return w_sol
end

function get_transpose_matrix(hf_dst::HDF5.Dataset)
    dims = size(hf_dst)
    out_dims = reverse(dims)
    T = eltype(hf_dst)

    # Allocate the final output matrix
    A_out = Matrix{T}(undef, out_dims...)
    
    chunk_cols = 1000
    
    # Pre-allocate a single buffer matrix outside the loop
    buffer = Matrix{T}(undef, dims[1], chunk_cols)

    for i in 1:chunk_cols:dims[2]
        idx_end = min(i + chunk_cols - 1, dims[2])
        current_cols = idx_end - i + 1
        
        # Create a contiguous view of the buffer for the current chunk size
        buffer_view = view(buffer, :, 1:current_cols)

        # Perform zero-allocation direct read from HDF5 to the buffer
        copyto!(buffer_view, hf_dst, :, i:idx_end)

        # Write the transposed chunk into the pre-allocated output matrix
        @views A_out[i:idx_end, :] .= transpose(buffer_view)
    end

    return A_out
end

function load_nnls_array(filename::String)
    h5open(filename, "r") do file
        return get_transpose_matrix(file["A"]), read(file["B"]), read(file["B_depth_start"])
    end
end

function load_fail_array(filename::String)
    h5open(filename, "r") do file
        return get_transpose_matrix(file["A_fail"])
    end
end

function nnls_solve_partial(X::AbstractMatrix{T}, y::AbstractVector{T},
                            idx::Vector{Int}, verbose::Bool=false;
                            gram_mem_gb::Real=1.0,
                            w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    BLAS.set_num_threads(1)
    datafit = NNLSDatafit{T}()
    penalty = NNLSPenalty()

    options = AndersonCDOptions{T}(
        fit_intercept=false,
        verbose=verbose,
        gram_mem_gb=gram_mem_gb
    )

    Xv = view(X, :, idx .+ 1)
    if w_init === nothing
        w_sol, _, stop_crit = solve(Xv, y, datafit, penalty, options)
    else
        length(w_init) == length(idx) || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected length(idx) = $(length(idx))"))
        w0 = convert(Vector{T}, w_init)
        w_sol, _, stop_crit = solve(Xv, y, datafit, penalty, options, w0, Xv * w0)
    end
    predict_suc_B = Xv * w_sol

    return w_sol, predict_suc_B
end


using Base.Threads
# Note: `gram_mem_gb` applies per subproblem; concurrent tasks can use up to
# n_tasks × gram_mem_gb in total.
function run_nnls_map(X::AbstractMatrix{T}, y::AbstractVector{T},
                      thread_data::Vector{Vector{Int}}, verbose::Bool=false;
                      gram_mem_gb::Real=1.0,
                      w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    if w_init === nothing
        tasks = [
            @spawn nnls_solve_partial(X, y, idx, verbose; gram_mem_gb=gram_mem_gb)
            for idx in thread_data
        ]
    else
        length(w_init) == size(X, 2) || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected full feature count = $(size(X, 2))"))
        tasks = [
            @spawn nnls_solve_partial(
                X, y, idx, verbose;
                gram_mem_gb=gram_mem_gb,
                w_init=T.(view(w_init, idx .+ 1))
            )
            for idx in thread_data
        ]
    end
    return [fetch(t) for t in tasks]
end
