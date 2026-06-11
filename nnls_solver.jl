# nnls_solver.jl
#
# Julia port of the Adelie BVLS (bounded-variable least squares) solver,
# following the same drop-in API as anderson_nnls.jl (nnls_solve /
# nnls_solve_partial / run_nnls_map + HDF5 loaders).
#
# Ported from the locally cloned Adelie sources:
#   adelie/adelie/src/include/adelie_core/solver/solver_bvls.hpp
#       coordinate_descent / solve_active / fit / kkt_screen / solve
#   adelie/adelie/src/include/adelie_core/state/state_bvls.hpp(.ipp)
#       StateBVLS fields + input validation
#   adelie/adelie/solver.py :: bvls()
#       state construction and defaults (weights = 1/n, kappa = min(n, p),
#       max_iters = 1e5, tol = 1e-7, bounds clipped at ±max_solver_value,
#       beta0 = box vertex closest to the origin)
#   adelie/adelie/src/include/adelie_core/matrix/matrix_naive_dense.ipp
#       cmul / ctmul / mul  →  column dot / axpy! / Xᵀ·v
#
# The problem solved is
#     min_β 0.5‖y - Xβ‖²_W   s.t.  lower ≤ β ≤ upper,
# with uniform observation weights W = I/n (adelie's default when
# `weights=None`; the algorithm is what `adelie.solver.bvls(A, B, lower,
# upper)` runs in test.py). NNLS is the special case lower=0, upper=∞.
#
# Algorithm sketch (active-set coordinate descent with violation screening):
#   outer `solve!`: repeat { fit on screen set; stop if loss stagnates or all
#   KKT conditions hold; otherwise screen the top-kappa KKT violators }.
#   `fit!`: full CD passes over the screen set; coordinates that move are
#   added to the active set, which is solved to convergence (`solve_active!`)
#   and pruned of entries that hit their bounds.
#
# ── Performance layer (profiling-driven; see git history for the 1:1 port) ──
# Profiling on the real matrices (26693×514 / 26266×10091 Float32) shows
# >80% of the time in `solve_active` — hundreds of CD sweeps over the same
# few-hundred/thousand active columns, each visit streaming an n-length
# column from DRAM twice, with 66–81% of visits actually moving — ~13% in
# the screen-set sweeps, ~3% in the KKT gradient. Optimizations, all
# preserving the adelie update order (Gauss–Seidel trajectory):
#
#  1. Gram-mode `solve_active` (`_solve_active_gram!`): gradients of the
#     ever-active coordinates are maintained in p-space via an incrementally
#     built Gram cache (BLAS-3 syrk/gemm over a gathered column buffer), so a
#     CD visit costs O(1) and a move costs one cache-resident axpy over the
#     Gram column instead of two O(n) DRAM streams. The residual is frozen
#     during the inner solve and trued up once at exit (resid -= XA·Δβ).
#     Mathematically identical updates with different FP rounding; the outer
#     KKT pass always re-checks the exact residual gradient. `use_gram=false`
#     restores the original bitwise-adelie-identical matrix-free path, which
#     also remains the fallback when the `gram_mem_gb` budget is exceeded.
#  2. Threaded setup/build/KKT: X_vars, the Gram build, the gradient
#     refresh / residual true-up gemvs and the KKT gradient (one-pass gemv
#     with alpha=1/n) use `thread` BLAS threads / chunked Julia tasks.
#
# The screen-set sweeps stay sequential by design: measured mover density is
# ~9–50% uniformly across sweeps (each sweep follows a converged active
# solve and finds fresh movers), so speculative batched gradients — computed
# against a frozen residual and discarded past the first mover — were
# measured a wash against the ~2× DRAM-bandwidth thread ceiling of n-length
# column dots. Their cost is streaming physics, not implementation.

using LinearAlgebra
using Printf, HDF5
using Base.Threads

# adelie_core Configs::max_solver_value_def: bounds are clipped to this range
# (±1e100). For Float32 this saturates to ±Inf32, which the solver handles.
const BVLS_MAX_SOLVER_VALUE = 1e100

#=
--------------------------------------------------------------------------------
                            THREADING HELPERS
--------------------------------------------------------------------------------
`n_threads` keeps anderson_nnls.jl's user-facing meaning (BLAS threads for the
big matrix products) and additionally drives the chunked Julia-thread loops
below. Both clamp to what is actually available.
=#

"Run `f(i0, i1)` over 1:n in `nt` chunks (`nt=1` or no Julia threads ⇒ inline)."
function _tforeach(f::F, n::Int, nt::Int) where F
    n <= 0 && return nothing
    nc = min(nt, Threads.nthreads(), n)
    if nc <= 1
        f(1, n)
        return nothing
    end
    sz = cld(n, nc)
    @sync for c in 1:nc
        i0 = (c - 1) * sz + 1
        i1 = min(c * sz, n)
        i0 <= i1 && Threads.@spawn f($i0, $i1)
    end
    return nothing
end

"Call `f()` with BLAS temporarily at `nt` threads (CD kernels need 1)."
@inline function _with_blas(f::F, nt::Int) where F
    nt <= 1 && return f()
    BLAS.set_num_threads(nt)
    r = f()
    BLAS.set_num_threads(1)
    return r
end

#=
--------------------------------------------------------------------------------
                                STATE
--------------------------------------------------------------------------------
Port of StateBVLS (state_bvls.hpp) + the setup in solver.py::bvls.
Uniform weights enter only as the scalar `inv_n` (= 1/n) factored out of the
weighted products, so all kernels run on raw columns of X.
=#

#=
Incrementally built Gram cache over the coordinates that ever became active.
`XA[:, 1:m]` gathers their columns (gemm-able dense buffer), `G[1:m, 1:m]`
holds (XAᵀXA)/n. Coordinates are appended on activation and never removed
(pruned coordinates re-enter for free); `g`/`dbeta` are per-`_solve_active_gram!`
scratch: the maintained gradient over the Gram coordinates and the net β
change used to true up the residual at exit.
=#
mutable struct GramCache{T <: AbstractFloat}
    cap::Int              # current column capacity of XA/G
    cap_max::Int          # capacity limit from gram_mem_gb
    m::Int                # number of cached coordinates
    set::Vector{Int}      # gram coordinates in admission order (length p, 1:m used)
    pos::Vector{Int}      # coordinate -> gram position (0 = absent)
    XA::Matrix{T}         # n × cap gathered columns
    G::Matrix{T}          # cap × cap, (XAᵀXA)·inv_n, full symmetric
    g::Vector{T}          # cap, maintained gradient (valid inside solve_active)
    dbeta::Vector{T}      # cap, net Δβ since entry (residual true-up at exit)
end

function GramCache{T}(n::Int, p::Int, gram_mem_gb::Float64) where T
    # sizeof(T)·(n·c + c²) ≤ gram_mem_gb·1e9  (XA + G; g/dbeta are O(c))
    budget = gram_mem_gb * 1e9 / sizeof(T)
    cap_max = floor(Int, (-n + sqrt(Float64(n)^2 + 4 * budget)) / 2)
    cap_max = min(cap_max, p)
    cap = min(256, cap_max)
    return GramCache{T}(
        cap, cap_max, 0, Vector{Int}(undef, p), zeros(Int, p),
        Matrix{T}(undef, n, cap), Matrix{T}(undef, cap, cap),
        Vector{T}(undef, cap), zeros(T, cap),
    )
end

mutable struct BVLSState{T <: AbstractFloat, M <: AbstractMatrix{T}}
    X::M

    # static states
    y_var::T              # Σ y²ᵢ wᵢ = ‖y‖²/n
    X_vars::Vector{T}     # Σ X²ᵢⱼ wᵢ = ‖X[:,j]‖²/n
    lower::Vector{T}
    upper::Vector{T}
    inv_n::T              # uniform observation weight 1/n

    # configurations
    kappa::Int            # violation batching size
    max_iters::Int
    tol::T
    # outer loss-stagnation exit: |Δloss| < loss_stag_rtol·y_var. adelie
    # hardcodes the 1e-6 literal in solver_bvls.hpp; 0 disables the exit so
    # the solver runs until the KKT screen passes (tol then sets the final
    # gradient precision: |g_k| ≲ √(tol·y_var·X_vars[k])).
    loss_stag_rtol::Float64
    n_threads::Int        # BLAS threads (gemv/syrk/gemm) + chunked Julia loops
    gram::Union{Nothing, GramCache{T}}   # nothing ⇒ matrix-free solve_active

    # dynamic states
    screen_set_size::Int
    screen_set::Vector{Int}
    is_screen::Vector{Bool}
    active_set_size::Int
    active_set::Vector{Int}
    is_active::Vector{Bool}
    beta::Vector{T}
    resid::Vector{T}      # y - Xβ
    grad::Vector{T}       # Xᵀ(resid .* w); also reused as the KKT violations
    loss::T               # 0.5 Σ residᵢ² wᵢ
    iters::Int
    n_kkt::Int
end

"""
Build a `BVLSState`, mirroring the state construction in
`adelie/solver.py::bvls` (weights=None ⇒ uniform 1/n).

`w_init` plays the role of adelie's `warm_start`: `beta` is the clamped
`w_init` and the screen/active sets start as the strictly-interior
coordinates — exactly what a previous solve's pruned active set contains.
(Unlike adelie's Python warm path, the sets are independent copies, not
aliased buffers.)
"""
function BVLSState(
    X::AbstractMatrix{T},
    y::AbstractVector{T},
    lower::AbstractVector,
    upper::AbstractVector;
    kappa::Union{Nothing, Int} = nothing,
    max_iters::Int = 100_000,
    tol::Real = 1e-7,
    loss_stag_rtol::Real = 1e-6,
    n_threads::Int = 1,
    w_init::Union{Nothing, AbstractVector} = nothing,
    use_gram::Bool = true,
    gram_mem_gb::Real = 1.0,
) where T <: AbstractFloat
    n, p = size(X)
    length(y) == n || throw(ArgumentError("y must be (n,) where X is (n, p)."))
    length(lower) == p || throw(ArgumentError("lower must be (p,) where X is (n, p)."))
    length(upper) == p || throw(ArgumentError("upper must be (p,) where X is (n, p)."))
    kap = kappa === nothing ? min(n, p) : kappa
    kap > 0 || throw(ArgumentError("kappa must be > 0."))
    tol >= 0 || throw(ArgumentError("tol must be >= 0."))
    nt = max(n_threads, 1)

    # The setup products are accumulated in Float64 and then cast to T,
    # matching adelie's Python wrapper where weights=np.full(n, 1/n) is
    # float64 and the numpy products promote (X_vars, y_var, loss).
    inv_n = one(T) / T(n)
    inv_n64 = 1.0 / n
    y_var = T(_dot64(y, y) * inv_n64)

    # X.sq_mul(weights, X_vars): one full pass over X
    X_vars = Vector{T}(undef, p)
    _tforeach(p, nt) do j0, j1
        @inbounds for j in j0:j1
            col = view(X, :, j)
            X_vars[j] = T(_dot64(col, col) * inv_n64)
        end
    end

    lo = Vector{T}(T.(max.(lower, -BVLS_MAX_SOLVER_VALUE)))
    up = Vector{T}(T.(min.(upper, BVLS_MAX_SOLVER_VALUE)))

    screen_set = Vector{Int}(undef, p)
    is_screen = zeros(Bool, p)
    active_set = Vector{Int}(undef, p)
    is_active = zeros(Bool, p)
    screen_set_size = 0
    active_set_size = 0

    if w_init === nothing
        # vertex of the box closest to the origin (0 for NNLS)
        beta = ifelse.(abs.(lo) .< abs.(up), lo, up)
    else
        length(w_init) == p || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected n_features = $p"))
        beta = clamp.(Vector{T}(w_init), lo, up)
        for k in 1:p
            if lo[k] < beta[k] < up[k]
                active_set_size += 1
                active_set[active_set_size] = k
                is_active[k] = true
            end
        end
        screen_set_size = active_set_size
        copyto!(screen_set, 1, active_set, 1, active_set_size)
        copyto!(is_screen, is_active)
    end

    resid = Vector{T}(undef, n)
    if all(iszero, beta)
        copyto!(resid, y)
    else
        mul!(resid, X, beta)
        @. resid = y - resid
    end
    grad = Vector{T}(undef, p)
    loss = T(0.5 * _dot64(resid, resid) * inv_n64)

    gram = use_gram ? GramCache{T}(n, p, Float64(gram_mem_gb)) : nothing

    return BVLSState{T, typeof(X)}(
        X, y_var, X_vars, lo, up, inv_n,
        kap, max_iters, T(tol), Float64(loss_stag_rtol), nt, gram,
        screen_set_size, screen_set, is_screen,
        active_set_size, active_set, is_active,
        beta, resid, grad, loss, 0, 0,
    )
end

#=
--------------------------------------------------------------------------------
                            MATRIX KERNELS
--------------------------------------------------------------------------------
MatrixNaiveDense ops with the uniform weight folded into a scalar. Column
views are contiguous (also for whole-column views of a column-subset view),
so dot/axpy! dispatch to BLAS for Float32/Float64.
=#

# Float64-accumulated dot, used only for the one-time setup products where
# adelie's Python wrapper works in promoted float64 numpy arithmetic.
_dot64(x::AbstractVector{Float64}, y::AbstractVector{Float64}) = dot(x, y)
function _dot64(x::AbstractVector{T}, y::AbstractVector{T}) where T <: AbstractFloat
    s = 0.0
    @inbounds @simd for i in eachindex(x)
        s += Float64(x[i]) * Float64(y[i])
    end
    return s
end

# cmul(j, v, w) = X[:,j]ᵀ(v .* w) = (X[:,j]ᵀ v)/n
@inline _cmul(X::AbstractMatrix{T}, j::Int, v::Vector{T}, inv_n::T) where T =
    dot(view(X, :, j), v) * inv_n

# ctmul(j, c, out): out .+= c .* X[:,j]
@inline _ctmul!(out::Vector{T}, X::AbstractMatrix{T}, j::Int, c::T) where T =
    axpy!(c, view(X, :, j), out)

# mul(v, w, out): out = Xᵀ(v .* w) = (Xᵀv)/n — the KKT-pass gradient.
# One-pass BLAS gemv (alpha = 1/n) for plain (strided) matrices; chunked
# per-column BLAS dots otherwise (e.g. the column-subset views used by
# nnls_solve_partial).
function _mul!(out::Vector{T}, X::StridedMatrix{T}, v::Vector{T}, inv_n::T,
               nt::Int) where T <: Union{Float32, Float64}
    _with_blas(nt) do
        mul!(out, transpose(X), v, inv_n, zero(T))
    end
    return out
end
function _mul!(out::Vector{T}, X::AbstractMatrix{T}, v::Vector{T}, inv_n::T,
               nt::Int) where T <: AbstractFloat
    _tforeach(length(out), nt) do j0, j1
        @inbounds for j in j0:j1
            out[j] = dot(view(X, :, j), v) * inv_n
        end
    end
    return out
end

#=
--------------------------------------------------------------------------------
                            SOLVER CORE
--------------------------------------------------------------------------------
Port of solver_bvls.hpp. `early_exit` is omitted: adelie's user-facing solve
path passes a constant-false functor.
=#

#=
Shared per-coordinate CD update given the gradient `gk` (the body of adelie's
coordinate_descent loop minus the gradient computation). Writes `beta[k]` and
returns the step `del` (0 ⇒ no move). All sweep variants funnel through this
so the arithmetic is identical in every path.
=#
@inline function _cd_step!(beta::Vector{T}, k::Int, vk::T, lk::T, uk::T,
                           gk::T) where T <: AbstractFloat
    bk_old = beta[k]
    step = vk <= 0 ? zero(T) : gk / vk
    bk = clamp(bk_old + step, lk, uk)
    beta[k] = bk
    return bk - bk_old
end

# C++: `loss -= del * gk - 0.5 * scaled_del_sq` — the 0.5 literal makes the
# right-hand side a double expression even for value_t = float.
@inline _loss_update(loss::T, del::T, gk::T, scaled_del_sq::T) where T =
    T(Float64(loss) - (Float64(del * gk) - 0.5 * Float64(scaled_del_sq)))

"""
One coordinate-descent sweep over the indices `ks` (`bvls::coordinate_descent`).
Updates `beta`, `resid` and `loss` in place and returns the convergence
measure max_k X_vars[k]·Δβ_k². With `track_active=true`, coordinates that
move are appended to the active set (the `add_active` step used by `fit!`).
"""
function _coordinate_descent!(
    state::BVLSState{T},
    ks::AbstractVector{Int},
    track_active::Bool,
) where T <: AbstractFloat
    X = state.X
    X_vars = state.X_vars
    lower = state.lower
    upper = state.upper
    beta = state.beta
    resid = state.resid
    inv_n = state.inv_n

    convg_measure = zero(T)
    loss = state.loss
    @inbounds for k in ks
        vk = X_vars[k]
        gk = _cmul(X, k, resid, inv_n)
        del = _cd_step!(beta, k, vk, lower[k], upper[k], gk)
        del == 0 && continue
        scaled_del_sq = vk * del * del
        convg_measure = max(convg_measure, scaled_del_sq)
        loss = _loss_update(loss, del, gk, scaled_del_sq)
        _ctmul!(resid, X, k, -del)
        if track_active && !state.is_active[k]
            state.active_set_size += 1
            state.active_set[state.active_set_size] = k
            state.is_active[k] = true
        end
    end
    state.loss = loss
    return convg_measure
end

"""
Drop active coordinates that sit on a bound (the `prune` lambda in
`bvls::fit`); they stay screened and only re-enter the active set if a later
sweep moves them off the bound.
"""
function _prune!(state::BVLSState)
    n_active = 0
    @inbounds for i in 1:state.active_set_size
        k = state.active_set[i]
        bk = state.beta[k]
        if bk <= state.lower[k] || bk >= state.upper[k]
            state.is_active[k] = false
            continue
        end
        n_active += 1
        state.active_set[n_active] = k
    end
    state.active_set_size = n_active
    return nothing
end

"CD on the active set until convergence (`bvls::solve_active`), matrix-free."
function _solve_active!(state::BVLSState{T}) where T <: AbstractFloat
    tol_scaled = state.tol * state.y_var
    while true
        state.iters += 1
        ks = view(state.active_set, 1:state.active_set_size)
        convg_measure = _coordinate_descent!(state, ks, false)
        state.iters >= state.max_iters && error("bvls: max iterations reached!")
        convg_measure <= tol_scaled && break
    end
    return nothing
end

#=
--------------------------------------------------------------------------------
                        GRAM-MODE ACTIVE SOLVE
--------------------------------------------------------------------------------
=#

"Grow XA/G to at least `need` columns (≤ cap_max), geometric to amortize."
function _gram_grow!(gc::GramCache{T}, need::Int) where T
    new_cap = min(max(2 * gc.cap, need), gc.cap_max)
    n = size(gc.XA, 1)
    XA = Matrix{T}(undef, n, new_cap)
    copyto!(XA, 1, gc.XA, 1, n * gc.m)
    G = Matrix{T}(undef, new_cap, new_cap)
    @inbounds for j in 1:gc.m
        copyto!(view(G, 1:gc.m, j), view(gc.G, 1:gc.m, j))
    end
    g = Vector{T}(undef, new_cap)
    dbeta = zeros(T, new_cap)
    gc.XA, gc.G, gc.g, gc.dbeta, gc.cap = XA, G, g, dbeta, new_cap
    return nothing
end

#=
Admit every active coordinate not yet cached: gather its column into XA and
extend G by the cross block (gemm) and the new diagonal block (syrk), both
scaled by 1/n and mirrored to keep G fully symmetric (the per-move axpy reads
whole columns). Returns false iff the active set cannot fit in cap_max — the
caller then falls back to the matrix-free path.
=#
function _gram_admit!(gc::GramCache{T}, X::AbstractMatrix{T},
                      active_set::Vector{Int}, active_set_size::Int,
                      inv_n::T, nt::Int) where T
    q = 0
    @inbounds for i in 1:active_set_size
        k = active_set[i]
        if gc.pos[k] == 0
            q += 1
            gc.set[gc.m + q] = k          # stage; pos assigned after fit check
        end
    end
    q == 0 && return true
    m_new = gc.m + q
    m_new > gc.cap_max && return false
    m_new > gc.cap && _gram_grow!(gc, m_new)

    m = gc.m
    @inbounds for i in 1:q
        gc.pos[gc.set[m + i]] = m + i
    end
    # gather new columns (chunked over the new columns; each copy is one
    # contiguous column even when X is a column-subset view)
    XA = gc.XA
    _tforeach(q, nt) do i0, i1
        @inbounds for i in i0:i1
            copyto!(view(XA, :, m + i), view(X, :, gc.set[m + i]))
        end
    end

    G = gc.G
    Xnew = view(XA, :, (m + 1):m_new)
    _with_blas(nt) do
        if m > 0
            mul!(view(G, 1:m, (m + 1):m_new), transpose(view(XA, :, 1:m)),
                 Xnew, inv_n, zero(T))
        end
        BLAS.syrk!('U', 'T', inv_n, Xnew,
                   zero(T), view(G, (m + 1):m_new, (m + 1):m_new))
    end
    # mirror: cross block and the syrk upper triangle
    @inbounds for j in 1:m
        for i in 1:q
            G[m + i, j] = G[j, m + i]
        end
    end
    @inbounds for j2 in 1:q, j1 in 1:(j2 - 1)
        G[m + j2, m + j1] = G[m + j1, m + j2]
    end
    gc.m = m_new
    return true
end

#=
CD on the active set until convergence in Gram mode. The residual is frozen:
gradients over the Gram coordinates are refreshed once from it at entry
(gemv), maintained per move through G columns, and the residual is trued up
once at exit from the net Δβ. Same sweep order, update formula, tolerance
and iteration accounting as `_solve_active!`.
=#
function _solve_active_gram!(state::BVLSState{T}, gc::GramCache{T}) where T <: AbstractFloat
    X_vars = state.X_vars
    lower = state.lower
    upper = state.upper
    beta = state.beta
    active_set = state.active_set
    pos = gc.pos
    G = gc.G
    g = gc.g
    dbeta = gc.dbeta
    m = gc.m
    nt = state.n_threads
    tol_scaled = state.tol * state.y_var

    # refresh maintained gradients from the (true) residual: g = XAᵀ resid / n
    _with_blas(nt) do
        mul!(view(g, 1:m), transpose(view(gc.XA, :, 1:m)), state.resid,
             state.inv_n, zero(T))
    end

    loss = state.loss
    moved_any = false
    GC.@preserve G g begin
        gptr = pointer(g)
        while true
            state.iters += 1
            convg_measure = zero(T)
            @inbounds for i in 1:state.active_set_size
                k = active_set[i]
                pk = pos[k]
                vk = X_vars[k]
                gk = g[pk]
                del = _cd_step!(beta, k, vk, lower[k], upper[k], gk)
                del == 0 && continue
                scaled_del_sq = vk * del * del
                convg_measure = max(convg_measure, scaled_del_sq)
                loss = _loss_update(loss, del, gk, scaled_del_sq)
                dbeta[pk] += del
                moved_any = true
                # g[1:m] .-= del .* G[1:m, pk]  (cache-resident Gram column)
                BLAS.axpy!(m, -del, pointer(G, (pk - 1) * size(G, 1) + 1), 1,
                           gptr, 1)
            end
            state.iters >= state.max_iters && error("bvls: max iterations reached!")
            convg_measure <= tol_scaled && break
        end
    end
    state.loss = loss

    # true up the residual: resid -= XA[:, 1:m] * dbeta[1:m]
    if moved_any
        _with_blas(nt) do
            mul!(state.resid, view(gc.XA, :, 1:m), view(dbeta, 1:m),
                 -one(T), one(T))
        end
        fill!(view(dbeta, 1:m), zero(T))
    end
    return nothing
end

#=
Dispatch for the inner active-set solve: Gram mode whenever enabled and the
active set fits the memory budget; otherwise the original matrix-free path
(also the `use_gram=false` path, bitwise identical to the 1:1 port).
=#
function _solve_active_dispatch!(state::BVLSState{T}) where T <: AbstractFloat
    gc = state.gram
    if gc !== nothing
        if _gram_admit!(gc, state.X, state.active_set, state.active_set_size,
                        state.inv_n, state.n_threads)
            _solve_active_gram!(state, gc)
            return nothing
        end
        state.gram = nothing   # over budget: matrix-free from here on
    end
    _solve_active!(state)
    return nothing
end

"""
Solve to convergence on the current screen set (`bvls::fit`): alternate one
sweep over the screen set (activating coordinates that move) with full
convergence on the active set, pruning bound-hitters after each round.
"""
function _fit!(state::BVLSState{T}) where T <: AbstractFloat
    tol_scaled = state.tol * state.y_var
    while true
        state.iters += 1
        ks = view(state.screen_set, 1:state.screen_set_size)
        convg_measure = _coordinate_descent!(state, ks, true)
        state.iters >= state.max_iters && error("bvls: max iterations reached!")
        if convg_measure <= tol_scaled
            _prune!(state)
            break
        end
        _solve_active_dispatch!(state)
        _prune!(state)
    end
    return nothing
end

"""
KKT check + screening (`bvls::kkt_screen`): compute the full gradient,
turn it into per-coordinate KKT violations (in place, as in adelie), and
screen up to `kappa` of the largest violators. Returns `true` iff no
unscreened coordinate violates the KKT conditions.
"""
function _kkt_screen!(state::BVLSState{T}, viols_order::Vector{Int}) where T <: AbstractFloat
    beta = state.beta
    lower = state.lower
    upper = state.upper
    grad = state.grad

    state.n_kkt += 1

    # gradient via gemv — the one full-matrix op of the pass. CD sweeps are
    # level-1 BLAS where extra threads only add fork/join overhead (measured
    # slower), so BLAS stays at 1 thread and is raised just for this call.
    _mul!(grad, state.X, state.resid, state.inv_n, state.n_threads)

    # compute violations (viols aliases grad):
    # viol_k = max(g_k,0)·1[β_k<u_k] - min(g_k,0)·1[β_k>l_k]
    @inbounds for k in eachindex(grad)
        g = grad[k]
        grad[k] = max(g, zero(T)) * (beta[k] < upper[k]) -
                  min(g, zero(T)) * (beta[k] > lower[k])
    end

    # sort violations in decreasing order
    sortperm!(viols_order, grad; rev = true)

    screen_set_size_old = state.screen_set_size
    kkt_passed = true

    # check KKT and screen
    @inbounds for j in eachindex(viols_order)
        k = viols_order[j]
        vk = grad[k]
        (state.is_screen[k] || vk <= 0) && continue
        kkt_passed = false
        # break if reached max screen capacity for this pass
        state.screen_set_size >= screen_set_size_old + state.kappa && break
        # otherwise add violator to the screen set
        state.screen_set_size += 1
        state.screen_set[state.screen_set_size] = k
        state.is_screen[k] = true
    end

    return kkt_passed
end

"""
Main loop (`bvls::solve`): fit on the screen set, stop when the loss
stagnates (|Δloss| < loss_stag_rtol·y_var, only after the first KKT pass;
adelie hardcodes 1e-6) or the KKT conditions hold everywhere; otherwise
screen more violators and refit.
"""
function solve!(state::BVLSState{T}; verbose::Bool = false) where T <: AbstractFloat
    p = length(state.grad)
    viols_order = collect(1:p)

    while true
        loss_prev = state.loss
        _fit!(state)
        if verbose
            gm = state.gram === nothing ? -1 : state.gram.m
            @printf("[bvls] pass %d: iters=%d  screen=%d  active=%d  gram=%d  loss=%.6e\n",
                    state.n_kkt + 1, state.iters, state.screen_set_size,
                    state.active_set_size, gm, state.loss)
        end
        if state.n_kkt > 0 &&
           abs(state.loss - loss_prev) < state.loss_stag_rtol * abs(state.y_var)
            verbose && println("[bvls] converged: loss stagnation")
            return state
        end
        if _kkt_screen!(state, viols_order)
            verbose && println("[bvls] converged: KKT satisfied")
            return state
        end
    end
end

#=
--------------------------------------------------------------------------------
                            PUBLIC API
--------------------------------------------------------------------------------
=#

"""
    bvls(X, y, lower, upper; kappa=nothing, max_iters=100_000, tol=1e-7,
         loss_stag_rtol=1e-6, n_threads=1, w_init=nothing, use_gram=true,
         gram_mem_gb=1.0, verbose=false) -> BVLSState

Julia counterpart of `adelie.solver.bvls(X, y, lower, upper)` with uniform
observation weights 1/n. Returns the solved state; the solution is
`state.beta`. Defaults match adelie (`kappa = min(n, p)`,
`loss_stag_rtol = 1e-6`); `loss_stag_rtol = 0` disables the early
loss-stagnation exit so the solve runs to the KKT screen (tighter, slower,
no longer adelie-identical).

`use_gram=true` runs the inner active-set solve in Gram mode (same update
order, FP rounding differs from adelie at machine precision; the KKT pass
re-checks the exact residual gradient). `use_gram=false` reproduces the 1:1
matrix-free port bit for bit. `gram_mem_gb` caps the Gram buffers (GB = 1e9
bytes, anderson_nnls convention); past the cap the solver falls back to
matrix-free. `n_threads` drives BLAS (Gram build, KKT gemv, residual
true-up) and the chunked Julia-thread loops (column gather, X_vars setup).
"""
function bvls(
    X::AbstractMatrix{T},
    y::AbstractVector{T},
    lower::AbstractVector,
    upper::AbstractVector;
    kappa::Union{Nothing, Int} = nothing,
    max_iters::Int = 100_000,
    tol::Real = 1e-7,
    loss_stag_rtol::Real = 1e-6,
    n_threads::Int = 1,
    w_init::Union{Nothing, AbstractVector} = nothing,
    use_gram::Bool = true,
    gram_mem_gb::Real = 1.0,
    verbose::Bool = false,
) where T <: AbstractFloat
    state = BVLSState(X, y, lower, upper;
                      kappa = kappa, max_iters = max_iters, tol = tol,
                      loss_stag_rtol = loss_stag_rtol,
                      n_threads = n_threads, w_init = w_init,
                      use_gram = use_gram, gram_mem_gb = gram_mem_gb)
    return solve!(state; verbose = verbose)
end

"""
Number of physical cores (Linux sysfs; falls back to logical count). BLAS is
never faster on SMT siblings, so `nnls_solve` clamps its thread argument to
this.
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
    nnls_solve(X, y, thread=1, verbose=false; max_iters=100_000, tol=1e-7,
               use_gram=true, gram_mem_gb=1.0, w_init=nothing)

Solve min ‖y - Xw‖² s.t. w ≥ 0 with the Adelie BVLS solver (lower=0,
upper=∞). Same call shape as `anderson_nnls.jl`'s `nnls_solve`. `thread`
(clamped to the physical core count) drives the Gram build / KKT gradient /
residual true-up BLAS calls and the chunked Julia-thread loops; the
sequential CD kernels always run on single-threaded BLAS.

`w_init` warm-starts from a previous solution (length n_features): its
positive coordinates seed the screen/active sets, mirroring adelie's
`warm_start`.
"""
function nnls_solve(X::AbstractMatrix{T}, y::AbstractVector{T}, thread::Int=1, verbose::Bool=false;
                    max_iters::Int=100_000, tol::Real=1e-7, loss_stag_rtol::Real=1e-6,
                    use_gram::Bool=true, gram_mem_gb::Real=1.0,
                    w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    nt = min(thread, _physical_cores())
    BLAS.set_num_threads(1)   # CD sweeps must run on single-threaded BLAS
    p = size(X, 2)
    state = bvls(X, y, zeros(T, p), fill(T(Inf), p);
                 max_iters = max_iters, tol = tol, loss_stag_rtol = loss_stag_rtol,
                 n_threads = nt, w_init = w_init,
                 use_gram = use_gram, gram_mem_gb = gram_mem_gb, verbose = verbose)
    # leave BLAS at `thread` like anderson_nnls.jl (downstream X*w predicts)
    BLAS.set_num_threads(nt)
    return state.beta
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
                            max_iters::Int=100_000, tol::Real=1e-7,
                            use_gram::Bool=true, gram_mem_gb::Real=1.0,
                            w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    BLAS.set_num_threads(1)
    Xv = view(X, :, idx .+ 1)
    p = length(idx)

    if w_init !== nothing
        length(w_init) == p || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected length(idx) = $p"))
    end
    state = bvls(Xv, y, zeros(T, p), fill(T(Inf), p);
                 max_iters = max_iters, tol = tol, w_init = w_init,
                 use_gram = use_gram, gram_mem_gb = gram_mem_gb,
                 verbose = verbose)
    w_sol = state.beta
    predict_suc_B = Xv * w_sol

    return w_sol, predict_suc_B
end

# Note: `gram_mem_gb` applies per subproblem; concurrent tasks can use up to
# n_tasks × gram_mem_gb in total.
# `w_init`: a single FULL-SIZE warm start over all n_features (e.g. a previous
# full solve); each group warm-starts from its own slice w_init[idx] (idx
# 0-based, as in nnls_solve_partial).
function run_nnls_map(X::AbstractMatrix{T}, y::AbstractVector{T},
                      thread_data::Vector{Vector{Int}}, verbose::Bool=false;
                      max_iters::Int=100_000, tol::Real=1e-7,
                      use_gram::Bool=true, gram_mem_gb::Real=1.0,
                      w_init::Union{Nothing, AbstractVector}=nothing) where T <: AbstractFloat
    if w_init !== nothing
        length(w_init) == size(X, 2) || throw(ArgumentError(
            "w_init has length $(length(w_init)), expected n_features = $(size(X, 2)) " *
            "(full-size; each group takes its w_init[idx] slice)"))
    end
    tasks = [@spawn nnls_solve_partial(X, y, idx, verbose;
                                       max_iters=max_iters, tol=tol,
                                       use_gram=use_gram, gram_mem_gb=gram_mem_gb,
                                       w_init=w_init === nothing ? nothing : w_init[idx .+ 1])
             for idx in thread_data]
    return [fetch(t) for t in tasks]
end
