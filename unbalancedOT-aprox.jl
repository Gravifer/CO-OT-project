module UnbalancedOptimalTransport

using LinearAlgebra: norm

export DiscreteMeasure,
    OT!, optimal_coupling!, sinkhorn_divergence!, unbalanced_sinkhorn!, cost_matrix

struct DiscreteMeasure{P,LP,S,T}
    density::P
    log_density::LP
    set::S
    dual_potential::Vector{T}
    cache::Vector{T}
end

function DiscreteMeasure(density::P, log_density::LP, set::S) where {P,LP,S}
    T = eltype(density)
    n = length(density)
    n == length(log_density) ||
    throw(ArgumentError("`density`, `log_density`, and `set` (if supplied) must have equal length"))
    set !== nothing &&
    length(set) != n &&
    throw(ArgumentError("`density`, `log_density`, and `set` must have equal length"))
    dual_potential = zeros(T, n)
    cache = similar(dual_potential)
    DiscreteMeasure{P,LP,S,T}(density, log_density, set, dual_potential, cache)
end

DiscreteMeasure(density, set = nothing) = DiscreteMeasure(density, log.(density), set)

Base.eltype(::DiscreteMeasure{P,LP,S,T}) where {P,LP,S,T} = T
Base.length(a::DiscreteMeasure) = length(a.density)

abstract type AbstractDivergence end
function φstar end
function aprox end


# The following formulas for `φstar` and `aprox` are found in Section 2.4 of SFVTP19.
struct Balanced <: AbstractDivergence end
φstar(::Balanced, q) = q
aprox(::Balanced, ϵ, x) = x

struct KL{ρ} <: AbstractDivergence end
KL(ρ::Number) = KL{ρ}()
KL() = KL(1.0)

φstar(::KL{ρ}, q) where {ρ} = ρ * (exp(q / ρ) - 1)
aprox(::KL{ρ}, ϵ, x) where {ρ} = inv(one(ρ) + ϵ / ρ) * x

struct RG{l,u} <: AbstractDivergence end

function RG(l::Number, u::Number)
    l <= u || throw(DomainError(u - l, "Need l <= u"))
    RG{l,u}()
end

φstar(::RG{a,b}, q) where {a,b} = max(a * q, b * q)

function aprox(::RG{a,b}, ϵ, x) where {a,b}
    if (s = x - ϵ * log(a)) < 0
        s
    elseif (t = x - ϵ * log(b)) > 0
        t
    else
        0
    end
end

struct TV{ρ} <: AbstractDivergence end
TV(ρ::Number) = TV{ρ}()
TV() = TV(1.0)

φstar(::TV{ρ}, q) where {ρ} =
    q <= ρ ? max(-ρ, q) : throw(DomainError(q, "Must have q <= ρ"))

function aprox(::TV{ρ}, ϵ, x) where {ρ}
    if -ρ <= x <= ρ
        x
    elseif x < -ρ
        -ρ
    else
        ρ
    end
end

function initialize_dual_potential!(::AbstractDivergence, a::DiscreteMeasure)
    a.dual_potential .= 0
end


function unbalanced_sinkhorn!(
    D::AbstractDivergence,
    C::AbstractMatrix,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    tol = 1e-5,
    max_iters = 10^5,
    warn::Bool = true,
)
    if D isa Balanced && warn && sum(a.density) ≉ sum(b.density)
        @warn "Should have `sum(a.density) ≈ sum(b.density)` for `D==Balanced()`"
    end

    initialize_dual_potential!(D, a)
    initialize_dual_potential!(D, b)

    f = a.dual_potential
    tmp_f = a.cache

    g = b.dual_potential
    tmp_g = b.cache

    max_residual = Inf
    iters = 0
    while iters < max_iters && max_residual > tol
        iters += 1
        max_residual = 0.0
        @inbounds for j in eachindex(g)
            for i in eachindex(a.log_density, f, tmp_f)
                tmp_f[i] = a.log_density[i] + (f[i] - C[i, j]) / ϵ
            end
            new_g = -ϵ * logsumexp!(tmp_f)
            new_g = -aprox(D, ϵ, -new_g)
            diff = abs(g[j] - new_g)
            if diff > max_residual
                max_residual = diff
            end
            g[j] = new_g
        end
        @inbounds for j in eachindex(f)
            for i in eachindex(b.log_density, g, tmp_g)
                tmp_g[i] = b.log_density[i] + (g[i] - C[j, i]) / ϵ
            end
            new_f = -ϵ * logsumexp!(tmp_g)
            new_f = -aprox(D, ϵ, -new_f)
            diff = abs(f[j] - new_f)
            if diff > max_residual
                max_residual = diff
            end
            f[j] = new_f
        end
    end

    if warn && iters == max_iters
        @warn "Maximum iterations ($max_iters) reached" max_residual
    end

    return (iters = iters, max_residual = max_residual)
end

function OT!(
    D::AbstractDivergence,
    C::AbstractMatrix,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    kwargs...,
)

    unbalanced_sinkhorn!(D, C, a, b, ϵ; kwargs...)

    f = a.dual_potential
    g = b.dual_potential

    T = promote_type(eltype(a), eltype(b))
    _nφstar = q -> -φstar(D, -q)

    t1 = fdot(_nφstar, a.density, f)
    t2 = fdot(_nφstar, b.density, g)
    t3 = zero(T)
    for i in eachindex(f), j in eachindex(g)
        t3 -= ϵ * a.density[i] * b.density[j] * (exp((f[i] + g[j] - C[i, j]) / ϵ) - one(T))
    end
    return t1 + t2 + t3
end

# Generic method
function generic_sinkhorn_divergence!(
    D::AbstractDivergence,
    C,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ;
    kwargs...,
)
    OT_aa = OT!(D, C, a, a, ϵ; kwargs...)
    OT_bb = OT!(D, C, b, b, ϵ; kwargs...)
    OT_ab = OT!(D, C, a, b, ϵ; kwargs...)

    OT_ab + (-OT_aa - OT_bb + ϵ * (sum(a.density) - sum(b.density))^2) / 2
end

sinkhorn_divergence!(
    D::AbstractDivergence,
    C,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    kwargs...,
) = generic_sinkhorn_divergence!(D, C, a, b, ϵ; kwargs...)

sinkhorn_divergence!(
    D::AbstractDivergence,
    C::AbstractMatrix,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    kwargs...,
) = throw(ArgumentError("Must pass a cost function `C`, not a cost matrix."))

function optimal_coupling!(
    D::AbstractDivergence,
    C::AbstractMatrix,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    dual_potentials_populated::Bool = false,
    kwargs...,
)
    if !dual_potentials_populated
        unbalanced_sinkhorn!(D, C, a, b, ϵ; kwargs...)
    end

    f = a.dual_potential
    g = b.dual_potential
    return [
        exp((f[i] + g[j] - C[i, j]) / ϵ) * a.density[i] * b.density[j]
        for i in eachindex(a.density), j in eachindex(b.density)
    ]
end


# Provide default choice of norm as cost function
for fun in (:unbalanced_sinkhorn!, :OT!, :sinkhorn_divergence!, :optimal_coupling!)
    @eval begin
        $fun(
            D::AbstractDivergence,
            a::DiscreteMeasure,
            b::DiscreteMeasure,
            args...;
            kwargs...,
        ) = $fun(D, (x, y) -> norm(x - y), a, b, args...; kwargs...)
    end
end

# If a function is provided, precompute the cost matrix. sinkhorn_divergence is omitted since it requires three different cost matrices.
for fun in (:unbalanced_sinkhorn!, :OT!, :optimal_coupling!)
    @eval begin
        $fun(
            D::AbstractDivergence,
            C, # Not restricting to Function to allow callable structs etc.
            a::DiscreteMeasure,
            b::DiscreteMeasure,
            args...;
            kwargs...,
        ) = $fun(D, cost_matrix(C, a, b), a, b, args...; kwargs...)
    end
end

# The formulas for `initialize_dual_potential!` are in Section 6.1.2 of SFVTP19.
function initialize_dual_potential!(::KL{ρ}, a::DiscreteMeasure) where {ρ}
    c = -ρ * log(length(a.log_density))
    a.dual_potential .= c
end

function initialize_dual_potential!(::TV{ρ}, a::DiscreteMeasure) where {ρ}
    c = -ρ * sign(log(length(a.log_density)))
    a.dual_potential .= c
end

# Specialized implementation for the KL-divergence
# Prop. 12 of SFVTP19.
function sinkhorn_divergence!(
    D::KL{ρ},
    C,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    kwargs...,
) where {ρ}
    f = a.dual_potential
    g = b.dual_potential
    scaled_exp = z -> exp(-z / ρ)
    unbalanced_sinkhorn!(D, C, a, a, ϵ; kwargs...)
    term_aa = fdot(scaled_exp, a.density, f)
    unbalanced_sinkhorn!(D, C, b, b, ϵ; kwargs...)
    term_bb = fdot(scaled_exp, b.density, g)
    unbalanced_sinkhorn!(D, C, a, b, ϵ; kwargs...)
    term_ab = fdot(scaled_exp, a.density, f) + fdot(scaled_exp, b.density, g)
    return -(ρ + ϵ / 2) * (term_ab - term_aa - term_bb)
end

# Needed to avoid ambiguity errors
sinkhorn_divergence!(
    D::KL{ρ},
    C::AbstractMatrix,
    a::DiscreteMeasure,
    b::DiscreteMeasure,
    ϵ = 1e-1;
    kwargs...,
) where {ρ} = throw(ArgumentError("Must pass a cost function `C`, not a cost matrix."))

################################
# Utilities
################################

# Numerically stable implementation of `w -> log(sum(exp, w))`.
# https://discourse.julialang.org/t/fast-logsumexp/22827/9
function logsumexp!(w)
    N = length(w)
    offset, maxind = findmax(w)
    w .= exp.(w .- offset)
    Σ = _sum_all_but(w, maxind)
    log1p(Σ) + offset
end

# Add all elements of vector `w` except for index `i`.
# The element at index `i` is assumed to have value 1
function _sum_all_but(w, i)
    w[i] -= 1
    s = sum(w)
    w[i] += 1
    s
end

"""
    fdot(f, u, v) -> Number

A generic, allocation-free implementation of `dot(u, f.(v))`. It may be faster
to provide a specialized method to dispatch to BLAS or so forth.
"""
function fdot(f, u, v)
    T = promote_type(eltype(u), eltype(v))
    s = zero(T)
    @inbounds for i in eachindex(u, v)
        s += conj(u[i]) * f(v[i])
    end
    s
end

"""
    cost_matrix([C,] a, b) -> Matrix

Precompute the cost matrix given a cost function `C`. If no function `C` is provided, the default is `C(x,y) = norm(x-y)`.
"""
cost_matrix(C, a, b) = [C(x, y) for x in a.set, y in b.set]
cost_matrix(a, b) = cost_matrix((x, y) -> norm(x - y), a, b)
end
