# SPDX-License-Identifier: GPL-3.0-or-later

# ForwardDiff support. The eigenvalue search (secant / Brent / bisection with
# deflation and restarts) is iterative with value-dependent control flow, so
# it is not run in dual arithmetic. Instead the search runs entirely on
# ForwardDiff *values*, and each converged eigenvalue is made dual-exact with
# a single Newton step of the (undeflated) characteristic function evaluated
# in full dual arithmetic — by the implicit function theorem this transfers
# the exact partials ∂x/∂θ = -(∂Δ/∂θ)/(∂Δ/∂x) regardless of how many
# iterations the value search took. Eigenvectors, normalization,
# perturbations, Richardson extrapolation and the mode summation are smooth
# arithmetic and simply re-run with dual numbers. See PORTING_NOTES §3.

import ForwardDiff

_value(x::Real) = x
_value(d::ForwardDiff.Dual) = _value(ForwardDiff.value(d))
_value(z::Complex) = complex(_value(real(z)), _value(imag(z)))
_valuetype(::Type{T}) where {T<:Real} = typeof(_value(one(T)))

_value_halfspace(hs::Halfspace) =
    Halfspace{Float64}(hs.bc, _value(hs.cp), _value(hs.cs), _value(hs.rho))

_value_medium(m::Medium) =
    Medium{Float64}(_value(m.z0), _value(m.z1), _value.(m.zn), _value.(m.cpn),
                    _value.(m.csn), _value.(m.rhon), m.spline, m.ng)

"Strip ForwardDiff partials from a `Problem` for the value-space search."
_value_problem(prob::Problem) =
    Problem{Float64}(_value(prob.freq), _value_medium.(prob.media),
                     _value.(prob.sigma), _value_halfspace(prob.hstop),
                     _value_halfspace(prob.hsbot), _value(prob.clow),
                     _value(prob.chigh), _value(prob.rmax))

"""
    dual_newton(w, prob, x0) -> x

One Newton step of the undeflated characteristic function in dual
arithmetic, starting from the value-converged root `x0`. ∂Δ/∂x is computed
by the same central-difference device `Normalize` (kraken(c).f90) uses for
the admittance derivatives.
"""
function dual_newton(w::MeshWorkspace{T,X}, prob::Problem{T}, x0::X) where {T,X}
    δ = 1.0e-7
    noroots = X[]
    d0, ip0, _ = funct(w, prob, x0, noroots, 0; deflate=false)
    x1 = (1 - δ) * x0
    x2 = (1 + δ) * x0
    d1, ip1, _ = funct(w, prob, x1, noroots, 0; deflate=false)
    d2, ip2, _ = funct(w, prob, x2, noroots, 0; deflate=false)
    f1 = d1 * 10.0^(ip1 - ip0)      # bring to a common power of ten
    f2 = d2 * 10.0^(ip2 - ip0)
    dfdx = (f2 - f1) / (x2 - x1)
    x0 - d0 / dfdx
end

function solve_modes(prob::Problem{T}, ngs::Vector{Int};
                     complex_solver::Bool=true, robust::Bool=false,
                     threads::Int=1) where {T<:ForwardDiff.Dual}
    # 1. run the eigenvalue search on values only
    core = _solve_modes_core(_value_problem(prob), ngs;
                             complex_solver, robust, threads)

    X = complex_solver ? Complex{T} : T
    nmedia = length(prob.media)
    w = MeshWorkspace{T,X}(nmedia)
    Mmax = size(core.evmat, 2)
    evmat = zeros(X, NSETS, Mmax)
    extrap = zeros(X, NSETS, Mmax)

    # 2. dual-correct the eigenvalues on each mesh used, replaying the
    # Richardson extrapolation of Solve (kraken(c).f90) in dual arithmetic
    for iset in 1:core.nsets_used
        for m in 1:nmedia
            w.n[m] = ngs[m] * NV[iset]
            w.h[m] = (prob.media[m].z1 - prob.media[m].z0) / w.n[m]
        end
        initialize!(w, prob)
        M = core.ms[iset]
        if threads <= 1 || M < 2 * threads
            for mode in 1:M
                evmat[iset, mode] = dual_newton(w, prob, X(core.evmat[iset, mode]))
            end
        else       # funct only reads the workspace, so modes can be chunked
            tasks = map(_chunks(M, threads)) do idxs
                Threads.@spawn for mode in idxs
                    evmat[iset, mode] = dual_newton(w, prob, X(core.evmat[iset, mode]))
                end
            end
            foreach(wait, tasks)
        end
        extrap[iset, 1:M] .= evmat[iset, 1:M]
        if iset > 1
            for j in iset-1:-1:1, mode in 1:M
                x1 = T(NV[j])^2
                x2 = T(NV[iset])^2
                f1 = extrap[j, mode]
                f2 = extrap[j+1, mode]
                extrap[j, mode] = f2 - (f1 - f2) / (x2 / x1 - 1)
            end
        end
    end

    # 3. eigenvectors / normalization / perturbations on the first mesh,
    # entirely in dual arithmetic
    for m in 1:nmedia
        w.n[m] = ngs[m] * NV[1]
        w.h[m] = (prob.media[m].z1 - prob.media[m].z0) / w.n[m]
    end
    initialize!(w, prob)
    M = core.M
    kpert = zeros(Complex{T}, Mmax)
    vg = zeros(T, Mmax)
    phi, zmesh = vector!(w, prob, evmat, kpert, vg, M, complex_solver, threads)

    k = [sqrt(Complex{T}(extrap[1, i]) + kpert[i]) for i in 1:M]
    if complex_solver
        k = [_value(imag(ki)) > 0 ? Complex{T}(real(ki)) : ki for ki in k]
    end
    ModeResult{T}(k, phi[:, 1:M], zmesh, vg[1:M])
end
