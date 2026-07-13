# SPDX-License-Identifier: GPL-3.0-or-later

# Eigenvalue search. Ports FUNCT / AcousticLayers / Solve / Solve1 / Solve2 /
# Bisection (kraken.f90, krakenc.f90), ZBRENTX (RootFinderBrent.f90) and
# ZSecantX / ZSecantCX (misc/RootFinderSecantMod.f90).

"""
    acoustic_layers(w, x, f, g, ipower) -> (f, g, ipower, modecount)

Shoot (towards the surface) through the acoustic layers (ports
`AcousticLayers`, kraken(c).f90). Generic over real (KRAKEN) / complex
(KRAKENC) `x`, `f`, `g`. `count_modes` enables the Sturm sign-change count of
the real path.
"""
function acoustic_layers(w::MeshWorkspace{T,X}, x, f, g, ipower::Int;
                         count_modes::Bool=false) where {T,X}
    modecount = 0
    w.first_acoustic == 0 && return f, g, ipower, modecount

    for medium in w.last_acoustic:-1:w.first_acoustic
        h2k2 = w.h[medium]^2 * x
        ii = w.loc[medium] + w.n[medium] + 1
        rho_medium = w.rho[w.loc[medium]+1]  # density at top of each medium

        p1 = -2 * g
        p2 = (w.b1[ii] - h2k2) * g - 2 * w.h[medium] * f * rho_medium
        p0 = zero(p2)

        # Shoot (towards surface) through a single medium
        for jj in w.loc[medium]+w.n[medium]:-1:w.loc[medium]+1
            p0 = p1
            p1 = p2
            p2 = (h2k2 - w.b1[jj]) * p1 - p0

            if count_modes && real(p0 * p1) <= 0
                modecount += 1
            end

            while abs(real(p2)) > ROOF          # scale if necessary
                p0 *= FLOOR_
                p1 *= FLOOR_
                p2 *= FLOOR_
                ipower -= IPOW_F
            end
        end

        # f = P'/rho and g = -P since f P + g P'/rho = 0
        f = -(p2 - p0) / (2 * w.h[medium]) / rho_medium
        g = -p1
    end
    f, g, ipower, modecount
end

"""
    funct(w, prob, x, roots, nroots; count_modes=false, deflate=true)
        -> (delta, ipower, modecount)

The dispersion relation: FUNCT(x) = 0 at the eigenvalues. Ports `FUNCT`
(krakenc.f90 for complex `x`, kraken.f90 for real `x`). `roots[1:nroots]` are
previously found roots for deflation (KRAKENC always deflates; KRAKEN only
when elastic layers are present — the caller controls via `deflate`).
"""
function funct(w::MeshWorkspace{T,X}, prob::Problem{T}, x,
               roots::AbstractVector, nroots::Int;
               count_modes::Bool=false, deflate::Bool=true) where {T,X}
    cnt = Ref(0)
    fc, gc, ipower = bc_impedance(w, prob, x, :bot, prob.hsbot, cnt)  # bottom
    if X <: Complex
        f, g = fc, gc
    else
        f, g = real(fc), real(gc)
    end
    f, g, ipower, modecount = acoustic_layers(w, x, f, g, ipower; count_modes)
    fbot, gbot, ipowerbot = bc_impedance(w, prob, x, :top, prob.hstop, cnt)  # top

    if X <: Complex
        delta = f * gbot - g * fbot
    else
        delta = real(f * gbot - g * fbot)
        g * delta > 0 && (modecount += 1)
    end
    ipower += ipowerbot
    modecount += cnt[]

    # Deflate previous roots
    if deflate && nroots > 0
        for j in 1:nroots
            delta = delta / (x - roots[j])
            while abs(real(delta)) < FLOOR_ && abs(delta) > 0
                delta *= ROOF
                ipower -= IPOW_R
            end
            while abs(real(delta)) > ROOF
                delta *= FLOOR_
                ipower -= IPOW_F
            end
        end
    end

    delta, ipower, modecount
end

"""
    secant(x2, tolerance, maxiteration, fun) -> (x2, errmsg)

Secant method with extended-range function values (ports `ZSecantX` /
`ZSecantCX`, misc/RootFinderSecantMod.f90). `fun(x) -> (f, ipower)`.
"""
function secant(x2::X, tolerance::Real, maxiteration::Int, fun) where {X}
    tolerance > 0 || return x2, "Non-positive tolerance specified"
    # ZSecantX uses 10·tolerance for the second point, ZSecantCX uses 100·
    x1 = x2 + (X <: Complex ? 100 : 10) * tolerance
    f1, ipower1 = fun(x1)
    for _ in 1:maxiteration
        x0 = x1
        f0 = f1
        ipower0 = ipower1
        x1 = x2
        f1, ipower1 = fun(x1)

        # block overflows by forcing the shift to be bounded
        cnum = f1 * (x1 - x0)
        cden = f1 - f0 * 10.0^(ipower0 - ipower1)
        shift = abs(cnum) >= abs(cden * x1) ? X(0.1 * tolerance) : cnum / cden
        x2 = x1 - shift
        if abs(x2 - x1) + abs(x2 - x0) < tolerance
            return x2, ""
        end
    end
    x2, "Failure to converge in RootFinderSecant"
end

"""
    zbrentx(a, b, t, fun) -> (x, errmsg)

Brent's method, extended-range version (ports `ZBRENTX`,
RootFinderBrent.f90). `fun(x) -> (f, ipower)` with f·10^ipower the true value.
"""
function zbrentx(a::T, b::T, t::T, fun) where {T<:Real}
    macheps = T(1.0e-16)
    fa, iexpa = fun(a)
    fb, iexpb = fun(b)
    if (fa > 0 && fb > 0) || (fa < 0 && fb < 0)
        return b, "Function sign is the same at the interval endpoints"
    end
    local c::T, fc::T, d::T, e::T, F1::T, F2::T
    local iexpc::Int
    restart = true
    while true
        if restart                        # label 2000: internal root
            c, fc, iexpc = a, fa, iexpa
            e = b - a
            d = e
            if iexpa < iexpb              # external root
                F1 = fc * T(10)^(iexpc - iexpb)
                F2 = fb
            else
                F1 = fc
                F2 = fb * T(10)^(iexpb - iexpc)
            end
            restart = false
        end
        # label 3000
        if abs(F1) < abs(F2)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
            iexpa, iexpb, iexpc = iexpb, iexpc, iexpb
        end
        tol = 2 * macheps * abs(b) + t
        m = (c - b) / 2
        (abs(m) > tol && fb != 0) || return b, ""

        # see if a bisection is forced
        if iexpa < iexpb
            F1 = fa * T(10)^(iexpa - iexpb)
            F2 = fb
        else
            F1 = fa
            F2 = fb * T(10)^(iexpb - iexpa)
        end
        if abs(e) < tol || abs(F1) <= abs(F2)
            e = m
            d = e
        else
            s = fb / fa * T(10)^(iexpb - iexpa)
            local p::T, q::T
            if a == c                     # linear interpolation
                p = 2 * m * s
                q = 1 - s
            else                          # inverse quadratic interpolation
                q = fa / fc * T(10)^(iexpa - iexpc)
                r = fb / fc * T(10)^(iexpb - iexpc)
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            end
            if p > 0
                q = -q
            else
                p = -p
            end
            s = e
            e = d
            if 2 * p < 3 * m * q - abs(tol * q) && p < abs(s * q / 2)
                d = p / q
            else
                e = m
                d = e
            end
        end
        a, fa, iexpa = b, fb, iexpb
        if abs(d) > tol
            b += d
        else
            b += m > 0 ? tol : -tol
        end
        fb, iexpb = fun(b)
        (fb > 0) == (fc > 0) && (restart = true)
    end
end

"""
    bisection(w, prob, xmin, xmax, M) -> (xl, xr)

Sturm-sequence bisection returning an isolating interval per eigenvalue in
[xmin, xmax] (ports `Bisection`, kraken.f90).
"""
function bisection(w::MeshWorkspace{T,T}, prob::Problem{T}, xmin::T, xmax::T,
                   M::Int) where {T}
    maxbisections = 50
    xl = fill(xmin, M + 1)
    xr = fill(xmax, M + 1)
    noroots = T[]

    _, _, nzer1 = funct(w, prob, xmax, noroots, 0; count_modes=true, deflate=false)
    M == 1 && return xl, xr

    for mode in 1:M-1
        if xl[mode] == xmin
            x2 = xr[mode]
            x1 = max(maximum(view(xl, mode+1:M)), xmin)
            for _ in 1:maxbisections
                x = x1 + (x2 - x1) / 2
                _, _, cnt = funct(w, prob, x, noroots, 0; count_modes=true,
                                  deflate=false)
                nzeros = cnt - nzer1
                if nzeros < mode          # not too many zeros: new right bdry
                    x2 = x
                    xr[mode] = x
                else                      # new left bdry
                    x1 = x
                    nzeros + 1 <= M + 1 && xr[nzeros+1] >= x && (xr[nzeros+1] = x)
                    nzeros >= 1 && nzeros <= M + 1 && xl[nzeros] <= x && (xl[nzeros] = x)
                end
                xl[mode] != xmin && break
            end
        end
    end
    xl, xr
end

"""
    solve1_real!(w, prob, evmat, iset) -> M

KRAKEN's bisection + Brent solver for the first two meshes (ports `Solve1`,
kraken.f90). Fills `evmat[iset, 1:M]` and returns M (also sizes `evmat` on
the first mesh via the returned count — caller allocates).
"""
function count_modes_real(w::MeshWorkspace{T,T}, prob::Problem{T}) where {T}
    ω² = (2π * prob.freq)^2
    xmin = T(1.00001) * ω² / w.chigh^2
    _, _, cnt = funct(w, prob, xmin, T[], 0; count_modes=true, deflate=false)
    cnt
end

function solve1_real!(w::MeshWorkspace{T,T}, prob::Problem{T},
                      evmat::Matrix{T}, iset::Int) where {T}
    ω² = (2π * prob.freq)^2
    xmin = T(1.00001) * ω² / w.chigh^2
    Mtotal = count_modes_real(w, prob)

    xmax = ω² / w.clow^2
    _, _, cnt = funct(w, prob, xmax, T[], 0; count_modes=true, deflate=false)
    M = Mtotal - cnt
    M == 0 && error("KrakenJL: No modes for given phase speed interval")
    M = min(M, size(evmat, 2))

    xl, xr = bisection(w, prob, xmin, xmax, M)

    # call ZBRENT to refine each eigenvalue in turn
    noroots = T[]
    for mode in 1:M
        x1 = xl[mode]
        x2 = xr[mode]
        eps = abs(x2) * T(10)^(2 - precision_digits(T))
        x, errmsg = zbrentx(x1, x2, eps, xx ->
            (dp = funct(w, prob, xx, noroots, 0; deflate=false);
             (dp[1], dp[2])))
        # errmsg is a warning only in the Fortran; the root is kept either way
        evmat[iset, mode] = x
    end
    M
end

# Fortran PRECISION(double) = 15
precision_digits(::Type{Float64}) = 15
precision_digits(::Type{T}) where {T<:Real} = floor(Int, -log10(eps(T))) - 1

"""
    solve2_real!(w, prob, evmat, hv, iset, M, deflate) -> M

KRAKEN's secant solver, used for meshes ≥ 3 (or throughout when elastic
layers are present). Ports `Solve2` (kraken.f90).
"""
function solve2_real!(w::MeshWorkspace{T,T}, prob::Problem{T},
                      evmat::Matrix{T}, hv::Vector{T}, iset::Int, M::Int,
                      deflate::Bool) where {T}
    ω² = (2π * prob.freq)^2
    maxiteration = 2000
    x = ω² / w.clow^2
    P = zeros(T, 10)

    for mode in 1:M
        x = T(1.00001) * x
        if iset >= 2
            for j in 1:iset-1
                P[j] = evmat[j, mode]
            end
            if iset >= 3
                for ii in 1:iset-2, j in 1:iset-ii-1
                    x1 = hv[j]^2
                    x2 = hv[j+ii]^2
                    P[j] = ((hv[iset]^2 - x2) * P[j] - (hv[iset]^2 - x1) * P[j+1]) /
                           (x1 - x2)
                end
                x = P[1]
            end
        end

        tolerance = abs(x) * length(w.b1) * T(10)^(1 - precision_digits(T))
        roots = view(evmat, iset, 1:mode-1)
        x, errmsg = secant(x, tolerance, maxiteration, xx ->
            (dp = funct(w, prob, xx, roots, mode - 1; deflate);
             (dp[1], dp[2])))
        isempty(errmsg) || (x = floatmin(T))    # make sure value is discarded

        evmat[iset, mode] = x

        if ω² / w.chigh^2 > x       # toss out modes outside the spectrum
            return mode - 1
        end
    end
    M
end

"""
    solve2_complex!(w, prob, evmat, hv, iset, M, robust, rng) -> M

KRAKENC's secant solver with deflation and random restarts (ports `Solve2`,
krakenc.f90). `rng` supplies the restart randomness — a deterministic
generator in this port (see PORTING_NOTES §3).
"""
function solve2_complex!(w::MeshWorkspace{T,Complex{T}}, prob::Problem{T},
                         evmat::Matrix{Complex{T}}, hv::Vector{T}, iset::Int,
                         M::Int, robust::Bool, rng) where {T}
    ω = 2π * prob.freq
    ω² = ω^2
    maxiteration = 1000
    x = Complex{T}(ω² / w.clow^2)
    P = zeros(Complex{T}, 10)

    maxtries = robust && iset <= 2 ? max(10, length(w.b1) ÷ 50) : 1
    itry = 1
    mode = 1
    while true
        x = T(1.00001) * x
        if iset >= 2
            for j in 1:iset-1
                P[j] = evmat[j, mode]
            end
            if iset >= 3
                for ii in 1:iset-2, j in 1:iset-ii-1
                    x1 = hv[j]^2
                    x2 = hv[j+ii]^2
                    P[j] = ((hv[iset]^2 - x2) * P[j] - (hv[iset]^2 - x1) * P[j+1]) /
                           (x1 - x2)
                end
                x = P[1]
            end
        end

        tolerance = abs(x) * length(w.b1) * T(10)^(1 - precision_digits(T))
        roots = view(evmat, iset, 1:mode-1)
        x, errmsg = secant(x, tolerance, maxiteration, xx ->
            (dp = funct(w, prob, xx, roots, mode - 1);
             (dp[1], dp[2])))

        if ω / w.chigh > real(sqrt(x)) || !isempty(errmsg)
            # mode outside the spectral limits: restart at a random point
            # (uniform distribution in vertical wavenumber, krakenc.f90)
            rvar1 = rand01(rng)
            rvar2 = rand01(rng)
            kzhigh = ω / w.clow
            kztry = rvar1 * kzhigh + im * T(0.01) * rvar2 * kzhigh
            x = ω² / w.clow^2 - kztry^2
            if itry < maxtries
                itry += 1
                continue
            else
                break
            end
        else
            evmat[iset, mode] = x
            mode += 1
            mode > M && break
        end
    end
    M = mode - 1
    M == 0 && error("KrakenJL: No modes for given phase speed interval")

    # order eigenvalues by real part (descending), Sort (misc/SortMod.f90)
    xs = evmat[iset, 1:M]
    sort!(xs; by=real, rev=true)
    evmat[iset, 1:M] .= xs
    M
end

# Deterministic substitute for the Fortran RANDOM_NUMBER in the (rare)
# restart path — a simple xorshift64*, seeded per solve for reproducibility.
mutable struct RestartRNG
    state::UInt64
end
RestartRNG() = RestartRNG(0x9E3779B97F4A7C15)
function rand01(rng::RestartRNG)
    x = rng.state
    x ⊻= x >> 12
    x ⊻= x << 25
    x ⊻= x >> 27
    rng.state = x
    Float64((x * 0x2545F4914F6CDD1D) >> 11) / Float64(1 << 53)
end

"""
    solve_modes(prob, ngs; complex_solver=true, robust=false, threads=1)
        -> ModeResult

Top-level driver: solve the eigenproblem over a sequence of up to `NSETS`
meshes with Richardson extrapolation (ports the main program + `Solve` of
kraken.f90 / krakenc.f90). `ngs` are the base mesh counts per medium.
"""
function solve_modes(prob::Problem{T}, ngs::Vector{Int};
                     complex_solver::Bool=true, robust::Bool=false,
                     threads::Int=1) where {T}
    X = complex_solver ? Complex{T} : T
    nmedia = length(prob.media)
    w = MeshWorkspace{T,X}(nmedia)
    ω = 2π * prob.freq
    ω² = ω^2
    hv = zeros(T, NSETS)
    rng = RestartRNG()

    evmat = Matrix{X}(undef, 0, 0)
    extrap = Matrix{X}(undef, 0, 0)
    kpert = Complex{T}[]
    vg = T[]
    phi = Matrix{Complex{T}}(undef, 0, 0)
    zmesh = T[]
    M = 0
    err = T(Inf)
    no_elastic_layers = true

    for iset in 1:NSETS
        for m in 1:nmedia
            w.n[m] = ngs[m] * NV[iset]
            w.h[m] = (prob.media[m].z1 - prob.media[m].z0) / w.n[m]
        end
        hv[iset] = w.h[1]

        initialize!(w, prob)
        nacoustic = w.last_acoustic - w.first_acoustic + 1
        no_elastic_layers = nmedia <= nacoustic

        if iset == 1
            # size the eigenvalue storage (kraken: mode count / 3000;
            # krakenc: MaxM)
            M0 = if complex_solver
                MAX_MODES
            elseif no_elastic_layers
                max(count_modes_real(w, prob), 1)
            else
                3000
            end
            evmat = zeros(X, NSETS, M0)
            extrap = zeros(X, NSETS, M0)
            kpert = zeros(Complex{T}, M0)
            vg = zeros(T, M0)
            M = M0
        end

        # Solve for the eigenvalues on this mesh
        if complex_solver
            M = solve2_complex!(w, prob, evmat, hv, iset, M, robust, rng)
        elseif iset <= 2 && no_elastic_layers
            M = solve1_real!(w, prob, evmat, iset)
        else
            M = solve2_real!(w, prob, evmat, hv, iset, M, !no_elastic_layers)
        end

        extrap[iset, 1:M] .= evmat[iset, 1:M]

        # Remove eigenvalues outside the spectral limits (they are sorted
        # descending by real part, so this is the MINLOC trim of the Fortran)
        lim = ω² / w.chigh^2
        Mtrim = 0
        for i in 1:M
            real(extrap[1, i]) > lim && (Mtrim = i)
        end
        M = Mtrim
        M == 0 && error("KrakenJL: No modes for given phase speed interval")

        if iset == 1     # compute the eigenvectors on the first mesh
            phi, zmesh = vector!(w, prob, evmat, kpert, vg, M, complex_solver,
                                 threads)
        end

        # Richardson extrapolation to improve the accuracy (Solve)
        err = T(1.0e10)
        key = 2 * M ÷ 3 + 1
        if iset > 1
            t1 = extrap[1, key]
            for j in iset-1:-1:1, mode in 1:M
                x1 = T(NV[j])^2
                x2 = T(NV[iset])^2
                f1 = extrap[j, mode]
                f2 = extrap[j+1, mode]
                extrap[j, mode] = f2 - (f1 - f2) / (x2 / x1 - 1)
            end
            t2 = extrap[1, key]
            err = abs(t2 - t1)
        end

        err * 1000 * prob.rmax < 1 && break   # convergence achieved
    end

    # discard modes with phase velocity above cHigh (final trim)
    lim = ω² / w.chigh^2
    Mtrim = 0
    for i in 1:M
        real(extrap[1, i]) > lim && (Mtrim = i)
    end
    M = Mtrim

    # k() contains scatter/attenuation perturbations; combine with the
    # extrapolated eigenvalues (kraken(c).f90 main)
    k = [sqrt(Complex{T}(extrap[1, i]) + kpert[i]) for i in 1:M]
    if complex_solver
        # zero out positive imaginary parts which would cause growth in range
        k = [imag(ki) > 0 ? Complex{T}(real(ki)) : ki for ki in k]
    end

    ModeResult{T}(k, phi[:, 1:M], zmesh, vg[1:M])
end
