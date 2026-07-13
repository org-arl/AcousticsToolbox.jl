# SPDX-License-Identifier: GPL-3.0-or-later

# Eigenvectors, normalization and perturbations. Ports Vector / Normalize /
# ScatterLoss (kraken.f90, krakenc.f90), InverseIterationMod.f90 and
# KupIng (Kraken/Scattering.f90).

"""
    inverse_iteration!(phi, d, e) -> iError

Inverse iteration for the eigenvector of a singular tridiagonal symmetric
matrix (ports `InverseIterationD`/`InverseIterationZ`,
InverseIterationMod.f90 — a modified EISPACK TINVIT). `d` is the diagonal,
`e[2:n]` the sub-diagonal (`e[1]`, `e[n+1]` arbitrary). Generic over
real/complex element type.
"""
function inverse_iteration!(phi::Vector{X}, d::Vector{X}, e::Vector{X}) where {X}
    n = length(d)
    maxiteration = 3
    rv1 = Vector{X}(undef, n)
    rv2 = Vector{X}(undef, n)
    rv3 = Vector{X}(undef, n)
    rv4 = Vector{X}(undef, n)

    # infinity norm of the matrix (L1 of real+imag parts in the Fortran)
    norm_ = sum(j -> abs(real(d[j])) + abs(imag(d[j])), 1:n) +
            sum(j -> abs(real(e[j])) + abs(imag(e[j])), 2:n)

    # factor of 100: some Scholte modes weren't getting amplified enough (mbp 8/2010)
    eps3 = 100 * eps(_valuetype(typeof(norm_))) * norm_
    uk = real(n)
    eps4 = uk * eps3
    uk = eps4 / sqrt(uk)

    # elimination with interchanges
    xu = one(X)
    u = d[1]
    v = e[2]
    for i in 2:n
        if abs(e[i]) >= abs(u)
            xu = u / e[i]
            rv4[i] = xu
            rv1[i-1] = e[i]
            rv2[i-1] = d[i]
            rv3[i-1] = e[i+1]
            u = v - xu * rv2[i-1]
            v = -xu * rv3[i-1]
        else
            xu = e[i] / u
            rv4[i] = xu
            rv1[i-1] = u
            rv2[i-1] = v
            rv3[i-1] = zero(X)
            u = d[i] - xu * v
            v = e[i+1]
        end
    end
    u == 0 && (u = X(eps3))
    rv3[n-1] = zero(X)
    rv1[n] = u
    rv2[n] = zero(X)
    rv3[n] = zero(X)

    # main loop of inverse iteration
    fill!(phi, X(uk))
    for _ in 1:maxiteration
        # back substitution
        for i in n:-1:1
            phi[i] = (phi[i] - u * rv2[i] - v * rv3[i]) / rv1[i]
            v = u
            u = phi[i]
        end
        # norm of vector; test for doneness
        norm_ = sum(x -> abs(real(x)) + abs(imag(x)), phi)
        norm_ >= 1 && return 0

        xu = X(eps4 / norm_)     # scale the vector down
        phi .*= xu
        # forward elimination
        for i in 2:n
            u = phi[i]
            if rv1[i-1] == e[i]  # rows switched during triangularization
                u = phi[i-1]
                phi[i-1] = phi[i]
            end
            phi[i] = u - rv4[i] * phi[i-1]
        end
    end
    -1   # failure to converge
end

"""
Kuperman–Ingenito interfacial-roughness perturbation (ports `KupIng`,
Kraken/Scattering.f90). `P` = pressure at the interface, `U` = P'/ρ.
"""
function kuping(sigma::T, eta1sq::Complex{T}, rho1::T, eta2sq::Complex{T},
                rho2::T, P::Complex{T}, U::Complex{T}) where {T}
    sigma == 0 && return zero(Complex{T})
    scatter_root(z) = real(z) >= 0 ? sqrt(z) : -im * sqrt(-z)
    eta1 = scatter_root(eta1sq)
    eta2 = scatter_root(eta2sq)
    del = rho1 * eta2 + rho2 * eta1
    del == 0 && return zero(Complex{T})
    a11 = (eta1sq - eta2sq) / 2 - (rho2 * eta1sq - rho1 * eta2sq) * (eta1 + eta2) / del
    a12 = im * (rho2 - rho1)^2 * eta1 * eta2 / del
    a21 = -im * (rho2 * eta1sq - rho1 * eta2sq)^2 / (rho1 * rho2 * del)
    a22 = (eta1sq - eta2sq) / 2 + (rho2 - rho1) * eta1 * eta2 * (eta1 + eta2) / del
    -sigma^2 * (-a21 * P^2 + (a11 - a22) * P * U + a12 * U^2)
end

"""
    scatter_loss(w, prob, phi, x) -> Complex

Interfacial scatter-loss perturbation to k² (ports `ScatterLoss`,
kraken(c).f90). `phi` is the normalized eigenvector on the acoustic mesh,
`x` the eigenvalue (real for KRAKEN, complex for KRAKENC).
"""
function scatter_loss(w::MeshWorkspace{T,X}, prob::Problem{T},
                      phi::AbstractVector{Complex{T}}, x) where {T,X}
    ω² = (2π * prob.freq)^2
    pert = zero(Complex{T})
    j = 1
    L = w.loc[w.first_acoustic]
    nmedia = length(prob.media)

    for medium in w.first_acoustic-1:w.last_acoustic
        local rho1::T, eta1sq::Complex{T}, U::Complex{T}
        if medium == w.first_acoustic - 1    # top properties
            if prob.hstop.bc === :acoustoelastic
                rho1 = prob.hstop.rho
                eta1sq = x - ω² / prob.hstop.cp^2
                root1 = X <: Complex ? pekeris_root(eta1sq) : sqrt(eta1sq)
                U = root1 * phi[1] / prob.hstop.rho
            elseif prob.hstop.bc === :vacuum
                rho1 = T(1.0e-9)
                eta1sq = one(Complex{T})
                rho_inside = w.rho[w.loc[w.first_acoustic]+1]
                U = phi[2] / w.h[w.first_acoustic] / rho_inside
            else                             # rigid
                rho1 = T(1.0e9)
                eta1sq = one(Complex{T})
                U = zero(Complex{T})
            end
        else
            h2 = w.h[medium]^2
            j += w.n[medium]
            L = w.loc[medium] + w.n[medium] + 1
            rho1 = w.rho[L]
            eta1sq = (2 + w.b1[L]) / h2 - x
            U = (-phi[j-1] - (w.b1[L] - h2 * x) * phi[j] / 2) / (w.h[medium] * rho1)
        end

        local rho2::T, eta2sq::Complex{T}
        if medium == w.last_acoustic         # bottom properties
            if prob.hsbot.bc === :acoustoelastic
                rho2 = prob.hsbot.rho
                eta2sq = ω² / prob.hsbot.cp^2 - x
            elseif prob.hsbot.bc === :vacuum
                rho2 = T(1.0e-9)
                eta2sq = one(Complex{T})
            else                             # rigid
                rho2 = T(1.0e9)
                eta2sq = one(Complex{T})
            end
        else
            rho2 = w.rho[L+1]
            eta2sq = (2 + w.b1[L+1]) / w.h[medium+1]^2 - x
        end

        pert += kuping(prob.sigma[medium+1], eta1sq, rho1, eta2sq, rho2,
                       phi[j], U)
    end
    pert
end

"""
    normalize!(phi, w, prob, x, iturningpoint, complex_solver)
        -> (vg, pert)

Normalize the eigenvector and compute the group speed and (real path) the
attenuation perturbation (ports `Normalize` of krakenc.f90 for the complex
path and kraken.f90 for the real path). Mutates `phi` in place. Returns the
group speed (0 for the real path — KRAKEN leaves VG unset) and the
perturbation to k² *excluding* interfacial scatter.
"""
function normalize!(phi::AbstractVector{Complex{T}}, w::MeshWorkspace{T,X},
                    prob::Problem{T}, x, iturningpoint::Int) where {T,X}
    ω = 2π * prob.freq
    ω² = ω^2
    ntotal1 = length(phi)
    sqnorm = zero(Complex{T})
    slow = zero(Complex{T})
    pert = zero(Complex{T})

    # contribution from the top halfspace
    if prob.hstop.bc === :acoustoelastic
        if X <: Complex
            slow += phi[1]^2 / (2 * sqrt(x - ω² / prob.hstop.cp^2)) /
                    (prob.hstop.rho * prob.hstop.cp^2)
        else
            del = im * imag(sqrt(complex(x - ω² / prob.hstop.cp^2)))
            pert -= del * phi[1]^2 / prob.hstop.rho
            slow += phi[1]^2 / (2 * sqrt(complex(x - real(ω² / prob.hstop.cp^2)))) /
                    (prob.hstop.rho * real(prob.hstop.cp)^2)
        end
    end

    # contribution from the volume
    L = w.loc[w.first_acoustic]
    j = 1
    for medium in w.first_acoustic:w.last_acoustic
        L += 1
        rho_medium = w.rho[L]
        rho_omega_h2 = rho_medium * ω² * w.h[medium]^2

        # top interface
        sqnorm += w.h[medium] * phi[j]^2 / rho_medium / 2
        slow += w.h[medium] * (w.b1[L] + 2) * phi[j]^2 / rho_omega_h2 / 2
        X <: Complex ||
            (pert += w.h[medium] * im * w.b1c[L] * phi[j]^2 / rho_medium / 2)

        # medium
        L1 = L + 1
        L += w.n[medium] - 1
        j1 = j + 1
        j += w.n[medium] - 1
        for (LL, jj) in zip(L1:L, j1:j)
            sqnorm += w.h[medium] * phi[jj]^2 / rho_medium
            slow += w.h[medium] * (w.b1[LL] + 2) * phi[jj]^2 / rho_omega_h2
            X <: Complex ||
                (pert += w.h[medium] * im * w.b1c[LL] * phi[jj]^2 / rho_medium)
        end

        # bottom interface
        L += 1
        j += 1
        sqnorm += w.h[medium] * phi[j]^2 / rho_medium / 2
        slow += w.h[medium] * (w.b1[L] + 2) * phi[j]^2 / rho_omega_h2 / 2
        X <: Complex ||
            (pert += w.h[medium] * im * w.b1c[L] * phi[j]^2 / rho_medium / 2)
    end

    # contribution from the bottom halfspace
    if prob.hsbot.bc === :acoustoelastic
        if X <: Complex
            slow += phi[j]^2 / (2 * sqrt(x - ω² / prob.hsbot.cp^2)) /
                    (prob.hsbot.rho * prob.hsbot.cp^2)
        else
            # loop-length formula that works for leaky modes (kraken.f90):
            # difference of the ComplexFlag=.TRUE./.FALSE. impedances — which
            # is nonzero only for an acoustic halfspace (elastic halfspaces
            # ignore ComplexFlag in BCImpedanceMod.f90)
            if real(prob.hsbot.cs) == 0
                cnt = Ref(0)
                f1, g1, _ = bc_impedance(w, prob, x, :bot, prob.hsbot, cnt)
                gammaP = sqrt(x - ω² / prob.hsbot.cp^2)
                f2, g2 = gammaP, Complex{T}(prob.hsbot.rho)
                del = f2 / g2 - f1 / g1
                pert -= del * phi[j]^2
            end
            # group velocity still uses the perturbation method
            slow += phi[j]^2 / (2 * sqrt(complex(x - real(ω² / prob.hsbot.cp^2)))) /
                    (prob.hsbot.rho * real(prob.hsbot.cp)^2)
        end
    end

    # derivative of the top and bottom admittances w.r.t. x
    x1 = T(0.9999999) * x
    x2 = T(1.0000001) * x
    cnt = Ref(0)
    ftop1, gtop1, _ = bc_impedance(w, prob, x1, :top, prob.hstop, cnt)
    ftop2, gtop2, _ = bc_impedance(w, prob, x2, :top, prob.hstop, cnt)
    drhodx = zero(Complex{T})
    gtop1 != 0 && (drhodx = (ftop2 / gtop2 - ftop1 / gtop1) / (x2 - x1))

    fbot1, gbot1, _ = bc_impedance(w, prob, x1, :bot, prob.hsbot, cnt)
    fbot2, gbot2, _ = bc_impedance(w, prob, x2, :bot, prob.hsbot, cnt)
    detadx = zero(Complex{T})
    gbot1 != 0 && (detadx = (fbot2 / gbot2 - fbot1 / gbot1) / (x2 - x1))

    # scale the mode
    if X <: Complex
        rn = sqnorm - drhodx * phi[1]^2 + detadx * phi[ntotal1]^2
        scalefactor = 1 / sqrt(rn)      # generally a complex number
        real(scalefactor * phi[iturningpoint]) < 0 && (scalefactor = -scalefactor)
    else
        rn = real(sqnorm) - real(drhodx) * real(phi[1])^2 +
             real(detadx) * real(phi[ntotal1])^2
        rn <= 0 && (rn = -rn)           # grid-too-coarse warning case
        scalefactor = Complex{T}(1 / sqrt(rn))
        real(phi[iturningpoint]) < 0 && (scalefactor = -scalefactor)
    end

    phi .*= scalefactor
    slow *= scalefactor^2 * ω / sqrt(Complex{T}(x))
    pert *= scalefactor^2

    vg = X <: Complex ? real(1 / slow) : zero(T)   # KRAKEN leaves VG unset
    vg, pert
end

"""
    vector!(w, prob, evmat, kpert, vg, M, complex_solver, threads)
        -> (phi, zmesh)

Inverse iteration for all `M` eigenvectors on the first mesh (ports `Vector`,
kraken(c).f90). Fills `kpert` (perturbation incl. scatter loss) and `vg`.
Modes are independent, so they are optionally chunked across threads (each
task with `local` workspaces; deterministic — each mode writes its own
column).
"""
function vector!(w::MeshWorkspace{T,X}, prob::Problem{T}, evmat::Matrix{X},
                 kpert::Vector{Complex{T}}, vg::Vector{T}, M::Int,
                 complex_solver::Bool, threads::Int) where {T,X}
    ntotal = sum(w.n[w.first_acoustic:w.last_acoustic])
    ntotal1 = ntotal + 1

    # tabulate z-coordinates and off-diagonals of the matrix
    z = zeros(T, ntotal1)
    e0 = zeros(X, ntotal1 + 1)
    j = 1
    z[1] = prob.media[w.first_acoustic].z0
    h_rho = zero(T)
    for medium in w.first_acoustic:w.last_acoustic
        h_rho = w.h[medium] * w.rho[w.loc[medium]+1]   # rho at top of layer
        for jj in 1:w.n[medium]
            e0[j+jj] = 1 / h_rho
            z[j+jj] = z[j] + w.h[medium] * jj
        end
        j += w.n[medium]
    end
    e0[ntotal1+1] = 1 / h_rho    # dummy value; never used

    phi = zeros(Complex{T}, ntotal1, M)

    solve1mode! = (mode, d, e, phim) -> begin
        x = evmat[1, mode]
        cnt = Ref(0)
        copyto!(e, e0)

        # corner element requires the top impedance
        ftop, gtop, _ = bc_impedance(w, prob, x, :top, prob.hstop, cnt)
        if gtop == 0
            d[1] = one(X)
            e[2] = zero(X)
        else
            L = w.loc[w.first_acoustic] + 1
            xh2 = x * w.h[w.first_acoustic]^2
            hr = w.h[w.first_acoustic] * w.rho[L]
            fg = X <: Complex ? ftop / gtop : real(ftop / gtop)
            d[1] = (w.b1[L] - xh2) / hr / 2 + fg
        end

        # set up the diagonal; iTurningPoint = index closest to the surface
        # where the mode is still oscillatory
        iturningpoint = ntotal
        jj = 1
        L = w.loc[w.first_acoustic] + 1
        for medium in w.first_acoustic:w.last_acoustic
            xh2 = x * w.h[medium]^2
            hr = w.h[medium] * w.rho[w.loc[medium]+1]
            if medium >= w.first_acoustic + 1
                L += 1
                d[jj] = (d[jj] + (w.b1[L] - xh2) / hr) / 2
            end
            for _ in 1:w.n[medium]
                jj += 1
                L += 1
                d[jj] = (w.b1[L] - xh2) / hr
                if real(w.b1[L] - xh2) + 2 > 0    # turning point nearest top
                    iturningpoint = min(jj, iturningpoint)
                end
            end
        end

        # corner element requires the bottom impedance
        fbot, gbot, _ = bc_impedance(w, prob, x, :bot, prob.hsbot, cnt)
        if gbot == 0
            d[ntotal1] = one(X)
            e[ntotal1] = zero(X)
        else
            fg = X <: Complex ? fbot / gbot : real(fbot / gbot)
            d[ntotal1] = d[ntotal1] / 2 - fg
        end

        phix = Vector{X}(undef, ntotal1)
        ierr = inverse_iteration!(phix, d, e)
        if ierr != 0
            fill!(phim, zero(Complex{T}))   # zero out the errant eigenvector
        else
            phim .= Complex{T}.(phix)
            vgm, pert = normalize!(phim, w, prob, x, iturningpoint)
            vg[mode] = vgm
            kpert[mode] = pert + scatter_loss(w, prob, phim, x)
        end
    end

    if threads <= 1 || M < 2 * threads
        d = Vector{X}(undef, ntotal1)
        e = Vector{X}(undef, ntotal1 + 1)
        for mode in 1:M
            solve1mode!(mode, d, e, view(phi, :, mode))
        end
    else
        tasks = map(_chunks(M, threads)) do idxs
            Threads.@spawn begin
                local td = Vector{X}(undef, ntotal1)
                local te = Vector{X}(undef, ntotal1 + 1)
                for mode in idxs
                    solve1mode!(mode, td, te, view(phi, :, mode))
                end
            end
        end
        foreach(wait, tasks)
    end

    phi, z
end

"Partition `1:n` into at most `k` contiguous, near-equal chunks."
_chunks(n, k) = [(1 + div((j - 1) * n, k)):div(j * n, k) for j in 1:min(k, n)]
