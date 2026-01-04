using LinearAlgebra
using Statistics


"""
    sliding_counts_xy(x, y; L, step)

Count particles in a sliding square window of side length `L` (nm) over the XY-plane.
The window moves by `step` (nm) in X and Y. Returns:

- `counts` :: Vector{Int} counts per window
- `nx, ny` :: number of windows in X and Y
"""
function sliding_counts_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                           L::Real, step::Real)
    xmax = maximum(x)
    ymax = maximum(y)

    nx = max(1, floor(Int, (xmax - L) / step) + 1)
    ny = max(1, floor(Int, (ymax - L) / step) + 1)

    counts = Vector{Int}(undef, nx * ny)
    k = 1

    for i in 0:(nx - 1)
        x1 = i * step
        x2 = x1 + L

        idx_x = findall(@. (x >= x1) & (x <= x2))
        y_sel = y[idx_x]

        for j in 0:(ny - 1)
            y1 = j * step
            y2 = y1 + L
            counts[k] = count(@. (y_sel >= y1) & (y_sel <= y2))
            k += 1
        end
    end

    return counts, nx, ny
end

"""
    kbi_G22_from_coords(coords; sizes, steps, shift_to_zero=true)

Compute KBI estimate G₂₂ using number fluctuations in sliding windows.

- `coords` must be `N×3` (x,y,z) in nm (already cropped/rotated if desired)
- `sizes`: vector of window side lengths L (nm)
- `steps`: vector of XY step sizes (nm)

Returns a NamedTuple with vectors: `L`, `step`, `G22`.
"""
function kbi_G22_from_coords(coords::AbstractMatrix{<:Real};
                             sizes::AbstractVector{<:Real},
                             steps::AbstractVector{<:Real},
                             shift_to_zero::Bool=true)

    x = Float64.(coords[:, 1])
    y = Float64.(coords[:, 2])
    z = Float64.(coords[:, 3])

    if shift_to_zero
        x .-= minimum(x)
        y .-= minimum(y)
        z .-= minimum(z)
    end

    z_height = maximum(z) - minimum(z)

    L_vals = Float64[]
    step_vals = Float64[]
    G22_vals = Float64[]

    for step in steps
        for L in sizes
            counts, nx, ny = sliding_counts_xy(x, y; L=L, step=step)

            N_mean = mean(counts)
            G22 = if N_mean <= 0
                NaN
            else
                N_sq_mean = mean(counts .^ 2)
                variance = N_sq_mean - N_mean^2
                volume = L^2 * z_height  # nm^3
                (variance / N_mean^2 - 1 / N_mean) * volume
            end

            push!(L_vals, Float64(L))
            push!(step_vals, Float64(step))
            push!(G22_vals, G22)
        end
    end

    return (L=L_vals, step=step_vals, G22=G22_vals)
end

"""
    kbi_fit_G22_vs_sv(res; height, sv_min, sv_max)

Linear fit of `G22` vs `S/V` for each `step` contained in `res`.

Model:  G22 = slope*(S/V) + intercept

Returns a Vector of NamedTuples:
(step, slope, intercept, r2, n)

`res` must provide fields: `L`, `step`, `G22` (same format as returned by `kbi_G22_from_coords`).
"""
function kbi_fit_G22_vs_sv(res;
                           height::Real,
                           sv_min::Real,
                           sv_max::Real)

    steps_unique = unique(res.step)
    out = NamedTuple[]

    for st in steps_unique
        idx = findall((res.step .== st) .& isfinite.(res.G22) .& isfinite.(res.L))

        sv = kbi_surface_to_volume.(res.L[idx], height)
        y  = res.G22[idx]

        mask = (sv .>= sv_min) .& (sv .<= sv_max) .& isfinite.(sv) .& isfinite.(y)
        sv_fit = sv[mask]
        y_fit  = y[mask]

        if length(sv_fit) < 2
            push!(out, (step=st, slope=NaN, intercept=NaN, r2=NaN, n=length(sv_fit)))
            continue
        end

        A = hcat(sv_fit, ones(length(sv_fit)))
        m, b = A \ y_fit

        y_pred = m .* sv_fit .+ b
        ss_res = sum((y_fit .- y_pred).^2)
        ss_tot = sum((y_fit .- mean(y_fit)).^2)
        r2 = ss_tot == 0 ? NaN : 1 - ss_res / ss_tot

        push!(out, (step=st, slope=m, intercept=b, r2=r2, n=length(sv_fit)))
    end

    return out
end

"""
Surface-to-volume ratio S/V [nm^-1] for a box L×L×H (H is the cropped height Lz).
S = 2(L^2 + 2LH), V = L^2 H => S/V = 2/H + 4/L
"""
kbi_surface_to_volume(L::Real, H::Real) = 2/H + 4/L

"""
R² score like sklearn.metrics.r2_score.
Returns NaN if undefined (e.g., zero variance).
"""
function kbi_r2_score(y::AbstractVector{<:Real}, yhat::AbstractVector{<:Real})
    n = length(y)
    n == 0 && return NaN
    μ = sum(y) / n
    ss_res = sum((y .- yhat).^2)
    ss_tot = sum((y .- μ).^2)
    ss_tot == 0 ? NaN : 1 - ss_res / ss_tot
end

"""
Linear fit y = m*x + b (like np.polyfit(x,y,1)).
"""
function kbi_linear_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    A = hcat(collect(x), ones(length(x)))
    m, b = A \ collect(y)
    return m, b
end

"""
Scan windows and choose the best R².

- i_range: candidate start indices (1-based)
- j_range: candidate stop indices (EXCLUSIVE, like Python slicing x[i:j])
- min_points: minimum number of points in a window

Returns a NamedTuple:
(best_i, best_j, slope, intercept, r2, n_used, used_idx, yhat_full)

where best_j is the *exclusive* stop index.
"""
function kbi_best_window_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
    i_range::AbstractVector{<:Integer}=collect(1:3),
    j_range::AbstractVector{<:Integer}=collect(10:30),
    min_points::Integer=6
)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have same length"))

    best_r2 = -Inf
    best = nothing

    # ensure sorted by x (S/V)
    ord = sortperm(x)
    xs = collect(x)[ord]
    ys = collect(y)[ord]

    for i in i_range
        for j in j_range
            # Python-like slice: i:(j-1)
            if i < 1 || j < 2
                continue
            end
            if i >= j
                continue
            end
            lo = i
            hi = j - 1
            if hi > n
                continue
            end
            idx = lo:hi
            if length(idx) < min_points
                continue
            end

            m, b = kbi_linear_fit(xs[idx], ys[idx])
            yhat_win = m .* xs[idx] .+ b
            r2 = kbi_r2_score(ys[idx], yhat_win)

            if isfinite(r2) && (r2 > best_r2)
                best_r2 = r2
                best = (i=i, j=j, m=m, b=b, r2=r2, idx=idx, ord=ord)
            end
        end
    end

    if best === nothing
        return (best_i=missing, best_j=missing, slope=NaN, intercept=NaN, r2=NaN,
                n_used=0, used_idx=Int[], yhat_full=fill(NaN, n), x_sorted=xs, y_sorted=ys)
    end

    yhat_full = best.m .* xs .+ best.b

    return (best_i=best.i,
            best_j=best.j,               # exclusive stop index
            slope=best.m,
            intercept=best.b,            # this is G22,0
            r2=best.r2,
            n_used=length(best.idx),
            used_idx=collect(best.idx),
            yhat_full=yhat_full,
            x_sorted=xs,
            y_sorted=ys)
end
