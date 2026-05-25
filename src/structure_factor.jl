using LinearAlgebra
using ProgressMeter
using DelimitedFiles
using Plots
using PlotlyBase

# Safe sinc: sin(x)/x with small-x handling (MATLAB sinc_function.m behavior)
sinc_safe(x::Real) = abs(x) < 1e-5 ? 1.0 : sin(x) / x

"""
    structure_factor_debye_histogram(coords_nm; qmin=5e-5, qmax=0.6, nq=500, stepA=1.0, show_progress=true)

Compute S(q) from particle positions using Debye formula and a histogram of pair distances.

- `coords_nm` : N×3 coordinates in **nm**
- q is returned in **Å⁻¹** (matches your MATLAB script)
- `stepA` : histogram bin width in **Å** (default 1 Å = 0.1 nm)
"""
function structure_factor_debye_histogram(coords_nm::AbstractMatrix{<:Real};
    qmin::Float64 = 0.00005,
    qmax::Float64 = 0.6,
    nq::Int = 500,
    stepA::Float64 = 1.0,
    show_progress::Bool = true
)
    N = size(coords_nm, 1)
    N < 2 && error("Need at least 2 particles for structure factor.")

    # Convert nm → Å (to match MATLAB script)
    coordsA = 10.0 .* Float64.(coords_nm)

    # Upper bound for max distance: diagonal of bounding box
    xmin, xmax = minimum(coordsA[:,1]), maximum(coordsA[:,1])
    ymin, ymax = minimum(coordsA[:,2]), maximum(coordsA[:,2])
    zmin, zmax = minimum(coordsA[:,3]), maximum(coordsA[:,3])
    maxdist = sqrt((xmax-xmin)^2 + (ymax-ymin)^2 + (zmax-zmin)^2)

    nbins = Int(ceil(maxdist / stepA)) + 1
    counts = zeros(Int, nbins)

    # Histogram pair distances without storing them all
    total_pairs = N*(N-1) ÷ 2
    p = show_progress ? Progress(total_pairs; dt=1.0) : nothing

    @inbounds for i in 1:(N-1)
        xi, yi, zi = coordsA[i,1], coordsA[i,2], coordsA[i,3]
        for j in (i+1):N
            dx = xi - coordsA[j,1]
            dy = yi - coordsA[j,2]
            dz = zi - coordsA[j,3]
            r = sqrt(dx*dx + dy*dy + dz*dz)

            b = Int(floor(r / stepA)) + 1
            b = clamp(b, 1, nbins)
            counts[b] += 1

            show_progress && next!(p)
        end
    end
    show_progress && finish!(p)

    # Bin centers (Å), matching MATLAB: edges_cnt(1:end-1) + step/2
    r_centers = (0:(nbins-1)) .* stepA .+ stepA/2

    q = collect(range(qmin, qmax; length=nq))
    Iq = zeros(Float64, nq)

    @inbounds for k in 1:nq
        qq = q[k]
        s = 0.0
        for b in 1:nbins
            c = counts[b]
            c == 0 && continue
            s += c * sinc_safe(qq * r_centers[b])
        end
        Iq[k] = N + 2.0*s
    end

    Sq = Iq ./ N
    return (q=q, Sq=Sq, r=r_centers, counts=counts)
end

function save_structure_factor(q, Sq, outpath::AbstractString)
    writedlm(outpath, hcat(q, Sq), '\t')
    return true
end

function plot_structure_factor(q, Sq; xlim_tuple=nothing, ylim_tuple=nothing, title_str="Structure factor S(q)")
    plotly()
    p = plot(q, Sq;
        label="",
        xlabel="q [Å⁻¹]",
        ylabel="S(q) [-]",
        title=title_str,
    )
    xlim_tuple !== nothing && xlims!(p, xlim_tuple...)
    ylim_tuple !== nothing && ylims!(p, ylim_tuple...)
    return p
end
