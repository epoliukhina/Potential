# --- Geometry helpers (same math, with Julia-style docstrings) ---

"""
    sphere_cut_volume(Rs, A, B)

Volume of the wedge left when two perpendicular planes at distances A and B
intersect a sphere of radius Rs (used for box/sphere intersection).
"""
function sphere_cut_volume(Rs::Float64, A::Float64, B::Float64)
    Root = sqrt(Rs^2 - A^2 - B^2)
    Vcut = (1/6) * Rs^3 * (pi - 2 * atan(A * B / (Rs * Root)))
    Vcut += (1/2) * (atan(A / Root) - pi / 2) * (Rs^2 * B - (1/3) * B^3)
    Vcut += (1/2) * (atan(B / Root) - pi / 2) * (Rs^2 * A - (1/3) * A^3)
    Vcut += (1/3) * A * B * Root
    return Vcut
end

"""
    octant_volume(Rs, xb, yb, zb)

Intersection volume between a sphere octant and a box corner in the first octant.
"""
function octant_volume(Rs::Float64, xb::Float64, yb::Float64, zb::Float64)
    if xb^2 + yb^2 + zb^2 < Rs^2
        return xb * yb * zb
    end

    V = (1/8) * (4/3) * pi * Rs^3

    for B in (xb, yb, zb)
        if B < Rs
            V -= (pi/4) * ((2/3) * Rs^3 - B * Rs^2 + (1/3) * B^3)
        end
    end

    for (a, b) in ((xb, yb), (xb, zb), (yb, zb))
        if a^2 + b^2 < Rs^2
            V += sphere_cut_volume(Rs, a, b)
        end
    end

    return V
end

"""
    sphere_volume(Rs, bounds)

Intersection volume of a sphere (centered at origin) with a rectangular box.
`bounds = [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax]` (distances from origin, can be ±).
"""
function sphere_volume(Rs::Float64, bounds::Vector{Float64})
    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = bounds
    V = 0.0
    for xb in (Xmin, Xmax), yb in (Ymin, Ymax), zb in (Zmin, Zmax)
        V += octant_volume(Rs, abs(xb), abs(yb), abs(zb))
    end
    return V
end

"""
    shell_volume(Rmin, Rmax, bounds)

Intersection volume of a spherical shell with a box.
"""
function shell_volume(Rmin::Float64, Rmax::Float64, bounds::Vector{Float64})
    Rmin = max(Rmin, 0.0)
    return sphere_volume(Rmax, bounds) - sphere_volume(Rmin, bounds)
end

# --- RDF core ---

# Find all bin indices k such that abs(r[k]-d) <= dr/2 (exactly your original rule).
@inline function _bin_range(r::AbstractVector{<:Real}, d::Float64, dr::Float64)
    lo = searchsortedfirst(r, d - dr/2)
    hi = searchsortedlast(r,  d + dr/2)
    return lo, hi
end

_thread_buffer_count() = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()

"""
    rdf_analytic_norm(particles, r, dr; threaded=true, show_progress=false)

Compute g(r) with analytic normalization using the intersection volume between spherical shells and the
tight bounding box of the particle set.
"""
function rdf_analytic_norm(particles::AbstractMatrix{<:Real}, r, dr::Real;
                           threaded::Bool=true, show_progress::Bool=false)

    size(particles, 2) == 3 || error("rdf_analytic_norm expects an N×3 matrix")
    dr = Float64(dr)

    nbins = length(r)
    r = collect(Float64, r)

    x = Float64.(view(particles, :, 1))
    y = Float64.(view(particles, :, 2))
    z = Float64.(view(particles, :, 3))

    bounds = [minimum(x), maximum(x),
              minimum(y), maximum(y),
              minimum(z), maximum(z)]

    Lx = bounds[2] - bounds[1]
    Ly = bounds[4] - bounds[3]
    Lz = bounds[6] - bounds[5]
    mean_density = length(x) / (Lx * Ly * Lz)

    n = length(x)
    nt = threaded ? _thread_buffer_count() : 1
    gr_threads = [zeros(Float64, nbins) for _ in 1:nt]
    shells_threads = [zeros(Float64, nbins) for _ in 1:nt]

    # ProgressMeter is not thread-safe; only show in single-thread mode.
    p = (show_progress && !threaded) ? Progress(n, dt=1.0) : nothing

    work = 1:n
    if threaded && Threads.nthreads() > 1
        Threads.@threads for central in work
            tid = Threads.threadid()
            local_gr = zeros(Float64, nbins)
            local_shells = zeros(Float64, nbins)

            for neigh in 1:n
                if neigh == central
                    continue
                end
                dx = x[central] - x[neigh]
                dy = y[central] - y[neigh]
                dz = z[central] - z[neigh]
                d = sqrt(dx*dx + dy*dy + dz*dz)

                lo, hi = _bin_range(r, d, dr)
                @inbounds for k in lo:hi
                    local_gr[k] += 1.0
                end
            end

            local_box = [bounds[1] - x[central], bounds[2] - x[central],
                         bounds[3] - y[central], bounds[4] - y[central],
                         bounds[5] - z[central], bounds[6] - z[central]]

            @inbounds for k in 1:nbins
                sv = shell_volume(r[k] - dr/2, r[k] + dr/2, local_box)
                if sv > 0.0
                    local_gr[k] /= sv
                    local_shells[k] += 1.0
                end
            end

            local_gr ./= mean_density

            thread_gr = gr_threads[tid]
            thread_shells = shells_threads[tid]
            @inbounds for k in 1:nbins
                thread_gr[k] += local_gr[k]
                thread_shells[k] += local_shells[k]
            end
        end
    else
        for central in work
            local_gr = zeros(Float64, nbins)
            local_shells = zeros(Float64, nbins)

            for neigh in 1:n
                if neigh == central
                    continue
                end
                dx = x[central] - x[neigh]
                dy = y[central] - y[neigh]
                dz = z[central] - z[neigh]
                d = sqrt(dx*dx + dy*dy + dz*dz)

                lo, hi = _bin_range(r, d, dr)
                @inbounds for k in lo:hi
                    local_gr[k] += 1.0
                end
            end

            local_box = [bounds[1] - x[central], bounds[2] - x[central],
                         bounds[3] - y[central], bounds[4] - y[central],
                         bounds[5] - z[central], bounds[6] - z[central]]

            @inbounds for k in 1:nbins
                sv = shell_volume(r[k] - dr/2, r[k] + dr/2, local_box)
                if sv > 0.0
                    local_gr[k] /= sv
                    local_shells[k] += 1.0
                end
            end

            local_gr ./= mean_density

            @inbounds for k in 1:nbins
                gr_threads[1][k] += local_gr[k]
                shells_threads[1][k] += local_shells[k]
            end

            p === nothing || next!(p)
        end
        p === nothing || finish!(p)
    end

    global_gr = zeros(Float64, nbins)
    nonempty = zeros(Float64, nbins)
    @inbounds for t in 1:nt, k in 1:nbins
        global_gr[k] += gr_threads[t][k]
        nonempty[k] += shells_threads[t][k]
    end

    @inbounds for k in 1:nbins
        if nonempty[k] > 0
            global_gr[k] /= nonempty[k]
        end
    end

    return global_gr
end

"""
    rdf_analytic_norm(particles1, particles2, r, dr; threaded=true, show_progress=false)

Cross-RDF: particles1 as centers, particles2 as neighbours. Normalization uses density of particles2.
"""
function rdf_analytic_norm(p1::AbstractMatrix{<:Real}, p2::AbstractMatrix{<:Real}, r, dr::Real;
                           threaded::Bool=true, show_progress::Bool=false)

    size(p1, 2) == 3 || error("rdf_analytic_norm expects an N×3 matrix for p1")
    size(p2, 2) == 3 || error("rdf_analytic_norm expects an N×3 matrix for p2")
    dr = Float64(dr)

    nbins = length(r)
    r = collect(Float64, r)

    x1 = Float64.(view(p1, :, 1)); y1 = Float64.(view(p1, :, 2)); z1 = Float64.(view(p1, :, 3))
    x2 = Float64.(view(p2, :, 1)); y2 = Float64.(view(p2, :, 2)); z2 = Float64.(view(p2, :, 3))

    bounds = [min(minimum(x1), minimum(x2)), max(maximum(x1), maximum(x2)),
              min(minimum(y1), minimum(y2)), max(maximum(y1), maximum(y2)),
              min(minimum(z1), minimum(z2)), max(maximum(z1), maximum(z2))]

    Lx = bounds[2] - bounds[1]
    Ly = bounds[4] - bounds[3]
    Lz = bounds[6] - bounds[5]
    mean_density = length(x2) / (Lx * Ly * Lz)

    n1 = length(x1)
    n2 = length(x2)

    nt = threaded ? _thread_buffer_count() : 1
    gr_threads = [zeros(Float64, nbins) for _ in 1:nt]
    shells_threads = [zeros(Float64, nbins) for _ in 1:nt]

    p = (show_progress && !threaded) ? Progress(n1, dt=1.0) : nothing

    if threaded && Threads.nthreads() > 1
        Threads.@threads for central in 1:n1
            tid = Threads.threadid()
            local_gr = zeros(Float64, nbins)
            local_shells = zeros(Float64, nbins)

            for neigh in 1:n2
                dx = x1[central] - x2[neigh]
                dy = y1[central] - y2[neigh]
                dz = z1[central] - z2[neigh]
                d = sqrt(dx*dx + dy*dy + dz*dz)

                lo, hi = _bin_range(r, d, dr)
                @inbounds for k in lo:hi
                    local_gr[k] += 1.0
                end
            end

            local_box = [bounds[1] - x1[central], bounds[2] - x1[central],
                         bounds[3] - y1[central], bounds[4] - y1[central],
                         bounds[5] - z1[central], bounds[6] - z1[central]]

            @inbounds for k in 1:nbins
                sv = shell_volume(r[k] - dr/2, r[k] + dr/2, local_box)
                if sv > 0.0
                    local_gr[k] /= sv
                    local_shells[k] += 1.0
                end
            end

            local_gr ./= mean_density

            thread_gr = gr_threads[tid]
            thread_shells = shells_threads[tid]
            @inbounds for k in 1:nbins
                thread_gr[k] += local_gr[k]
                thread_shells[k] += local_shells[k]
            end
        end
    else
        for central in 1:n1
            local_gr = zeros(Float64, nbins)
            local_shells = zeros(Float64, nbins)

            for neigh in 1:n2
                dx = x1[central] - x2[neigh]
                dy = y1[central] - y2[neigh]
                dz = z1[central] - z2[neigh]
                d = sqrt(dx*dx + dy*dy + dz*dz)

                lo, hi = _bin_range(r, d, dr)
                @inbounds for k in lo:hi
                    local_gr[k] += 1.0
                end
            end

            local_box = [bounds[1] - x1[central], bounds[2] - x1[central],
                         bounds[3] - y1[central], bounds[4] - y1[central],
                         bounds[5] - z1[central], bounds[6] - z1[central]]

            @inbounds for k in 1:nbins
                sv = shell_volume(r[k] - dr/2, r[k] + dr/2, local_box)
                if sv > 0.0
                    local_gr[k] /= sv
                    local_shells[k] += 1.0
                end
            end

            local_gr ./= mean_density

            @inbounds for k in 1:nbins
                gr_threads[1][k] += local_gr[k]
                shells_threads[1][k] += local_shells[k]
            end

            p === nothing || next!(p)
        end
        p === nothing || finish!(p)
    end

    global_gr = zeros(Float64, nbins)
    nonempty = zeros(Float64, nbins)
    @inbounds for t in 1:nt, k in 1:nbins
        global_gr[k] += gr_threads[t][k]
        nonempty[k] += shells_threads[t][k]
    end

    @inbounds for k in 1:nbins
        if nonempty[k] > 0
            global_gr[k] /= nonempty[k]
        end
    end

    return global_gr
end

"""
    process_particle_coords(particles, r, dr, out, sim_name)

Convenience wrapper: compute g(r), save txt, save pdf.
Returns `(gr, txt_path, pdf_path)`.
"""
function process_particle_coords(particles, r, dr, out, sim_name)
    gr = rdf_analytic_norm(particles, r, dr)
    txt = save_rdf_data(r, gr, dr, sim_name, out)
    pdf = plot_rdf(r, gr, sim_name, out)
    return gr, txt, pdf
end
