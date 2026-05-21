"""
    read_particle_coords(path; factor=1.0)

Read a coordinate file with 3 columns (x y z) into an `N×3 Matrix{Float64}`.
Accepts comma-separated, tab-separated, or whitespace-separated files.
Skips empty lines, lines starting with `#`, and one non-numeric header row.
Applies an optional scaling factor to all coordinates.
"""
function read_particle_coords(path::AbstractString; factor::Real = 1.0)
    coords = Vector{NTuple{3, Float64}}()

    open(path, "r") do io
        for (lineno, line) in enumerate(eachline(io))
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue

            parts = occursin(',', s) ? split(s, ',') : split(s)
            parts = strip.(parts)
            if length(parts) < 3
                error("Expected 3 columns at $path:$lineno, got: '$s'")
            end

            parsed = tryparse.(Float64, parts[1:3])
            if any(isnothing, parsed)
                if isempty(coords) && all(isnothing, parsed)
                    continue
                end
                error("Expected numeric coordinates at $path:$lineno, got: '$s'")
            end

            x, y, z = Float64.(parsed)
            push!(coords, (x, y, z))
        end
    end

    n = length(coords)
    particles = Matrix{Float64}(undef, n, 3)
    @inbounds for i in 1:n
        particles[i, 1] = coords[i][1]
        particles[i, 2] = coords[i][2]
        particles[i, 3] = coords[i][3]
    end

    return particles .* Float64(factor)
end

# Accept either an output directory OR a file path and derive an output directory.
_outdir(path::AbstractString) = isdir(path) ? path : dirname(path)

"""
    save_rdf_data(r, gr, dr, sim_name, out; source_path=nothing)

Write RDF data to `<outdir>/<sim_name>_RDF.txt`.
`out` may be a directory or a filepath (directory will be inferred).
"""
function save_rdf_data(r, gr, dr, sim_name::AbstractString, out::AbstractString; source_path=nothing)
    outdir = _outdir(out)
    mkpath(outdir)

    outpath = joinpath(outdir, "$(sim_name)_RDF.txt")
    open(outpath, "w") do io
        write(io, "Radial Distribution Function for Experiment: ", sim_name, "\n")
        if source_path !== nothing
            write(io, "Original Path: ", string(source_path), "\n")
        end
        write(io, "Radial Bin Width: ", string(dr), "\n")
        write(io, "# r [nm]\t g(r)\n")
        @inbounds for k in eachindex(r)
            write(io, string(r[k]), "\t", string(gr[k]), "\n")
        end
    end

    return outpath
end

"""
    plot_rdf(r, gr, sim_name, out)

Save `<outdir>/<sim_name>_Gr.pdf`.
"""
function plot_rdf(r, gr, sim_name::AbstractString, out::AbstractString)
    outdir = _outdir(out)
    mkpath(outdir)
    outpath = joinpath(outdir, "$(sim_name)_Gr.pdf")

    plot(r, gr, label="", xlabel="r", ylabel="g(r)")
    plot!(r, ones(length(r)), linestyle=:dash, linewidth=1, label="")
    savefig(outpath)

    return outpath
end
