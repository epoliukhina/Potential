### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 53f83c0f-2d6d-4aa9-8f77-84f0903ebebd
begin
    using Plots, ProgressMeter, PlutoUI, DelimitedFiles, Statistics, Pkg

    const TITLE_FS = 11
    const GUIDE_FS = 11
    const TICK_FS  = 9

    # Fixed marker sizes (pixels)
    const MS_RAW_3D = 1
    const MS_RAW_2D = 1
    const MS_CUT_3D = 1
    const MS_CUT_2D = 1

    # Use GR as the default backend (PDF saving works reliably).
    gr()

    # small helper: ensure directory string ends with path separator
    dirprefix(d::AbstractString) = endswith(d, string(Base.Filesystem.path_separator)) ? d :
                                  d * string(Base.Filesystem.path_separator)

    # --- Hide code by default; button toggles showing inputs ---
    PlutoUI.HTML("""
    <style>
      body.hide-code pluto-input { display: none !important; }
      #toggle-code-btn{
        position: fixed;
        top: 12px;
        right: 12px;
        z-index: 999999;
        padding: 8px 10px;
        border-radius: 10px;
        border: 1px solid rgba(0,0,0,0.18);
        background: rgba(255,255,255,0.92);
        cursor: pointer;
        font-size: 12px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.08);
      }
      #toggle-code-btn:hover{ background: rgba(255,255,255,1.0); }
    </style>

    <button id="toggle-code-btn" type="button">Show/Hide code</button>

    <script>
      document.body.classList.add("hide-code");
      const btn = document.getElementById("toggle-code-btn");
      btn.addEventListener("click", (ev) => {
        ev.preventDefault();
        ev.stopPropagation();
        document.body.classList.toggle("hide-code");
      }, true);
    </script>
    """)

    TableOfContents()
end

# ╔═╡ 8e7fa67f-85ae-4e66-aaa3-cb27df2aaa2d
begin
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

# ╔═╡ 9fd66a2c-dc50-4dfe-bc64-76b8a57e051c
using Potential

# ╔═╡ 9a0b0e4c-ef77-49e4-bc7f-4e97a27c7e6e
begin
    # ===== Structure factor utilities are loaded from an external file (no function definitions in Pluto) =====
    const SF_CANDIDATES = [
        joinpath(@__DIR__, "structure_factor.jl"),
        joinpath(@__DIR__, "..", "structure_factor.jl"),
        joinpath(@__DIR__, "..", "src", "structure_factor.jl"),
        joinpath(@__DIR__, "..", "Potential", "src", "structure_factor.jl"),
    ]

    sf_idx = findfirst(isfile, SF_CANDIDATES)
    sf_idx === nothing && error("structure_factor.jl not found. Place it next to the notebook or in ../src/.")

    include(SF_CANDIDATES[sf_idx])

    # Aliases for readability (match RDF style)
    sf_debye_hist = structure_factor_debye_histogram
    save_sf_data  = save_structure_factor
    plot_sf       = plot_structure_factor

    nothing
end

# ╔═╡ dc937899-81ef-4326-8acd-94a09e8cee0d
begin
    # Aliases for functions (keeps notebook code readable)
    read_particle_coords = Potential.readPartCoor
    rdf_analytic_norm    = Potential.RDF_AnalyticNorm
    save_rdf_data        = Potential.saveRDFdata
    plot_rdf             = Potential.plotRDF

    rotate_particles     = Potential.rotate_particles
    linear_fit           = Potential.linear_fit

    has_kbi = isdefined(Potential, :kbi_G22_from_coords)
    has_sv  = isdefined(Potential, :kbi_surface_to_volume)

    nothing
end

# ╔═╡ 9266fa4c-6757-4d8d-837c-06ffb8f7fd6c
begin
    plot_layout = @layout [a b c]
    nothing
end

# ╔═╡ 0f8dcadf-e874-4c0f-bca5-b6063d538134
md"""
# cryo-POT: RDF, PMF, G₂₂, S(q)

## 1) Load 3D coordinates
Input file (txt/tsv/csv): $(@bind initial_data_raw FilePicker())
"""

# ╔═╡ 13cc1e0c-07f9-4d73-97a3-8cfcfa7d902a
begin
    PlutoUI.HTML("""
    <div style="margin-top: 0.35rem;">
      <b>Selected file:</b>
      <div style="
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
        word-break: break-all;
        white-space: normal;
        line-height: 1.25;">
        $(get(initial_data_raw, "name", ""))
      </div>
    </div>
    """)
end

# ╔═╡ bff7c13d-c197-4a0a-88d7-4d174d2c5f6c
md"""
Output directory (base): $(@bind outdir_base_raw TextField((80, 1), default=joinpath(homedir(), "git", "Calculations")))
"""

# ╔═╡ 2cf2a7c5-0b7c-40df-a1f5-8d2cbd0ddc8d
begin
    outdir = joinpath(outdir_base_raw, splitext(initial_data_raw["name"])[1])
    mkpath(outdir)

    particle_file = joinpath(outdir, initial_data_raw["name"])
    particle_file_raw = joinpath(outdir, "Raw_" * initial_data_raw["name"])

    initial_data_text = String(copy(initial_data_raw["data"]))
    open(particle_file_raw, "w") do io
        write(io, initial_data_text)
    end

    nothing
end

# ╔═╡ 4b8cc26f-d7d6-4c8c-8ad5-4d3fbbf605f1
begin
    sim_name   = splitext(basename(particle_file))[1]
    out_prefix = dirprefix(outdir)

    rdf_txt_path = joinpath(outdir, "$(sim_name)_RDF.txt")
    rdf_pdf_path = joinpath(outdir, "$(sim_name)_Gr.pdf")

    sq_txt_path  = joinpath(outdir, "$(sim_name)_Sq.txt")
    sq_pdf_path  = joinpath(outdir, "$(sim_name)_Sq.pdf")

    kbi_txt_path       = joinpath(outdir, "$(sim_name)_KBI.txt")
    kbi_pdf_path       = joinpath(outdir, "$(sim_name)_KBI.pdf")
    kbi_svfit_pdf_path = joinpath(outdir, "$(sim_name)_KBI_SVfit.pdf")

    nothing
end

# ╔═╡ 5b3b0edb-3d84-46e9-9c7c-0d7b0b0e5b35
md"""
## 2) Basic parameters
"""

# ╔═╡ 1f9ad9f0-4f5e-4a35-8b8f-f5c03d7b1c67
function basic_params_inputs(vars::Vector{String}; mw_default="66500", factor_default="1")
    PlutoUI.combine() do Child
        md"""
Molecular weight (Da) — used to calculate the concentration in mg/mL:  
$(Child(vars[1], TextField(default=mw_default)))

Scale factor (units per pixel) — px → nm conversion:  
$(Child(vars[2], TextField(default=factor_default)))

**Note:** Use **nm/px** (typically pixel spacing × binning).  
If coordinates are already in **nm**, set **factor = 1**.
"""
    end
end

# ╔═╡ 33a2df48-7e7d-4f49-8d59-28c17c7c529e
let prev = @isdefined(basic_params) ? basic_params : ["66500", "1"]
    @bind basic_params confirm(basic_params_inputs(
        ["mw_str", "factor_str"];
        mw_default=string(prev[1]),
        factor_default=string(prev[2]),
    ))
end

# ╔═╡ 3c6bd3c2-45a4-493b-8f2e-4fa6de7fcbce
md"""
## 3) Rotation
"""

# ╔═╡ 4d6714b1-2d73-4b4b-a901-1d2c83cfe63f
begin
    function read_coords_any(path::AbstractString)
        ext = lowercase(splitext(path)[2])
        if ext == ".csv"
            try
                mat = readdlm(path, ',', Float64)
                return mat[:, 1:3]
            catch
                mat, _ = readdlm(path, ','; header=true)
                mat = Float64.(mat)
                return mat[:, 1:3]
            end
        elseif ext == ".tsv" || ext == ".tab"
            mat = Float64.(readdlm(path, '\t'))
            return mat[:, 1:3]
        else
            mat = Float64.(readdlm(path))  # spaces + tabs
            return mat[:, 1:3]
        end
    end

    factor = parse(Float64, basic_params[2])
    data_matrix = factor .* read_coords_any(particle_file_raw)
    shift = minimum(data_matrix, dims=1)
    data_matrix .-= shift

    x, y, z = data_matrix[:, 1], data_matrix[:, 2], data_matrix[:, 3]

    # 3D: try Plotly for interactivity; fallback to GR; then restore GR
    local p
    try
        plotly()
        p = scatter3d(
            x, y, z;
            markersize=MS_RAW_3D,
            label="",
            title="Raw coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
    catch
        gr()
        p = scatter3d(
            x, y, z;
            markersize=MS_RAW_3D,
            label="",
            title="Raw coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
    end
    gr()
    p
end

# ╔═╡ 7a80e2a6-2c96-4cf3-bc55-43c45f36ab06
begin
    points_3d = hcat(x, y, z)

    points_3d_rot = rotate_particles(points_3d)
    x_rot = points_3d_rot[:, 1]
    y_rot = points_3d_rot[:, 2]
    z_rot = points_3d_rot[:, 3]

    kx,  bx  = linear_fit(hcat(x, z)); fx(xv)  = kx  * xv + bx
    ky,  by  = linear_fit(hcat(y, z)); fy(yv)  = ky  * yv + by
    kxy, bxy = linear_fit(hcat(x, y)); fxy(xv) = kxy * xv + bxy

    gr()
    p1 = scatter(x, z; label="data", xlabel="X [nm]", ylabel="Z [nm]", markersize=MS_RAW_2D)
    plot!(x, fx.(x); label="fit", lw=3)

    p2 = scatter(y, z; label="data", xlabel="Y [nm]", ylabel="Z [nm]", markersize=MS_RAW_2D)
    plot!(y, fy.(y); label="fit", lw=3)

    p3 = scatter(x, y; label="data", xlabel="X [nm]", ylabel="Y [nm]", markersize=MS_RAW_2D)
    plot!(x, fxy.(x); label="fit", lw=3)

    plot(
        p1, p2, p3;
        layout=plot_layout,
        size=(1450, 450),
        plot_title="Projections + linear fits (used for rotation)",
        titlefontsize=TITLE_FS,
        guidefontsize=GUIDE_FS,
        tickfontsize=TICK_FS,
        left_margin=12Plots.mm, bottom_margin=12Plots.mm, top_margin=4Plots.mm, right_margin=6Plots.mm,
    )
end

# ╔═╡ 6a516c61-4b34-4b9b-a4be-ec12e0efc1b7
begin
    local p
    try
        plotly()
        p = scatter3d(
            x, y, z;
            markersize=MS_RAW_3D,
            label="raw",
            title="Raw vs rotated coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
        scatter3d!(x_rot, y_rot, z_rot; markersize=MS_RAW_3D, label="rotated")
    catch
        gr()
        p = scatter3d(
            x, y, z;
            markersize=MS_RAW_3D,
            label="raw",
            title="Raw vs rotated coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
        scatter3d!(x_rot, y_rot, z_rot; markersize=MS_RAW_3D, label="rotated")
    end
    gr()
    p
end

# ╔═╡ 0c865bd7-9e54-4e56-9327-e5d7239a2c2b
function crop_inputs(vars::Vector{String}; cut_z_default="50", cut_xy_default="350")
    PlutoUI.combine() do Child
        md"""
## 4) Cropping

Cut half-length in **Z** (nm):  
$(Child(vars[1], TextField(default=cut_z_default)))

Cut half-length in **X and Y** (nm):  
$(Child(vars[2], TextField(default=cut_xy_default)))
"""
    end
end

# ╔═╡ 2a7f6200-1a16-43f3-8a90-fb992fa0d8a2
let prev = @isdefined(crop_params) ? crop_params : ["50", "350"]
    @bind crop_params confirm(crop_inputs(
        ["cut_z_str", "cut_xy_str"];
        cut_z_default=string(prev[1]),
        cut_xy_default=string(prev[2]),
    ))
end

# ╔═╡ 7b0109f6-2cf8-4cf5-acde-8427df12b545
begin
    cut_length_z  = parse(Float64, crop_params[1])
    cut_length_xy = parse(Float64, crop_params[2])

    mx, my, mz = mean(x_rot), mean(y_rot), mean(z_rot)

    x_lo, x_hi = mx - cut_length_xy, mx + cut_length_xy
    y_lo, y_hi = my - cut_length_xy, my + cut_length_xy
    z_lo, z_hi = mz - cut_length_z,  mz + cut_length_z

    keep = findall(
        (z_rot .>= z_lo) .& (z_rot .<= z_hi) .&
        (x_rot .>= x_lo) .& (x_rot .<= x_hi) .&
        (y_rot .>= y_lo) .& (y_rot .<= y_hi)
    )

    points_3d_rot_cut = points_3d_rot[keep, :]
    x_rot_cut = points_3d_rot_cut[:, 1]
    y_rot_cut = points_3d_rot_cut[:, 2]
    z_rot_cut = points_3d_rot_cut[:, 3]

    local p
    try
        plotly()
        p = scatter3d(
            x_rot_cut, y_rot_cut, z_rot_cut;
            markersize=MS_CUT_3D,
            label="",
            title="Cropped coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
    catch
        gr()
        p = scatter3d(
            x_rot_cut, y_rot_cut, z_rot_cut;
            markersize=MS_CUT_3D,
            label="",
            title="Cropped coordinates",
            titlefontsize=TITLE_FS,
            xlabel="X [nm]", ylabel="Y [nm]", zlabel="Z [nm]",
            guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
        )
    end
    gr()
    p
end

# ╔═╡ 59ec80a7-689f-4258-b73b-a25aec2a9880
begin
    gr()

    # X–Z
    p1c = scatter(
        x_rot_cut, z_rot_cut;
        label="",
        xlabel="X [nm]", ylabel="Z [nm]",
        markersize=MS_CUT_2D,
        titlefontsize=TITLE_FS, guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
    )

    # Y–Z
    p2c = scatter(
        y_rot_cut, z_rot_cut;
        label="",
        xlabel="Y [nm]", ylabel="Z [nm]",
        markersize=MS_CUT_2D,
        titlefontsize=TITLE_FS, guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
    )

    # X–Y
    p3c = scatter(
        x_rot_cut, y_rot_cut;
        label="",
        xlabel="X [nm]", ylabel="Y [nm]",
        markersize=MS_CUT_2D,
        titlefontsize=TITLE_FS, guidefontsize=GUIDE_FS, tickfontsize=TICK_FS
    )

    plot(
        p1c, p2c, p3c;
        layout=plot_layout,
        size=(1450, 450),
        plot_title="Cropped projections",
        titlefontsize=TITLE_FS,
        guidefontsize=GUIDE_FS,
        tickfontsize=TICK_FS,
        left_margin=12Plots.mm, bottom_margin=12Plots.mm, top_margin=4Plots.mm, right_margin=6Plots.mm,
    )
end

# ╔═╡ 5c375e8a-8c35-4a3a-9f88-56b5f9a01b69
begin
    Lz_kbi = maximum(points_3d_rot_cut[:, 3]) - minimum(points_3d_rot_cut[:, 3])
    nothing
end

# ╔═╡ 64db4a58-c750-4a75-84ad-323f9f5ec593
begin
    mw = parse(Float64, basic_params[1])
    n_particles = size(points_3d_rot_cut, 1)

    Lx = maximum(x_rot_cut) - minimum(x_rot_cut)
    Ly = maximum(y_rot_cut) - minimum(y_rot_cut)
    Lz = maximum(z_rot_cut) - minimum(z_rot_cut)

    V_nm3 = Lx * Ly * Lz
    number_density = n_particles / V_nm3 * 1e21

    NA = 6.02214076e23
    concentration_mg_per_ml = number_density / NA * mw * 1e3

    md"""
## 5) Box summary

| parameter | value |
|---|---:|
| Number of particles | **$(n_particles)** |
| Number density (particles/cm³) | **$(round(number_density, sigdigits=3))**|
| Lx (nm) | **$(round(Lx; digits=1))** |
| Ly (nm) | **$(round(Ly; digits=1))** |
| Lz (nm) | **$(round(Lz; digits=1))** |
| Volume (nm³) | **$(round(V_nm3, sigdigits=3))** |
| Concentration (mg/mL) | **$(round(concentration_mg_per_ml; digits=1))** |
"""
end

# ╔═╡ 0f73b3bc-25c0-4c33-b082-6b4d2d3c4518
begin
    writedlm(particle_file, points_3d_rot_cut, ' ')
    nothing
end

# ╔═╡ 6eaa7bb3-2320-43ad-9c3e-b1c2e66fbd42
function rdf_inputs(vars::Vector; RMax_default="50", NBins_default="100", DeltaR_default="0.3")
    PlutoUI.combine() do Child
        md"""
## 6) RDF g(r) and PMF W(r)

Max radial distance **RMax** (nm): $(Child(vars[1], TextField(default=RMax_default)))

Number of bins **NBins**: $(Child(vars[2], TextField(default=NBins_default)))

Bin width **Δr** (nm): $(Child(vars[3], TextField(default=DeltaR_default)))
"""
    end
end

# ╔═╡ 2b1b55a9-06a8-486f-8b1b-79f11b8bc496
let prev = @isdefined(RDF_variables) ? RDF_variables : ["50", "100", "0.3"]
    @bind RDF_variables confirm(rdf_inputs(
        ["RMax_str", "NBins_str", "DeltaR_str"];
        RMax_default=string(prev[1]),
        NBins_default=string(prev[2]),
        DeltaR_default=string(prev[3]),
    ))
end

# ╔═╡ 2e21cccf-9860-4f52-9a0d-9f7c5f229bea
begin
    RMax   = parse(Float64, RDF_variables[1])
    NBins  = parse(Int,     RDF_variables[2])
    DeltaR = parse(Float64, RDF_variables[3])

    rList = range(0, stop=RMax, length=NBins)
    Gr = rdf_analytic_norm(points_3d_rot_cut, rList, DeltaR)

    _ = save_rdf_data(rList, Gr, DeltaR, sim_name, out_prefix)

    # Ensure GR backend for PDF saving inside plotRDF
    gr()
    _ = plot_rdf(rList, Gr, sim_name, out_prefix)

    nothing
end

# ╔═╡ 88d4a3ea-2f6f-4e48-a7f0-08c7e50a31b0
begin
    gr()
    plot(rList, ones(length(rList)); color="grey", linestyle=:dash, linewidth=2, label="")
    plot!(rList, Gr; color="green", linewidth=2, label="", xlabel="r [nm]", ylabel="g(r)",
          title="Radial distribution function", titlefontsize=TITLE_FS, guidefontsize=GUIDE_FS, tickfontsize=TICK_FS)
end

# ╔═╡ 947c3043-5c0d-4af4-8a62-a9ee23f5c1f9
begin
    gr()
    plot(rList, zeros(length(rList)); color="grey", linestyle=:dash, linewidth=2, label="")

    plot!(rList, -log.(Gr);
          color="green", linewidth=2, label="",
          xlabel="r [nm]",
          ylabel="W [kT]",
          title="Potential of mean force",
          titlefontsize=TITLE_FS, guidefontsize=GUIDE_FS, tickfontsize=TICK_FS)
end

# ╔═╡ 5e8d1b43-fc7c-4c53-9da6-1877c50fd7cc
function sf_inputs(
    vars::Vector{String};
    qmin_default="0.00005",
    qmax_default="0.6",
    nq_default="500",
    stepA_default="1.0",
    use_window_default=true,
    save_sf_default=true
)
    PlutoUI.combine() do Child
        md"""
## 7) Structure factor S(q)

**qmin** (Å⁻¹): $(Child(vars[1], TextField(default=qmin_default)))

**qmax** (Å⁻¹): $(Child(vars[2], TextField(default=qmax_default)))

Number of q points **N**: $(Child(vars[3], TextField(default=nq_default)))

Histogram bin width **Δr** (Å): $(Child(vars[4], TextField(default=stepA_default)))

Use SAXS view window (xlim 0.019–0.2 Å⁻¹, ylim 0–1.5): $(Child(vars[5], CheckBox(default=use_window_default)))

Save outputs (`*_Sq.txt`, `*_Sq.pdf`): $(Child(vars[6], CheckBox(default=save_sf_default)))
"""
    end
end

# ╔═╡ 7a77a3f2-0b8c-47b9-87a5-1111e0a8c0b6
let prev = @isdefined(SF_params) ? SF_params : ["0.00005","0.6","500","1.0", true, true]
    @bind SF_params confirm(sf_inputs(
        ["qmin_str","qmax_str","nq_str","stepA_str","use_window","save_sf"];
        qmin_default=string(prev[1]),
        qmax_default=string(prev[2]),
        nq_default=string(prev[3]),
        stepA_default=string(prev[4]),
        use_window_default=Bool(prev[5]),
        save_sf_default=Bool(prev[6]),
    ))
end

# ╔═╡ 24f0e9ea-95f6-4f3c-87d2-0c2d9e7b8b79
begin
    qmin  = parse(Float64, SF_params[1])
    qmax  = parse(Float64, SF_params[2])
    nq    = parse(Int,     SF_params[3])
    stepA = parse(Float64, SF_params[4])

    use_window = Bool(SF_params[5])
    save_sf    = Bool(SF_params[6])

    res_sq = sf_debye_hist(points_3d_rot_cut; qmin=qmin, qmax=qmax, nq=nq, stepA=stepA)

    xlim_t = use_window ? (0.019, 0.2) : nothing
    ylim_t = use_window ? (0.0, 1.5)   : nothing

    # IMPORTANT: plot_sf must be called only with keywords it supports
    p_sq = plot_sf(res_sq.q, res_sq.Sq; xlim_tuple=xlim_t, ylim_tuple=ylim_t,
                   title_str="Structure factor S(q)")

    # Restyle the curve that plot_sf already created
    p_sq.series_list[1].plotattributes[:linewidth]   = 3
    p_sq.series_list[1].plotattributes[:linecolor]   = :green
    p_sq.series_list[1].plotattributes[:seriescolor] = :green

    # Reference line S(q)=1
    plot!(p_sq, res_sq.q, ones(length(res_sq.q));
          color=:grey, linestyle=:dash, linewidth=2, label="")

    plot!(p_sq;
          titlefontsize=TITLE_FS,
          guidefontsize=GUIDE_FS,
          tickfontsize=TICK_FS,
          size=(600, 500),
          left_margin=12Plots.mm, bottom_margin=12Plots.mm, top_margin=8Plots.mm, right_margin=6Plots.mm,
          label=""
    )

    if save_sf
        save_sf_data(res_sq.q, res_sq.Sq, sq_txt_path)
        gr()
        savefig(p_sq, sq_pdf_path)
    end

    p_sq
end

# ╔═╡ 4ec7d1a8-2d54-4ed1-92fb-4a5be6b7d2fd
function kbi_inputs(vars::Vector; sizes_default="vcat(10:1:99, 100:10:600)", steps_default="3:1:3", save_kbi_default=true)
    PlutoUI.combine() do Child
        md"""
## 8) KBI G₂₂

Box sizes **sizes** (nm, Julia expression):  
$(Child(vars[1], TextField(default=sizes_default)))

Step(s) **steps** (nm, Julia expression):  
$(Child(vars[2], TextField(default=steps_default)))

**Note**: Build a vector of values: 10–99 with step 1, then 100–600 with step 10 using `vcat(10:1:99, 100:10:600)`.
		
Save KBI outputs (`*_KBI.txt`, `*_KBI.pdf`, `*_KBI_SVfit.pdf`):  
$(Child(vars[3], CheckBox(default=save_kbi_default)))
"""
    end
end

# ╔═╡ 5e6e2bfa-c2d3-4b76-a2c4-4a7f81d86f4a
let prev = @isdefined(KBI_params) ? KBI_params : ["vcat(10:1:99, 100:10:600)", "3:1:3", true]
    @bind KBI_params confirm(kbi_inputs(
        ["sizes_expr", "steps_expr", "save_kbi"];
        sizes_default=string(prev[1]),
        steps_default=string(prev[2]),
        save_kbi_default=Bool(prev[3]),
    ))
end

# ╔═╡ 1c80f764-9d2d-4dff-9a9a-1771169f788f
begin
    function eval_vec_expr(s::AbstractString)
        ex = Meta.parse(strip(s))
        v  = Base.eval(Main, ex)
        collect(v)
    end
    function eval_int_vec_expr(s::AbstractString)
        Int.(eval_vec_expr(s))
    end
    nothing
end

# ╔═╡ 1166aef4-3c77-498b-8eb3-a2066f2a1fe5
begin
    r2_score(y, yhat) = begin
        μ = mean(y)
        ss_res = sum((y .- yhat).^2)
        ss_tot = sum((y .- μ).^2)
        ss_tot == 0 ? NaN : 1 - ss_res/ss_tot
    end

    linfit(x, y) = begin
        A = hcat(collect(x), ones(length(x)))
        m, b = A \ collect(y)
        m, b
    end

    function best_window_fit(x, y; i_range::Vector{Int}, j_range::Vector{Int}, min_points::Int=8, from_end::Bool=true)
        ord = sortperm(x)
        xs = collect(x)[ord]
        ys = collect(y)[ord]
        n  = length(xs)

        xscan = from_end ? reverse(xs) : xs
        yscan = from_end ? reverse(ys) : ys

        best_r2 = -Inf
        best = nothing

        for i in i_range, j in j_range
            (i < 1 || j < 2 || i >= j) && continue
            hi = j - 1
            hi > n && continue
            idx = i:hi
            length(idx) < min_points && continue

            m, b = linfit(xscan[idx], yscan[idx])
            yhat = m .* xscan[idx] .+ b
            r2   = r2_score(yscan[idx], yhat)

            if isfinite(r2) && r2 > best_r2
                best_r2 = r2
                best = (i=i, j=j, m=m, b=b, r2=r2, idx=collect(idx))
            end
        end

        if best === nothing
            return (best_i=missing, best_j=missing, slope=NaN, intercept=NaN, r2=NaN,
                    n_used=0, used_idx=Int[], x_sorted=xs, y_sorted=ys, yhat_full=fill(NaN, n), from_end=from_end)
        end

        yhat_sorted = best.m .* xs .+ best.b
        used_sorted = from_end ? sort(n .- best.idx .+ 1) : sort(best.idx)

        (best_i=best.i, best_j=best.j, slope=best.m, intercept=best.b, r2=best.r2,
         n_used=length(best.idx), used_idx=used_sorted, x_sorted=xs, y_sorted=ys,
         yhat_full=yhat_sorted, from_end=from_end)
    end

    nothing
end

# ╔═╡ 4b75e7f9-3c76-4b04-9d80-3657f0c7d6ee
begin
    if !has_kbi
        md"""**KBI not available:** `Potential.kbi_G22_from_coords` not found."""
    else
        sizes = Float64.(eval_vec_expr(KBI_params[1]))
        steps = Float64.(eval_vec_expr(KBI_params[2]))
        save_kbi = Bool(KBI_params[3])

        kbi_res = Potential.kbi_G22_from_coords(points_3d_rot_cut; sizes=sizes, steps=steps, shift_to_zero=true)

        L_vals    = kbi_res.L
        step_vals = kbi_res.step
        G22_vals  = kbi_res.G22

        sv_func = (has_sv && isdefined(Potential, :kbi_surface_to_volume)) ?
            (L -> Potential.kbi_surface_to_volume(L, Lz_kbi)) :
            (L -> 2/Lz_kbi + 4/L)

        gr()
        kbi_plot_sv = plot(
            title="Kirkwood-Buff integral",
            xlabel="S/V [nm⁻¹]",
            ylabel="G₂₂ [nm³]",
            titlefontsize=14, guidefontsize=14, tickfontsize=11,
            legend=:bottomleft,
            size=(900, 650),
            left_margin=12Plots.mm, bottom_margin=12Plots.mm, top_margin=8Plots.mm, right_margin=6Plots.mm,
        )

        for st in sort(unique(step_vals))
            idx = findall(step_vals .== st)
            sv  = sv_func.(L_vals[idx])
            g   = G22_vals[idx]
            ord = sortperm(sv)

            plot!(kbi_plot_sv, sv[ord], g[ord]; label="step = $(st) nm")
            scatter!(kbi_plot_sv, sv[ord], g[ord]; label="", markersize=3)
        end

        if save_kbi
            _ = writedlm(kbi_txt_path, hcat(L_vals, step_vals, G22_vals))
            gr()
            _ = savefig(kbi_plot_sv, kbi_pdf_path)
        end

        kbi_plot_sv
    end
end

# ╔═╡ 1a10e4d8-6f52-4ff1-8b35-87c5d92f40b6
begin
    function kbi_fit_inputs(
        vars::Vector{String},
        steps_available;
        i_range_default="1:10",
        j_range_default="30:40",
        min_pts_default="8",
        from_end_default=true,
        step_pick_default=nothing
    )
        opts = string.(steps_available)
        if step_pick_default === nothing || !(step_pick_default in opts)
            step_pick_default = isempty(opts) ? "" : first(opts)
        end

        PlutoUI.combine() do Child
            md"""
**Auto-fit window**

Start index range **i_range**:  
$(Child(vars[1], TextField(default=i_range_default)))

Stop index range **j_range** (exclusive):  
$(Child(vars[2], TextField(default=j_range_default)))

Minimum points in window:  
$(Child(vars[3], TextField(default=min_pts_default)))

Count indices from the end:  
$(Child(vars[4], CheckBox(default=from_end_default)))

Plot step (nm):  
$(Child(vars[5], Select(opts; default=step_pick_default)))
"""
        end
    end

    steps_available = (@isdefined(step_vals) && !isempty(step_vals)) ? sort(unique(step_vals)) : Float64[]

    if isempty(steps_available)
        md"Run KBI first to enable fit options."
    else
        let prev = @isdefined(kbi_fit_params) ? kbi_fit_params :
                   ["1:10", "30:40", "8", true, string(first(steps_available))]

            let opts = string.(steps_available)
                step_default_str = string(prev[5])
                step_default_str = (step_default_str in opts) ? step_default_str : first(opts)

                @bind kbi_fit_params confirm(kbi_fit_inputs(
                    ["i_range_expr", "j_range_expr", "min_pts_str", "from_end_cb", "step_pick_str"],
                    steps_available;
                    i_range_default=string(prev[1]),
                    j_range_default=string(prev[2]),
                    min_pts_default=string(prev[3]),
                    from_end_default=Bool(prev[4]),
                    step_pick_default=step_default_str
                ))
            end
        end
    end
end

# ╔═╡ 1f215fbb-7b63-42aa-9108-12f487c62b83
begin
    kbi_fit_summary_md = md""

    if !(@isdefined(kbi_fit_params)) || !(@isdefined(L_vals)) || !(@isdefined(step_vals)) || !(@isdefined(G22_vals))
        kbi_fit_summary_md = md"Run KBI first."
        md"Run KBI first."
    else
        i_range   = eval_int_vec_expr(kbi_fit_params[1])
        j_range   = eval_int_vec_expr(kbi_fit_params[2])
        min_pts   = parse(Int, kbi_fit_params[3])
        from_end  = Bool(kbi_fit_params[4])
        step_pick = parse(Float64, kbi_fit_params[5])

        idx = findall(step_vals .== step_pick)
        Ls  = L_vals[idx]
        Gs  = G22_vals[idx]

        sv = (has_sv && isdefined(Potential, :kbi_surface_to_volume)) ?
            Potential.kbi_surface_to_volume.(Ls, Lz_kbi) :
            (2/Lz_kbi .+ 4 ./ Ls)

        fit = best_window_fit(sv, Gs; i_range=i_range, j_range=j_range, min_points=min_pts, from_end=from_end)

        gr()
        pfit = plot(
            title="KBI fit (maximize R²) — step=$(step_pick) nm",
            xlabel="S/V [nm⁻¹]",
            ylabel="G₂₂ [nm³]",
            titlefontsize=TITLE_FS,
            guidefontsize=GUIDE_FS,
            tickfontsize=TICK_FS,
            legend=false,
            size=(900, 650),
            left_margin=12Plots.mm, bottom_margin=12Plots.mm, top_margin=8Plots.mm, right_margin=6Plots.mm,
        )

        scatter!(pfit, fit.x_sorted, fit.y_sorted; markersize=3)
        if fit.n_used > 0
            scatter!(pfit, fit.x_sorted[fit.used_idx], fit.y_sorted[fit.used_idx]; markersize=6)
        end
        if isfinite(fit.slope) && isfinite(fit.intercept)
            plot!(pfit, fit.x_sorted, fit.yhat_full; lw=2)
        end

        kbi_fit_summary_md = md"""
**Fit summary**

| parameter | value |
|---|---:|
| step (nm) | **$(step_pick)** |
| from end | **$(from_end)** |
| best i | **$(fit.best_i)** |
| best j (exclusive) | **$(fit.best_j)** |
| points used | **$(fit.n_used)** |
| R² | **$(round(fit.r2, sigdigits=4))** |
| **G₂₂,0 (intercept)** | **$(round(fit.intercept, digits=1)) nm³** |
| slope | **$(round(fit.slope, digits=1))** |
"""

        if Bool(KBI_params[3])
            gr()
            _ = savefig(pfit, kbi_svfit_pdf_path)
        end

        pfit
    end
end

# ╔═╡ 1b4da5d1-9b20-4d88-ae20-83fc1a2a6d2b
kbi_fit_summary_md

# ╔═╡ 08c6f524-9a35-44d8-91d7-3d1fa7a6a1e7
md"""
## 7) Output

- Raw uploaded file copy: `Raw_*`
- Cropped coordinates: saved as the original filename in the output folder
- RDF: `*_RDF.txt` and `*_Gr.pdf`
- Structure factor: `*_Sq.txt` and `*_Sq.pdf`
- KBI: `*_KBI.txt`, `*_KBI.pdf`, and `*_KBI_SVfit.pdf`
"""

# ╔═╡ Cell order:
# ╟─53f83c0f-2d6d-4aa9-8f77-84f0903ebebd
# ╟─8e7fa67f-85ae-4e66-aaa3-cb27df2aaa2d
# ╟─9fd66a2c-dc50-4dfe-bc64-76b8a57e051c
# ╟─dc937899-81ef-4326-8acd-94a09e8cee0d
# ╟─9266fa4c-6757-4d8d-837c-06ffb8f7fd6c
# ╟─0f8dcadf-e874-4c0f-bca5-b6063d538134
# ╟─13cc1e0c-07f9-4d73-97a3-8cfcfa7d902a
# ╟─bff7c13d-c197-4a0a-88d7-4d174d2c5f6c
# ╟─2cf2a7c5-0b7c-40df-a1f5-8d2cbd0ddc8d
# ╟─4b8cc26f-d7d6-4c8c-8ad5-4d3fbbf605f1
# ╟─5b3b0edb-3d84-46e9-9c7c-0d7b0b0e5b35
# ╟─1f9ad9f0-4f5e-4a35-8b8f-f5c03d7b1c67
# ╟─33a2df48-7e7d-4f49-8d59-28c17c7c529e
# ╟─3c6bd3c2-45a4-493b-8f2e-4fa6de7fcbce
# ╟─4d6714b1-2d73-4b4b-a901-1d2c83cfe63f
# ╟─7a80e2a6-2c96-4cf3-bc55-43c45f36ab06
# ╟─6a516c61-4b34-4b9b-a4be-ec12e0efc1b7
# ╟─0c865bd7-9e54-4e56-9327-e5d7239a2c2b
# ╟─2a7f6200-1a16-43f3-8a90-fb992fa0d8a2
# ╟─7b0109f6-2cf8-4cf5-acde-8427df12b545
# ╟─59ec80a7-689f-4258-b73b-a25aec2a9880
# ╟─5c375e8a-8c35-4a3a-9f88-56b5f9a01b69
# ╟─64db4a58-c750-4a75-84ad-323f9f5ec593
# ╟─0f73b3bc-25c0-4c33-b082-6b4d2d3c4518
# ╟─6eaa7bb3-2320-43ad-9c3e-b1c2e66fbd42
# ╟─2b1b55a9-06a8-486f-8b1b-79f11b8bc496
# ╟─2e21cccf-9860-4f52-9a0d-9f7c5f229bea
# ╟─88d4a3ea-2f6f-4e48-a7f0-08c7e50a31b0
# ╟─947c3043-5c0d-4af4-8a62-a9ee23f5c1f9
# ╟─9a0b0e4c-ef77-49e4-bc7f-4e97a27c7e6e
# ╟─5e8d1b43-fc7c-4c53-9da6-1877c50fd7cc
# ╟─7a77a3f2-0b8c-47b9-87a5-1111e0a8c0b6
# ╟─24f0e9ea-95f6-4f3c-87d2-0c2d9e7b8b79
# ╟─4ec7d1a8-2d54-4ed1-92fb-4a5be6b7d2fd
# ╟─5e6e2bfa-c2d3-4b76-a2c4-4a7f81d86f4a
# ╟─1c80f764-9d2d-4dff-9a9a-1771169f788f
# ╟─1166aef4-3c77-498b-8eb3-a2066f2a1fe5
# ╟─4b75e7f9-3c76-4b04-9d80-3657f0c7d6ee
# ╟─1a10e4d8-6f52-4ff1-8b35-87c5d92f40b6
# ╟─1f215fbb-7b63-42aa-9108-12f487c62b83
# ╟─1b4da5d1-9b20-4d88-ae20-83fc1a2a6d2b
# ╟─08c6f524-9a35-44d8-91d7-3d1fa7a6a1e7
