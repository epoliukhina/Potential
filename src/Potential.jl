module Potential

using LinearAlgebra
using ProgressMeter
using Plots
using Base.Threads
using Statistics

include("io.jl")
include("coordinates.jl")
include("rdf.jl")
include("kbi.jl")
include("structure_factor.jl")
export structure_factor_debye_histogram, save_structure_factor, plot_structure_factor

export kbi_G22_from_coords, sliding_counts_xy, kbi_surface_to_volume, kbi_fit_G22_vs_sv, kbi_best_window_fit
export read_particle_coords, save_rdf_data, plot_rdf
export rotate_particles, linear_fit
export rdf_analytic_norm, process_particle_coords

# Backwards-compatible aliases (so your notebook keeps working)
const readPartCoor     = read_particle_coords
const saveRDFdata      = save_rdf_data
const plotRDF          = plot_rdf
const RDF_AnalyticNorm = rdf_analytic_norm
const ProcessPartCoor  = process_particle_coords

end # module
