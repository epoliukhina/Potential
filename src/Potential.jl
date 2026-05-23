module Potential

using LinearAlgebra
using Logging
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

const _TESTED_JULIA_VERSIONS = (v"1.12",)
const _DID_WARN_UNTESTED_JULIA = Ref(false)

function _julia_minor_version(version::VersionNumber)
    return VersionNumber(version.major, version.minor)
end

function _julia_minor_version_text(version::VersionNumber)
    return "$(version.major).$(version.minor)"
end

function _is_tested_julia_version(version::VersionNumber = VERSION)
    minor_version = _julia_minor_version(version)
    return any(==(minor_version), _TESTED_JULIA_VERSIONS)
end

function _tested_julia_versions_text()
    return join(_julia_minor_version_text.(_TESTED_JULIA_VERSIONS), ", ")
end

function _untested_julia_warning_message(version::VersionNumber = VERSION)
    current = _julia_minor_version_text(version)
    tested = _tested_julia_versions_text()
    fallback = _julia_minor_version_text(first(_TESTED_JULIA_VERSIONS))
    return "Potential.jl has not been tested on Julia $(current). CI tested on Julia $(tested). If you hit install, precompile, or runtime issues, try Julia $(fallback) or open an issue with your Julia version and platform."
end

function _warn_if_untested_julia_version(version::VersionNumber = VERSION)
    _is_tested_julia_version(version) && return nothing
    if version == VERSION && _DID_WARN_UNTESTED_JULIA[]
        return nothing
    end
    @warn _untested_julia_warning_message(version)
    version == VERSION && (_DID_WARN_UNTESTED_JULIA[] = true)
    return nothing
end

_warn_if_untested_julia_version()

function __init__()
    _warn_if_untested_julia_version()
    return nothing
end

end # module
