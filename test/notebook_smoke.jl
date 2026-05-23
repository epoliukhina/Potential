ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

notebook_path = normpath(joinpath(@__DIR__, "..", "example", "cryoPOT_pluto_notebook.jl"))

@info "Running Pluto notebook smoke test" notebook_path
include(notebook_path)
@info "Pluto notebook smoke test completed"
