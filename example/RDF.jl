using Plots
using Potential

"""
This script calculates and plots g(r) using the analytic
normalization procedure
"""

#DASHBOARD
RMax = 50    # maximal radial distance
DeltaR = 0.3 # bin width
NBins = 100  # number of equally spaced bins

factor = 1.0 # ? 
particle_size  = 6.0 # ?

dir_name    = "./data/"
init_data =  "Ferritin_PBS_30mgml_position1.dat"
sim_name = init_data[1:end-4]

println("Current Simulation Name: ", sim_name)

## Read the particle coordinates
Particles = readPartCoor(dir_name*init_data; factor=factor)

## Perform rotation of data 
Particles_rot = rotate_particles(Particles)

##create the list with radial distances
rList = range(0, stop=RMax, length=NBins)

println("Computing RDF:")
ProcessPartCoor(Particles, rList, DeltaR, dir_name, sim_name)
 
println("Program Finished")





