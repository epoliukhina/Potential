using Plots
using Potentials


"""
This script calculates and plots g(r) using the analytic
normalization procedure
"""

#DASHBOARD
RMax = 50 # maximal radial distance
DeltaR = 0.3 # bin width
NBins = 100 # number of equally spaced bins

dir_name = "data/"

InitailData =  dir_name * "Ferritin_PBS_30mgml_position1.dat"
SimuName = InitailData[1:end-4]

Particles = readPartCoor(dir_name*InitailData)


# # Perform rotation of data 


# #create the list with radial distances
# rList = range(0, stop=RMax, length=NBins)

# println("Computing RDF:")
# ProcessPartCoor(ParticleFile, rList, DeltaR, SimuName)


# println("Current FilePath: ", PartFile)



# # calculate the RDF
# Gr = RDF_AnalyticNorm(Particles, rList, DeltaR)

# # save the RDF data to a file
# saveRDFdata(rList, Gr, DeltaR, SimuName, PartFile)

# # plot RDF
# plotRDF(rList, Gr, SimuName)


# println("Program Finished")





