import sys
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

from scipy.spatial.distance import pdist

import plotly.express as px
import plotly.graph_objs as go




def linear_fit(points):
    # construct the matrix A and the vector b
    A = np.column_stack((points[:, 0], np.ones(len(points))))
    b = points[:, 1]
    
    # use the least squares method to find the line parameters
    m, c = np.linalg.lstsq(A, b, rcond=None)[0]
    
    return m, c

#Initial values
# "blobs_expBSA_19_Position_1_SIRT20i_1n.inv.csv"
# "test.dat"

try:
    path_to_data, path_to_save = sys.argv[1].split(' ')
    dir_to_save = path_to_save.split('/')[0]
except:
    raise(Exception("Specify path for data"))
#print(path_to_data)
#data_frame = pd.read_csv(path_to_data)
#x = 0.3 * data_frame['axis-2'].values
#y = 0.3 * data_frame['axis-1'].values
#z = 0.3 * data_frame['axis-0'].values


data_frame = np.loadtxt(path_to_data, delimiter = '\t')
x = 0.3 * data_frame[:, 0]
y = 0.3 * data_frame[:, 1]
z = 0.3 * data_frame[:, 2]


points_3d = np.array([x,y,z]).T

fig = px.scatter_3d(x=x,y=y,z=z)
fig.update_traces(marker_size = 2)
fig.show()

points_xz = np.array([x,z]).T
kx, bx = linear_fit(points_xz)
phi_xz= np.arctan(kx)

fx = lambda x : kx*x + bx

points_yz = np.array([y,z]).T
ky, by = linear_fit(points_yz)
phi_yz= - np.arctan(ky)

fy = lambda x : ky*x + by

points_xy = np.array([x,y]).T
kxy, bxy = linear_fit(points_xy)
phi_xy= -np.arctan(kxy)

fxy = lambda x : kxy*x + bxy

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

# Plot the first graph
axes[0].scatter(x, z)
axes[0].scatter(x, fx(x))
axes[0].set_xlabel('X [nm]')
axes[0].set_ylabel('Z [nm]')
axes[0].set_title('Projection on OXZ')

# Plot the second graph
axes[1].scatter(y, z)
axes[1].scatter(y, fy(y))
axes[1].set_xlabel('Y [nm]')
axes[1].set_ylabel('Z [nm]')
axes[1].set_title('Projection on OYZ')

# Plot the third graph
axes[2].scatter(x, y)
axes[2].scatter(x, fxy(x))
axes[2].set_xlabel('X [nm]')
axes[2].set_ylabel('Y [nm]')
axes[2].set_title('Projection on OXY')

# Adjust the spacing between subplots
plt.tight_layout()
plt.savefig(dir_to_save+"/projection_before_rotation.pdf")


U_x = np.array([[1.0, 0.0,             0.0],
                [0.0, np.cos(phi_yz), -np.sin(phi_yz)],
                [0.0, np.sin(phi_yz), np.cos(phi_yz)]])

U_y = np.array([[np.cos(phi_xz), 0.0, np.sin(phi_xz)],
                [0.0,            1.0, 0.0           ],
                [-np.sin(phi_xz),0.0, np.cos(phi_xz)]])
print(f'Theta_XZ = {round(phi_xz/ np.pi * 180, 1)} degr')
print(f'Theta_YZ = {round(phi_yz/ np.pi * 180, 1)} degr')
print(f'Theta_XY = {round(phi_xy/ np.pi * 180, 1)} degr')

points_3d_new = np.dot(np.dot(U_x,U_y),points_3d.T).T

distances = pdist(points_3d)
distances_new = pdist(points_3d_new)

# Print the pairwise distances difference
# print(distances-distances_new)
print(f"Change in pairwise distances after rotation= {sum(distances-distances_new )/len(points_3d)} nm")

x_new = points_3d_new.T[0]
y_new = points_3d_new.T[1]
z_new = points_3d_new.T[2]

fig = px.scatter_3d(x=x_new,y=y_new,z=z_new)
fig.update_traces(marker_size = 2)
fig.show()

# create trace for first scatter plot
trace1 = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=2,
        color='blue',
        opacity=0.5
    ),
    name='Init'
)

# create trace for second scatter plot
trace2 = go.Scatter3d(
    x=x_new,
    y=y_new,
    z=z_new,
    mode='markers',
    marker=dict(
        size=2,
        color='red',
        opacity=0.5
    ),
    name='Rotated'
)

# create layout
layout = go.Layout(
    title='Two Scatter Plots',
    scene=dict(
        xaxis=dict(title='X Axis'),
        yaxis=dict(title='Y Axis'),
        zaxis=dict(title='Z Axis')
    )
)

# create figure and add traces to it
fig = go.Figure(data=[trace1, trace2], layout=layout)

# show the plot
fig.show()

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

# Plot the first graph
axes[0].scatter(x, z)
axes[0].scatter(x_new, z_new)
axes[0].set_xlabel('X [nm]')
axes[0].set_ylabel('Z [nm]')
axes[0].set_title('Projection on OXZ')

# Plot the second graph
axes[1].scatter(y, z)
axes[1].scatter(y_new, z_new)
axes[1].set_xlabel('Y [nm]')
axes[1].set_ylabel('Z [nm]')
axes[1].set_title('Projection on OYZ')

# Plot the third graph
axes[2].scatter(x, y)
axes[2].scatter(x_new, y_new)
axes[2].set_xlabel('X [nm]')
axes[2].set_ylabel('Y [nm]')
axes[2].set_title('Projection on OXY')

# Adjust the spacing between subplots
plt.tight_layout()
plt.savefig(dir_to_save+"/projection_comparison_after_rotation.pdf")

print(f"delta LenX = {round(abs((np.max(x_new) - np.min(x_new))-(np.max(x) - np.min(x))), 1)} nm")
print(f"delta LenY = {round(abs((np.max(y_new) - np.min(y_new))-(np.max(y) - np.min(y))), 1)} nm")
print(f"delta LenZ = {round(abs((np.max(z_new) - np.min(z_new))-(np.max(z) - np.min(z))), 1)} nm")

cut_length = 25
low_lim = np.mean(z_new) - cut_length
high_lim = np.mean(z_new) + cut_length
indices_to_keep = np.where((z_new >= low_lim) & (z_new <= high_lim))

points_3d_new_cut = points_3d_new[indices_to_keep]

x_cut = points_3d_new_cut.T[0]
y_cut = points_3d_new_cut.T[1]
z_cut = points_3d_new_cut.T[2]
fig = px.scatter_3d(x=x_cut,y=y_cut,z=z_cut)
fig.update_traces(marker_size = 2)
fig.show()
                           
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

# Plot the first graph
axes[0].scatter(x_cut, z_cut)
axes[0].set_xlabel('X [nm]')
axes[0].set_ylabel('Z [nm]')
axes[0].set_title('Projection on OXZ')

# Plot the second graph
axes[1].scatter(y_cut, z_cut)
axes[1].set_xlabel('Y [nm]')
axes[1].set_ylabel('Z [nm]')
axes[1].set_title('Projection on OYZ')

# Plot the third graph
axes[2].scatter(x_cut, y_cut)
axes[2].set_xlabel('X [nm]')
axes[2].set_ylabel('Y [nm]')
axes[2].set_title('Projection on OXY')

# Adjust the spacing between subplots
plt.tight_layout()

plt.savefig(dir_to_save+"/projection_after_cutting.pdf")

points_to_save = points_3d_new_cut
V = (np.max(x_cut) - np.min(x_cut))*(np.max(y_cut) - np.min(y_cut))*(np.max(z_cut) - np.min(z_cut))

print(f'Final Number of particles = {len(points_to_save)}')
print(f'Final Number density = {len(points_to_save) / V} nm^(-3)')

np.savetxt(path_to_save,points_to_save)   