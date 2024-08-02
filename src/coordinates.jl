#Functions for rotation and processing of coordinates

# import Pkg;Pkg.add("Plotly")
# particle_size = parse(Float64, size_str)
# factor = parse(Float64, factor_str) 
# data_matrix = factor .* readdlm(ParticleFileRaw)
# x,y,z = data_matrix[:,1],data_matrix[:,2],data_matrix[:,3]
# plotlyjs() 
# Plots.scatter3d(x,y,z, markersize = particle_size/4)


function linear_fit(points)
    # construct the matrix A and the vector b
    A = [points[:, 1] ones(size(points, 1))]
    b = points[:, 2]

    # use the least squares method to find the line parameters
    m, c = A \ b

    return m, c
end

function rotate_particles(particles)
    x,y,z = particles[:,1],particles[:,2],particles[:,3]

    points_3d = hcat(x, y, z)
    points = hcat(x, z)
    kx, bx = linear_fit(points)
    phi_xz = atan(kx)

    fx(x) = kx * x + bx

    points = hcat(y, z)
    ky, by = linear_fit(points)
    phi_yz = -atan(ky)

    fy(y) = ky * y + by

    points = hcat(x,y)
    kxy, bxy = linear_fit(points)
    phi_xy= -atan(kxy)

    fxy(x) = kxy*x + bxy	

    U_x = [1.0 0.0 0.0;
            0.0 cos(phi_yz) -sin(phi_yz);
            0.0 sin(phi_yz) cos(phi_yz)]

    U_y = [cos(phi_xz) 0.0 sin(phi_xz);
            0.0 1.0 0.0;
            -sin(phi_xz) 0.0 cos(phi_xz)]

    points_3d_new = (U_x * U_y * points_3d')'

    # x_new = points_3d_new[:, 1]
    # y_new = points_3d_new[:, 2]
    # z_new = points_3d_new[:, 3]

    return points_3d_new
end

