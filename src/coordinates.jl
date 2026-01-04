"""
    linear_fit(x, y) -> (m, c)

Least-squares fit of `y ≈ m*x + c`.
"""
function linear_fit(x::AbstractVector, y::AbstractVector)
    A = hcat(x, ones(length(x)))
    m, c = A \ y
    return m, c
end

# Backward-compatible helper: Nx2 matrix input
linear_fit(points::AbstractMatrix) = linear_fit(view(points, :, 1), view(points, :, 2))

"""
    rotate_particles(particles; return_meta=false)

Rotate an `N×3` coordinate set by estimating tilt from linear fits of (x,z) and (y,z).
Returns rotated particles. If `return_meta=true`, also returns angles and rotation matrices.
"""
function rotate_particles(particles::AbstractMatrix{<:Real}; return_meta::Bool=false)
    size(particles, 2) == 3 || error("rotate_particles expects an N×3 matrix")

    x = Float64.(view(particles, :, 1))
    y = Float64.(view(particles, :, 2))
    z = Float64.(view(particles, :, 3))

    kx, bx = linear_fit(x, z)
    ky, by = linear_fit(y, z)

    phi_xz = atan(kx)
    phi_yz = -atan(ky)

    Ux = [1.0 0.0 0.0;
          0.0 cos(phi_yz) -sin(phi_yz);
          0.0 sin(phi_yz)  cos(phi_yz)]

    Uy = [ cos(phi_xz) 0.0 sin(phi_xz);
           0.0        1.0 0.0;
          -sin(phi_xz) 0.0 cos(phi_xz)]

    pts = hcat(x, y, z)
    rotated = (Ux * Uy * pts')'

    if return_meta
        return (particles = rotated,
                phi_xz = phi_xz, phi_yz = phi_yz,
                kx = kx, bx = bx, ky = ky, by = by,
                Ux = Ux, Uy = Uy)
    else
        return rotated
    end
end
