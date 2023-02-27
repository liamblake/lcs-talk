using LinearAlgebra

"""
	star_grid(x, δx)
Construct a star grid around a point x, for calculating a finite difference
approximation of a spatial derivative. The number of points in the grid is
determined from the dimension of x; if the length of x is n, then 2n points are
calculated.
"""
function star_grid(x, δx)
    n = length(x)

    # Matrix of coordinate shifts
    # ⌈1  0  0  0  ⋯ 0⌉
    # |-1 0  0  0  ⋯ 0|
    # |0  1  0  0  ⋯ 0|
    # |0  -1 0  0  ⋯ 0|
    # |⋮     ⋱       ⋮|
    # |0  ⋯     0  ⋯ 1|
    # ⌊0  ⋯     0  ⋯ 0⌋
    A = zeros(2 * n, n)
    A[1:(2*n+2):(2*n^2)] .= 1
    A[2:(2*n+2):(2*n^2)] .= -1

    return repeat(x', 2 * n) + δx * A
end

"""
    ∇F(star_values, n, δx)
Approximate the flow map gradient with a centered finite-difference approximation, given a star grid
of values.
"""
function ∇F_fd(star_values, n, δx)
    return 1 / (2 * δx) * (star_values[1:2:(2*n), :] - star_values[2:2:(2*n), :])'
end