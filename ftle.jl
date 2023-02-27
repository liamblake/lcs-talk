using LinearAlgebra

using DifferentialEquations


function ftle(u, x₀, t₀, T, δx; kwargs...)
    n = size(x₀)[1]
    d = size(x₀)[2]

    # Form a star grid around each point
    # TODO: Vectorise this
    # Coordinate shifts for computing the star grid
    A = zeros(2 * d, d)
    A[1:(2*d+2):(2*d^2)] .= 1
    A[2:(2*d+2):(2*d^2)] .= -1

    # Index as (initial condition star point, dim)
    grid = zeros(n, 2 * d, d)
    for (i, x) in enumerate(eachrow(x₀))
        grid[i, :, :] = repeat(x', 2 * d) + δx * A
    end

    # Solve forward as an Ensemble ODEProblem
    # Advect these points forwards to the initial time
    # Mapping between grid and cartesian index i: sol[i] ⟺ grid[cld(i, 2d), mod1(i, 2d), :]
    advected = Array{Float64}(undef, n, 2 * d, d)
    prob = ODEProblem(u, grid[1, 1, :], (t₀, T))
    # Hack the output function to directly place the solution in advected. Can ignore the output of solve.
    ensemble = EnsembleProblem(
        prob;
        prob_func=(prob, i, _) -> remake(prob; u0=grid[cld(i, 2 * d), mod1(i, 2 * d), :]),
        output_func=function (sol, i)
            advected[cld(i, 2 * d), mod1(i, 2 * d), :] = sol[:, 2]
            return (nothing, false)
        end
    )
    _ = solve(ensemble, EnsembleThreads(); save_everystep=false, trajectories=2 * d * n)

    # Compute each flow map gradient, at each time step
    ∇F = Vector{Matrix}(undef, n)
    for i = 1:n
        ∇F[i] = 1 / (2 * δx) * ones(d, d)
        for l = 1:d
            for m = 1:d
                ∇F[i][l, m] *= advected[i, 2*m-1, l] - advected[i, 2*m, l]
            end
        end
    end

    # Compute FTLE
    return 1 / T * log.(opnorm.(∇F))

end