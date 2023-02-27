using Interpolations
using NCDatasets

include("earth.jl")


function read_nc_vel(file)
    # println("Processing $(file)...")

    nc = Dataset(file)
    lat = nc["latitude"][:, :]
    lon = nc["longitude"][:, :]
    time = nc["time"][:, :]
    v = nc["vo"][:, :]
    u = nc["uo"][:, :]

    # Drop depth dimension
    v = v[:, :, 1, :]
    u = u[:, :, 1, :]

    # Convert velocities to deg/day
    vsize = size(v)
    v = v .* repeat(arc_to_parallel(WGS84, lat)'; outer=(vsize[1], 1, vsize[3])) * 86400
    u = u .* repeat(arc_to_meridonal(WGS84, lat)'; outer=(vsize[1], 1, vsize[3])) * 86400

    # Replace missing data (land) with zero
    land_u = ismissing.(u)
    land_v = ismissing.(u)
    land = prod(land_u .|| land_v, dims=3)[:, :, 1]
    u[land_u] .= 0.0
    v[land_v] .= 0.0

    # Rescale time to start from 0
    mil_to_days = s -> s.value / (1000 * 60 * 60 * 24)
    days = mil_to_days.(time .- time[1])

    # Construct interpolations
    u_interp = linear_interpolation((lon, lat, days), u, extrapolation_bc=0.0)
    v_interp = linear_interpolation((lon, lat, days), v, extrapolation_bc=0.0)

    return lat, lon, days, u_interp, v_interp, land
end
