using Interpolations
using NCDatasets

include("earth.jl")

DATA_DIR = "data/"
ALL_DATA = [
    "cmems_mod_glo_phy_my_0.083_P1D-m_1674715488827.nc",
    "cmems_mod_glo_phy_my_0.083_P1D-m_1674716401167.nc",
    "cmems_mod_glo_phy_my_0.083_P1D-m_1674717969963.nc",
]

function read_nc_vel(file)
    println("Processing $(file)...")

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

    return lat, lon, time, v, u
end

function read_all_files(; kwargs...)

    # Read the first file, use as a baseline
    lat, lon, time, v, u = read_nc_vel(DATA_DIR * ALL_DATA[1])

    for dfile in ALL_DATA[2:end]
        nlat, nlon, ntime, nv, nu = read_nc_vel(DATA_DIR * dfile)

        # Ensure that the new data is on the same spatial grid
        @assert nlat == lat
        @assert nlon == lon

        time = cat(time, ntime, dims=1)
        u = cat(u, nu, dims=3)
        v = cat(v, nv, dims=3)
    end

    # Replace missing data (land) with zero
    land_u = prod(ismissing.(u), dims=3)
    land_v = prod(ismissing.(u), dims=3)
    land = land_u .* land_v
    u[land] = 0.0
    v[land] = 0.0

    # Rescale time to start from 0
    mil_to_days = s -> s.value / (1000 * 60 * 60 * 24)
    days = mil_to_days.(time .- time[1])

    # Construct interpolations
    u_interp = linear_interpolation((lon, lat, days), u; kwargs...)
    v_interp = linear_interpolation((lon, lat, days), v; kwargs...)

    return lat, lon, days, u_interp, v_interp, land
end


function load_vel_interps()
    return read_all_files(; extrapolate_bc=Line())
end
