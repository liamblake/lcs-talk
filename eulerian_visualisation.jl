using CairoMakie
using LaTeXStrings

function eulerian_animation(glon, glat, time, u, v, outdir; kwargs...)

    lon = glon[1, :]
    lat = glat[:, 1]

    fig, ax, hm = heatmap(lat, lon, sqrt.(u.(glon, glat, 0.0) .^ 2 + v.(glon, glat, 0.0) .^ 2)', axis=(; title=L"u\,\left(x,\,0\right)"); kwargs...)
    hidedecorations!(ax)

    # Animation itself
    record(fig, outdir, time[2:end]; framerate=15) do t
        hm[3] = sqrt.(u.(glon, glat, t) .^ 2 + v.(glon, glat, t) .^ 2)'
        ax.title = L"u\,\left(x,\,%$(t)\right)"
    end
end