struct EarthModel
    a::Float64
    inv_f::Float64
end

function e(self::EarthModel)
    f = 1.0 / self.inv_f
    return 2 * f - f^2
end

WGS84 = EarthModel(6378187.0, 298.257223563);

function arc_to_meridonal(self::EarthModel, φ)
    e² = e(self)^2
    return @. 180 / pi * (1 - e²) / self.a * (1 - e² * sind(φ)^2)^(3 / 2)

end

function arc_to_parallel(self::EarthModel, φ)
    e² = e(self)^2
    return @. 180 / pi * (1 - e²) / self.a * (1 - e² * sind.(φ)^2)^(1 / 2) / cosd.(φ)

end