
function gauss(N::Int, a, b)
    λ, Q = eigen(SymTridiagonal(zeros(N), [n / sqrt(4n^2 - 1) for n ∈ 1:N-1]))
    @. (λ + 1) * (b - a) / 2 + a, [2Q[1, i]^2 for i ∈ 1:N] * (b - a) / 2
end


function integrate(N::Int, lower, upper)
    x, w = gauss(N, lower, upper)
    function itg(func, args...; fields)
        ret = map(a -> func(a, args...), x)
        map(f -> sum(map(r -> getproperty(r, f), ret) .* w), fields)
    end
end


solar_declination(doy::Int) = -0.4093 * cos(0.0172 * (doy + 10))


solar_constant(doy::Int) = 1370 * (1 + 0.033 * cos(0.0172 * (doy - 10)))


et_solar_radiation(sc, inc) = max(0.0, sc * cos(inc))


svp_fn(ta) = 6.1078 * exp(17.269 * ta / (ta + 237.3))


saturated_vapor_pressure(ta) = svp_fn(ta)


vapor_pressure(ta, tdew) = svp_fn(min(ta, tdew))


relative_humidity(svp, vp) = 100 * vp / svp


function _ab(lat, dec)
    a = sin(lat) * sin(dec)
    b = cos(lat) * cos(dec)
    (; a, b)
end


function twilight_hours(lat, dec)
    ab = _ab(lat, dec)
    tss = 12 + (12 / π) * acos(-ab.a / ab.b)  # sunset
    tsr = 24 - tss  # sunrise
    tsr, tss
end


function solar_position(th, lat, dec)
    ab = _ab(lat, dec)
    ha = π / 12 * (th - 12)  # hour angle

    inc = min(0.5 * π, acos(ab.a + ab.b * cos(ha)))  # solar inclination
    hgt = 0.5 * π - inc  # solar height/elevation

    # azimuth (angle from North in a clockwise direction):
    a = cos(dec) * (cos(lat) * tan(dec) + sin(lat) * cos(ha)) / cos(hgt)
    acosa = acos(max(-1, min(1, a)))
    azi = th <= 12 ? acosa : (2π - acosa)

    inc, hgt, azi
end


function day_et_solar_radiation(doy::Int, lat, dec)
    ab = _ab(lat, dec)
    aob = ab.a / ab.b
    sc = solar_constant(doy)
    0.027501974 * sc * (ab.a * acos(-aob) + ab.b * sqrt(1 - aob^2))
end


function air_temperature(th, tmin, tmax, tsr, tss, lag=1.5)
    dl = tss - tsr   # day length
    if (tsr + lag) <= th <= tss
        n1 = π * (th - tsr - lag) / dl
        ta = tmin + (tmax - tmin) * sin(n1)
    else
        tset = tmin + (tmax - tmin) * sin(π * (dl - lag) / dl)
        t = th < (tsr + lag) ? tsr : -tss
        ta = tset + ((tmin - tset) * (th + t)) / ((tsr + lag) + tsr)
    end
    ta
end


function solar_radiation(rh, etrad)
    kt = 1.1595 - 0.0106 * rh  # sky clearness index = total radiation / ET radiation
    it = etrad * kt    # total
    # diffuse partitioning based on Khatib et al. (2012) for 28 Malaysian sites:
    ratio = 0.9505 + 0.91634 * kt - 4.851 * kt^2 + 3.2353 * kt^3
    idf = it * ratio    # diffuse
    idr = it - idf      # direct
    it, idr, idf
end


function daily_solar_radiation(lat)
    function radfn(doy, tmin, tmax, tsr, tss, dec, lat)
        solarcon = solar_constant(doy)
        dew_temp = 23.0

        function hour_wthr(th)  # th = local solar time (hours)
            air_temp = air_temperature(th, tmin, tmax, tsr, tss)
            svp = saturated_vapor_pressure(air_temp)
            vp = vapor_pressure(air_temp, dew_temp)
            rh = relative_humidity(svp, vp)
            solarinc, _, _ = solar_position(th, lat, dec)
            etrad = et_solar_radiation(solarcon, solarinc)
            totrad, drrad, dfrad = solar_radiation(rh, etrad)
            (; totrad, drrad, dfrad)
        end
    end

    function op(doy, tmin, tmax)
        dec = solar_declination(doy)
        tsr, tss = twilight_hours(lat, dec)
        fn = radfn(doy, tmin, tmax, tsr, tss, dec, lat)
        itg = integrate(7, tsr, tss)
        res = itg(fn; fields=[:totrad, :drrad, :dfrad])
        res .* 3600 / 10^6  # in MJ/m2/day
    end
end


function add_solar_radiation!(df::AbstractDataFrame, lat)
    !all(["tmin", "tmax"] .∈ Ref(names(df))) && return
    op = daily_solar_radiation(lat)
    args = [:doy, :tmin, :tmax]
    transform!(df, args => ByRow(op) => [:totrad, :drrad, :dfrad])
end
