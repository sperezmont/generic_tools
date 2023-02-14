# =============================
#     Program: generic_ins.jl
#     Aim: This program contains functions to calculate insolation
#     Author: Sergio Pérez-Montero (github/sperezmont)
# =============================
@doc """
    calc_artificial_insolation:
        Compute daily average insolation through linear composition
"""
function calc_artificial_insolation(time::Real ;ins_min=425.0, ins_max=565.0, P_obl=0.7, P_pre=0.2, P_exc=0.1, tau_obl=41e3, tau_pre=21e3, tau_exc=100e3)
    # Calculate reference and amplitude
    ins_ref = (ins_max + ins_min) / 2
    A_ins = (ins_max - ins_min) / 2

    # Return artificial insolation
    return ins_ref + A_ins * ( P_obl * cos(2.0 * pi * time / tau_obl) + P_pre * cos(2.0 * pi * time / tau_pre) + P_exc * cos(2.0 * pi * time / tau_exc))
end

@doc """
    calc_solar_longitude:
        Compute solar longitude given day of the year and orbital parameters
        Adapted from The Climate Laboratory of Brian E. J. Rose (https://brian-rose.github.io/ClimateLaboratoryBook/home.html)
"""
function calc_solar_longitude(day, long_peri, ecc; days_per_year=365.2422)
    long_peri_rad = deg2rad(long_peri)
    delta_omega = (day - 80.0) * 2 * pi / days_per_year
    beta = sqrt(1 - ecc^2)

    omega_long_m = -2 * ((ecc / 2 + (ecc^3) / 8) * (1 + beta) * sin(-long_peri_rad) -
                         (ecc^2) / 4 * (1 / 2 + beta) * sin(-2 * long_peri_rad) +
                         (ecc^3) / 8 * (1 / 3 + beta) * sin(-3 * long_peri_rad)) + delta_omega

    omega_long = omega_long_m + (2 * ecc - (ecc^3) / 4) * sin(omega_long_m - long_peri_rad) +
                 (5 / 4) * (ecc^2) * sin(2 * (omega_long_m - long_peri_rad)) +
                 (13 / 12) * (ecc^3) * sin(3 * (omega_long_m - long_peri_rad))
    return omega_long
end

@doc """
    calc_real_insolation: 
        Compute daily average insolation given latitude and day of the year
        Adapted from The Climate Laboratory of Brian E. J. Rose (https://brian-rose.github.io/ClimateLaboratoryBook/home.html)
"""
function calc_real_insolation(t::Real; lat=65.0, day=170.0, S0=1365.2, day_type=1, days_per_year=365.2422, orbParams=Dict())
    # First, calculate the orbital parameters at t (years) since J2000 epoch
    if orbParams == Dict()
    	# this option needs Insolation.jl (rad, rad, --), only valid for -5 Myr to 5 Myr:
    	# Pkg.add("Insolation.jl")
    	# using Insolation
    	long_peri, obliquity, ecc = orbital_params(t) 
    else
    	long_peri, obliquity, ecc = orbParams["long_peri"], orbParams["obliquity"], orbParams["ecc"]  
    end
    
    long_peri, obliquity = rad2deg(long_peri), rad2deg(obliquity)
    phi = deg2rad(lat)

    # -- calculate solar longitude
    if day_type == 1
        omega_long = calc_solar_longitude(day, long_peri, ecc, days_per_year=days_per_year)
    elseif day_type == 2
        omega_long = deg2rad(day)
    end

    # -- declination angle of the Sun
    delta = asin(sin(deg2rad(obliquity)) * sin(omega_long))

    # -- hour angle at sunrise / sunset
    if abs(delta) - pi / 2 + abs(phi) < 0.0 # there is sunset/sunrise
        Ho = acos(-tan(phi) * tan(delta))
    else # otherwise figure out if it's all night or all day
        if phi * delta > 0.0
            Ho = pi
        else
            Ho = 0.0
        end
    end

    # -- cos zenith angle
    coszen = Ho * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(Ho)

    # Compute insolation as in Berger 1978 (equation 10)
    Fsw = S0 / pi * ((1 + ecc * cos(omega_long - deg2rad(long_peri)))^2 / (1 - ecc^2)^2 * coszen)
    return Fsw
end



