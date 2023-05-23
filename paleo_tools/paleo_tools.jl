"""
	convert_d18O_to_temp(d180; a)
converts d18O data to temperature anomaly through
            d18O = aT + b       (based on Cuffey et al. (1994, 1995))
## Attributes
`d18O` Vector or value to convert
`a` slope

## Default values
Default value of a is selected on a calibration of Barker et al., 2011 d18O values to -20ºC anomaly in LGM

"""
function convert_d18O_to_temp(d18O::Any; a::Real=0.3)
    b = 31.47 * a - 34.87
    T = (d18O .- b) ./ a 
    return T
end
