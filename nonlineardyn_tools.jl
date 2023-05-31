"""
    calc_recurrence(d; tol)
calculates the recurrence matrix of a time series `d`
"""
function calc_recurrence(d::Vector; tol::Real=0.1)
    R = Matrix{Any}(undef, length(d), length(d))
    for i in eachindex(d), j in eachindex(d)
        (abs(d[i] - d[j]) <= tol) ? (R[i, j] = 1.0) : (R[i, j] = 0.0) 
    end
    return R
end
