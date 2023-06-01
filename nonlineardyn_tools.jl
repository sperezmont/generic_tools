using LinearAlgebra

"""
    calc_recurrence(u; tol)
calculates the recurrence matrix of the vector `u` that represents 
the time series of the state of a dynamic process
"""
function calc_recurrence(u::Vector; tol::Real=0.1)
    R = Matrix{Any}(undef, length(u), length(u))
    for i in eachindex(u), j in eachindex(u)
        (norm(u[i] - u[j]) <= tol) ? (R[i, j] = 1.0) : (R[i, j] = 0.0) 
    end
    return R
end
