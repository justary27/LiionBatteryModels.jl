"""
Simultaneous ODE solver based on the Runge-Kutta 4ᵗʰ order method.
"""
function rk4(df, X0::Vector{Float64}, t::Vector{Float64}, tspan, args...)

    # Time series length
    N = length(t)

    # Number of equations in f
    e_count = length(X0)

    # Initialising the output data vector
    Ypred = zeros(e_count, N)

    # Preparing for calculation
    Ypred[:, 1] = X0
    h = tspan

    # Constructing result
    for i in 1 : N-1
        args = (args..., i)
        k1 = h*df(t[i], Ypred[:, i], args...)
        k2 = h*df(t[i] + h/2, Ypred[:, i] + k1/2, args...)
        k3 = h*df(t[i] + h/2, Ypred[:, i] + k2/2, args...)
        k4 = h*df(t[i] + h, Ypred[:, i] + k3, args...)

        Ypred[:, i+1] = Ypred[:, i] + 1/6*(k1 + 2*k2 + 2*k3 + k4)

    end

    return Ypred
end