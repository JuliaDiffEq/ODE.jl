# ODE_MS Fixed-step, fixed-order multi-step numerical method
#   with Adams-Bashforth-Moulton coefficients
function ode_ms{Ty,T}(F, x0::Ty, tspan::AbstractVector{T}, order::Integer; kargs...)

    if !isleaftype(T)
        error("The output times have to be of a concrete type.")
    elseif !(T <:AbstractFloat)
        error("The time variable should be a floating point number.")
    end

    if !isleaftype(Ty) & !isleaftype(eltype(Ty))
        error("The initial data has to be of a concrete type (or an array)")
    end

    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    if 1 <= order <= 4
        b = ms_coefficients4
    else
        b = zeros(order, order)
        b[1:4, 1:4] = ms_coefficients4
        for s = 5:order
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * polyval(p_int, 1))
            end
        end
    end

    # TODO: use a better data structure here (should be an order-element circ buffer)
    xdot = similar(x)
    for i = 1:length(tspan)-1
        # Need to run the first several steps at reduced order
        steporder = min(i, order)
        xdot[i] = F(tspan[i], x[i])

        x[i+1] = x[i]
        for j = 1:steporder
            x[i+1] += h[i]*b[steporder, j]*xdot[i-(steporder-1) + (j-1)]
        end
    end
    return vcat(tspan), x
end

# Use order 4 by default
ode4ms(F, x0, tspan; kargs...) = ode_ms(F, x0, tspan, 4; kargs...)
ode5ms(F, x0, tspan; kargs...) = ODE.ode_ms(F, x0, tspan, 5; kargs...)
