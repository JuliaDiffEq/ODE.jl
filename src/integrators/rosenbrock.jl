#ODEROSENBROCK Solve stiff differential equations, Rosenbrock method
#   with provided coefficients.
function oderosenbrock{Ty,T}(F, x0::Ty, tspan::AbstractVector{T},
                             gamma, a, b, c;
                             jacobian = forward_jacobian(F,x0),
                             kargs...)

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

    solstep = 1
    while solstep < length(tspan)
        ts = tspan[solstep]
        hs = h[solstep]
        xs = x[solstep]
        dFdx = jacobian(ts, xs)
        # FIXME
        if size(dFdx,1) == 1
            jac = 1/gamma/hs - dFdx[1]
        else
            jac = eye(dFdx)/gamma/hs - dFdx
        end

        g = Array(typeof(x0), size(a,1))
        g[1] = (jac \ F(ts + b[1]*hs, xs))
        x[solstep+1] = x[solstep] + b[1]*g[1]

        for i = 2:size(a,1)
            dx = zero(x0)
            dF = zero(x0/hs)
            for j = 1:i-1
                dx += a[i,j]*g[j]
                dF += c[i,j]*g[j]
            end
            g[i] = (jac \ (F(ts + b[i]*hs, xs + dx) + dF/hs))
            x[solstep+1] += b[i]*g[i]
        end
        solstep += 1
    end
    return vcat(tspan), x
end


# Kaps-Rentrop coefficients
const kr4_coefficients = (0.231,
                          [0              0             0 0
                           2              0             0 0
                           4.452470820736 4.16352878860 0 0
                           4.452470820736 4.16352878860 0 0],
                          [3.95750374663  4.62489238836 0.617477263873 1.28261294568],
                          [ 0               0                0        0
                           -5.07167533877   0                0        0
                            6.02015272865   0.1597500684673  0        0
                           -1.856343618677 -8.50538085819   -2.08407513602 0],)

ode4s_kr(F, x0, tspan; kargs...) = oderosenbrock(F, x0, tspan, kr4_coefficients...; kargs...)

# Shampine coefficients
const s4_coefficients = (0.5,
                         [ 0    0    0 0
                           2    0    0 0
                          48/25 6/25 0 0
                          48/25 6/25 0 0],
                         [19/9 1/2 25/108 125/108],
                         [   0       0      0   0
                            -8       0      0   0
                           372/25   12/5    0   0
                          -112/125 -54/125 -2/5 0],)

ode4s_s(F, x0, tspan; kargs...) = oderosenbrock(F, x0, tspan, s4_coefficients...; kargs...)

# Use Shampine coefficients by default (matching Numerical Recipes)
const ode4s = ode4s_s

const ms_coefficients4 = [ 1      0      0     0
                          -1/2    3/2    0     0
                          5/12  -4/3  23/12 0
                          -9/24   37/24 -59/24 55/24]
