## TODO: reactivate when introducing events/rootfinding
# """

# A simple bisection algorithm for finding a root of a solution f(x)=0
# starting within the range x∈rng, the result is a point x₀ which is
# located within the distance eps from the true root of f(x)=0.  For
# this algorithm to work we need f(rng[1]) to have a different sign then
# f(rng[2]).

# """
# function findroot(f,rng,eps)
#     xl, xr = rng
#     fl, fr = f(xl), f(xr)

#     if fl*fr > 0 || xl > xr
#         error("Inconsistent bracket")
#     end

#     while xr-xl > eps
#         xm = (xl+xr)/2
#         fm = f(xm)

#         if fm*fr > 0
#             xr = xm
#             fr = fm
#         else
#             xl = xm
#             fl = fm
#         end
#     end

#     return (xr+xl)/2
# end

# generate a jacobian using ForwardDiff
function forward_jacobian(F,y0::AbstractArray)
    (t,y)->ForwardDiff.jacobian(y->F(t,y),y)
end

function forward_jacobian(F,y0)
    (t,y)->ForwardDiff.derivative(y->F(t,y),y)
end

function forward_jacobian!(F!,tmp)
    jac!(t,y,J)=ForwardDiff.jacobian!(J,(dy,y)->F!(t,y,dy),tmp,y)
    return jac!
end

function forward_jacobian_implicit!(F!,tmp)
    error("Not implemented yet")
end
