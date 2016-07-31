"""

We assume that the initial data y0 is given at tspan[1], and that
tspan[end] is the last integration time.

"""
function ode{T,Y,M<:AbstractSolver}(F, y0::Y,
                                    tout::AbstractVector{T};
                                    solver::Type{M} = RKIntegratorAdaptive{:rk45},
                                    points = :all,
                                    kargs...)

    t0 = tout[1]

    # construct a Problem
    equation  = explicit_ineff(t0,y0,F;kargs...)
    if points == :all
        prob = solve(equation, solver;
                     tout = tout,
                     kargs...)
    elseif points == :specified
        prob = solve(equation,
                     DenseOutput{solver};
                     tout = tout,
                     kargs...)
    else
        error("Unsupported points value (should be :all or :specified)")
    end

    # determine if we have to unpack y
    extract = Y <: Number

    to = Array(T,0)
    yo = Array(Y,0)
    for (t,y) in prob
        push!(to,t)
        push!(yo, extract ? y[1] : copy(y))
    end

    return (to,yo)
end

"""
    ODE.odeXX(F,y0,t0;kargs...)

Solves an ODE `y'=F(t,y)` with initial conditions `y0` and `t0`.
"""

ode23s(F,y0,t0;kargs...)        = ode_conv(F,y0,t0;solver = ModifiedRosenbrockIntegrator, kargs...)
ode1(F,y0,t0;kargs...)          = ode_conv(F,y0,t0;solver = RKIntegratorFixed{:feuler}, kargs...)
ode2_midpoint(F,y0,t0;kargs...) = ode_conv(F,y0,t0;solver = RKIntegratorFixed{:midpoint}, kargs...)
ode2_heun(F,y0,t0;kargs...)     = ode_conv(F,y0,t0;solver = RKIntegratorFixed{:heun}, kargs...)
ode4(F,y0,t0;kargs...)          = ode_conv(F,y0,t0;solver = RKIntegratorFixed{:rk4}, kargs...)
ode21(F,y0,t0;kargs...)         = ode_conv(F,y0,t0;solver = RKIntegratorAdaptive{:rk21}, kargs...)
ode23(F,y0,t0;kargs...)         = ode_conv(F,y0,t0;solver = RKIntegratorAdaptive{:rk23}, kargs...)
ode45_fe(F,y0,t0;kargs...)      = ode_conv(F,y0,t0;solver = RKIntegratorAdaptive{:rk45}, kargs...)
ode45_dp(F,y0,t0;kargs...)      = ode_conv(F,y0,t0;solver = RKIntegratorAdaptive{:dopri5}, kargs...)
const ode45 = ode45_dp
ode78(F,y0,t0;kargs...)         = ode_conv(F,y0,t0;solver = RKIntegratorAdaptive{:feh78}, kargs...)


function ode_conv{Ty,T}(F,y0::Ty,t0::AbstractVector{T};kargs...)

    if !isleaftype(T)
        error("The output times have to be of a concrete type.")
    elseif !(T <:AbstractFloat)
        error("The time variable should be a floating point number.")
    end

    if !isleaftype(Ty) & !isleaftype(eltype(Ty))
        error("The initial data has to be of a concrete type (or an array)")
    end

    ode(F,y0,t0;kargs...)

end



"""

Convert a out-of-place explicitly defined ODE function to
ExplicitODE.  As the name suggests, the result is not going to be very
efficient.

"""
function explicit_ineff{T,Y}(t0::T, y0::AbstractArray{Y}, F::Function; kargs...)
    F!(t,y,dy) =copy!(dy,F(t,y))
    return ExplicitODE(t0,y0,F!; kargs...)
end

# A temporary solution for handling scalars, should be faster then the
# previous implementation.  Should be used only at the top level
# interface.  This function cheats by converting scalar functions F
# and jac to vector functions F! and jac!.  Still, solving this ODE
# will result in a vector of length one result, so additional external
# conversion is necessary.
function explicit_ineff{T,Y}(t0::T, y0::Y, F::Function; kargs...)
    F!(t,y,dy) =(dy[1]=F(t,y[1]))
    new_y0 = Array(Y,1)
    new_y0[1] = y0
    return ExplicitODE(t0,new_y0,F!; kargs...)
end
