# This file contains the implementation of explicit Runkge-Kutta
# solver from (Hairer & Wanner 1992 p.134, p.165-169).

###########################################
# Tableaus for explicit Runge-Kutta methods
###########################################
immutable TableauRKExplicit{T} <: Tableau{T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    isFSAL::Bool
    s::Int
    name::String
    function TableauRKExplicit(name,order,a,b,c)
        s = length(c)
        @assert c[1]==0
        @assert istril(a)
        @assert s==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        @assert norm(sum(a,2)-c'',Inf)<T(1e-10) # consistency.
        isFSAL = (a[end,:]==b[1,:] && c[end]==1)
        new(order,a,b,c,isFSAL,s,name)
    end
end

lengthks(tab::TableauRKExplicit) = length(tab.c)

Base.convert{Tnew<:Real,T}(::Type{TableauRKExplicit{Tnew}}, tab::TableauRKExplicit{T}) =
    TableauRKExplicit{Tnew}(tab.name, tab.order,
                            convert(Matrix{Tnew},tab.a),
                            convert(Matrix{Tnew},tab.b),
                            convert(Vector{Tnew},tab.c))

isexplicit(b::TableauRKExplicit) = istril(b.a) # Test whether it's an explicit method
isadaptive(b::TableauRKExplicit) = size(b.b, 1)==2

#######################
# Solver implementation
#######################
"""

A general Runge-Kutta integrator (it can represent either, a fixed step
or an adaptive step algorithm).

"""
immutable RKIntegrator{Kind,Name,T,OP<:Options} <: AbstractIntegrator{T}
    tableau::TableauRKExplicit{T}
    opts::OP
end

typealias RKIntegratorFixed    RKIntegrator{:fixed}
typealias RKIntegratorAdaptive RKIntegrator{:adaptive}

"""

A constructor for an explicit Runge-Kutta method.  Works only for
explicit differential equations.

Notes:

- `Kind` is either `:adaptive` or `:fixed`, corresponding to adaptive
  step size method or a fixed step size method

- `Name` is the name of a Butcher tableau based on which the method is
  constructed.  The kind (adaptive or fixed) of the Butcher tableau
  has to correspond to `Kind` (`:adaptive` or `:fixed`).

Input:

- `ode::ExplicitODE`

- `opts` options for the method, supports the same basic options as
  other adaptive steppers (see `AdaptiveOptions` for the complete
  list).

Output:

- `::RKIntegrator{Kind,Name}`

"""
@compat function (::Type{RKIntegrator{Kind,Name}}){Kind,Name,T}(ode::ExplicitODE{T};opts...)
    tab = convert(TableauRKExplicit{T},tableaus_rk_explicit[Name])
    if Kind == :fixed
        opts = FixedOptions{T}(;opts...)
        if isadaptive(tab)
            error("Cannot construct a fixed step method from an adaptive step tableau")
        end
    elseif Kind == :adaptive
        opts = AdaptiveOptions{T}(;opts...)
        if !isadaptive(tab)
            error("Cannot construct an adaptive step method from an fixed step tableau")
        end
    end
    RKIntegrator{Kind,Name,T,typeof(opts)}(tab,opts)
end


order(integ::RKIntegrator) = minimum(order(integ.tableau))

name(integ::RKIntegrator) = integ.tableau.name

tdir(ode::ExplicitODE, integ::RKIntegrator) = sign(integ.opts.tstop - ode.t0)

# lower level interface

# explicit RK integrator

"""

Pre allocated arrays to store temporary data.  Used only by
Runge-Kutta integrator.

"""
type RKWorkArrays{Y}
    y   ::Y
    ynew::Y
    yerr::Y
    ks  ::Vector{Y}
end


"""
State for the Runge-Kutta integrator.
"""
type RKState{T,Y} <: AbstractState{T,Y}
    step    ::Step{T,Y}
    dt      ::T
    newdt   ::T
    work    ::RKWorkArrays{Y}
    timeout ::Int
    iters   ::Int # iters<maxiters TODO: not implemented
end

output(st::RKState) = st.step.t, st.step.y, st.work.ks[1]

function show(io::IO, state::RKState)
    show(state.step)
    println("dt      = $(state.dt)")
    println("timeout = $(state.timeout)")
    println("work    = $(state.work)")
end


function init(ode::ExplicitODE,integ::RKIntegrator)
    t0, dt0, y0 = ode.t0, integ.opts.initstep, ode.y0

    # clip the dt0 if t0+dt0 exceeds tstop
    dt0 = tdir(ode,integ)*min(abs(dt0),abs(integ.opts.tstop-t0))

    lk = lengthks(integ.tableau)
    work = RKWorkArrays(zero(y0), # y
                        zero(y0), # ynew
                        zero(y0), # yerr
                        Array(typeof(y0), lk)) # ks

    # we have to allocate each component separately
    for i = 1:lk
        work.ks[i]=zero(y0)
    end

    # pre-initialize work.ks[1]
    ode.F!(t0,y0,work.ks[1])

    step = Step(t0,copy(y0),copy(work.ks[1]))

    timeout = 0 # for step control
    return RKState(step,dt0,dt0,work,timeout,0)
end


#####################
# Fixed step method #
#####################


function onestep!(ode::ExplicitODE, integ::RKIntegratorFixed, state::RKState)
    step = state.step
    work = state.work

    td = tdir(ode,integ)

    if td*step.t >= td*integ.opts.tstop
        # nothing left to integrate
        return finish
    end

    dof  = length(step.y)
    b    = integ.tableau.b
    dt   = td*min(abs(state.dt),abs(integ.opts.tstop-step.t))

    copy!(work.ynew,step.y)

    for k=1:length(b)
        calc_next_k!(work, k, ode, integ.tableau, step, dt)
        for d=1:dof
            work.ynew[d] += dt * b[k]*work.ks[k][d]
        end
    end
    step.t += dt
    copy!(step.y,work.ynew)
    return cont
end


########################
# Adaptive step method #
########################


const timeout_const = 5

# `trialstep!` ends with a step computed for the stepsize `state.dt`
# and stores it in `work.y`, so `work.y` contains a candidate for
# `y(t+dt)` with `dt=state.dt`.
function trialstep!(ode::ExplicitODE, integ::RKIntegratorAdaptive, state::RKState)
    work    = state.work
    step    = state.step
    tableau = integ.tableau
    opts = integ.opts

    td = tdir(ode,integ)

    # use the proposed step size to perform the computations
    state.dt = state.newdt
    dt = state.dt

    if td*step.t >= td*opts.tstop
        # nothing left to integrate
        return finish
    end

    if abs(dt) < opts.minstep
        # TODO: use some sort of logging system
        warn("Minimum step size reached")
        return abort
    end

    # work.y and work.yerr and work.ks are updated after this step
    rk_embedded_step!(work, ode, tableau, step, dt)

    return cont
end

# computes the error for the candidate solution `y(t+dt)` with
# `dt=state.dt` and proposes a new time step
function errorcontrol!(ode::ExplicitODE,
                       integ::RKIntegratorAdaptive,
                       state::RKState)
    work = state.work
    step = state.step
    tableau = integ.tableau
    timeout = state.timeout
    opts = integ.opts
    err, state.newdt, state.timeout =
        stepsize_hw92!(work, step, tableau, state.dt, state.timeout, opts)

    td = tdir(ode,integ)

    # trim in case newdt > dt
    state.newdt = td*min(abs(state.newdt), abs(opts.tstop-(state.step.t+state.dt)))

    if err > 1
        # The error is too large, the step will be rejected.  We reset
        # the timeout and set the new stepsize
        state.timeout = timeout_const
    end

    return err, cont
end

# Here we assume that trialstep! and errorcontrol! have already been
# called, that is `work.y` holds `y(t+dt)` with `dt=state.dt`, and
# error was small enough for us to keep `y(t+dt)` as the next step.
function accept!(ode::ExplicitODE,
                 integ::RKIntegratorAdaptive,
                 state::RKState)
    work    = state.work
    step    = state.step
    tableau = integ.tableau

    # preload ks[1] for the next step
    if tableau.isFSAL
        copy!(work.ks[1],work.ks[end])
    else
        ode.F!(step.t+state.dt, work.ynew, work.ks[1])
    end

    # Swap bindings of y and ytrial, avoids one copy
    step.y, work.ynew = work.ynew, step.y
    # state.dt holds the size of the last successful step
    step.t += state.dt

    return cont
end


##########################
# Lower level algorithms #
##########################

function rk_embedded_step!(work      ::RKWorkArrays,
                           ode       ::ExplicitODE,
                           tableau   ::Tableau,
                           last_step ::Step,
                           dt)
    # Does one embedded R-K step updating work.ynew, work.yerr and work.ks.
    # Assumes that work.ks[:,1] is already calculated!
    # Modifies work.y, work.ynew and work.yerr only

    y      = last_step.y
    dof    = length(y)
    b      = tableau.b

    fill!(work.ynew, zero(eltype(y)))
    fill!(work.yerr, zero(eltype(y)))

    for s=1:lengthks(tableau)
        # we skip the first step beacause we assume that work.ks[1] is
        # already computed
        if s > 1
            calc_next_k!(work, s, ode, tableau, last_step, dt)
        end
        for d=1:dof
            work.ynew[d] += b[1,s]*work.ks[s][d]
            work.yerr[d] += b[2,s]*work.ks[s][d]
        end
    end

    for d=1:dof
        work.yerr[d] = dt*(work.ynew[d]-work.yerr[d])
        work.ynew[d] = y[d] + dt*work.ynew[d]
    end

    return nothing
end


function stepsize_hw92!{T}(work,
                           last_step ::Step,
                           tableau   ::TableauRKExplicit,
                           dt        ::T,
                           timeout   ::Int,
                           opts   ::Options{T})
    # Estimates the error and a new step size following Hairer &
    # Wanner 1992, p167 (with some modifications)
    #
    # If timeout>0 no step size increase is allowed, timeout is
    # decremented in here.
    #
    # Returns the error, newdt and the number of timeout-steps
    #
    # TODO:
    # - allow component-wise reltol and abstol?
    # - allow other norms

    ord = minimum(order(tableau))
    timout_after_nan = 5
    # fac = T[0.8, 0.9, (0.25)^(1/(ord+1)), (0.38)^(1/(ord+1))][1]
    fac = T(8//10)
    facmax = T(5) # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?
    dof = length(last_step.y)

    # in-place calculate yerr./tol
    for d=1:dof

        # TODO: for some reason calling opts.isoutofdomain
        # generates a lot of allocations

        # if opts.isoutofdomain(work.y[d])::Bool
        if isnan(work.y[d])
            return T(10), dt*facmin, timout_after_nan
        end

        y0 = last_step.y[d] # TODO: is this supposed to be the last successful step?
        y1 = work.ynew[d]    # the approximation to the next step
        sci = (opts.abstol + opts.reltol*max(norm(y0),norm(y1)))
        work.yerr[d] ./= sci # Eq 4.10
    end

    # TOOD: should we use opts.norm here as well?
    err   = opts.norm(work.yerr) # Eq. 4.11
    newdt = sign(dt)*min(opts.maxstep, abs(dt)*clamp(fac*(1/err)^(1/(ord+1)),facmin,facmax)) # Eq 4.13 modified

    if timeout > 0
        newdt = sign(dt)*min(abs(newdt), abs(dt))
        timeout -= 1
    end

    return err, newdt, timeout
end


# For clarity we pass the RKWorkArrays part of the state separately,
# this is the only part of state that can be changed here
function calc_next_k!(work      ::RKWorkArrays,
                      i         ::Int,
                      ode       ::ExplicitODE,
                      tableau   ::Tableau,
                      last_step ::Step,
                      dt)
    dof = length(last_step.y)
    t, a, c = last_step.t, tableau.a, tableau.c

    copy!(work.y,last_step.y)
    for j=1:i-1
        for d=1:dof
            work.y[d] += dt * work.ks[j][d] * a[i,j]
        end
    end
    ode.F!(t + c[i]*dt, work.y, work.ks[i])
    return nothing
end

###################################
## Tableaus for explicit RK methods
###################################

# Fixed step:
const tableaus_rk_explicit = Dict{Symbol,ODE.TableauRKExplicit{Rational{Int}}}()

tableaus_rk_explicit[:feuler] =
    TableauRKExplicit{Rational{Int64}}("Forward Euler",(1,),
                                       zeros(Int,1,1),
                                       [1]',
                                       [0]
                                       )

tableaus_rk_explicit[:midpoint] =
    TableauRKExplicit{Rational{Int64}}("Midpoint",(2,),
                                       [0  0
                                        1//2  0],
                                       [0, 1]',
                                       [0, 1//2]
                                       )

tableaus_rk_explicit[:heun] =
    TableauRKExplicit{Rational{Int64}}("Heun",(2,),
                                       [0  0
                                        1  0],
                                       [1//2, 1//2]',
                                       [0, 1])

tableaus_rk_explicit[:rk4] =
    TableauRKExplicit{Rational{Int64}}("Runge-Kutta(4)",(4,),
                                       [0    0    0 0
                                        1//2 0    0 0
                                        0    1//2 0 0
                                        0    0    1 0],
                                       [1//6, 1//3, 1//3, 1//6]',
                                       [0, 1//2, 1//2, 1])

# Adaptive step:
# Heun Euler https://en.wikipedia.org/wiki/Runge–Kutta_methods
tableaus_rk_explicit[:rk21] =
    TableauRKExplicit{Rational{Int64}}("Heun Euler",(2,1),
                                       [0     0
                                        1     0],
                                       [1//2  1//2
                                        1     0],
                                       [0,    1])

# Bogacki–Shampine coefficients
tableaus_rk_explicit[:rk23] =
    TableauRKExplicit{Rational{Int64}}("Bogacki-Shampine",(2,3),
                                       [0           0      0      0
                                        1/2         0      0      0
                                        0         3/4      0      0
                                        2/9       1/3     4/9     0],
                                       [7/24 1/4 1/3 1/8
                                        2/9 1/3 4/9 0],
                                       [0, 1//2, 3//4, 1]
                                       )

# Fehlberg https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
tableaus_rk_explicit[:rk45] =
    TableauRKExplicit{Rational{Int64}}("Fehlberg",(4,5),
                      [  0           0           0            0         0     0
                         1//4        0           0            0         0     0
                         3//32       9//32       0            0         0     0
                         1932//2197 -7200//2197  7296//2197      0         0     0
                         439//216     -8        3680//513    -845//4104   0     0
                         -8//27       2       -3544//2565   1859//4104 -11//40 0 ],
[25//216      0        1408//2565   2197//4104  -1//5  0
 16//135      0        6656//12825 28561//56430 -9//50 2//55],
[0,          1//4,       3//8,       12//13,    1,    1//2])

# Dormand-Prince https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
tableaus_rk_explicit[:dopri5] =
    TableauRKExplicit{Rational{Int64}}("Dormand-Prince", (5,4),
                     [0           0            0                   0            0      0 0
                      1//5        0            0                   0            0      0 0
                      3//40       9//40        0                   0            0      0 0
                      44//45      -56//15      32//9               0            0      0 0
                      19372//6561 -25360//2187 64448//6561 -212//729            0      0 0
                      9017//3168  -355//33     46732//5247   49//176 -5103//18656      0 0
                      35//384     0            500//1113    125//192  -2187//6784 11//84 0],
                     [35//384         0     500//1113       125//192      -2187//6784         11//84      0
                      5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
                     [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
                     )

# Fehlberg 7(8) coefficients
# Values from pag. 65, Fehlberg, Erwin. "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control".
# National Aeronautics and Space Administration.
tableaus_rk_explicit[:feh78] =
    TableauRKExplicit{Rational{Int64}}("Fehlberg(7,8)", (7,8),
                            [     0       0      0       0         0          0        0        0      0       0     0 0 0
                                  2//27   0      0       0         0          0        0        0      0       0     0 0 0
                                  1//36   1//12  0       0         0          0        0        0      0       0     0 0 0
                                  1//24   0      1//8    0         0          0        0        0      0       0     0 0 0
                                  5//12   0    -25//16  25//16     0          0        0        0      0       0     0 0 0
                                  1//20   0      0       1//4      1//5       0        0        0      0       0     0 0 0
                                -25//108  0      0     125//108  -65//27    125//54    0        0      0       0     0 0 0
                                 31//300  0      0       0       61//225    -2//9     13//900   0      0       0     0 0 0
                                  2       0      0     -53//6    704//45   -107//9    67//90    3      0       0     0 0 0
                                -91//108  0      0      23//108 -976//135   311//54  -19//60   17//6  -1//12   0     0 0 0
                               2383//4100 0      0    -341//164 4496//1025 -301//82 2133//4100 45//82 45//164 18//41 0 0 0
                                  3//205  0      0       0        0        -6//41     -3//205  -3//41  3//41   6//41 0 0 0
                              -1777//4100 0      0    -341//164 4496//1025 -289//82 2193//4100 51//82 33//164 12//41 0 1 0],
                              [41//840 0 0 0 0 34//105 9//35 9//35 9//280 9//280 41//840 0 0
                               0 0 0 0 0 34//105 9//35 9//35 9//280 9//280 0     41//840 41//840],
                               [0,    2//27, 1//9, 1//6 , 5//12, 1//2 , 5//6 , 1//6 , 2//3 , 1//3 , 1 , 0, 1]
                            )
