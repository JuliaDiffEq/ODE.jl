Below you will find the simplest implementation of a reasonable
generic solver that finds solutions to explicit ODE equations.

```julia
type EulerIntegrator <: AbstractIntegrator
    # nothing in here yet
end

type EulerState
    t
    y
    dy
end

output(state::EulerState) = state.t, state.y, state.dy

# see we don't need a fancy constructor
function solver(ode::ExplicitODE,
                ::Type{EulerIntegrator};
                options...)
    return Problem(ode,EulerIntegrator())
end

function init(ode::ExplicitODE, integ::EulerIntegrator)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState(t0,y0,dy0)
end

function onestep!(ode::ExplicitODE, integ::EulerIntegrator, state::EulerState)
    t, y, dy = output(state)
    dt = 0.1
    y += dt*dy
    t += dt
    ode.F!(t0,y0,dy0)           # update the derivative
    return cont
end
```

There are several problems with the above implementation.  First of
all, it has a constant prescribed step size.  This could easily be
fixed by changing the type definition and the `solver` to


```julia
type EulerIntegrator <: AbstractIntegrator
    initstep
end

function solver(ode::ExplicitODE,
                ::Type{EulerIntegrator};
                initstep = 0.1,
                options...)
    return Problem(ode,EulerIntegrator(initstep))
end
```

we should also change the line `dt = 0.1` in the `onestep!` function
to `dt = stepper.initstep`.  Now we can run our integrator with a
custom step size!

```julia
sol = solver(ode,EulerIntegrator,initstep = 0.01)
for (t,y) in sol
    if t > 1
        print(t,y)
        break
    end
end
```

Another issue is type stability, to make `EulerIntegrator` perform
better we should explicitly annotate the fields in both
`EulerIntegrator` and `EulerState` like this

```julia
type EulerIntegrator{T,Y} <: AbstractIntegrator{T,Y}
    initstep::T
end

type EulerState{T,Y}
    t::T
    y::Y
    dy::Y
end

function solver(ode::ExplicitODE{T,Y},
                ::Type{EulerIntegrator};
                initstep::T = T(0.1),
                options...)
    return Problem(ode,EulerIntegrator{T,Y}(initstep))
end
```

But the `EulerState{T,Y}` is exactly the same as `Step` from
`base.jl`, so we can simplify it a bit more

```julia
type EulerState{T,Y}
    step::Step{T,Y}
end
```

Once we do that, in the `init` we should change
`EulerState(t0,y0,dy0)` to `EulerState(Step(t0,y0,dy0))` and redefine
`output` to `output(state::EulerState)=output(state.step)`
(`output(::Step)` is already defined in `base.jl`).

One could even replace `EulerState` with `Step` completely, but this
won't allow us to extend the state with some additional variables and
storage space in the future.

The last thing is that our stepper will continue the integration
forever: it doesn't have a stopping condition.  We could add one as an
option.

```julia
type EulerIntegrator{T,Y} <: AbstractIntegrator{T,Y}
    initstep::T
    tstop::T
end

function solver(ode::ExplicitODE{T,Y},
                ::Type{EulerIntegrator};
                initstep::T = T(0.1),
                tstop::T = T(Inf)
                options...)
    return Problem(ode,EulerIntegrator{T,Y}(initstep,tstop))
end

function onestep!(ode::ExplicitODE, integ::EulerIntegrator, state::EulerState)
    t, y, dy = output(state)

    if t > integ.tstop
        return finished
    end

    dt = integ.initstep
    y += dt*dy
    t += dt
    ode.F!(t0,y0,dy0)           # update the derivative
    return cont
end
```

As a final improvement, we can (although this is not necessary) use a
structure `FixedOptions` from `options.jl` to keep our options in one
structure.  A corresponding options type for adaptive solver is
`AdaptiveOptions`.  This way we can use the standarized defaults for
most options and keep our solver in line with the standard
keywords.  Naturally, we have to update `onestep!` to use the subtype.

```julia
type EulerIntegrator{T,Y} <: AbstractIntegrator{T,Y}
    options::FixedOptions{T,Y}
end

function solver(ode::ExplicitODE{T,Y},
                ::Type{EulerIntegrator};
                options...)
    options = FixedOptions{T}(;options...)
    return Problem(ode,EulerIntegrator{T,Y}(options))
end
```
