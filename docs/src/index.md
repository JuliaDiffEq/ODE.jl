```@contents
Pages = [
    "tutorials/euler_integrator.md",
    "man/basics.md",
    "man/base.md"
    ]
```

# ODE.jl

## Top level interface

If you are looking to getting a solution without the additional hussle
of handling an iterator we provide the wrappers `ODE.odeXX`.  They
provide a simplistic way of handling explicit differential equations
of the form `y'=F(t,y)` with `y` being either a `Number` or an
`AbstractArray` of any dimension.  Below we solve a simple initial
value problem given by `y'=y` with initial data at `t0=0.0` and
`y0=1.0` on the interval `[t0,1]`.

```@example ode
using ODE
tspan  = [0.0,1.0]
y0     = 1.0
F(t,y) = y
(t,y)  = ODE.ode(F,y0,tspan)
```

The vectors `t` and `y` store the time and solution values at the
corresponding times.

You might find the basic interface limiting.  First of all, it stores
all the results, so if you are only interested in the final value of
`y` it still stores all the intermediate steps.  Secondly, you cannot
process the results on the fly (e.g. plot the current state of a
solution).  If you need more control you should consider using the
iterator interface.

## Iterator interface

To offeset the limitations of the `ODE.ode` interface we implemented a
general.  First we define an initial value problem, in our case this is
an explicit differential equation `y'=y` with inital data `y0=[1.0]`
given at the time `t0=0.0`.

```@example iterator
using ODE
t0  = 0.0
y0  = [1.0]
F!(t,y,dy) = dy[1]=y[1]
ode = ODE.ExplicitODE(t0,y0,F!)
```

Note that unlike in `ODE.ode` we now have to supply an in place
function `F!` instead of an explicit function `F`.  Now we are ready
to produce the iterator that solvese to our problem.

```@example iterator
sol = ODE.solve(ode)
for (t,y) in sol
    @show (t,y)
    if t > 1
        break
    end
end
```

Note that we had to interrupt the loop because `sol` would be
producing solutions ad infinitum (in theory, in practice we will get
to the point where the solver won't be able to produce reasonable
solution anymore).  To set the final integration time and other
parameters of the integrator `integ` we can pass optional arguments to
`ODE.solver`.

```@example iterator
sol = ODE.solve(ode; tstop = 1)
for (t,y) in sol
    @show (t,y)
end
```

This approach has the added benefit of the solution never exceeding
the final time.  Apart from the time and value `(t,y)` the `ODE.solve`
returns also the derivative, you can retrive it as the third argument
in the returned tuple.  In the following example we use it compute the
absolute error.

```@example iterator
sol = ODE.solve(ode; tstop = 1)
for (t,y,dy) in sol
    err = norm(y-dy)
    @show err
end
```

With `tstop` specified we can also get all results at once using
`collect`.

```@example iterator
res = collect(sol)
```

Note that `collect` returns a vector of triples `(t,y,dy)`.

## Options

Both `ODE.ode` and `ODE.solve` accept the following keyword arguments.

- `integ`: the type of integrator to use, defaults to a adaptive
  Runge-Kutta method of order 4/5.  To see the list of available
  integrators see [`Integrators`](@ref).

- `initstep`: The initial step size, defaults to `eps(T)^(1/3)`.

- `tstop`: The final integration time, never exceeded by the
  integrator.  In case of `ODE.ode(F,y0,tspan)` this option defaults
  to the last element of `tspan` if it is a vector.  In `ODE.solve`
  the default is `tstop=Inf`.  If `tstop` is smaller then `t0` the
  integration runs backwards in time.

Apart from these general options, each integrator has its own keyword
arguments.  In particular all integrators with adaptive step size
can be cotrolled with

- `reltol`, `abstol`: The relative and absolute error tolerances.  The
  solution guarantees that at each step we have
  `norm((y-yc)*reltol.+abstol)<=1`, where `yc` is a true solution to
  and IVP.  Defaults are `reltol=eps(T)^(1/3)/10`,
  `abstol=eps(T)^(1/2)/10`.

- `norm`: The norm used to measure error in the formula above,
  defaults to `y->Base.vecnorm(y,Inf)`.  You can specify it to assign
  different weights to different components of `y`.

- `minstep`, `maxstep`: Minimal and maximal stepsize for the
  integrator.  If at any point the stepsize exceeds these limits the
  integrator will yield an error and cease producing
  results.  Deafaults are `minstep=10*eps(T)` and `maxstep=1/minstep`.

- `maxiters`: The number of iterations before the integrator ceases to
  work, defaults to `Inf`.  Useful as a safeguard from iterator
  continuing ad infinitum.

- `isoutofdomain`: Applied to each component of `y`, if
  `isoutofdomain(y[i])==true` the integrator stops.  Defaults to
  `Base.isnan`.

Apart from these, each integrator may support additional options.

## Integrators

### Explicit Runge-Kutta integrators

### Rosenbrock methods

### Backwards differential formula (BDF) methods

### ???
