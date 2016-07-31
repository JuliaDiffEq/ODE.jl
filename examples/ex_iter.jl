using ODE

# pick your solver
integ = [ODE.RKIntegratorAdaptive{:rk45},
           ODE.ModifiedRosenbrockIntegrator][2]

# Define IVP-instance which holds the mathematical problem definition:
t0 = 0.0
y0 = [1.0]
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])


### Forward time integration

# options for the solver
opts = Dict(:initstep=>0.1,
            :tstop=>1.0,
            :reltol=>1e-5,
            :abstol=>1e-5)

# create a Problem instance
prob = ODE.solve(ode,integ;opts...)

 # iterate over the solution
println("t, y, err")
for (t,y) in prob
    println((t,y[1],abs(y[1]-e.^t)))
end

# or collect it
println(collect(prob))

### Reverse time integration, rest as above
t0 = 1.0
y0 = [1.0]
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
opts = Dict(:initstep=>0.1,
            :tstop=>0.0,
            :reltol=>1e-5,
            :abstol=>1e-5)

prob = ODE.solve(ode,integ;opts...)

println("t, y, err")
for (t,y) in prob # iterate over the solution
    println((t,y[1],abs(y[1]-e.^(t-1))))
end

println(collect(prob))


### Dense output
opts = Dict(:initstep=>0.1,
            :tstop=>0.0,
            :reltol=>1e-5,
            :abstol=>1e-5,
            :tout=>[t0:-0.1:0;])

prob = ODE.solve(ode,ODE.DenseOutput{integ};opts...)
println("t, y, err")
for (t,y) in prob # iterate over the solution
    println((t,y[1],abs(y[1]-e.^(t-1))))
end

println(collect(prob))
