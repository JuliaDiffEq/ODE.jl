include("../src/ODE.jl")

module Test

using ODE

T = Float64
Y = Vector{T}
t0 = zero(T)
y0 = T[one(T)]

solvers = [# ODE.RKIntegratorAdaptive{:rk45},
            # ODE.RKIntegratorFixed{:feuler},
            ODE.DenseOutput]

for st in solvers
    ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
    opts = Dict(:initstep=>0.1,
                :tspan=>[0.,0.5,1.],
                :points=>:specified,
                :reltol=>1e-5,
                :abstol=>1e-5)

    prob = ODE.solve(ode,st;opts...)

    println("Raw iterator")
    for (t,y) in prob
        println((t,y,norm(y-[exp(t)])))
    end

    println(collect(prob))
end

end
