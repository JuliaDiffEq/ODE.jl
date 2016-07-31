using ODE
using Base.Test

const tol = 1e-2

solvers = [
           ## Non-stiff
           # fixed step
           ODE.ode1,
           ODE.ode2_midpoint,
           ODE.ode2_heun,
           ODE.ode4,
           ODE.ode4ms,
           ODE.ode5ms,
           # adaptive
           ODE.ode21,
           ODE.ode23,
           ODE.ode45_dp,
           ODE.ode45_fe,
           ODE.ode78,

           ## Stiff
           # fixed-step
           ODE.ode4s_s,
           ODE.ode4s_kr,
           # adaptive
           ODE.ode23s]

for solver in solvers
    println("using $solver")

    # dy
    # -- = 6 ==> y = 6t
    # dt
    # we need to fix initstep for the fixed-step methods
    t,y=solver((t,y)->6.0, 0., [0:.1:1;], initstep=.1)
    @test maximum(abs(y-6t)) < tol
    tj,yj=solver((t,y)->6.0, 0., [0:.1:1;], initstep=.1, J! = (t,y,dy)->dy[1]=0.0)
    @test maximum(abs(yj-6tj)) < tol
    @test norm(yj-y,Inf)<eps(1.)

    # dy
    # -- = 2t ==> y = t.^2
    # dt
    t,y  =solver((t,y)->2t, 0., [0:.001:1;], initstep=0.001)
    @test maximum(abs(y-t.^2)) < tol
    tj,yj=solver((t,y)->2t, 0., [0:.001:1;], initstep=0.001, J! = (t,y,dy)->dy[1]=0.0)
    @test maximum(abs(yj-tj.^2)) < tol
    @test norm(yj-y,Inf)<eps(1.)

    # dy
    # -- = y ==> y = y0*e.^t
    # dt
    t,y=solver((t,y)->y, 1., [0:.001:1;], initstep=0.001)
    @test maximum(abs(y-e.^t)) < tol
    tj,yj=solver((t,y)->y, 1., [0:.001:1;], initstep=0.001, J! = (t,y,dy)->dy[1]=1.0)
    @test maximum(abs(yj-e.^tj)) < tol
    @test norm(yj-y,Inf)<eps(1.)

    # reverse time integration
    t,y=solver((t,y)->y, 1., [1:-.001:0;], initstep=0.001)
    @test maximum(abs(y-e.^(t-1))) < tol
    tj,yj=solver((t,y)->y, 1., [1:-.001:0;], initstep=0.001, J! = (t,y,dy)->dy[1]=1.0)
    @test maximum(abs(yj-e.^(tj-1))) < tol
    @test norm(yj-y,Inf)<eps(1.)

    # dv       dw
    # -- = -w, -- = v ==> v = v0*cos(t) - w0*sin(t), w = w0*cos(t) + v0*sin(t)
    # dt       dt
    #
    # y = [v, w]
    t,y=solver((t,y)->[-y[2]; y[1]], [1., 2.], [0:.001:2*pi;], initstep=0.001)
    ys = hcat(y...).'   # convert Vector{Vector{Float}} to Matrix{Float}
    @test maximum(abs(ys-[cos(t)-2*sin(t) 2*cos(t)+sin(t)])) < tol
    tj,yj=solver((t,y)->[-y[2]; y[1]], [1., 2.], [0:.001:2*pi;], initstep=0.001, J! = (t,y,dy)->copy!(dy,Float64[[0,1] [-1,0]]))
    ysj = hcat(yj...).'   # convert Vector{Vector{Float}} to Matrix{Float}
    @test maximum(abs(ysj-[cos(tj)-2*sin(tj) 2*cos(tj)+sin(tj)])) < tol
    @test norm(map(norm,yj-y),Inf)<eps(1.)

    # test typeof(tspan)==Vector{Int} does not throw
    @test_throws ErrorException t,y=solver((t,y)->2y, 0., [0,1])
    # test typeof(y0)==Vector{Int} does not throw
    @test_throws ErrorException t,y=solver((t,y)->[2y], [0], [0,1])
    # test typeof(y0)==Int does not throw
    @test_throws ErrorException t,y=solver((t,y)->2y, 0, [0,1])
    # test if we can deal with a mixed case
    @test_throws ErrorException t,y=solver((t,y)->2y, Number[1,1.1,BigInt(1)], Rational[0,1])
end

# Test negative starting times ODE.ode23s
@test length(ODE.ode23s((t,y)->[-y[2]; y[1]], [1., 2.], [-5., 0])[1]) > 1

# rober testcase from http://www.unige.ch/~hairer/testset/testset.html
let
    println("ROBER test case")
    function f(t, y)
        ydot = similar(y)
        ydot[1] = -0.04*y[1] + 1.0e4*y[2]*y[3]
        ydot[3] = 3.0e7*y[2]*y[2]
        ydot[2] = -ydot[1] - ydot[3]
        ydot
    end
    t = [0., 1e11]
    t,y = ODE.ode23s(f, [1.0, 0.0, 0.0], t; abstol=1e-8, reltol=1e-8,
                     maxstep=1e11/10, minstep=1e11/1e18)

    refsol = [0.2083340149701255e-07,
              0.8333360770334713e-13,
              0.9999999791665050] # reference solution at tspan[2]
    @test norm(refsol-y[end], Inf) < 2.1e-10
end

include("interface-tests.jl")
include("iterators.jl")

println("All looks OK")
