const testsets = [
            Dict(
                 :F!    => (t,y,dy)->dy[1]=6.0,
                 :y0    => [0.],
                 :tout => [0:0.1:1;],
                 :jac   => (t,y,dy)->dy[1]=0.0,
                 :sol   => t->[6t],
                 :isscalar => true,
                 :name  => "y'=6",
                 :initstep => 0.1),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=2t,
                 :y0    => [0.],
                 :tout => [0:0.001:1;],
                 :jac   => (t,y,dy)->dy[1]=0.0,
                 :sol   => t->[t^2],
                 :isscalar => true,
                 :name  => "y'=2t",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=y[1],
                 :y0    => [1.0],
                 :tout => [0:0.001:1;],
                 :jac   => (t,y,dy)->dy[1]=1.0,
                 :sol   => t->[exp(t)],
                 :isscalar => true,
                 :name  => "y'=y",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=y[1],
                 :y0    => [1.0],
                 :tout => [1:-0.001:0;],
                 :jac   => (t,y,dy)->dy[1]=1.0,
                 :sol   => t->[exp(t-1)],
                 :isscalar => true,
                 :name  => "y'=y backwards",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->(dy[1]=-y[2];dy[2]=y[1]),
                 :y0    => [1.0,2.0],
                 :tout => [0:.1:1;],
                 :jac   => (t,y,dy)->copy!(dy,Float64[[0,1] [-1,0]]),
                 :sol   => t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)],
                 :isscalar => false,
                 :name  => "pendulum",
                 :initstep => 0.001)
            ]


# Testing function ode
const integrators = [ODE.RKIntegratorFixed{:feuler},
                  ODE.RKIntegratorFixed{:midpoint},
                  ODE.RKIntegratorFixed{:heun},
                  ODE.RKIntegratorFixed{:rk4},
                  ODE.RKIntegratorAdaptive{:rk21},
                  ODE.RKIntegratorAdaptive{:rk23},
                  ODE.RKIntegratorAdaptive{:rk45},
                  ODE.RKIntegratorAdaptive{:dopri5},
                  ODE.RKIntegratorAdaptive{:feh78},
                  ODE.ModifiedRosenbrockIntegrator
                  ]

function test_ode()
    tol = 0.002

    for integ in integrators
        println("Testing $integ")
        for ts in testsets
            println("Testing problem $(ts[:name])")

            tout, h0, stepper = ts[:tout], ts[:initstep], integ

            y0, F!, jac!, sol = ts[:y0], ts[:F!], ts[:jac], ts[:sol]

            F(t,y) = (dy = similar(y); F!(t,y,dy); return dy)

            for points = [:specified, :all]
                if ts[:isscalar]
                    # test the ODE.odeXX scalar interface (if the equation is scalar)
                    Fscal = (t,y)->F(t,[y])[1]
                    y0scal = y0[1]
                    # with jacobian
                    tj,yj = ODE.ode(Fscal,y0scal,tout,
                                    solver = stepper,
                                    points = points,
                                    initstep = h0,
                                    J! = jac!)
                    @test_approx_eq_eps yj map(x->sol(x)[1],tj) tol
                    # without jacobian
                    t,y   = ODE.ode(Fscal,y0scal,tout,
                                    solver = stepper,
                                    points = points,
                                    initstep = h0)
                    @test_approx_eq_eps y  map(x->sol(x)[1],tj) tol

                    # results with and without jacobian should be exactly the same
                    @test_approx_eq yj y

                    if points == :specified
                        # test if we covered the whole timespan
                        @test length(tout) == length(t) == length(tj)
                        @test_approx_eq tout t
                        @test_approx_eq tout tj
                    end
                end

                # ODE.odeXX vector interface
                # with jacobian
                tj,yj = ODE.ode(F,y0,tout,
                                solver = stepper,
                                points = points,
                                initstep = h0,
                                J! = jac!)
                @test_approx_eq_eps hcat(yj...) hcat(map(sol,tj)...) tol
                # without jacobian
                t,y   = ODE.ode(F,y0,tout,
                                solver = stepper,
                                points = points,
                                initstep = h0)
                @test_approx_eq_eps hcat(y...)  hcat(map(sol,t)...) tol

                @test_approx_eq hcat(yj...) hcat(y...)

                # TODO: tests for `y::AbstractArray`
                # # ODE.odeXX array interface for arrays
                # # with jacobian
                # tj,yj = ODE.ode(F,reshape(y0,length(y0),1,1),tout,
                #                 solver = stepper,
                #                 points = points,
                #                 initstep = h0,
                #                 J! = jac!)
                # @test_approx_eq_eps hcat(yj...) hcat(map(sol,tj)...) tol
                # # without jacobian
                # t,y   = ODE.ode(F,reshape(y0,length(y0),1,1),tout,
                #                 solver = stepper,
                #                 points = points,
                #                 initstep = h0)
                # @test_approx_eq_eps hcat(y...)  hcat(map(sol,t)...) tol

                # @test_approx_eq hcat(yj...) hcat(y...)

                if points == :specified
                    # test if we covered the whole timespan
                    @test length(tout) == length(t) == length(tj)
                    @test_approx_eq tout t
                    @test_approx_eq tout tj
                end

                # test the iterator interface
                equation = ODE.ExplicitODE(tout[1],y0,F!)
                opts     = Dict(:tout => tout,
                                :initstep => h0,
                                :points => points)

                solver   = ODE.solve(equation,stepper;opts...)

                for (t,y) in solver
                    @test_approx_eq_eps y sol(t) tol
                end

                for (t,y) in ODE.solve(equation,ODE.DenseOutput{stepper}; opts...)
                    @test_approx_eq_eps y sol(t) tol
                end

            end
        end
    end
end

test_ode()
