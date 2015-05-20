# Here are tests which test what interface the solvers require.

################################################################################
# This is to test a scalar-like state variable
# (due to @acroy: https://gist.github.com/acroy/28be4f2384d01f38e577)
const delta0 = 0.
const V0 = 1.
const g0 = 0.

# define custom type ...
immutable CompSol
  rho::Matrix{Complex128}
  x::Float64
  p::Float64
 
  CompSol(r,x,p) = new(copy(r),x,p)
end
 
# ... which has to support the following operations
# to work with odeX
Base.norm(y::CompSol, p::Float64) = maximum([Base.norm(y.rho, p) abs(y.x) abs(y.p)])
Base.norm(y::CompSol) = norm(y::CompSol, 2.0)

+(y1::CompSol, y2::CompSol) = CompSol(y1.rho+y2.rho, y1.x+y2.x, y1.p+y2.p)
-(y1::CompSol, y2::CompSol) = CompSol(y1.rho-y2.rho, y1.x-y2.x, y1.p-y2.p)
*(y1::CompSol, s::Real) = CompSol(y1.rho*s, y1.x*s, y1.p*s)
*(s::Real, y1::CompSol) = y1*s
/(y1::CompSol, s::Real) = CompSol(y1.rho/s, y1.x/s, y1.p/s)
################################################################################
 
# define RHSs of differential equations
# delta, V and g are parameters
function rhs(t, y, delta, V, g)
  H = [[-delta/2 V]; [V delta/2]]
 
  rho_dot = -im*H*y.rho + im*y.rho*H
  x_dot = y.p
  p_dot = -y.x
 
  return CompSol( rho_dot, x_dot, p_dot)
end
 
# inital conditons
rho0 = zeros(2,2);
rho0[1,1]=1.;
y0 = CompSol(complex(rho0), 2., 1.)
 
# solve ODEs
endt = 2*pi;

t,y1 = ODE.ode45((t,y)->rhs(t, y, delta0, V0, g0), y0, [0., endt]);
print("Testing interface for scalar-like state... ")
for solver in solvers
    # these only work with some Array-like interface defined:
    if solver in [ODE.ode23s, ODE.ode4s_s, ODE.ode4s_kr]
        continue
    end
    t,y2 = solver((t,y)->rhs(t, y, delta0, V0, g0), y0, linspace(0., endt, 500));
    @test norm(y1[end]-y2[end])<0.1
end
println("ok.")

################################################################################
# TODO: test a vector-like state variable, i.e. one which can be indexed.
