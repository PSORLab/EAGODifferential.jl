using EAGO, EAGO_Differential

opt = EAGO.Optimizer()

f(x,x0,p,t) = (x[9]-0.35)^2
function h(H,x,p,t)
    H[1] = p[1]*x[1]*(1-x[1])
end
function hj(J,x,p,t)
    J[1,1] = p[1]*(1.0-2.0*x[1])
end
x0(p) = [0.1]

np = 1
nx = 1
nt = 10
s = 2

t_start = 0.0
t_end = 3.0
method = :AM

pL = [0.25]; pU = [1.0]
xL = [0.1]; xU = [1.0]

out1, out2 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt)

var_II200, opt_II200 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, 1, :AM, opt)
save_history("ii200", var_II200, opt_II200, "Implicit Kinetic Example, 200 Time-Steps, Implicit Euler")

var_AM200, opt_AM200 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, 2, :AM, opt)
save_history("am200", var_AM200, opt_AM200, "Implicit Kinetic Example, 200 Time-Steps, Trapezoidal Method")

var_BDF200, opt_BDF200 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, 2, :BDF, opt)
save_history("bdf200", var_BDF200, opt_BDF200, "Implicit Kinetic Example, 200 Time-Steps, Backwards Difference Formula")
