using EAGO, EAGO_Differential

opt = EAGO.Optimizer()

f(x,x0,p,t) = x[9]
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
t_end = 1.0
method = :BDF

pL = [1.4]; pU = [1.6]
xL = [0.1]; xU = [1.0]

out1, out2 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt)
