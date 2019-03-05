# Sahlodin 2011

using EAGO_Differential, EAGO, DataFrames, MathOptInterface

const MOI = MathOptInterface

f_3(x, x0, p, t) = 3.0*x[1]
g_3(x, x0, p, t) = [2.0*x[2]; 3.0*x[2]]

function h_3!(H,x,p,t)
    H[1] = -x[1]^2 + p[1]
end
function hj_3!(J,x,p,t)
    J[1,1] = -2.0*x[1]
end

np = 1
nx = 1
nt = 20
s = 2
t_start = 0.0
t_end = 1.0
x0(p) = [9.0]
method = :AM

pL = [-1.0]
pU = [1.0]
xL = [0.1]
xU = [9.0]

lower_eval = EAGO_Differential.save_bounds("/sahlodin2011_bnd_20_am2", "sahlodin2011", f_3, h_3!, np, nx, nt, s,
                              t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj_3!)

@time EAGO_Differential.save_bounds("/sahlodin2011_bnd_20_am2", "sahlodin2011", f_3, h_3!, np, nx, nt, s,
                              t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj_3!)
timer = 0.028465/10^3
