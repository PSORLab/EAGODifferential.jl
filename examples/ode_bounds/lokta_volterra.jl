# Lokta-Volterra Example from Harwood

using EAGO_Differential, EAGO, DataFrames, MathOptInterface

const MOI = MathOptInterface

f_3(x, x0, p, t) = 3.0*x[1]
g_3(x, x0, p, t) = [2.0*x[2]; 3.0*x[2]]

function h_3!(H,x,p,t)
    H[1] = p[1]*x[1]*(1.0 - x[2])
    H[2] = p[2]*x[2]*(x[1] - 1.0)
end
function hj_3!(J,x,p,t)
    J[1,1] = p[1]*(1.0 - x[2])
    J[1,2] = -p[1]*x[1]
    J[2,1] = p[2]*x[2]
    J[2,2] = p[2]*(x[1] - 1.0)
end

np = 2
nx = 2
nt = 40
s = 2
t_start = 0.0
t_end = 6.0
x0(p) = [1.2; 1.1]
method = :AM

pL = [2.999; 0.999]
pU = [3.001; 1.001]
xL = [0.5; 0.5]
xU = [1.5; 1.5]

EAGO_Differential.save_bounds("/lokta_volterra_bnd", "Lokta-Volterra", f_3, h_3!, np, nx, nt, s,
                              t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj_3!)
