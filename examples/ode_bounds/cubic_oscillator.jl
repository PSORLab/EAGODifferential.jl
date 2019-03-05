# Houska 2015 Case Study:

using EAGO_Differential, EAGO, DataFrames, MathOptInterface

const MOI = MathOptInterface

f_3(x, x0, p, t) = 3.0*x[1]
g_3(x, x0, p, t) = [2.0*x[2]; 3.0*x[2]]

function h_3!(H,x,p,t)
    H[1] = x[2] + 0.1*(x[1]-x[1]^3 - x[1]*x[2]^2)
    H[2] = -x[1] + 0.1*(x[2]-x[2]*x[1]^2 - x[2]^3) - 0.2*x[2]
end
function hj_3!(J,x,p,t)
    J[1,1] = 0.1*(1.0 - 3.0*x[1]^2 - x[2]^2)
    J[1,2] = -0.2*x[1]*x[2]
    J[2,1] = -1.0 - 0.2*(x[2]*x[1])
    J[2,2] = 0.1*(1.0 - x[1]^2 - 3.0*x[2]^2) - 0.2
end

np = 2
nx = 2
nt = 480
s = 2
t_start = 0.0
t_end = 20.0
x0(p) = [p[1]; p[2]]
method = :AM

pL = [2.24; -0.01]
pU = [2.26; 0.01]
xL = [-1.5; -1.75]
xU = [3.5; 1.5]

EAGO_Differential.save_bounds("/cubic_oscillator_bnd", "cubic_oscillator", f_3, h_3!, np, nx, nt, s,
                              t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj_3!)
