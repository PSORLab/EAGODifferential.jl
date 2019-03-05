# Shin Stirred-Tank Example (Singularity Error)

using EAGO_Differential, EAGO, DataFrames, MathOptInterface
const MOI = MathOptInterface

tau = 20.0
k2 = 0.4

pL = [0.9; 0.8; 10.0]
pU = [1.1; 1.0; 50.0]

# Need to determine upper bounds on xB = X, xD = Y
#
X = 0.5
Y = 0.5
xL = [0.0; 0.0; 0.0; 0.0]
xU = [0.2; X; 0.5; Y]

f(out, x, x0, p, t) = ()
g(out, x, x0, p, t) = ()

function h(out, x, p, t)
    out[1] = -p[3]*x[1]*x[2] - k2*x[1]*x[3] + (p[1] - 2.0*x[1])/tau
    out[2] = -p[3]*x[1]*x[2] + (p[2] - 2.0*x[2])/tau
    out[3] = p[3]*x[1]*x[2] - k2*x[1]*x[3] - (2.0/tau)*x[3]
    out[4] = k2*x[1]*x[3]  - (2.0/tau)*x[4]
end

function hj(out, x, p, t)

    fill!(out, zero(p[1]))

    out[1,1] = -p[3]*x[2] - k2*x[3] - 2.0/tau
    out[1,2] = -p[3]*x[1]
    out[1,3] =  -k2*x[1]

    out[2,1] = -p[3]*x[2]
    out[2,2] = -p[3]*x[1] - 2.0/tau

    out[3,1] = p[3]*x[2] - k2*x[3]
    out[3,2] = p[3]*x[1]
    out[3,3] =  -k2*x[1] - (2.0/tau)

    out[4,1] = k2*x[3]
    out[4,3] = k2*x[1]
    out[4,4] = -(2.0/tau)
end

x0(p) = [0.0; 0.0; 0.0; 0.0]

time_start = 0.0
time_end = 15.0
nx = 4
np = 3
nt = 21

x0(p) = zeros(4)
method = :AM
s = 2

EAGO_Differential.save_bounds("/stirred_tank_bnd", "Stirred-Tank", f, h, np, nx, nt, s,
                              t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)
