using EAGO_Differential, EAGO, DataFrames, CSV

opt = EAGO.Optimizer(absolute_tolerance = 0.01, relative_tolerance = 0.00001,
                     iteration_limit = Int(1E6),
                     verbosity = 1,
                     output_iterations = 100,
                     header_iterations = 2000)

T = 273.0
K2 = 46.0*exp(6500.0/T - 18.0)
K3 = 2.0*K2
k_1 = 53.0
k_1s = k_1*10^(-6)
k_5 = 0.0012
cO2 = 0.002

pL = [10.0; 10.0; 0.001]
pU = [1200.0; 1200.0; 40.0]

xL = [0.0; 0.0; 0.0; 0.0; 0.0]
xU = [140.0; 140.0; 140.0; 0.4; 140.0]

f(x,x0,p,t) = x[1]

x0(p) = [0.0; 0.0; 0.0; 0.4; 140]

function h(out,x,p,t)

    fill!(out, zero(p[1]))

    out[1] = k_1*x[4]*x[5] - cO2*(p[1] + p[2])*x[1] + (p[1]/K2)*x[3] + (p[2]/K3)*x[2] - k_5*x[1]^2
    out[2] = p[2]*x[1]*cO2 - (p[2]/K3 + p[3])*x[2]
    out[3] = p[1]*x[1]*cO2 - (p[1]/K2)*x[3]
    out[4] = -k_1s*x[4]*x[5]
    out[5] = -k_1*x[4]*x[5]

end

function hj(out,x,p,t)

    fill!(out, zero(p[1]))

    out[1,1] = -cO2*(p[1] + p[2]) - 2.0*k_5*x[1]
    out[1,2] = p[2]/K3
    out[1,3] = p[1]/K2
    out[1,4] = k_1*x[5]
    out[1,5] = k_1*x[4]

    out[2,1] = p[2]*cO2
    out[2,2] = -(p[2]/K3 + p[3])

    out[3,1] = p[1]*cO2
    out[3,3] = -p[1]/K2

    out[4,4] = -k_1s*x[5]
    out[4,5] = -k_1s*x[4]

    out[5,4] = -k_1*x[5]
    out[5,5] = -k_1*x[4]
end

np = 3
nx = 5
nt = 201
s = 1

t_start = 0.0
t_end = 2.0
method = :AM

out1, out2 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt)
