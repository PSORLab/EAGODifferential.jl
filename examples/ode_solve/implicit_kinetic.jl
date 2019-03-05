# Good Until Result

using EAGO_Differential, EAGO, DataFrames, CSV

opt = EAGO.Optimizer(absolute_tolerance = 1.00, relative_tolerance = 0.2,
                     iteration_limit = Int(1E6),
                     verbosity = 1,
                     output_iterations = 20,
                     header_iterations = 200)

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

data = CSV.read("C:/Users/wilhe/.julia/dev/EAGO_Differential/examples/ode_solve/kinetic_intensity_data.csv")

function f200(x,x0,p,t)
    temp = zero(p[1])
    for i in 1:200
        temp += (x[1,i] + (2.0/21.0)*x[2,i] + (2.0/21.0)*x[3,i] - data.intensity[i])^2
    end
    return temp
end

function f40(x,x0,p,t)
    temp = zero(p[1])
    for i in 1:40
        temp += (x[1,i] + (2.0/21.0)*x[2,i] + (2.0/21.0)*x[3,i] - data.intensity[5*i])^2
    end
    return temp
end

t_start = 0.0
t_end = 2.0
method = :AM


run_style = :AM40
if run_style === :II200
    var_II200, opt_II200 = solve_ode(f200, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 201, 1, :AM, opt)
    EAGO_Differential.save_history("ii200", var_II200, opt_II200, "Implicit Kinetic Example, 200 Time-Steps, Implicit Euler")
elseif run_style === :AM200
    var_AM200, opt_AM200 = solve_ode(f200, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 201, 2, :AM, opt)
    EAGO_Differential.save_history("am200", var_AM200, opt_AM200, "Implicit Kinetic Example, 200 Time-Steps, Trapezoidal Method")
elseif run_style == :BDF200
    var_BDF200, opt_BDF200 = solve_ode(f200, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 201, 2, :BDF, opt)
    save_history("bdf200", var_BDF200, opt_BDF200, "Implicit Kinetic Example, 200 Time-Steps, Backwards Difference Formula")
elseif run_style === :II40
    var_II40, opt_II40 = solve_ode(f40, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 41, 1, :AM, opt)
    EAGO_Differential.save_history("ii40", var_II40, opt_II40, "Implicit Kinetic Example, 40 Time-Steps, Implicit Euler")
elseif run_style === :AM40
    var_AM40, opt_AM40 = solve_ode(f40, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 41, 2, :AM, opt)
    EAGO_Differential.save_history("am40", var_AM40, opt_AM40, "Implicit Kinetic Example, 40 Time-Steps, Trapezoidal Method")
elseif run_style == :BDF40
    var_BDF40, opt_BDF40 = solve_ode(f40, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, 41, 2, :AM, opt)
    EAGO_Differential.save_history("bdf40", var_BDF40, opt_BDF40, "Implicit Kinetic Example, 40 Time-Steps, Backwards Difference Formula")
end
