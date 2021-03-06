using EAGO_Differential, EAGO, DataFrames, CSV

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

# soft build, then evaluate
lower_eval = ImplicitODELowerEvaluator{3}()

# build the basic evaluator (w/o inequality constraints)
EAGO_Differential.build_evaluator!(lower_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)

lower_vars = zeros(nx*(nt-1))
upper_vars = zeros(nx*(nt-1))

xL = [0.0; 0.0; 0.0; 0.0; 0.0]
xU = [140.0; 140.0; 140.0; 0.4; 140.0]
indx = 1
for i in 1:(nt-1)
    lower_vars[(5*(i-1)+1):(5*i)] = zeros(5)
    upper_vars[(5*(i-1)+1)] = 140.0
    upper_vars[(5*(i-1)+2)] = 140.0
    upper_vars[(5*(i-1)+3)] = 140.0
    upper_vars[(5*(i-1)+4)] = 0.4
    upper_vars[(5*(i-1)+5)] = 140.0
end
append!(lower_vars, pL)
append!(upper_vars, pU)
n = NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)


# build the basic evaluator (w/o inequality constraints)
upper_eval = ImplicitODEUpperEvaluator()
EAGO_Differential.build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)
y = [xL; (pL + pU)/2.0]
EAGO.set_current_node!(upper_eval, n)
EAGO_Differential.relax_ode_implicit!(upper_eval)

#=
state_relax = upper_eval.state_relax_n
pdi_kinetics_bounds = DataFrame(lower_kinetic_ip_1 = lo.(state_relax[1,:]),
                                upper_kinetic_ip_1 = hi.(state_relax[1,:]),
                                lower_kinetic_ip_2 = lo.(state_relax[2,:]),
                                upper_kinetic_ip_2 = hi.(state_relax[2,:]),
                                lower_kinetic_ip_3 = lo.(state_relax[3,:]),
                                upper_kinetic_ip_3 = hi.(state_relax[3,:]),
                                lower_kinetic_ip_4 = lo.(state_relax[4,:]),
                                upper_kinetic_ip_4 = hi.(state_relax[4,:]),
                                lower_kinetic_ip_5 = lo.(state_relax[5,:]),
                                upper_kinetic_ip_5 = hi.(state_relax[5,:]))

CSV.write("pi_kinetic_bounds_ii200.csv", pdi_kinetics_bounds)
=#

y = copy(upper_vars)
y[1001:1003] = (pL + pU)/2.0
EAGO.set_current_node!(lower_eval, n)
EAGO_Differential.relax_ode_implicit!(lower_eval, y)

y1 = (0.25*pL + 0.75*pU)/2.0

#EAGO_Differential.relax_ode_implicit!(lower_eval, y1)
