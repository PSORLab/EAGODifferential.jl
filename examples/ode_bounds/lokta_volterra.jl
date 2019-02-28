# Lokta-Volterra Example from Harwood

using EAGO_Differential, EAGO, DataFrames, MathOptInterface

const MOI = MathOptInterface

lower_eval = ImplicitODELowerEvaluator{2}()

f_3(x,p,t) = 3.0*x[1]
g_3(x,p,t) = [2.0*x[2]; 3.0*x[2]]

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
nt = 30
s = 1
t_start = 0.0
t_end = 10.0
x0(p) = [1.2; 1.1]
method = :AM

pL = [2.99; 0.99]
pU = [3.01; 1.01]
xL = [0.5; 0.5]
xU = [1.5; 1.5]

build_evaluator!(lower_eval, f_3, h_3!, np, nx, nt, s, t_start, t_end,
                 method, pL, pU, xL, xU, x0; g = g_3, hj = hj_3!)

lower_vars = Float64[]
upper_vars = Float64[]
for i in 1:(nt-1)
    append!(lower_vars, xL)
    append!(upper_vars, xU)
end
append!(lower_vars, pL)
append!(upper_vars, pU)

# update the node info for the evaluator to initial node
n1 = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)
EAGO.set_current_node!(lower_eval, n1)

y = [3.0; 1.0]
fval = MOI.eval_objective(lower_eval, y)

lower_eval1 = ImplicitODELowerEvaluator{2}()
build_evaluator!(lower_eval1, f_3, h_3!, np, nx, nt, s, t_start, t_end,
                 method, pL, pU, xL, xU, x0; g = g_3, hj = hj_3!)
