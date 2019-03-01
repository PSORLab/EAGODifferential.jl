using EAGO_Differential, EAGO, DataFrames, CSV

opt = EAGO.Optimizer(absolute_tolerance = 0.1, relative_tolerance = 0.01,
                     iteration_limit = Int(1E6),
                     verbosity = 1,
                     output_iterations = 20,
                     header_iterations = 200)

pL = [0.9; 0.3]
pU = [1.0; 0.4]

nodes = 21
tau = 10.0
dx = 1.0/(nodes-1)
xL = zeros(nodes)
xU = ones(nodes)

function x0(p)
    out = zeros(typeof(p[1]),nodes)
    out[1] = p[1]
    return out
end

# Residual for method of lines
function h(out,x,p,t)
    fill!(out, zero(p[1]))
    out[1] = x[1] - 1;
    for i = 2:nodes
        out[i] = x[i-1]/dx - (1.0/dx + tau*p[2])*x[i]
    end
end

function hj(out,x,p,t)
    fill!(out, zero(p[1]))
    out[1,1] = one(p[1])
    for i = 2:nodes
        out[i,i] = (-1.0/dx)*x[i] - p[1]
        out[i,i-1] = 1.0/dx
    end
end

f(x,x0,p,t) = x[1000]

np = 2
nx = nodes
nt = 51
s = 1

t_start = 0.0
t_end = 1.0
method = :AM

# calculate aprior lower bounds via ODE solution

# calculate aprior upper bounds via ODE solution
#

out1, out2 = solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt)
