CNH40=0.00277
# initialize dimensions of pfr tank
tankLength = 60
width = 10
height = 25
# area is perpendicular to flow aka cross section of tank
area = width*height
vol = tankLength*area

# initialize flow rate, velocity of bulk
# Q = flow rate of bulk = 36 L/day = 1.5L/hr = 1500 cm3/hr
Q = 1500
v0=Q/area
tau = tankLength/v0 # retention/residence time aka amount of time fluid spends in volume. units is hrs
k = 0.35 # reaction rate constant in hrs^-1 (first order assumption)

nodes = 20
dx = (1/(nodes-1))

# initial condition
function x(p)
    out = zeros(nodes)
    out[1] = 1.0
end

# state variable
function h!(out,x,p,t)
    fill!(0.0, out);
    out[1] = x[1] - 1;   # no change in inlet concentration
    for i = 2:nodes
        out[i] = -(x[i] - x[i-1])/(dx) - Da*x[i]
    end
end

# jacobian
function hj!(out,x,p,t)
    fill!(0.0, out);
    out[1,1] = one(p[1]);   # no change in inlet concentration
    for i = 2:nodes
        out[i,i] = -one(p[1])/dx - Da
        out[i,i-1] = one(p[1])/dx
    end
end
