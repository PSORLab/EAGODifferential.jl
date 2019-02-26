@testset "Parametric Interval Preprocessing" begin
    interval_preprocess_ode!(x::EAGO.Optimizer, y::EAGO.NodeBB)
end

@testset "Create MidPoint Node" begin
    ymid = create_mid_node(y::NodeBB, nx::Int, np::Int, nt::Int)
end

@testset "Interval MidPoint Upper Bound" begin
    midpoint_upper_bnd_ode!(x::EAGO.Optimizer, y::NodeBB)
end
