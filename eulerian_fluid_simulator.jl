using GLMakie

mutable struct Fluid
    ρ::Float64 # density
    Δx::Float64
    Δt::Float64
    nx::Int64 # number of cells by x
    ny::Int64 # number of cells by y
    v::Array{Float64, 2} # velocity field by x
    vn::Array{Float64, 2} # velocity field by x new
    u::Array{Float64, 2} # velocity field by y
    un::Array{Float64, 2} # velocity field by y new
    p::Array{Float64, 2} # static pressure field
    o::Array{Float64, 2} # solid occupancy field

    function Fluid(ρ::Float64, Δx::Float64, Δt::Float64, nx::Int64, ny::Int64)
        v  = zeros(Float64, nx, ny)
        vn = zeros(Float64, nx, ny)
        u  = zeros(Float64, nx, ny)
        un = zeros(Float64, nx, ny)
        p  = zeros(Float64, nx, ny)
        o  = zeros(Float64, nx, ny)
        return new(ρ, Δx, Δt, nx, ny, v, vn, u, un, p, o)
    end
end