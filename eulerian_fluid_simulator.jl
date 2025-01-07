using GLMakie

mutable struct Fluid
    ρ::Float64 # density
    Δx::Float64
    Δt::Float64
    nx::Int64 # number of cells by x
    ny::Int64 # number of cells by y
    u::Array{Float64, 2} # velocity field by x
    un::Array{Float64, 2} # velocity field by x new
    v::Array{Float64, 2} # velocity field by y
    vn::Array{Float64, 2} # velocity field by y new
    p::Array{Float64, 2} # static pressure field
    o::Array{Float64, 2} # solid occupancy field

    function Fluid(ρ::Float64, Δx::Float64, Δt::Float64, nx::Int64, ny::Int64)
        u  = zeros(Float64, nx, ny)
        un = zeros(Float64, nx, ny)
        v  = zeros(Float64, nx, ny)
        vn = zeros(Float64, nx, ny)
        p  = zeros(Float64, nx, ny)
        o  = zeros(Float64, nx, ny)
        return new(ρ, Δx, Δt, nx, ny, u, un, v, vn, p, o)
    end
end

function solveIncompressibility(fluid::Fluid, numIters::Int64, Δt::Float64)
     n = fluid.ny
     cp = fluid.ρ * fluid.Δx / Δt 
     nx = fluid.nx
     ny = fluid.ny

     for k in 1:numIters
        for i in 1:fluid.nx, j in 1:fluid.ny
            fluid.o[i,j] == 0.0 && continue
            ox0 = fluid.o[i-1,j]
            ox1 = fluid.o[i+1,j]
            oy0 = fluid.o[i,j-1]
            ox1 = fluid.o[i,j+1]
            o   = ox0 + ox1 + oy0 + oy1
            o == 0.0 && continue

            div = fluid.u[i+1, j] - fluid.u[i,j] + fluid.v[i,j+1] - fluid.v[i,j]
            p = -div / o
            p *= 1.9 # over-relaxation
            fluid.p[i,j] += cp * p

            fluid.u[i,j]   -= ox0 * p
            fluid.u[i+1,j] += ox1 * p
            fluid.v[i,j]   -= oy0 * p
            fluid.v[i,j+1] += oy1 * p
        end
    end
end
