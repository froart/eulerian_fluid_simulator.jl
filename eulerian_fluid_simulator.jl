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
            fluid.o[j,i] == 0.0 && continue
            ox0 = fluid.o[j,i-1]
            ox1 = fluid.o[j,i+1]
            oy0 = fluid.o[j-1,i]
            ox1 = fluid.o[j+1,i]
            o   = ox0 + ox1 + oy0 + oy1
            o == 0.0 && continue

            div = fluid.u[j,i+1] - fluid.u[j,i] + fluid.v[j+1,i] - fluid.v[j,i]
            p = -div / o
            p *= 1.9 # over-relaxation
            fluid.p[j,i] += cp * p

            fluid.u[j,i]   -= ox0 * p
            fluid.u[j,i+1] += ox1 * p
            fluid.v[j,i]   -= oy0 * p
            fluid.v[j+1,i] += oy1 * p
        end
    end
end

function extrapolate(fluid::Fluid)
   for i in 1:fluid.nx
       fluid.u[0,i] = fluid.u[1,i]
       fluid.u[fluid.ny,i] = fluid.u[fluid.ny-1,i]
   end 
   for j in 1:fluid.ny
       fluid.v[j,0] = fluid.v[j,1]
       fluid.v[j,fluid.nx] = fluid.v[j,fluid.ny-1]
   end
end

function sample_field(fluid::Fluid, x::Float64, y::Float64, field)
    Δx = fluid.Δx
    h₁ = 1.0 / Δx
    h₂ = 0.5 * Δx

    x = max(min(x, fluid.nx * Δx), Δx)
    y = max(min(y, fluid.ny * Δx), Δx)

    dx = 0.0
    dy = 0.0

    f
    field == "U_FIELD" ? (f = fluid.u; dy = h₂) : 
    field == "V_FIELD" ? (f = fluid.v; dx = h₂) : nothing

    x₀ = min(floor((x - dx) * h₁), fluid.nx)
    tx = ((x - dx) - x₀ * Δx) * h₁
    x₁ = min(x₀ + 1, fluid.nx)

    y₀ = min(floor((y - dy) * h₁), fluid.ny)
    tx = ((y - dy) - y₀ * Δx) * h₁
    y₁ = min(y₀ + 1, fluid.ny)

    sx = 1.0 - tx
    sy = 1.0 - ty

    val = sx * sy * f[y₀,x₀]
        + tx * sy * f[y₀,x₁]
        + tx * ty * f[y₁,x₁]
        + sx * ty * f[y₁,x₀]

    return val 
end