using GLMakie, Observables

mutable struct Fluid
    ρ::Float64 # density
    dl::Float64
    nx::Int64 # number of cells by x
    ny::Int64 # number of cells by y
    u::Array{Float64, 2} # velocity field by x
    un::Array{Float64, 2} # velocity field by x new
    v::Array{Float64, 2} # velocity field by y
    vn::Array{Float64, 2} # velocity field by y new
    p::Array{Float64, 2} # static pressure field
    o::Array{Float64, 2} # fluid occupancy field

    function Fluid(ρ::Float64, dl::Float64, nx::Int64, ny::Int64)
        u  = zeros(Float64, nx, ny)
        un = zeros(Float64, nx, ny)
        v  = zeros(Float64, nx, ny)
        vn = zeros(Float64, nx, ny)
        p  = zeros(Float64, nx, ny)
        o  = ones(Float64, nx, ny)
        for i in 1:nx, j in 1:ny # defining border
            i == 1  && (o[j,i] = 0.0)
            i == nx && (o[j,i] = 0.0)
            j == 1  && (o[j,i] = 0.0)
            j == ny && (o[j,i] = 0.0)
        end
        return new(ρ, dl, nx, ny, u, un, v, vn, p, o)
    end
end

function solveIncompressibility(fluid::Fluid, numIters::Int64, dt::Float64)
     n = fluid.ny
     cp = fluid.ρ * fluid.dl / dt
     nx = fluid.nx
     ny = fluid.ny

     for k in 1:numIters
        for i in 1:fluid.nx, j in 1:fluid.ny
            fluid.o[j,i] == 0.0 && continue
            ox0 = fluid.o[j,i-1]
            ox1 = fluid.o[j,i+1]
            oy0 = fluid.o[j-1,i]
            oy1 = fluid.o[j+1,i]
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
       fluid.u[1,i] = fluid.u[2,i]
       fluid.u[fluid.ny,i] = fluid.u[fluid.ny-1,i]
   end 
   for j in 1:fluid.ny
       fluid.v[j,1] = fluid.v[j,2]
       fluid.v[j,fluid.nx] = fluid.v[j,fluid.ny-1]
   end
end

function sample_field(fluid::Fluid, x::Float64, y::Float64, field::String)
    dl = fluid.dl
    h₁ = 1.0 / dl 
    h₂ = 0.5 * dl 

    x = max(min(x, fluid.nx * dl), dl)
    y = max(min(y, fluid.ny * dl), dl)

    dx = 0.0
    dy = 0.0

    f = similar(fluid.u)
    field == "U_FIELD" ? (f = fluid.u; dy = h₂) : 
    field == "V_FIELD" ? (f = fluid.v; dx = h₂) : nothing

    x₀ = Int(min(floor((x - dx) * h₁), fluid.nx))
    tx = ((x - dx) - x₀ * dl) * h₁
    x₁ = Int(min(x₀ + 1, fluid.nx))

    y₀ = Int(min(floor((y - dy) * h₁), fluid.ny))
    ty = ((y - dy) - y₀ * dl) * h₁
    y₁ = Int(min(y₀ + 1, fluid.ny))

    sx = 1.0 - tx
    sy = 1.0 - ty

    val = sx * sy * f[y₀,x₀]
        + tx * sy * f[y₀,x₁]
        + tx * ty * f[y₁,x₁]
        + sx * ty * f[y₁,x₀]

    return val 
end

@inline function avgU(fluid::Fluid, j::Int64, i::Int64)
    return (fluid.u[j-1,i] + fluid.u[j,i] + fluid.u[j-1,i+1] + fluid.u[j,i+1]) * 0.25
end

@inline function avgV(fluid::Fluid, j::Int64, i::Int64)
    return (fluid.v[j,i-1] + fluid.v[j,i] + fluid.v[j+1,i-1] + fluid.v[j+1,i]) * 0.25
end

function advect_velocity(fluid::Fluid, dt::Float64)
    fluid.un, fluid.u = fluid.u, fluid.un
    fluid.vn, fluid.v = fluid.v, fluid.vn

    dl = fluid.dl
    h₂ = 0.5 * dl 

    for i in 1:fluid.nx, j in 1:fluid.ny
        if fluid.o[j,i] != 0.0 && fluid.o[j,i-1] != 0.0 && j < fluid.ny
            x = i * dl 
            y = j * dl + h₂
            u = fluid.u[j,i]
            v = avgV(fluid,j,i)
            x = x - dt * u
            y = y - dt * v
            u = sample_field(fluid,x,y,"U_FIELD")
            fluid.un[j,i] = u
        end
        if fluid.o[j,i] != 0.0 && fluid.o[j-1,i] != 0.0 && i < fluid.nx
            x = i * dl + h₂
            y = j * dl 
            u = avgU(fluid,j,i)
            v = fluid.v[j,i]
            x = x - dt * u
            y = y - dt * v
            v = sample_field(fluid,x,y,"V_FIELD")
            fluid.un[j,i] = v
        end
    end

    fluid.u, fluid.un = fluid.un, fluid.u
    fluid.v, fluid.vn = fluid.vn, fluid.v

end

function simulate(fluid::Fluid, dt::Float64, numIters::Int64)
    fill!(fluid.p, 0.0)
    solveIncompressibility(fluid, numIters, dt)
    extrapolate(fluid)
    advect_velocity(fluid, dt)
end

function to_pixels(value_in_cm::Float64)
   dpi = sqrt(1920^2 + 1080^2) / 15.6 # should be fixed with each new monitor
   return (value_in_cm / 2.54) * dpi
end

window_width    = 25.0 # (cm)
window_height   = 25.0 # (cm)
axis_width      = 20.0
axis_height     = 20.0
ρ = 1000.0; dl = 0.001; nx = 100; ny = 100

fluid = Fluid(ρ, dl, nx, ny)
data = Observable(similar(fluid.p))
fill!(data[], 0.0)

set_theme!(theme_dark())
fig = Figure(size = (to_pixels(window_width), to_pixels(window_height)))
ax  = Axis(fig[1,1][1,1], title = "Pressure field", aspect = 1) 
heatmap!(ax, data, colormap = :viridis, colorrange = (-10.0, 10.0))
display(fig)

# Configuring UI
quit = Observable(false)
on(events(fig).keyboardbutton) do event
   if event.action == Keyboard.press && event.key == Keyboard.escape # ESC key
      quit[] = true
   end
end

fps = 60

# initial speeds
fluid.u[Int(ny/2),Int(nx/2)] = -10.0
fluid.v[Int(nx/2),Int(ny/2)] = -10.0

println("The simulation is running...")
while !quit[]
    # run simulation and update figure data
    simulate(fluid, 0.1, 5)
    data[] = copy(fluid.p')
    notify(data)
    sleep(1/fps)
    display(fig)
end
println("The simulation stopped...")
GLMakie.closeall()