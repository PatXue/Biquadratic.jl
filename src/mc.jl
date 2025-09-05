# Note: Using temperature in units of energy (k_B = 1)
# All energy units are in terms of J2a (best to set J2a = 1)
struct MC{AlgType} <: AbstractMC
    T::Float64     # Temperature
    J1::Float64    # Nearest neighbor coupling energy
    J2a::Float64   # Next-nearest neighbor coupling energy (NE-SW direction)
    J2b::Float64   # Next-nearest neighbor coupling energy (NW-SE direction)
    K::Float64     # Biquadratic coupling energy

    outdir::String # Output directory for local spin current plots
    savefreq::Int  # No. of sweeps between saving local spin current

    spins::PeriodicMatrix{SpinVector}
end

function MC{AlgType}(; T=0.5, J1=0.1, J2a=1.0, J2b=-1.0, K=0.1, Lx::Int=20,
                     Ly::Int=20, outdir=".", savefreq=0) where {AlgType}
    MC{AlgType}(T, J1, J2a, J2b, K, outdir, savefreq,
                fill(zeros(SpinVector), (Lx, Ly)))
end

function MC{AlgType}(params::AbstractDict) where {AlgType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    J1 = params[:J1]
    J2a = params[:J2a]
    J2b = params[:J2b]
    K = params[:K]

    outdir = haskey(params, :outdir) ? params[:outdir] : "."
    savefreq = haskey(params, :savefreq) ? params[:savefreq] : 0

    return MC{AlgType}(; T, J1, J2a, J2b, K, Lx, Ly, outdir, savefreq)
end

function Carlo.init!(mc::MC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = params[:init_type]
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :orth
        init_orth!(mc.spins)
    elseif init_type == :eag
        init_eag!(mc.spins)
    elseif init_type == :dir
        dir::String = params[:dir] # Directory to read spin dump file from
        h5open("$dir/run0001.dump.h5") do file
            Carlo.read_checkpoint!(mc, file["simulation"])
        end
    else
        rand!(ctx.rng, mc.spins)
    end
    return nothing
end

function Carlo.sweep!(mc::MC, ctx::Carlo.MCContext)
    Carlo.sweep!(mc, ctx.rng)
end

# Returns if spin current should be saved on this sweep (assuming thermalized)
function is_save_sweep(mc::MC, ctx::Carlo.MCContext)
    measure_sweeps = ctx.sweeps - ctx.thermalization_sweeps
    return mc.savefreq > 0 && measure_sweeps % mc.savefreq == 0
end

function save_spin(mc::MC, ctx::Carlo.MCContext)
    jldopen("$(mc.outdir)/spins.jld2", "a+") do file
        sweep_name = "sweep$(ctx.sweeps - ctx.thermalization_sweeps)"
        file[sweep_name] = Matrix(mc.spins)
    end
end

function Carlo.measure!(mc::MC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Magnetization per lattice site
    mag = norm(sum(mc.spins)) / N
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Mag4, mag^4)

    # Energy per lattice site
    energy = 0.0
    # Averaged adjacent dot products
    Dx0 = Dy0 = 0.0
    Dxπ = Dyπ = 0.0
    # Spin current
    spin_curr = zeros(3)
    x_hat = SVector(1, 0, 0)

    for y in 1:Ly
        for x in 1:Lx
            energy += half_energy(mc, x, y)

            s = mc.spins[x, y]
            sx = mc.spins[x+1, y]
            sy = mc.spins[x, y+1]
            x_dot = s ⋅ sx
            y_dot = s ⋅ sy
            Dx0 += x_dot
            Dy0 += y_dot
            Dxπ += x_dot * (-1)^(x+y)
            Dyπ += y_dot * (-1)^(x+y)

            spin_curr += s × sx
        end
    end
    energy /= N
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)
    Dx0 /= N
    Dy0 /= N
    measure!(ctx, :Dx0, Dx0)
    measure!(ctx, :Dy0, Dy0)
    Dxπ /= N
    Dyπ /= N
    measure!(ctx, :Dxπ, abs(Dxπ))
    measure!(ctx, :Dyπ, abs(Dyπ))
    spin_curr /= N
    measure!(ctx, :J_s, norm(spin_curr))
    measure!(ctx, :P, norm(x_hat × spin_curr))

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{MC{AlgType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType}
    T = params[:T]
    J2a = params[:J2a]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N * J2a/T * (mag2 - mag^2)
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end

function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["spins"] = mc.spins
    return nothing
end
function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    raw_spins = read(in, "spins")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_spins)
    return nothing
end
