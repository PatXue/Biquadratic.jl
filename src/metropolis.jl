# Calculate the energy at a lattice site (x, y) if it had spin s
function energy(mc::MC{:Metropolis}, s::SpinVector, x, y)
    nn = nn_sum(mc.spins, x, y)
    nnna = nnna_sum(mc.spins, x, y)
    nnnb = nnnb_sum(mc.spins, x, y)
    H0 = s ⋅ (mc.J1 * nn + mc.J2a * nnna + mc.J2b * nnnb)
    biquad = mc.K * ((s ⋅ mc.spins[x-1, y])^2 + (s ⋅ mc.spins[x+1, y])^2)
    return H0 + biquad
end

function Carlo.sweep!(mc::MC{:Metropolis}, rng::AbstractRNG=default_rng())
    Lx, Ly = size(mc.spins)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        x = rand(rng, 1:Lx)
        y = rand(rng, 1:Ly)
        old_s = mc.spins[x, y]
        # Propose new spin vector
        new_s = rand(rng, SpinVector)
        ΔE = energy(mc, new_s, x, y) - energy(mc, old_s, x, y)

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / mc.T)
        if prob >= 1.0 || rand(rng) < prob
            mc.spins[x, y] = new_s
        end
    end
    return nothing
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::MC{:Metropolis}, x, y)
    s = mc.spins[x, y]
    nn = mc.spins[x+1, y] + mc.spins[x, y+1]
    nnna = mc.spins[x+1, y+1]
    nnnb = mc.spins[x+1, y-1]
    H0 = s ⋅ (mc.J1 * nn + mc.J2a * nnna + mc.J2b * nnnb)
    biquad = mc.K * (s ⋅ mc.spins[x+1, y])^2
    return H0 + biquad
end
