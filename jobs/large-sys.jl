import Pkg
Pkg.activate("..")

using Biquadratic
using Carlo
using Carlo.JobTools

tm = TaskMaker()

tm.Lx = tm.Ly = 80
tm.sweeps = 150000
tm.binsize = 500
# tm.savefreq = 5000

tm.J2a = 1.0
tm.J2b = -1.0
tm.J1 = 0.1
# Ks = (-0.05, -0.02, -0.01, -0.005, -0.003, -0.001, 0.001, 0.003, 0.005, 0.01, 0.05)
Ks = (-0.005, 0.005)
Ts = 0.0:0.05:0.7
for K in Ks
    tm.init_type = K < 0 ? :eag : :orth
    tm.K = K
    for T in Ts
        tm.thermalization = (0.25 ≤ T ≤ 0.45) ? 200000 : 100000
        tm.T = max(0.01, T)
        task(tm)
    end
end

tm.sweeps = 200000
tm.binsize = 1000
tm.Lx = tm.Ly = 120
for K in Ks
    tm.init_type = K < 0 ? :eag : :orth
    tm.K = K
    for T in Ts
        tm.thermalization = (0.25 ≤ T ≤ 0.45) ? 400000 : 200000
        tm.T = max(0.01, T)
        task(tm)
    end
end

job = JobInfo("large-sys", Biquadratic.MC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)