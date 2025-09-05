import Pkg
Pkg.activate("..")

using Biquadratic
using Carlo
using Carlo.JobTools

tm = TaskMaker()

tm.J2a = 1.0
tm.J2b = -1.0
tm.J1 = 0.1
# Ks = (-0.05, -0.02, -0.01, -0.005, -0.003, -0.001, 0.001, 0.003, 0.005, 0.01, 0.05)
Ks = (-0.005, 0.005)
Ts = 0.05:0.05:0.7

tm.dir = "."

tm.sweeps = 500000
tm.thermalization = 500000
tm.binsize = 1000

tm.Lx = tm.Ly = 120
for K in Ks
    tm.K = K
    tm.init_type = K < 0 ? :eag : :orth
    for T in Ts
        tm.T = max(0.01, T)
        spins_dir = "large-sys.data/$(current_task_name(tm))"
        task(tm)
        tm.init_type = :dir
        tm.dir = spins_dir
    end
end

job = JobInfo("large-sys", Biquadratic.MC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
