import Pkg
Pkg.activate("..")

using Biquadratic
using Carlo
using Carlo.JobTools

tm = TaskMaker()
tm.init_type = :rand

L = 40
tm.Lx = tm.Ly = L
tm.sweeps = 20000
tm.thermalization = 40000
tm.binsize = 100

tm.J2a = 1.0
tm.J2b = -1.0
tm.J1 = 0.1
# Ks = (-0.05, -0.02, -0.01, -0.005, -0.003, -0.001, 0.001, 0.003, 0.005, 0.01, 0.05)
Ks = (-0.005, 0.005)
Ts = 0.0:0.05:0.7
for K in Ks
    tm.K = K
    for T in Ts
        tm.T = max(0.01, T)
        task(tm)
    end
end

job = JobInfo("temp-sweep", Biquadratic.MC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)