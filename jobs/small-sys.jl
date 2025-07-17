import Pkg
Pkg.activate("..")

using Biquadratic
using Carlo
using Carlo.JobTools

tm = TaskMaker()
tm.init_type = :rand

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 20000
tm.thermalization = 0
tm.binsize = 100

tm.T = 0.02
tm.J2a = 1.0
tm.J2b = -1.0
Ks = (-0.2, 0.2)
J1s = -2.0:0.1:2.0
for K in Ks
    tm.K = K
    for J1 in J1s
        tm.J1 = J1
        task(tm)
    end
end

job = JobInfo("small-sys", Biquadratic.MC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)