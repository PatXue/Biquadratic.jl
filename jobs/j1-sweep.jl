import Pkg
Pkg.activate("..")

using Biquadratic
using Carlo
using Carlo.JobTools

tm = TaskMaker()
tm.rand_init = true

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 5000
tm.thermalization = 1000
tm.binsize = 50

tm.T = 0.02
tm.J2a = 1.0
tm.J2b = -1.0
J1s = -2.0:0.1:2.0
for J1 in J1s
    tm.J1 = J1
    tm.K = 0.2
    task(tm)
    tm.K = -0.2
    task(tm)
end

job = JobInfo("j1-sweep", Biquadratic.MC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)