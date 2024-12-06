run("../../startup.m")
run("examples/runMICP3DCase.m")
assert (isfile("statesCO2afterMICP-00030.vtu"), "The runMICP3DCase.m failed")
