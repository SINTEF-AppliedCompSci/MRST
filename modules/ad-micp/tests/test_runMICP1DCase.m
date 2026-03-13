run("../../startup.m")
run("examples/runMICP1DCase.m")
assert (isfile("states-00520.vtu"), "The runMICP1DCase.m failed")
