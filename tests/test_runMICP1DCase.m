run("../../startup.m")
run("examples/runMICP1DCase.m")
assert (isfile("examples/vtk_1DCase/states-00520.vtu"), "The runMICP1DCase.m failed")
