run ../../startup.m
addpath(fullfile(ROOTDIR, 'mex', 'mex_opmpressure'))
addpath(fullfile(ROOTDIR, 'solvers', 'tpfa', 'compressible'))
[x1, x2, x3] = test_compTPFA2;
