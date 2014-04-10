
clear

mrstModule add ad-refactor ad-fi
N = 10;
rhs = (1:N).';
A = speye(N);

eq1 = ADI(rhs          , {A, A});
eq2 = ADI(rhs(end:-1:1), {A, A});
eqs = {eq1, eq2};
%%
P = linearProblem(eqs, {'cell', 'cell'}, {'forward', 'reverse'});
solver = mldivideSolverAD();


sol = solver.solveLinearProblem(P);

%%
