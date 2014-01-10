ecl_dir = @(m, kr, g) ...
   ['/work/test-data/ecl-compare/ECL_runs_initsat/', ...
    'ecl_', m, '_', kr, '_', g];

models = { 'vertSlize' };
kr     = { 'LIN', 'NL' };
grav   = { 'NOG', 'GR' };

for m = models, for k = kr, for g = grav,

d = ecl_dir(m{1}, k{1}, g{1});

clc
[w, x, gg, fluid, deck, report] = runAndCompareComp(d, false);

fprintf('\n\n\n***** %s *****\n\n\n', d);

pause(5)
end, end, end
