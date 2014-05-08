function make_hybrid_table

clear all;
load rho_huge;

% Doubling span in the hybrid table
Ph = double_span(P);
Th = double_span(T);

% sampling full set of values from coarsest table
pp = zeros(Ph.num, 1);
for i = 1:Ph.num
    pp(i) = Ph.span(1) + (i-1) * Ph.stepsize;
end
%pp = linspace(Ph.span(1), Ph.span(2), Ph.num);

tt = zeros(Th.num, 1);
for i = 1:Th.num
    tt(i) = Th.span(1) + (i-1) * Th.stepsize;
end
%tt = linspace(Th.span(1), Th.span(2), Th.num);

[tgrid, pgrid] = meshgrid(tt, pp);

CO2 = CO2props('rho_small', '');
rho = reshape(CO2.rho(pgrid(:), tgrid(:)), Ph.num, Th.num);
CO2.dispose();

% sampling values from medium table, wherever there is coverage
CO2 = CO2props('rho_big_trunc', '');
rho_tmp = reshape(CO2.rho(pgrid(:), tgrid(:)), Ph.num, Th.num);
CO2.dispose();

ix = find(~isnan(rho_tmp(:)));
rho(ix) = rho_tmp(ix);

% inserting finely sampled values where applicable
rho(1:2000, 1:2000) = vals; % defined when loading 'rho_huge'

% saving result as new table
P = Ph;
T = Th;
vals = rho;

save rho_hybrid P T vals;

end
%% END MAIN FUNCTION

function V = double_span(Vorig)

    V.num = Vorig.num * 2 - 1;
    V.span = [Vorig.span(1), Vorig.span(2) + diff(Vorig.span)];
    V.stepsize = Vorig.stepsize;
end
