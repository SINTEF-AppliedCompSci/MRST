%% Parallel scaling of AD assembly
mrstModule add ad-core deckformat ad-blackoil vista ad-props compositional libgeometry
if ~exist('n_strong', 'var')
    n_strong = 2e6;
end

if ~exist('doSave', 'var')
    doSave = false;
end

if ~exist('max_threads', 'var')
    maxNumCompThreads('automatic');
    max_threads = maxNumCompThreads();
end

if ~exist('modeltype', 'var')
    modeltype = 'immiscible';
end
threads = 1:max_threads;
%% Test strong scaling
stored_data = sprintf('scaling_strong_%d.mat', n_strong);

if exist(stored_data, 'file')
    load(stored_data);
else
    G = makeGrid(n_strong);
    strong = cell(max_threads, 1);
    for i = threads
        fprintf('Testing strong scaling with %d threads\n', i);
        maxNumCompThreads(i);
        n = ceil((n_strong)^(1/3));
        strong{i} = assemblyBenchmarkAD(n, 'diagonal-rowmajor', modeltype, 'none', 10);
    end

    if doSave
        save(sprintf('scaling_strong_%d', n_strong), 'strong', 'max_threads', 'modeltype', 'n_strong');
    end
end

%% Define plots
threads = 1:max_threads;
pos = [   601   724   765   226];

f = 'time_avg';
%     f = 'time_eqs'
getData = @(f) cellfun(@(x) x.results{1}.(f), strong);
t_total = getData('time_avg');
t_eqs = getData('time_eqs');
t_init = getData('time_state');

t_strong = cellfun(@(x) x.results{1}.(f), strong);

nt = (1:max_threads)';

%% Strong scaling
p1 = 0.9;
p2 = 0.95;

colors = lines(5);
lw = 1.5;

t = t_total;
figure('position', pos); hold on;

if 0
    l = {'Actual'};
    plot(nt, t(1)./t, '.', 'MarkerSize', 18, 'color', colors(1, :));
else
    l = {'Init', 'Equations', 'Total'};
    plot(nt, t_init(1)./t_init, 'x', 'MarkerSize', 6, 'linewidth', 2 , 'color', colors(1, :));
    plot(nt, t_eqs(1)./t_eqs, '+', 'MarkerSize', 6, 'linewidth', 2, 'color', colors(4, :));
    plot(nt, t_total(1)./t_total, '.', 'MarkerSize', 18, 'color', colors(5, :));

end
plot(nt, nt, '-', 'linewidth', lw, 'color', [.9, .2, .2]);
yl = ylim();
amdahl = @(F_s) 1./(F_s  + (1-F_s)./nt);

amdahl1 = amdahl(1-p1);
amdahl2 = amdahl(1-p2);

plot(nt, amdahl1, 'color', colors(2, :), 'linewidth', lw);
plot(nt, amdahl2, 'color', colors(3, :), 'linewidth', lw);
ylim([1, max_threads])
xlim([1, max_threads])
s = 'Amdahl''s law, F_s = %1.2f';
legend(l{:}, 'Ideal', sprintf(s, 1-p1), sprintf(s, 1-p2), 'location', 'northwest')
xlabel('Number of threads');
ylabel('Speedup')
%%
function G = makeGrid(n)
    dims = ceil(n.^(1/3))*[1, 1, 1];
    G = cartGrid(dims);
    G = computeGeometry(G);
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
