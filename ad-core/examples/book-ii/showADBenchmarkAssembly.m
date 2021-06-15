mrstModule add ad-core ad-props ad-blackoil compositional libgeometry deckformat
%% Define physics, iterations and backend choice
niter = 10;
physics = {'1ph', 'immiscible', 'blackoil', 'overall'};
backends = {'diagonal-rowmajor-deferred'};
nphys = numel(physics);
%% Sizes used in book chapter
if ~exist('sizes', 'var')
    sizes = [20, 25, 30, 35];
    % Uncomment if you want to run bigger stuff
    % sizes = [20, 25, 30, 35, 50, 75, 100, 126];
end
%% Benchmark with and without wells/facility
well_out = assemblyBenchmarkAD(sizes, backends, physics, 'wells', 10);
nowell_out = assemblyBenchmarkAD(sizes, backends, physics, 'none', 10);

%% Name the different physics
physics = well_out.physics;

phys_names = cell(nphys, 1);
for i = 1:nphys
    switch physics{i}
        case '1ph'
            nn = 'single-phase';
        case 'immiscible'
            nn = '3ph, immiscible';
        case 'blackoil'
            nn = '3ph, blackoil';
        case 'overall'
            nn = '6c, compositional';
        case 'natural'
            nn = '3ph, natural compositional';
        otherwise
            error('Unknown physics %s', physics{i})
    end
    phys_names{i} = nn;
end
% Accomodate loading results from file
[backends, physics, results_w, results] = deal(well_out.backends, well_out.physics, well_out.results, nowell_out.results);
ncell = cellfun(@(x) x.ncell, well_out.results(1, :, 1));
sizes = round(ncell.^(1/3));

nb = numel(backends);
nz = numel(sizes);
np = numel(physics);
%% Show scaling plots
figure('position', [680   652   698   326]);
hold on
backend_styles = {'-o', '--o'};
colors = lines(np);

lh = [];
lnames = {};
alldata = cell(2, nb);
passed = false;
for sNo = 1:2
    if sNo == 1
        r = results;
        style = '-o';
        sn = '';
        passed = true;
    else
        r = results_w;
        if passed
            style = '--';
            sn = ' (wells)';
        else
            style = 'o-';
            sn = '';
        end
    end
    for bNo = 1:nb
        data = nan(nz, np);
        for zNo = 1:nz
            for pNo = 1:np
                s = ncell(zNo);
                s = 1;
                data(zNo, pNo) = r{bNo, zNo, pNo}.time_avg/s;
            end
        end
        alldata{sNo, bNo} = data;
        for pNo = 1:np
            h = plot(ncell, data(:, pNo), style, 'color', colors(pNo, :), 'markerFaceColor', colors(pNo, :), 'linewidth', 1.5);
            lh = [lh; h];
            lnames{end+1} = [phys_names{pNo}, sn];
        end
    end
end
upper = 0.5*ncell/1e4;
lower = 2*ncell/1e7;
lh(end+1) = plot(ncell, upper, 'k--', 'linewidth', 1);
lh(end+1) = plot(ncell, lower, 'k--', 'linewidth', 1);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Number of cells');
legend(lh, [lnames, '0.5 s / million cells', '50 s / million cells'], 'location', 'northwest')
axis tight
ylabel('Assembly time [s]')

% xlim([1e4, 3.5e6])
set(gca, 'XTick', [10000      100000     1000000 2e6]);
set(gca, 'XTickLabel', {'1\times10^4', '1\times10^5', '1\times10^6', '2\times10^6'});
ylim([0, 1000])
%% Plot the time taken
figure(1); clf
scale = (sizes.^3)';

d = alldata{1};
dw = alldata{2};

w = max(dw - d, 0);
pos = find(sizes == max(sizes));
% yyaxis left
bdata = [d(pos, :); w(pos, :)]';
top = sum(bdata, 1);
bar(bdata)

grid on
grid minor

legend('Reservoir equations', 'Wells', 'Location', 'northwest')
set(gca, 'XTicklabel', phys_names)
ylabel('Assembly time [s]');
hA = gca;
hA.YRuler.MinorTick = 'on';
%% Print latex table
clc
for i = 1:size(w, 1)
    fprintf('$%d$ & %10d ', sizes(i), sizes(i)^3)
    for j = 1:size(d, 2)
        base = d(i, j);
        well = dw(i, j);
        delta = 100*max(well - base, 0)/base;
        fprintf('& %.2f s & %.2f s', base, well)
    end
    fprintf('\\\\\n');
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
