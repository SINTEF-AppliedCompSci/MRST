mrstModule add ad-core
%% Conceptual example
group = StateFunctionGrouping('mystorage'); % Store in state.mystorage
% Set the numbers
group = group.setStateFunction('A', TutorialNumber('a'));
group = group.setStateFunction('B', TutorialNumber('b'));
group = group.setStateFunction('X', TutorialNumber('x'));
group = group.setStateFunction('Y', TutorialNumber('y'));
% Create products
group = group.setStateFunction('AB', TutorialProduct('A', 'B'));
group = group.setStateFunction('XY', TutorialProduct('X', 'Y'));
% Create the final sum
group = group.setStateFunction('G', TutorialAdd('AB', 'XY'));
disp(group)
%%
if ~exist('n', 'var')
    n = 1e6;
end
if ~exist('its', 'var')
    its = 100;
end
useMex = true;
subset = 1:100:n; % Index into 1% of the values, at non-sequential positions
subset_seq = 1:numel(subset);
a = rand(n, 1);
b = rand(n, 1);
x = rand(n, 1);
y = rand(n, 1);

% types = {'double', 'Sparse AD', 'Diagonal AD', 'Diagonal AD, row major'};
types = {'double', 'Sparse AD', 'Diagonal AD'};

nt = numel(types);
[t_assign, t_compute, t_subset, t_subcompute] = deal(zeros(nt, 1));
opt = struct('useMex', useMex, 'types', []);
for i = 1:nt
    if strcmpi(types{i}, 'double')
        % We can use the values directly
        state0 = struct('a', a, 'b', b, 'x', x, 'y', y); % Initial state
    else
        switch types{i}
            case 'Sparse AD'
                [aAD, bAD, xAD, yAD] = initVariablesADI(a, b, x, y);
            case 'Diagonal AD'
                [aAD, bAD, xAD, yAD] = initVariablesAD_diagonal(a, b, x, y);
            case 'Diagonal AD, row major'
                [aAD, bAD, xAD, yAD] = initVariablesAD_diagonalRowMajor(a, b, x, y, opt);
            otherwise
                error('Bad name')
        end
        state0 = struct('a', aAD, 'b', bAD, 'x', xAD, 'y', yAD); % Initial state
    end
    for it = 1:its
        state = group.initStateFunctionContainer(state0); % Enable caching
        tic()
        G = group.get([], state, 'G');  
        t_compute(i) = t_compute(i) + toc();

        tic();
        G_sub = G(subset);
        t_subset(i) = t_subset(i) + toc();

        tic()
        v = 2*G_sub;
        t_subcompute(i) = t_subcompute(i) + toc();

        tic();
        G(subset) = v;
        t_assign(i) = t_assign(i) + toc();
    end
end
data = [t_compute, t_subset, t_subcompute, t_assign]/its;
opnames = {'Compute', 'Subset', 'Subcompute', 'Insertion'};
%%
figure(1); clf
bar(data)
set(gca, 'XTickLabel', types)
legend(opnames)
% set(gca, 'yscale', 'log')
%%
v = data(2:end, :)./data(1, :);

figure(2); clf; hold on
h = bar(v);
set(gca, 'XTickLabel', types(2:end))
set(gca, 'XTick', 1:size(v, 2))

for i = 1:size(v, 2)
    for j = 1:size(v, 1)
        vij = v(j, i);
        x = h(i).XEndPoints(j);
        y = h(i).YEndPoints(j);
        text(x, y, sprintf('%2.1f', vij), 'HorizontalAlignment', 'center', ...
                                          'VerticalAlignment',   'bottom');
    end
end
legend(opnames)
yl = ylim();
ylim([0, yl(2)])
ylabel('Time relative to double')
% set(gca, 'yscale', 'log')
%%
v = data./data(1, :);

figure('Position', [838.3333  745.0000  717.3334  179.3333  ]); clf; hold on
h = bar(data);
set(gca, 'XTickLabel', types)
set(gca, 'XTick', 1:size(v, 2))

for i = 1:size(v, 2)
    for j = 1:size(v, 1)
        vij = v(j, i);
        x = h(i).XEndPoints(j);
        y = h(i).YEndPoints(j);
        text(x, y, sprintf('%2.1f', vij), 'HorizontalAlignment', 'center', ...
                                          'VerticalAlignment',   'bottom');
    end
end
legend(opnames)
yl = ylim();
ylim([0, yl(2)])
ylabel('Time [s]')
% set(gca, 'yscale', 'log')
%%
v = (1:(5*10))';
[mod(v-1, 5)+1, v]

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
