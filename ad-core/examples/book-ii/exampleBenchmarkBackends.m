%% Benchmarking performance of different AD backends
% Switching between backends can yield significant performance
% improvements. In this example, we demonstrate how to systematically
% compare the performance for different operations.
mrstModule add ad-core
%%
if ~exist('dim', 'var')
    dim = [100, 100, 100];
end

if ~exist('it', 'var')
    it = 10;
end

if ~exist('bz', 'var')
    bz = 5;
end

mrstVerbose off
% Sparse, multiple blocks
sparseBlocks = SparseAutoDiffBackend();
% Sparse, single large b lock
sparseSingleBlock = SparseAutoDiffBackend('useBlocks', false);
% Diagonal, column major 
diagCol = DiagonalAutoDiffBackend('useMex', false, 'rowMajor', false);
% Diagonal, column major (with C++ acceleration)
diagColMex = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', false);
% Diagonal, row major
diagRow = DiagonalAutoDiffBackend('useMex', false, 'rowMajor', true);
% Diagonal, row major (with C++ acceleration)
diagRowMex = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
%% Perform benchmark
% It is possible to uncomment backends to add them to the benchmark
models = {};

backends = {};
backends{end+1} = sparseBlocks;
% backends{end+1} = sparseSingleBlock;
% backends{end+1} = diagCol;
backends{end+1} = diagColMex;
% backends{end+1} = diagRow;
backends{end+1} = diagRowMex;

[results, models, info] = benchmarkAutoDiffBackends(dim, backends, ....
                    'iterations', it, 'models', models, 'verbose', true, 'block_size', bz);

%% Plot a selected subset of tests
ops = {};
ops{end+1} = 'cell_xy';
ops{end+1} = 'face_xy';

% ops{end+1} = 'subset_small';
ops{end+1} = 'subset_large';
ops{end+1} = 'subasgn_large';

ops{end+1} = 'faceavg';
ops{end+1} = 'Grad';
ops{end+1} = 'upw';


ops{end+1} = 'Div';
ops{end+1} = 'AccDiv';

no = numel(ops);
nb = numel(backends);
data = nan(no, nb);
names = cell(no, 1);
for i = 1:no
    tmp = ops{i};
    switch tmp
        case 'cell_xy'
            tmp = 'x.*y (cell)';
        case {'subset_small', 'subset_large'}
            tmp = 'x(subs)';
        case {'subasgn_large'}
            tmp = 'x(subs)=y';
        case 'faceavg'
            tmp = 'favg(x)';
        case 'upw'
            tmp = 'upw(x)';
        case 'Grad'
            tmp = 'grad(x)';
        case 'face_xy'
            tmp = 'x.*y (on face)';
        case 'Div'
            tmp = 'div(x)';
        case 'AccDiv'
            tmp = 'x + div(y)';
        otherwise
            fprintf('Did not map %s\n', tmp);
    end
    names{i} = tmp;
end
for b = 1:nb
    r = results{b};
    for o = 1:no
        op = ops{o};
        t = r.(op).t_wall;
        data(o, b) = mean(t);
    end
end
intp = 'none';
figure(1); clf
set(gcf, 'position', [680   701   972   277]);
bar(data)
set(gca, 'XTickLabel', names, 'TickLabelInterpreter', intp)
ylabel('Time [s]')
names = cellfun(@(x) x.getBackendDescription, backends, 'UniformOutput', false);
legend(names);

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
