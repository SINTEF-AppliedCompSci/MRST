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