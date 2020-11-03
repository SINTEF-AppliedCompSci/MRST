%% Additive vs. multiplicative nonlinear domain decompositioning
% In this example, we compare additive and multiplicative nonlinear domain
% decomposition (NLDD), and show how we can obtain optimal ordering for
% multiplicative NLDD based on topological sorting of the intercell flux
% graph
mrstModule add ad-core ad-props ad-blackoil sequential example-suite ...
    mrst-gui compositional
mrstVerbose on

%% Get example
% We consider a slightly modified version of an example from Klemetsdal et.
% al, SPE RSC 2019, doi: 10.2118/193934-ms. The setup consists of a
% quarter five-spot example with CO2, Metahne, and nDecane.
example0 = MRSTExample('qfs_compositional', 'nsteps', 60); example = example0;
% We will use NLDD for the transport subproblem, so we construct a
% sequential pressure-transport model
mrstModule add linearsolvers
psolve = AMGCLSolverAD('tolerance', 1e-4);
example.model = getSequentialModelFromFI(example0.model, 'pressureLinearSolver', psolve);

%% Simulate reference problem
problem = example.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true);
simulatePackedProblem(problem);

%% Inspect results
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);

%%  Set up NLDD
% We construct disk-segment partition to define the subdomains
mrstModule add domain-decomposition coarsegrid agglom
x = example.model.G.cells.centroids;
r = linspace(150,sqrt(2)*1000,7);
partition1 = sum(sum(x.^2,2) < r.^2,2);
t = linspace(0,pi/2,6);
partition2 = sum(x(:,2)./x(:,1) < tan(t),2);
[~, ~, partition] = unique([partition1, partition2], 'rows');
partition(partition1 == 7) = max(partition) + 1;
partition = mergeBlocksByConnections(example.model.G, partition, example0.model.operators.T, 10);
% Randomize order (for illustrational purposes)
rng(2020), perm = randperm(max(partition))'; partition = perm(partition);
example.plot(partition);
% To ensure convergence of the full problem, it is smart to impose a
% slightly stricter convergence tolerance for the subdomain problems. We
% do this with the optional input argument ''subdomainTolerances'', where
% the tolerance names and corresponding reduction factors are given in a
% cell array
subtol = {'nonlinearTolerance', 0.5};
% Subdomain solves can also give verbose outputs by passing 1 (some) or 2
% (everything) to the optional input argument ''verboseSubmodel''
ddargs = {'subdomainTol', subtol, 'verboseSubmodel', 1};

%% Adaptive NLDD
% In adaptive NLDD, we solve each subdomain keeping all other subdomains
% fixed at the previous solution
exampleADD = example;
exampleADD.model.transportModel = DomainDecompositionModel(exampleADD.model.transportModel, partition, ddargs{:});
exampleADD.name = [exampleADD.name, '-add'];

%% Simulate
problemADD = exampleADD.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemADD, 'prompt', true);
simulatePackedProblem(problemADD);

%% Multiplicative NLDD
% In multiplicative NLDD, we solve the subdomains in the order defined by
% the partition vector, with already solved subdomains fixed at their
% updated solution
exampleMDD = example;
exampleMDD.model.transportModel = DomainDecompositionModel(exampleMDD.model.transportModel, partition, ...
                                                              ddargs{:}, 'strategy', 'multiplicative');
exampleMDD.name = [exampleMDD.name, '-mdd'];

%% Simulate
problemMDD = exampleMDD.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemMDD, 'prompt', true);
simulatePackedProblem(problemMDD);

%% Dynamic partitions
% The effectiveness of multiplicative NLDD will obiously depend on the
% order in which we solve the subdomains. The optimal order will generally
% depend on the reservoir state, and may therefore change from one timestep
% to the next. Partitions are implementet in a dedicated Partition class,
% which is added to the DomainDecompositionModel
dummyPartition = Partition();
disp(dummyPartition);
% During the call to model.prepareTimestep, the partition is computed (or
% simply returned, if it is static). The the class has a property
% ''numBlocks'', which is the target number of blocks in the partition. In
% addition it has a property ''wellPadding'', saying how many cells must be
% around each well (wells must be at least one cell a boundary between
% subdomains),. The class also implemnets routines for padding the wells.
% The actual partition computation is implemented in the mehtod ''compute''

%% Multiplicative NLDD with topological ordering
% We can significantly reduce the number of iterations for multiplicative
% NLDD by solving the subdomains in the direction of flow, as demonstrated
% for this example in Klemetsdal et. al, SPE RSC 2019. In fact, due to the
% upstream discretization, the transport equations will be converged after
% one sweep through the subdomains in absence of capillary and buoyancy
% effects. We define a dynamic TopologicalFluxPartition, which at each
% timestep gets a topological ordering of the grid cells based on the
% intercell fluxes from the presure step, and partitions the grid into a
% number of topologically sorted blocks
mrstModule add matlab_bgl
% We use the same number of subdomains (blocks) as in the static partition
topoPartition = TopologicalFluxPartition('numBlocks', max(partition));
exampleMDD_topo   = example;
exampleMDD_topo.model.transportModel = DomainDecompositionModel(exampleMDD_topo.model.transportModel, topoPartition, ...
                                                            ddargs{:}, 'strategy', 'multiplicative');
exampleMDD_topo.name = [exampleMDD_topo.name, '-mdd-topo'];

%% Simulate
problemMDD_topo = exampleMDD_topo.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemMDD_topo, 'prompt', true);
simulatePackedProblem(problemMDD_topo);

%% Multiplicative NLDD with coarse topological ordering
% We can also choose to topologically sort the original partition by using
% UpscaledPartition, which wraps around a Partition, and additionally takes
% in a base partition and the model. When we construct a partition,
% UpscaledPartition first upscales the state according to the base
% partition, reorders the coarse blocks of the base partition, and maps
% this partition onto the fine grid
topoPartition.numBlocks = inf; % Use all blocks in the partition
topoPartition.wellPadding = 0; % No need to pad around wells now!
blockTopoPartition = UpscaledPartition(topoPartition, example.model.pressureModel.parentModel, partition);
exampleMDD_block_topo   = example;
exampleMDD_block_topo.model.transportModel = DomainDecompositionModel(exampleMDD_block_topo.model.transportModel, blockTopoPartition, ...
                                                                      ddargs{:}, 'strategy', 'multiplicative');
exampleMDD_block_topo.name = [exampleMDD_block_topo.name, '-mdd-block-topo'];

%% Simulate
problemMDD_block_topo = exampleMDD_block_topo.getPackedSimulationProblem();
clearPackedSimulatorOutput(problemMDD_block_topo, 'prompt', true);
simulatePackedProblem(problemMDD_block_topo);

%% Inspect solutions
[allWellSols, allStates, allReports] ...
    = getMultiplePackedSimulatorOutputs({problem, problemADD, problemMDD, problemMDD_topo, problemMDD_block_topo}, ...
                                        'readFromDisk', false, 'readReportsFromDisk', true, 'readWellSolsFromDisk', true);
solver = 4; % Choose solver
example.plot(allStates{solver});

%% Compare nonlinear transport iterations
its = cellfun(@(r) getReportOutput(r, 'solver', 'TransportSolver'), allReports, 'UniformOutput', false);
itsTot = cellfun(@(it) sum(it.total), its);
figure('Position', [0,0,800,400]);
bar(itsTot); h = line([0, 6], repmat(numel(example.schedule.step.val), 1, 2), 'color', 'r', 'linew', 2);
set(gca, 'XTickLabel', {'Global', 'Additive', 'Multiplicative', 'Mult topo', 'Mult block topo'}, ...
         'XTickLabelRotation', -25, ...
         'FontSize', 13);
legend(h, 'Timesteps'); title('Nonlinear transport iterations')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
