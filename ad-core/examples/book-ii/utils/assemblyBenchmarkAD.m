function [out, models] = assemblyBenchmarkAD(sizes, backends, physics, force_type, niter)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

if nargin < 5
    niter = 10;
end
if nargin < 4 || isempty(force_type)
    force_type = 'wells';
else
    assert(ischar(force_type));
end

if nargin < 3 || isempty(physics)
    physics = {'immiscible', 'blackoil', 'overall', 'natural'};
end
physics = wrapcell(physics);

if nargin < 2 || isempty(backends)
    backends = {'sparse', 'diagonal'};
end
backends = wrapcell(backends);

if nargin < 1 || isempty(sizes)
    sizes = [10, 20, 40, 50];
end

require ad-core ad-props ad-blackoil compositional libgeometry deckformat

nb = numel(backends);
nz = numel(sizes);
np = numel(physics);

[models, results] = deal(cell(nb, nz, np));

eos = getBenchmarkMixture('spe5');
pth = getDatasetPath('spe9');
fn  = fullfile(pth, 'BENCH_SPE9.DATA');
deck = readEclipseDeck(fn);
fluid_spe = initDeckADIFluid(deck, 'useMex', true);

fluid = initSimpleADIFluid('c', [1, 1, 1]*1e-5/barsa);
dt = 1*year;
% results = repmat


for zNo = nz:-1:1
    n = sizes(zNo);
    for bNo = 1:nb
        backend = backends{bNo};
        for pNo = 1:np
            phys = physics{pNo};
            [G, rock, W, state0] = getProblem(n, eos, fluid_spe, lower(phys));
            if strcmpi(force_type, 'wells')
                
            elseif strcmpi(force_type, 'none')
                disp('Disabling wells!')
                W = [];
            else
                error('Bad force %s', force_type);
            end
            forces = struct('W', W, 'bc', [], 'src', []);
            if isempty(models{bNo, zNo, pNo})
                gravity reset on
                switch lower(phys)
                    case '1ph'
                        model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', false', 'gas', false);
                    case 'immiscible'
                        gravity reset off
                        model = GenericBlackOilModel(G, rock, fluid);
                    case 'blackoil'
                        model = GenericBlackOilModel(G, rock, fluid_spe, 'disgas', true);
                    case 'overall'
                        model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
                    case 'natural'
                        model = GenericNaturalVariablesModel(G, rock, fluid, eos, 'water', false);
                    otherwise
                        error('Unknown physics %s', phys);
                end
                clear G
%                 model.G.cells = rmfield(model.G.cells, {'facePos', 'indexMap', 'faces', 'volumes'});
%                 model.G = rmfield(model.G, 'nodes');
%                 model.G.faces = rmfield(model.G.faces, {'centroids', 'neighbors', 'areas', 'nodePos', 'nodes', 'tag', 'normals'});
                switch lower(backend)
                    case 'sparse'
                        bend = AutoDiffBackend();
                    case 'diagonal'
                        bend = DiagonalAutoDiffBackend('useMex', true);
                    case 'diagonal-rowmajor'
                        bend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
                    case 'diagonal-rowmajor-deferred'
                        bend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true, 'deferredAssembly', true);
                    otherwise
                        error('Unknown backend %s', backend);
                end
                model.AutoDiffBackend = bend;
                if pNo > 1
                    model.operators = ops;
                end
                model = model.validateModel(forces);
                if pNo == 1
                    ops = model.operators;
                end
                if nargout > 1
                    models{bNo, zNo, pNo} = model;
                end
            else
                model = models{bNo, zNo, pNo};
            end
            nc = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            model.FacilityModel.outputFluxes = false;
            model.outputFluxes = false;
            % Validate state and make it so that the previous state has
            % total mass etc computed as if it was from an earlier
            % time-step
            state0 = model.validateState(state0);
            state = state0;
            s0 = getState0_cached(model, state, state0);
            time_eqs = zeros(1, niter);
            time_state = zeros(1, niter);
            for it = 1:niter
                fprintf('%d', it);
                [t_state, t_eq] = bench(model, bend, s0, state, dt, forces);
                time_state(it) = t_state;
                time_eqs(it) = t_eq;
            end
            fprintf('\n');
            time_state = median(time_state);
            time_eqs = median(time_eqs);
            timing = time_state + time_eqs;
            ncell = model.G.cells.num;
            ndof = numel(forces.W)*(nph+1) + nc*ncell;
            
            s = struct('ndof', ndof, 'ncell', ncell, 'time_avg', timing, 'time_eqs', time_eqs, 'time_state', time_state);
            
            results{bNo, zNo, pNo} = s;
            fprintf('%s (%s) with %d cells:\n\t%f us per cell\n\t%f us per dof\n\t%f s total (%1.2f%% init, %1.2f%% equations)\n', ...
                phys, backend, ncell, timing/(micro*ncell), timing/(micro*ndof), timing, 100*time_state/timing, 100*time_eqs/timing);
            clear s0 model
        end
        clear ops
    end
end

out.results = results;
out.backends = backends;
out.physics = physics;
end

function x = wrapcell(x)
    if ~iscell(x)
        x = {x};
    end
end

function [G, rock, W, state0] = getProblem(n, eos, fluid_spe, name)
    ncomp = numel(eos.names);
    G = cartGrid([n, n, n], [1000, 1000, 10]);
    G = mcomputeGeometry(G);
    
    rock = makeRock(G, 1*darcy, 1);
    
    W = [];
    W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0, 0]);
    W = verticalWell(W, G, rock, 1, n, [], 'comp_i', [1, 0, 0]);
    W = verticalWell(W, G, rock, n, 1, [], 'comp_i', [1, 0, 0]);
    W = verticalWell(W, G, rock, n, n, [], 'comp_i', [1, 0, 0]);
    N = getNeighbourship(G);
    T = ones(size(N, 1), 1);
    G.cells.eMap = ':';
    G = simGridTPFA(G, rock, 'neighbors', N, 'porv', ones(G.cells.num, 1), 'depth', G.cells.centroids(:, 3), 'trans', T);
    
    state0 = initResSol(G, 50*rand(G.cells.num, 1) + 50*barsa, [0.5, 1, 0]);
    state0 = rmfield(state0, 'flux');
    twoPhase = 1:floor(0.1*n^3);
    state0.s(twoPhase, :) = repmat([0.2, 0.4, 0.4], numel(twoPhase), 1);

    if strcmpi(name, 'blackoil')
        sat = fluid_spe.rsSat(state0.pressure);
        state0.rs = 0.5*sat;
        state0.rs(twoPhase) = sat(twoPhase);
    elseif strcmpi(name, 'overall') || strcmpi(name, 'natural')
        state0.s = state0.s(:, 2:3)./sum(state0.s(:, 2:3), 2);
        state0.components = ones(G.cells.num, ncomp)/ncomp;
        state0.x = state0.components;
        state0.y = 0.5*state0.components;
        state0.K = ones(G.cells.num, ncomp);
        state0.L = state0.s(:, 2);
        state0.T = repmat(350, G.cells.num, 1);
        state0.Z_L = ones(G.cells.num, 1);
        state0.Z_V = state0.Z_L;
        for i = 1:numel(W)
            W(i).components = ones(1, ncomp)/ncomp;
            W(i).comp_i = [0, 1];
        end
    end
    state0.s = bsxfun(@rdivide, state0.s, sum(state0.s, 2));
end

function [t_state, t_eq] = bench(model, bend, s0, s, dt, forces)
    t = tic();
    [s, primaryVars] = model.getStateAD(s);
    t_state = toc(t);
    fprintf('.');
    t = tic();
    [eqs, names, types] = model.getModelEquations(s0, s, dt, forces);
    t_eq = toc(t);
    if isa(bend, 'DiagonalAutoDiffBackend') && bend.deferredAssembly
        % Other solvers ignore the problem
        s.FluxProps = [];
        s.FlowProps = [];
        s.PVTProps = [];
        problem = LinearizedProblem(eqs, types, names, primaryVars, s);
        t = tic();
        % Need to assemble this as well
        [A, b] = bend.getBlockSystemCSR(problem);
        t_eq = t_eq + toc(t);
    end
end

function s0 = getState0_cached(model, state, state0)
    s0 = model.getStateAD(state0, false);
    % model.getProps(s0, 'ComponentTotalMass', 'PhaseFlux');
    model.getProps(s0, 'ComponentTotalMass');
    s0 = model.reduceState(s0, false);
    s0 = model.updateAfterConvergence(state, s0, 1, struct());
end
