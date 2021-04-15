function [W, wellHistory, history] = ...
      optimizeWellPlacementDiagnostics(G, W, rock, objective, ...
                                       targetWells, D, wlimit, ...
                                       state0, fluid_ad, pv, T, s, varargin)
%Optimize the placement of wells using flow diagnostics
%
% SYNOPSIS:
%  W = optimizeWellPlacementDiagnostics(G, W, rock, objective, targetWells,...
%                                         D, wlimit, state0, fluid_ad, pv, T, s)
%
% DESCRIPTION:
%  This function optimizes the controls of some set of wells using
%  fast-to-evalute diagnostic proxy models.
%
% REQUIRED PARAMETERS:
%   G       - Grid structure.
%
%   W       - Well structured as defined by addWell.
%
%   rock    - Rock with valied perm field. Used to calculate well indices
%             for new well positions.
%
%   objective - Objective function as defined by getObjectiveDiagnostics.
%
%   targetWells - Indices into W indicating which wells are to be optimized.
%               If the rates *and* placement are to be optimized, all
%               targets must be of the same kind (i.e. all injectors or
%               producers) as the algorithm assumes that the total injected
%               fluid is constant.
%
%   D           - Struct of diagnostics quantities. Unused.
%
%   wlimit      - Well limits.- Passed onto optimizeTOF if optimizeSubsteps
%                 is enabled.
%
%   state0      - Initial state as defined by for example initResSol.
%
%   fluid_ad    - Valid ad_fluid. See initSimpleADIFluid.
%
%   pv          - Pore volume.
%
%   T           - Half-transmissibilities as produced by computeTrans.
%
%   s           - AD system. See initADISystem.
%
% OPTIONAL PARAMETERS:
%
%  maxIter  - Maximum outer iterations. An outer iteration does a single
%             pass over all wells, moving them until they improve the
%             solution or they reach the wellSteps constraint.
%
%  optimizeSubsteps - Optimize the rates of the wells after each placement.
%                     This may increase the complexity *substatially*!
%
%  searchRadius     - How large logical region should be used to evaluate
%                     candidate wells. The larger the region, the faster
%                     the convergence for simple cases, but at the same
%                     time the optimization assumptions may be less valid
%                     the further away from current wells a ghost well is
%                     moved. Should be carefully adjusted.
%
%  searchIncrement  - How many ghost wells should be placed? As adding
%                     wells is a bit costly, this allows coarser grained
%                     ghost well placement inside the searchRadius. If the
%                     search radius is 10 in each ij-direction and
%                     searchIncrement is set to 5, 4*4 = 16 ghost wells
%                     will be equally spaced.
%
% RETURNS:
%
%  W     - Optimized wells.
%
%  wellHistory - numel(W) long array containing the visitation history of
%               each well for plotting or debugging.
%
%  history  - Optimization history (if optimizeSubsteps is active). One
%               entry per optimization leading to improvement.
%


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


    opt = struct('iterThreshold', inf, ...
                 'maxIter', inf, ...
                 'optimizeSubsteps', false, ...
                 'searchRadius', 1, ...
                 'searchIncrement', 1, ...
                 'wellSteps', inf, ...
                 'linsolve', @mldivide, ...
                 'plotProgress', true);
    opt = merge_options(opt, varargin{:});

    ijk = gridLogicalIndices(G);
    wellHistory = {};
    history = [];
    hist = [];
    targets = targetWells;

    psolve = @(W, state) ...
       solveStationaryPressure(G, state, s, W, fluid_ad, pv, T, ...
                               'objective', objective, 'computeTracer', ...
                               false, 'linsolve', opt.linsolve);

    optim = @(W, state) optimizeTOF(G,W, fluid_ad, pv, T, s,...
                         state, wlimit, objective, ...
                         'targets', targets, 'linsolve', opt.linsolve);

    state = state0;
    state.wellSol = initWellSol(W, 0);
    [state, D0, grad0] = psolve(W, state);                      %#ok<ASGLU>
    bestobj = grad0.objective.val;

    outIter = 0;
    while true,
        if outIter >= opt.maxIter
            break;
        end
        changed = false;
        h = struct('objective', cell([numel(targetWells), 1]), ...
                   'visited'  , cell([numel(targetWells), 1]));

        for i = 1:numel(targetWells)
            tw = targetWells(i);

            fprintf('Optimizing well %d...\n', tw);
            state0.wellSol = initWellSol(W, 0);

            if opt.optimizeSubsteps
                [D, W, hist] = optim(W, state0);                %#ok<ASGLU>
                bestobj = min(hist.gradient(end).objective.val, bestobj);
            end

            visited = [];
            objlist = [];

            iter = 1;
            while true
                state0.wellSol = initWellSol(W, 0);
                [W, success, bestobj, D] = ...
                   optimizeGhost(G, rock, state0, W, ijk, tw, ...
                                 bestobj, optim, psolve, opt);

                changed = changed || success;
                vc = W(tw).cells(end);
                if any(cellfun(@(x) any(x(i).visited == vc), wellHistory)) ||...
                        any(visited == vc)
                    break
                end

                visited = [visited; vc];
                objlist = [objlist; bestobj];
                if iter > opt.iterThreshold || ~success
                    disp(['Unable to improve ' num2str(W(tw).name) '...'])
                    break
                end
                iter = iter + 1;

                if opt.plotProgress
                    if ishandle(2)
                        set(0, 'CurrentFigure', 2); clf;
                    else
                        figure(2); clf;
                    end
                    ftof = log10(D.tof(:,1));
                    plotCellData(G, ftof, 'EdgeAlpha', 0.125);
                    if G.griddim == 3
                        plotWell(G, W, 'height', 1, 'color', 'w');
                    else
                        plotGrid(G, vertcat(W.cells), 'FaceColor', 'red');
                    end
                    axis tight off
                    title(['Objective value: ', num2str(bestobj)])
                    drawnow
                end
                if iter > opt.wellSteps
                    break
                end
            end
            h(i).objective = objlist;
            h(i).visited = visited;
        end

        wellHistory = [wellHistory; h];
        history = [history; hist];
        outIter = outIter + 1;

        if ~changed;
            disp('No changes in wells during last step, all done'), break;
        else
            disp('Restarting loop')
        end
    end
end

function W_ghost = addGhostWells(G, W_target, rock, ijk, W_all, opt)
    r = opt.searchRadius;

    wc_all = vertcat(W_all.cells);
    existing = [ijk{1}(wc_all), ijk{2}(wc_all)];

    wc = W_target.cells;
    ix = unique(ijk{1}(wc));
    jx = unique(ijk{2}(wc));
    if 0
        assert(numel(ix) == 1);
        assert(numel(jx) == 1);
    else
        ix = ix(1);
        jx = jx(1);
    end

    W_ghost = [];
    for i = max(ix-r, 1):opt.searchIncrement:min(ix+r, G.cartDims(1))
        for j = max(jx-r, 1):opt.searchIncrement:min(jx+r, G.cartDims(2))
            if i == ix && j == jx || any(existing(:, 1) == i & existing(:, 2) == j)
                continue
            end
            w = verticalWell([], G, rock, double(i), double(j), [], ...
                             'Type', 'rate', 'val', 0, 'Name', ...
                             'MrGhost', 'Sign', 1, 'refDepth', min(G.cells.centroids(:, 3)));
            if ~isempty(w.cells)
                W_ghost = [W_ghost; w];                         %#ok<AGROW>
            end
        end
    end
end

function v = getObj(history)
    v = history.gradient(end).objective.val;
end

function [W, success, obj, D] = ...
      optimizeGhost(G, rock, state0, W0, ijk, target, ...
                    obj0, optimize, psolve, opt)
        success = false;

        W = reshape(W0, [], 1);

        ind = false(numel(W),1);
        ind(target) = true;

        W_static = W(~ind);
        W_target = W(ind);

        % Add ghost wells
        W_ghost = addGhostWells(G, W_target, rock, ijk, W_static, opt);
        ghostInd = numel(W_static) + numel(W_target) + 1;

        W_all = vertcat(W_static, W_target, W_ghost);

        state = state0;
        state.wellSol = initWellSol(W_all, 0);
        [state, D, grad] = psolve(W_all, state);                %#ok<ASGLU>

        if isinf(obj0)
            obj0 = grad.objective.val;
        end

        % Pick the ghost well with the largest gradient.
        if 0
            % Largest gradient in absolute value...
            ghostGrad = abs(grad.well(ghostInd:end));
            [v gi] = max(ghostGrad);
        else
            % Smallest gradient in value (i.e. maximum negative)
            ghostGrad = (grad.well(ghostInd:end));
            [v, gi] = min(ghostGrad);                           %#ok<ASGLU>
        end

        W_ghost(gi).val = W_target.val;
        W_ghost(gi).name = W_target.name;

        W(target) = W_ghost(gi);

        state0.wellSol = initWellSol(W, 0);

        if opt.optimizeSubsteps
            % Optimize the rates of the new well
            [D, W_all, history] = optimize(W, state0);          %#ok<ASGLU>
            obj = getObj(history);
        else
            [state, D, grad] = psolve(W, state0);               %#ok<ASGLU>
            obj = grad.objective.val;
        end


        if obj <= obj0
            fprintf('Improved value, %2.3f to %2.3f\n', obj0, obj);
            success = true;
        else
            fprintf('Did not improve value, %2.3f to %2.3f\n', obj0, obj);
            obj = obj0;
            W = W0;
        end

end
