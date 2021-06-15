function [state, eqs, meta] = solveLocalEqs(state0, state, dt, G, W, s, fluid, system, meta, varargin)
% Solves the equations only for cells where the residual is large.

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

    opt = struct('gmresflag', 0, 'gmresits', []);
    opt = merge_options(opt, varargin{:});


    % Get the equations given current states

    [eqs, state, history] = eqsfiBlackOil(state0, state, dt, G, W, s, fluid, system, 'history', meta.history);
    [converged, CNV, MB] = getConvergence(state, eqs, fluid, system, dt);
    [residuals, residuals] = getResiduals(meta, eqs, system, 1);       %#ok

    fprintf('Residuals before starting local solver\n');
    eqnnames = {'Oil', 'Water', 'Gas', 'Dis.Gas', 'qWs', 'qOs', 'qGs', 'pBHP'};
    printResidual(residuals, opt.gmresits, eqnnames, meta.iteration, CNV, MB);

    fprintf('Global iteration number: %d\n', meta.iteration);

    % We solve local equations for cells with highest residual.
    % We skip first some iterations,

    its=0;
    maxit = system.stepOptions.maxitLocal;
    system2 = system;
    system2.stepOptions.solveWellEqs = false;
    [eqs, state] = eqsfiBlackOil(state0, state, dt, G, W, s, fluid, system2);

    % For each residual, we select the 50 cells with the largest residuals.
    nc = numel(state.pressure);
    indxs = cell(3,1);
    % pv = system.s.pv;
    % pvsum = sum(pv);

    % BW = fluid.BW(state.pressure);
    % BW_avg = sum(BW)/nc;
    % BO = fluid.BO(state.pressure, state.rs, state.s(:,3)>0);
    % BO_avg = sum(BO)/nc;
    % BG = fluid.BG(state.pressure);
    % BG_avg = sum(BG)/nc;

    res = computeResidual(eqs, dt);

    indxs{1} = sortrows([res.mass(:,1), (1:nc)'], -1);
    indxs{2} = sortrows([res.mass(:,2), (1:nc)'], -1);
    indxs{3} = sortrows([res.mass(:,3), (1:nc)'], -1);

    M = 10;
    indx = indxs{1}(1:M,2);
    for i = 2:3
        indx = union(indx, indxs{i}(1:M,2));
    end
    % We convert to logical (maybe not necessary!)
    indx_log = false*ones(nc, 1);
    indx_log(indx) = true;
    indx = logical(indx_log);
    indx = or(indx, meta.ref_cells);
    meta.ref_cells = indx;


    % % We select cells using a threshold.
    % max_res = cellfun(@(x) norm(x.val, 'inf'), eqs);
    % indx =  abs(double(eqs{1})) > 0.1*max_res(1);
    % for i = 2 : 4
    %     indx = or(abs(double(eqs{i})) > 0.1*max_res(i), indx);
    % end

    fprintf('Number of cells in local iterations: %d\n', nnz(indx));


    % We compute the local residual equations
    system2.stepOptions.solveWellEqs = false;
    system2.stepOptions.dsMax = 0.1;
    [local_eqs, state] = eqsfiBlackOilLocal(indx, state0, state, dt, G, W, s, fluid, system2);
    res1 = computeLocalResidual(local_eqs, dt);
    state1 = state;

    % res_diff = inf;
    res.norm = inf;

    while (res.norm > 1e-7) && (its < maxit)

        local_eqs0 = local_eqs;

        % We solve the local equations.
        dx_indx = SolveEqsADI(local_eqs, []);

        % We set up the global dx.
        for i = 1 : 4
            dx{i} = zeros(nc, 1);
            if i > 1
                dx{i}(indx) = dx_indx{i - 1};
            end
        end
        for i = 5 : 8
            dx{i} = dx_indx{i - 1};
        end

        % We update the system.
        [state, nInc] = updateState(W, state, dx, fluid, system2);

        % We compute the local residual equations.
        [local_eqs, state] = eqsfiBlackOilLocal(indx, state0, state, dt, G, W, s, fluid, ...
                                                system2);

        % We estimate the improvement done in this iteration.
        res_diff = norm(double(vertcat(local_eqs{:}) - vertcat(local_eqs0{:})), inf);

        %  rs changes so we have to update BO.
        % BO = fluid.BO(state.pressure, state.rs, state.s(:,3)>0);
        % BO_avg = sum(BO)/nc;
        res = computeLocalResidual(local_eqs, dt);

        fprintf('Local interation number %d, res_diff = %d, value of residual = %d\n', its, ...
                res_diff, res.norm);
        its = its + 1;

    end

    fprintf('Number of local iterations: %d\n', its);

    if res1.norm < res.norm
        state = state1;
        fprintf('Local iterations did not give any improvement! We drop them\n');
    end

    [eqs, state, history] = eqsfiBlackOil(state0, state, dt, G, W, s, fluid, system, 'history', meta.history);
    [converged, CNV, MB] = getConvergence(state, eqs, fluid, system, dt);
    [residuals, residuals] = getResiduals(meta, eqs, system, 1);       %#ok

    fprintf('Residuals after running local solver\n');
    eqnnames = {'Oil', 'Water', 'Gas', 'Dis.Gas', 'qWs', 'qOs', 'qGs', 'pBHP'};
    printResidual(residuals, opt.gmresits, eqnnames, meta.iteration, CNV, MB);

end

function res = computeResidual(eqs, dt)
    res.mass = [dt*abs(eqs{1}.val), dt*abs(eqs{2}.val), dt*abs(eqs{3}.val), dt*abs(eqs{4}.val)];
    res.wells = vertcat(eqs{5:end});
    res.wells = res.wells.val;
    res.norm = max(norm(res.mass, inf), norm(res.wells, inf));
end

function res = computeLocalResidual(local_eqs, dt)
    res.mass = [dt*abs(local_eqs{1}.val), dt*abs(local_eqs{2}.val), dt*abs(local_eqs{3}.val)];
    res.wells = vertcat(local_eqs{4:end});
    res.wells = res.wells.val;
    res.norm = max(norm(res.mass, inf), norm(res.wells, inf));
end

