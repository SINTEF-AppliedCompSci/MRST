function obj = matchObservedTwinModelsOW_saved(model_1, states_1, schedule_1, model_2, states_2, schedule_2,observed, bd_cells,bd_faces,well_cells,well_data,varargin)
% Compute mismatch-function 

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
opt     = struct('WaterRateWeight',     [] , ...
                 'OilRateWeight',       [] , ...
                 'GasRateWeight',       [] , ...
                 'BHPWeight',           [] , ...
                 'ComputePartials',     false, ...
                 'tStep' ,              [], ...
                 'state1',               [],...
                 'state2',               [],...
                 'from_states',         false,...% can be false for generic models
                 'matchOnlyProducers',  false, ...
                 'mismatchSum',         true, ...
                 'EnergyMinimization', true, ...
                 'accumulateWells',       [], ...
                 'accumulateTypes',       []);
             
opt     = merge_options(opt, varargin{:});

dts   = schedule_1.step.val;
totTime = sum(dts);
nc = model_1.G.cells.num;
tSteps = opt.tStep;
if isempty(tSteps) %do all
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    sol_obs = observed{tSteps(step)};
    nw      = numel(sol_obs.wellSol);
    if opt.matchOnlyProducers
        matchCases = (vertcat(sol.sign) < 0);
    else
        matchCases = true(nw,1);
    end
    qWs_obs = vertcatIfPresent(sol_obs.wellSol, 'qWs', nw);
    qOs_obs = vertcatIfPresent(sol_obs.wellSol, 'qOs', nw);
    qGs_obs = vertcatIfPresent(sol_obs.wellSol, 'qGs', nw);
    bhp_obs = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
    status_obs = vertcat(sol_obs.wellSol.status);

    [ww, wo, wg, wp] = getWeights(qWs_obs, qOs_obs, qGs_obs, bhp_obs, opt);

    if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state_1 = model_1.getStateAD( states_1{tSteps(step)}, init);
            state_2 = model_2.getStateAD( states_2{tSteps(step)}, init);

        else
            state_1 = opt.state1;
            state_2 = opt.state2;

        end

        % Extract well and field variables for state_1
        if isa(state_1, 'ADI')
            qWs_1 = model_1.FacilityModel.getProp(state_1, 'qWs');
            qOs_1 = model_1.FacilityModel.getProp(state_1, 'qOs');
            qGs_1 = model_1.FacilityModel.getProp(state_1, 'qGs');
            bhp_1 = model_1.FacilityModel.getProp(state_1, 'bhp');

            pressure_1 = model_1.FacilityModel.getProp(state_1, 'pressure');
            saturation_1 = model_1.FacilityModel.getProp(state_1, 's');

            assert(isnumeric(qWs_1));
            status_1 = vertcat(state_1.wellSol.status);
        else
            state_1 = states_1{tSteps(step)};
            [qWs_1, qOs_1, qGs_1, bhp_1] = deal( ...
                vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
                vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
                vertcatIfPresent(state_1.wellSol, 'qGs', nw), ...
                vertcatIfPresent(state_1.wellSol, 'bhp', nw) );

            [pressure_1, saturation_1] = model_1.getProps(state_1,'pressure','s');
            
            assert(isnumeric(qWs_1));
            status_1 = vertcat(state_1.wellSol.status);
        end

        % Extract well and field variables for state_2
        if isa(state_2, 'ADI')
            qWs_2 = model_2.FacilityModel.getProp(state_2, 'qWs');
            qOs_2 = model_2.FacilityModel.getProp(state_2, 'qOs');
            qGs_2 = model_2.FacilityModel.getProp(state_2, 'qGs');
            bhp_2 = model_2.FacilityModel.getProp(state_2, 'bhp');

            pressure_2 = model_2.FacilityModel.getProp(state_2, 'pressure');
            saturation_2 = model_2.FacilityModel.getProp(state_2, 's');

            assert(isnumeric(qWs_2));
            status_2 = vertcat(state_2.wellSol.status);
        else
            state_2 = states_2{tSteps(step)};
            [qWs_2, qOs_2, qGs_2, bhp_2] = deal( ...
                vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
                vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
                vertcatIfPresent(state_2.wellSol, 'qGs', nw), ...
                vertcatIfPresent(state_2.wellSol, 'bhp', nw) );


            [pressure_2, saturation_2] = model_2.getProps(state_2,'pressure','s');

            assert(isnumeric(qWs_2));
            status_2 = vertcat(state_2.wellSol.status);
        end        
        


    else
        % Extract well and field variables for state_1 (non-ADI case)
        state_1 = states_1{tSteps(step)};
        [qWs_1, qOs_1, qGs_1, bhp_1] = deal( ...
            vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
            vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
            vertcatIfPresent(state_1.wellSol, 'qGs', nw), ...
            vertcatIfPresent(state_1.wellSol, 'bhp', nw));

        [pressure_1, saturation_1] = model_1.getProps(state_1,'pressure','s');

        assert(isnumeric(qWs_1));
        status_1 = vertcat(state_1.wellSol.status);

        % Extract well and field variables for state_2 (non-ADI case)
        state_2 = states_2{tSteps(step)};
        [qWs_2, qOs_2, qGs_2, bhp_2] = deal( ...
            vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
            vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
            vertcatIfPresent(state_2.wellSol, 'qGs', nw), ...
            vertcatIfPresent(state_2.wellSol, 'bhp', nw));


        [pressure_2, saturation_2] = model_2.getProps(state_2,'pressure','s');

        assert(isnumeric(qWs_2));
        status_2 = vertcat(state_2.wellSol.status);

    end
     
    if ~all(status_1) || ~all(status_2) ||~all(status_obs)
        [bhp_1, bhp_obs] = expandToFull(bhp_1, bhp_obs, status_1, status_obs, true);
        [qWs_1, qWs_obs] = expandToFull(qWs_1, qWs_obs, status_1, status_obs, false);
        [qOs_1, qOs_obs] = expandToFull(qOs_1, qOs_obs, status_1, status_obs, false);

        [bhp_2, bhp_obs] = expandToFull(bhp_2, bhp_obs, status_2, status_obs, true);
        [qWs_2, qWs_obs] = expandToFull(qWs_2, qWs_obs, status_2, status_obs, false);
        [qOs_2, qOs_obs] = expandToFull(qOs_2, qOs_obs, status_2, status_obs, false);
    end
    if opt.EnergyMinimization
        %% compute an energy norm in the global domain to add to the objective function
        cTotFlux1 = states_1{step}.FacilityFluxProps.ComponentTotalFlux;
        cTotFlux2 = states_2{step}.FacilityFluxProps.ComponentTotalFlux;
        cTotMass1 = states_2{step}.FlowProps.ComponentTotalMass;
        cTotMass2 = states_2{step}.FlowProps.ComponentTotalMass;
        Div = model_1.operators.Div;
        for i = 1:length(cTotMass1)
            L2_error = dts(step)./ (totTime).*(sum(Div(cTotFlux1{i}-cTotFlux2{i})).^2+ sum((cTotMass1{i}-cTotMass2{i})./dts(step)).^2);
        end
    end


    dt = dts(step);

    if opt.mismatchSum
      dtFactor  = dt / (totTime * nnz(matchCases));  
     obj{step} = dtFactor*sum( ...
                        (ww*matchCases.*(qWs_1-qWs_2)).^2 + ...
                        (wo*matchCases.*(qOs_1-qOs_2)).^2 + ...
                        (wp*matchCases.*(bhp_1-bhp_2)).^2 );
    else
        % output summands f_i^2 
        fac = dt/(totTime*nnz(matchCases));
        mm  = {fac*(ww*matchCases.*(qWs_1-qWs_2)).^2, ...
               fac*(wo*matchCases.*(qOs_1-qOs_2)).^2, ...true
               fac*(wp*matchCases.*(bhp_1-bhp_2)).^2};
        
        if isempty(opt.accumulateTypes)
            tmp = mm;
        else
            % sum squares of qWs/qOs/bhp
            pt = opt.accumulateTypes;
            tmp = num2cell(zeros(1, max(pt)));
            for k = 1:3
                if pt(k)>0
                    tmp{pt(k)} = tmp{pt(k)} + mm{k};
                end
            end
        end
        if ~isempty(opt.accumulateWells)
            % sum squares of values for wells (use sparse mult)
            pw  = opt.accumulateWells;
            M   = sparse(pw(pw>0), find(pw), 1);
            tmp = applyFunction(@(x)M*x, tmp);
        end
        obj{step} = vertcat(tmp{:});
    end
end
end

%--------------------------------------------------------------------------

function v = vertcatIfPresent(sol, fn, nw)
    if isfield(sol, fn)
        v = vertcat(sol.(fn));
        % Ensure v has the correct size
        if numel(v) < nw
            v = [v; zeros(nw - numel(v), 1)];
        end
        % Apply status filtering
        v = v(vertcat(sol.status));
    else
        v = zeros(nnz(vertcat(sol.status)), 1);
    end
    % Final assertion to ensure correctness
    assert(numel(v) == nnz(vertcat(sol.status)), 'Mismatch in well count');
end


%--------------------------------------------------------------------------

function [v, v_obs] = expandToFull(v, v_obs, status, status_obs, setToZero)
tmp = zeros(size(status));
if isa(v, 'ADI')
    tmp = double2ADI(tmp, v);
end
tmp(status) = v;
v = tmp;
%
tmp = zeros(size(status));
tmp(status_obs) = v_obs;
v_obs = tmp;
if setToZero
    ix = status ~= status_obs;
    v(ix)     = 0;
    v_obs(ix) = 0;
end

end
%--------------------------------------------------------------------------

function [ww, wo, wg, wp] = getWeights(qWs, qOs, qGs, bhp, opt)
    ww = opt.WaterRateWeight;
    wo = opt.OilRateWeight;
    wg = opt.GasRateWeight;
    wp = opt.BHPWeight;

    rw = sum(abs(qWs)) + sum(abs(qOs)) + sum(abs(qGs));

    if isempty(ww)
        % Set to zero if all water rates are zero
        if sum(abs(qWs)) == 0
            ww = 0;
        else
            ww = 1 / rw;
        end
    end

    if isempty(wo)
        % Set to zero if all oil rates are zero
        if sum(abs(qOs)) == 0
            wo = 0;
        else
            wo = 1 / rw;
        end
    end

    if isempty(wg)
        % Set to zero if all gas rates are zero
        if sum(abs(qGs)) == 0
            wg = 0;
        else
            wg = 1 / rw;
        end
    end

    if isempty(wp)
        % Set to zero if all BHP values are the same
        dp = max(bhp) - min(bhp);
        if dp == 0
            wp = 0;
        else
            wp = 1 / dp;
        end
    end
end


