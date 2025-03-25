function obj = matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2,observed, bd_cells,bd_faces,well_cells,well_data,varargin)
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
    'BHPWeight',           [] , ...
    'ComputePartials',     false, ...
    'tStep' ,              [], ...
    'state1',               [],...
    'state2',               [],...
    'from_states',         false,...% can be false for generic models
    'matchOnlyProducers',  false, ...
    'mismatchSum',         true, ...
    'accumulateWells',       [], ...
    'accumulateTypes',       []);

opt     = merge_options(opt, varargin{:});

dts   = schedule_1.step.val;
totTime = sum(dts);

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
    bhp_obs = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
    status_obs = vertcat(sol_obs.wellSol.status);

 [ww, wo, wp] = getWeights(qWs_obs, qOs_obs, bhp_obs, opt);
    if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state_1 = model_1.getStateAD( states_1{tSteps(step)}, init);
            state_2 = model_2.getStateAD( states_2{tSteps(step)}, init);

        else
            state_1 = opt.state1;
            state_2 = opt.state2;

        end

        if isa(state_1.pressure, 'ADI')

            qWs_1 = model_1.FacilityModel.getProp(state_1,'qWs');
            qOs_1 = model_1.FacilityModel.getProp(state_1,'qOs');
            bhp_1 = model_1.FacilityModel.getProp(state_1,'bhp');

            pressure_1 = model_1.getProp(state_1, 'pressure');
            saturation_1 = model_1.getProp(state_1, 's');

            assert(not(isnumeric(qWs_1)));
        else
            %             state_1 = states_1{tSteps(step)};
            [qWs_1, qOs_1, bhp_1] = deal( vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
                vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
                vertcatIfPresent(state_1.wellSol, 'bhp', nw) );

            [pressure_1, saturation_1] = deal(state_1.pressure,state_1.s);
        end


        if isa(state_2.pressure, 'ADI')

            qWs_2 = model_2.FacilityModel.getProp(state_2,'qWs');
            qOs_2 = model_2.FacilityModel.getProp(state_2,'qOs');
            bhp_2 = model_2.FacilityModel.getProp(state_2,'bhp');


            pressure_2 = model_2.getProp(state_2, 'pressure');
            saturation_2 = model_2.getProp(state_2, 's');

            assert(not(isnumeric(qWs_2)));

        else
            %             state_2 = states_2{tSteps(step)};
            [qWs_2, qOs_2, bhp_2] = deal( vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
                vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
                vertcatIfPresent(state_2.wellSol, 'bhp', nw) );

            [pressure_2, saturation_2] = deal(state_2.pressure,state_2.s);
        end

        status_1 = vertcat(state_1.wellSol.status);
        status_2 = vertcat(state_2.wellSol.status);

    else
        state_1 = states_1{tSteps(step)};
        [qWs_1, qOs_1, bhp_1] = deal( vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
            vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
            vertcatIfPresent(state_1.wellSol, 'bhp', nw) );
        [pressure_1, saturation_1] = deal(state_1.pressure,state_1.s);
        assert(isnumeric(qWs_1));
        status_1 = vertcat(state_1.wellSol.status);

        state_2 = states_2{tSteps(step)};
        [qWs_2, qOs_2, bhp_2] = deal( vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
            vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
            vertcatIfPresent(state_2.wellSol, 'bhp', nw) );
        [pressure_2, saturation_2] = deal(state_2.pressure,state_2.s);
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
    dt = dts(step);

    if ~iscell(pressure_1)
        pressure_1 = {pressure_1};
        saturation_1 = {saturation_1};
    end

    if ~iscell(pressure_2)
        pressure_2 = {pressure_2};
        saturation_2 = {saturation_2};
    end
%     if ~isa(qWs_1,'ADI')
%         [ww, wo, wp] = getWeights(qWs_1, qOs_1, bhp_1, opt);
%     else
%         [ww, wo, wp] = getWeights(value(qWs_1), value(qOs_2), value(bhp_1), opt);
%     end
% [ww, wo, wp] = [, max(qOs_1,qOs_2,eps),];
% ww = max(value(abs(qWs_1)),value(abs(qWs_2)));
% ww = 1./max(ww,eps);
% wo = max(value(abs(qOs_1)),value(abs(qOs_2)));
% wo = 1./max(wo,eps);
% wp =  1./max(value(bhp_1),value(bhp_2));



% Initialize error metrics
energyNormError = 0;
H1NormPressureError = 0;
L2NormPressureWellsError = 0;
interWellEnergyError = 0;

% Regularization parameters
regEnergy = 0.01;
regH1 = 0.001;

% Compute Energy Norm Error (Global Domain)
if true
    totalFlux1 = state_1.FluxDisc.ComponentTotalFlux;
    totalFlux2 = state_2.FluxDisc.ComponentTotalFlux;
    totalMass1 = state_1.FlowProps.ComponentTotalMass;
    totalMass2 = state_2.FlowProps.ComponentTotalMass;

    Div = model_1.operators.Div;

    for i = 1:length(totalMass1)
        massDiff = (totalMass1{i} - totalMass2{i}) ./ dts(step);
        fluxDivDiff = Div(totalFlux1{i} - totalFlux2{i});

        energyNormError = energyNormError + ...
            (dts(step) / totTime) .* (sum(fluxDivDiff.^2) + sum(massDiff.^2));
    end
end

% Compute H1 Norm for Pressure
if true
    Grad = model_1.operators.Grad;

    pressureDiff = pressure_1{1} - pressure_2{1};
    gradPressure = Grad(pressure_2{1});

    H1NormPressureError = H1NormPressureError + ...
        (dts(step) / totTime) * (sum(Grad(pressureDiff).^2) ./ sum(gradPressure.^2));
end

% Compute L2 Norm for Pressure Difference at Wells
if true
    wellCells = reshape([schedule_1.control(1).W.cells], [], 1);
    wellPressureDiff = (pressure_1{1}(wellCells) - pressure_2{1}(wellCells)) ./ pressure_2{1}(wellCells);

    L2NormPressureWellsError = L2NormPressureWellsError + ...
        (dts(step) / totTime) .* sum(wellPressureDiff.^2);
end

% Compute Energy Norm for Inter-Well Region
if true
    nw = numel(schedule_1.control(1).W);

    for wew = 1:nw

        facilityFlux1 = state_1.wellSol(wew).ComponentTotalFlux;
        facilityFlux2 = state_2.wellSol(wew).ComponentTotalFlux;
        wc = schedule_1.control(1).W(wew).cells;
        wellPressureDiff = (pressure_1{1}(wc) - pressure_2{1}(wc))./ pressure_2{1}(wc);

        for i = 1:length(totalMass1)
            interWellEnergyError = interWellEnergyError + ...
                (dts(step) / totTime) .* sum(abs((facilityFlux1(:,i) - facilityFlux2(:,i)) .* wellPressureDiff));
        end
    end
end


% Stabilize optimization using meaningful scaling factors
if opt.mismatchSum
    obj{step} = (dt/(totTime*nnz(matchCases)))*sum( ...
        (ww.*matchCases.*(qWs_1-qWs_2)).^2 + ...
            (wo.*matchCases.*(qOs_1-qOs_2)).^2 + ...
            (wp.*matchCases.*(bhp_1-bhp_2)).^2 ) + regEnergy.*(interWellEnergyError);

    else
     % output summands f_i^2
    fac = dt/(totTime*nnz(matchCases));
    mm  = {fac*(ww.*matchCases.*(qWs_1 - qWs_2)).^2, ...
        fac*(wo.*matchCases.*(qOs_1 - qOs_2)).^2, ...
        fac*(wp.*matchCases.*(bhp_1 - bhp_2)).^2, ...
         0.01.*interWellEnergyError, ...
         0.01.*L2NormPressureWellsError, ...
         0.01*H1NormPressureError, ...
        0.001*energyNormError};
    

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
    assert(numel(v)==nw);
    v = v(vertcat(sol.status));
else
    v = zeros(nnz(sol.status),1);
end
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

function  [ww, wo, wp] = getWeights(qWs, qOs, bhp, opt)
ww = opt.WaterRateWeight;
wo = opt.OilRateWeight;
wp = opt.BHPWeight;

rw = sum(abs(qWs)+abs(qOs));

if isempty(ww)
    % set to zero if all are zero
    if sum(abs(qWs))==0
        ww = 0;
    else
        ww = 1/rw;
    end
end

if isempty(wo)
    % set to zero if all are zero
    if sum(abs(qOs))==0
        wo = 0;
    else
        wo = 1/rw;
    end
end

if isempty(wp)
    % set to zero all are same
    dp = max(bhp)-min(bhp);
    if dp == 0
        wp = 0;
    else
        wp = 1/dp;
    end
end
end

