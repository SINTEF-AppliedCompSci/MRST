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
            assert(not(isnumeric(qWs_1)));
        else
            state_1 = states_1{tSteps(step)};        
            [qWs_1, qOs_1, bhp_1] = deal( vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state_1.wellSol, 'bhp', nw) );
        end


        if isa(state_2.pressure, 'ADI')
        
            qWs_2 = model_2.FacilityModel.getProp(state_2,'qWs'); 
            qOs_2 = model_2.FacilityModel.getProp(state_2,'qOs');
            bhp_2 = model_2.FacilityModel.getProp(state_2,'bhp'); 
            assert(not(isnumeric(qWs_2))); 
           
        else
            state_2 = states_2{tSteps(step)};        
            [qWs_2, qOs_2, bhp_2] = deal( vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state_2.wellSol, 'bhp', nw) );
        end

        status_1 = vertcat(state_1.wellSol.status);
        status_2 = vertcat(state_2.wellSol.status);

     else
        state_1 = states_1{tSteps(step)};
        [qWs_1, qOs_1, bhp_1] = deal( vertcatIfPresent(state_1.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state_1.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state_1.wellSol, 'bhp', nw) );
       assert(isnumeric(qWs_1));
       status_1 = vertcat(state_1.wellSol.status);

       state_2 = states_2{tSteps(step)};
       [qWs_2, qOs_2, bhp_2] = deal( vertcatIfPresent(state_2.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state_2.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state_2.wellSol, 'bhp', nw) );
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
    % one can satbilize optimization also using pressure  in connection cells
     dp =[];
     ds  =[];
     dp_cells = [];
     ds_cells = [];
     cells_data = 1:1:model_1.G.cells.num;
    if isa(state_1.pressure, 'ADI')&&isa(state_2.pressure, 'ADI')
    for w = 1:length(schedule_1.control(1).W)
        cells_w_step_schedule1 = schedule_1.control(schedule_1.step.control(step)).W(w).cells;
        cells_w_step_schedule2 = schedule_2.control(schedule_2.step.control(step)).W(w).cells;
        nbc = length(cells_w_step_schedule1);
        dp =[dp;sum(state_1.pressure(cells_w_step_schedule1)-state_2.pressure(cells_w_step_schedule2))./nbc];
        ds =[ds;sum(state_1.s{2}(cells_w_step_schedule1)-state_2.s{2}(cells_w_step_schedule2))./nbc]; 
    end
       dp_cells = model_2.operators.Grad(state_1.pressure - state_2.pressure)./(sum(model_1.operators.Grad(value(state_1.pressure)).^2)).^0.5;    
       ds_cells = (state_1.s{2}(cells_data) - state_2.s{2}(cells_data))./sum(value(state_1.s(cells_data,2)).^2).^0.5;

    end
       
    if isa(state_1.pressure, 'ADI')&&~isa(state_2.pressure, 'ADI')
    for w = 1:length(schedule_1.control(1).W)
        cells_w_step_schedule1 = schedule_1.control(schedule_1.step.control(step)).W(w).cells;
        cells_w_step_schedule2 = schedule_2.control(schedule_2.step.control(step)).W(w).cells;
        nbc = length(cells_w_step_schedule1);
        dp =[dp;sum((state_1.pressure(cells_w_step_schedule1)-state_2.pressure(cells_w_step_schedule2)))./nbc];
        ds =[ds;sum(state_1.s{2}(cells_w_step_schedule1)-state_2.s(cells_w_step_schedule2,2))./nbc];        
    end

    dp_cells = model_1.operators.Grad(state_1.pressure-state_2.pressure)./(sum(model_1.operators.Grad(state_2.pressure).^2)).^0.5;
    ds_cells = (state_1.s{2}(cells_data)-state_2.s(cells_data,2))./sum(state_2.s(cells_data,2).^2).^0.5;
    end

       
    if ~isa(state_1.pressure, 'ADI')&&~isa(state_2.pressure, 'ADI')
    for w = 1:length(schedule_1.control(1).W)
        cells_w_step_schedule1 = schedule_1.control(schedule_1.step.control(step)).W(w).cells;
        cells_w_step_schedule2 = schedule_2.control(schedule_2.step.control(step)).W(w).cells;
        nbc = length(cells_w_step_schedule1);
        dp =[dp;sum((state_1.pressure(cells_w_step_schedule1)-state_2.pressure(cells_w_step_schedule2)))./nbc];
        ds =[ds;sum(state_1.s(cells_w_step_schedule1,2)-state_2.s(cells_w_step_schedule2,2))./nbc];        
    end

        dp_cells =model_1.operators.Grad(state_1.pressure-state_2.pressure)./norm(model_1.operators.Grad(state_1.pressure));
        ds_cells =(state_1.s(cells_data,2)-state_2.s(cells_data,2))./norm(state_1.s(cells_data,2));
    end

    if ~isa(state_1.pressure, 'ADI')&&isa(state_2.pressure, 'ADI')
    for w = 1:length(schedule_1.control(1).W)

        cells_w_step_schedule1 = schedule_1.control(schedule_1.step.control(step)).W(w).cells;
        cells_w_step_schedule2 = schedule_2.control(schedule_2.step.control(step)).W(w).cells;
        nbc = length(cells_w_step_schedule1);
        dp =[dp;sum((state_1.pressure(cells_w_step_schedule1)-state_2.pressure(cells_w_step_schedule2)))./nbc];
        ds =[ds;sum(state_1.s(cells_w_step_schedule1,2)-state_2.s{2}(cells_w_step_schedule2))./nbc];        
    end

        dp_cells =(model_2.operators.Grad(state_1.pressure-state_2.pressure))./norm(model_2.operators.Grad(state_1.pressure));
        ds_cells =(state_1.s(cells_data,2)-state_2.s{2}(cells_data))./norm((state_1.s(cells_data,2)));
    end
    epsPT= 1.0e-3;
    epsPv= 1.0e-3;
    dist_bd_pressure = (states_1{step}.pressure(bd_cells)-states_2{step}.pressure(bd_cells))./(states_1{step}.pressure(bd_cells));
    dist_bd_flux = states_1{step}.flux(bd_faces)-states_2{step}.flux(bd_faces);
    if opt.mismatchSum
        obj{step} = 0.5.*(dt/(totTime*nnz(matchCases)))*sum( ...
                         (ww*matchCases.*(qWs_1-(qWs_2))).^2 + ...
                         (wo*matchCases.*(qOs_1-(qOs_2))).^2 + ...
                         (wp.*matchCases.*(bhp_1-(bhp_2))).^2+...
                         0.01.*matchCases.*ds.^2+1.*wp.*wp.*matchCases.*dp.^ 2)+...
                         2000.5.*(dt/(totTime))*sum(dp_cells.^2)+200.5.*(dt/(totTime)).*sum(ds_cells.^2)+...
                         epsPv.*wp.*(dt/(totTime))*(sum((model_2.operators.Grad(model_1.operators.pv).^2)))+...
               epsPT.*wp.*(dt/(totTime))*(sum((abs(model_2.operators.Grad(model_1.operators.pv)))))+...
               +0.1.*wo.*(dt/(totTime*nnz(matchCases))).*norm(model_2.operators.pv(well_cells)-well_data)./norm(well_cells)+...
               0.5.*10.*(dt/(totTime*nnz(matchCases)))*sum(abs(abs(dist_bd_pressure.*dist_bd_flux)));
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

