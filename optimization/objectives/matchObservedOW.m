function obj = matchObservedOW(model, states, schedule, observed, varargin)
% Compute mismatch-function 

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
                 'state',               [],...
                 'from_states',         true,...% can be false for generic models
                 'matchOnlyProducers',  false, ...
                 'mismatchSum',         true, ...
                 'accumulateWells',       [], ...
                 'accumulateTypes',       []);
             
opt     = merge_options(opt, varargin{:});

dts   = schedule.step.val;
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
            state = model.getStateAD( states{tSteps(step)}, init);
        else
            state = opt.state;
        end
        qWs = model.FacilityModel.getProp(state,'qWs');
        qOs = model.FacilityModel.getProp(state,'qOs');
        bhp = model.FacilityModel.getProp(state,'bhp');
        assert(not(isnumeric(qWs))); 
        status = vertcat(state.wellSol.status);
     else
        state = states{tSteps(step)};
        [qWs, qOs, bhp] = deal( vertcatIfPresent(state.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state.wellSol, 'bhp', nw) );
       assert(isnumeric(qWs));
       status = vertcat(state.wellSol.status);
    end
     
    if ~all(status) || ~all(status_obs) 
        [bhp, bhp_obs] = expandToFull(bhp, bhp_obs, status, status_obs, true);
        [qWs, qWs_obs] = expandToFull(qWs, qWs_obs, status, status_obs, false);
        [qOs, qOs_obs] = expandToFull(qOs, qOs_obs, status, status_obs, false);
    end
    dt = dts(step);
    if opt.mismatchSum
        obj{step} = (dt/(totTime*nnz(matchCases)))*sum( ...
                        (ww*matchCases.*(qWs-qWs_obs)).^2 + ...
                        (wo*matchCases.*(qOs-qOs_obs)).^2 + ...
                        (wp*matchCases.*(bhp-bhp_obs)).^2 );
    else
        % output summands f_i^2 
        fac = dt/(totTime*nnz(matchCases));
        mm  = {fac*(ww*matchCases.*(qWs-qWs_obs)).^2, ...
               fac*(wo*matchCases.*(qOs-qOs_obs)).^2, ...
               fac*(wp*matchCases.*(bhp-bhp_obs)).^2};
        
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

