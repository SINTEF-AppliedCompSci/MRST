function obj = matchObservedOW(model, states, schedule, observed, varargin)
% Compute mismatch-function 

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
opt     = struct('WaterRateWeight',     [] , ...
                 'OilRateWeight',       [] , ...
                 'BHPWeight',           [] , ...
                 'ComputePartials',     false, ...
                 'tStep' ,              [], ...
                 'state',               [],...
                 'from_states',         true,...% can be false for generic models
                 'matchOnlyProducers',  false);
             
opt     = merge_options(opt, varargin{:});


% pressure and saturaton vectors just used for place-holding
%p  = zeros(G.cells.num, 1);
%sW = zeros(G.cells.num, 1);

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
    %sol = wellSols{tSteps(step)};
    %[qWs, qOs, bhp] = deal(vertcat(sol.qWs), vertcat(sol.qOs), vertcat(sol.bhp) );
    sol_obs = observed{tSteps(step)};
    nw = numel([sol_obs.wellSol.bhp]);
    [qWs_obs, qOs_obs, bhp_obs] = deal( vertcatIfPresent(sol_obs.wellSol, 'qWs', nw), ...
                                        vertcatIfPresent(sol_obs.wellSol, 'qOs', nw), ...
                                        vertcatIfPresent(sol_obs.wellSol, 'bhp', nw) );

    if opt.matchOnlyProducers
        matchCases = (vertcat(sol.sign) < 0);
    else
        matchCases = true(nw,1);
    end
    
    [ww, wo, wp] = getWeights(qWs_obs, qOs_obs, bhp_obs, opt);
                                   
     if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state = model.getStateAD( states{tSteps(step)}, init);
        else
            state = opt.state;
        end
        qWs=model.FacilityModel.getProp(state,'qWs');
        qOs=model.FacilityModel.getProp(state,'qOs');
        bhp=model.FacilityModel.getProp(state,'bhp');
        assert(not(isnumeric(qWs))); 
     else
        state = states{tSteps(step)};
        [qWs, qOs, bhp] = deal( vertcatIfPresent(state.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state.wellSol, 'bhp', nw) );
       assert(isnumeric(qWs));                     
     end
    
    dt = dts(step);
    obj{step} = (dt/(totTime*nnz(matchCases)))*sum( ...
                                (ww*matchCases.*(qWs-qWs_obs)).^2 + ...
                                (wo*matchCases.*(qOs-qOs_obs)).^2 + ...
                                (wp*matchCases.*(bhp-bhp_obs)).^2 );
end
end

%--------------------------------------------------------------------------

function v = vertcatIfPresent(sol, fn, nw)
if isfield(sol, fn)
    v = vertcat(sol.(fn));
    assert(numel(v)==nw);
else
    v = zeros(nw,1);
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

