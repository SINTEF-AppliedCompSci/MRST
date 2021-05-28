function obj = NPVBlackOil(model, states, schedule, varargin)
% Compute net present value of a schedule with well solutions
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
opt     = struct('OilPrice',             1.0 , ...
                 'GasPrice',             0.1 , ...
                 'GasInjectionCost',     0.1 , ...
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               []);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
riw  = opt.WaterInjectionCost / stb;
rg  = opt.GasPrice   / stb;
rig  = opt.GasInjectionCost   / stb;


d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
%p  = zeros(G.cells.num, 1);
%sW = zeros(G.cells.num, 1);
%x  = zeros(G.cells.num, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    %sol = wellSols{tSteps(step)};
    state = states{tSteps(step)};
    nw = numel([state.wellSol.bhp]);
    if opt.ComputePartials
         %[~, ~, qWs, qOs, bhp] = ...
        %  initVariablesADI(p, sW, qWs, qOs, bhp);
        init=true;
        state = model.getStateAD( states{tSteps(step)}, init);
        qWs=model.FacilityModel.getProp(state,'qWs');
        qOs=model.FacilityModel.getProp(state,'qOs');
        qGs=model.FacilityModel.getProp(state,'qGs');
        bhp=model.FacilityModel.getProp(state,'bhp');
     else
        state = states{tSteps(step)};
        [qWs, qOs, qGs, bhp] = deal( vertcatIfPresent(state.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qGs', nw),...
                                vertcatIfPresent(state.wellSol, 'bhp', nw) );        
        injInx  = (vertcat(sol.sign) > 0);
        status = vertcat(sol.status);
        qWs = qWs(status);
        qOs = qOs(status);
        qGs = qGs(status);
        injInx = injInx(status);
    end
    %nW  = numel(qWs);
    %pBHP = zeros(nW, 1); %place-holder
  
    

    %if opt.ComputePartials
    %    [qWs, qWs, qWs, qWs, qOs, qGs, ignore] = ...
    %       initVariablesADI(p, sW, x, qWs, qOs, qGs, pBHP);       %#ok
    %
    %    clear ignore
    %end

    dt = dts(step);
    time = time + dt;

    prodInx = ~injInx;
    obj{step} = ( dt*(1+d)^(-time/year) )*...
                spones(ones(1, nW))*( (-ro*prodInx).*qOs +....
                              (-rg*prodInx - rig*injInx).*qGs ...
                             +(rw*prodInx - riw*injInx).*qWs );
end
