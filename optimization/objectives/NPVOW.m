function obj = NPVOW(G, wellSols, schedule, varargin)
% Compute net present value of a schedule with well solutions

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               [],   ...
                 'signChangePenaltyFactor', 0);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
p  = zeros(G.cells.num, 1);
sW = zeros(G.cells.num, 1);

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
    sol = wellSols{tSteps(step)};
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    
    % Remove closed well.
    status = vertcat(sol.status);
    qWs = qWs(status);
    qOs = qOs(status);
    
    injectors = (vertcat(sol.sign) > 0);
    injectors = injectors(status);
    injecting = (qWs + qOs) > 0;
    sgnCh     = (injectors & ~injecting) | (~injectors & injecting);
    
    nW  = numel(qWs);
    pBHP = zeros(nW, 1); %place-holder

    if opt.ComputePartials
        [qWs, qWs, qWs, qOs, ignore] = ...
           initVariablesADI(p, sW, qWs, qOs, pBHP);                    %#ok

        clear ignore
    end

    dt = dts(step);
    time = time + dt;

    penalty = opt.signChangePenaltyFactor;
    
    producing = ~injecting;
    producers = ~injectors;
    if penalty == 0 || any(~sgnCh)
        % just compute assuming potential oil injection has cost ro
        obj{step} = ( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( -ro*qOs + (rw*producing - ri*injecting).*qWs );
    else
        % penalize signchange by penalty*ro
        ii = injecting & injectors;
        pp = producing & producers;
        obj{step} = ( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( -(ro*pp).*qOs + (rw*pp - ri*ii).*qWs - (ro*penalty*sgnCh).*abs(qWs+qOs));
    end
end
end
