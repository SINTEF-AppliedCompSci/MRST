function obj = NPVOWPolymer(G, wellSols, schedule, varargin)
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
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'PolymerInjectionCost', 0.1,...
                 'tStep' ,               []);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
rp  = opt.PolymerInjectionCost * (kilo*gram/meter^3);
d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
p  = zeros(G.cells.num, 1);
sW = zeros(G.cells.num, 1);
c = zeros(G.cells.num, 1);

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
    nW  = numel(sol);
    pBHP = zeros(nW, 1); %place-holder
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    if isfield(sol, 'poly')
        poly = vertcat(sol.poly);
    else
        poly = vertcat(sol.cWPoly);
    end
    injPoly = [sol.qWs] > 0 & [sol.sign] == 1;


    if opt.ComputePartials
        [qWs, qWs, qWs, qWs, qOs, ignore] = ...
           initVariablesADI(p, sW, c, qWs, qOs, pBHP);           %#ok

        clear ignore
    end
    poly = poly(injPoly);

    dt = dts(step);
    time = time + dt;

    injInx  = (vertcat(sol.sign) > 0);
    prodInx = ~injInx;
    obj{step} = ( dt*(1+d)^(-time/year) )*...
                spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
                             +(rw*prodInx - ri*injInx).*qWs );

    obj{step} =  obj{step} - ( dt*(1-d)^(time/year) )*spones(ones(1, sum(injPoly)))*...
                             (rp*poly.*qWs(injPoly) );
end
