function eq = setupWellControlEquations(sol, pBH, q_s, status, mix_s, model)
%Setup well controll (residual) equations 
%
% SYNOPSIS:
%   eq = setupWellControlEquations(sol, pBH, q_s, status, mix_s, model)
%
% PARAMETERS:
%   sol         - List of current well solution structures containing
%                 control type for each well 
%   pBH         - Vector of well bhps
%   q_s         - List of vectors of well component volume-rates 
%                 (surface conds) 
%   status      - Logic vector of well statuses
%   mix_s       - List of vectors containing volumetric mixture of components 
%                 in wellbroe at connections (surface conds).
%   model       - Simulation model.
%
% RETURNS:
%   eq          - Well control equations
%
% SEE ALSO:
%   WellModel, computeWellContributionsNew.

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
[iw, io, ig] = deal([]);
if model.water
    iw = 1;
end
if model.oil
    io = 2;
end
if model.gas
    ig = 3;
end


type = {sol.type}';
val  = vertcat(sol.val);
bhp  = vertcat(sol.bhp);

qt_s = q_s{1};
for ph = 2:numel(q_s)
    qt_s = qt_s + q_s{ph};
end

setToZeroRate = and(val ==0, ~cellfun(@(x)strcmp('bhp',x), type));
disabledWells  = ~status;

eq = pBH; %just to initialize to whatever class pBH is
% bhp (injector or producer)
inx = getIndicesPerControl(type, status, 'bhp');
if ~isempty(inx)
    eq(inx) = pBH(inx) - val(inx);
end

%rate (injector
inx = getIndicesPerControl(type, status, 'rate');
if ~isempty(inx)
    eq(inx) = qt_s(inx) - val(inx);
end

%orat (producer)
inx = getIndicesPerControl(type, status, 'orat');
if ~isempty(inx)
    eq(inx) = q_s{io}(inx)-val(inx);
    prob    = mix_s(inx,io)==0;
    setToZeroRate(inx(prob)) = true;
end

%wrat (producer)
inx = getIndicesPerControl(type, status, 'wrat');
if ~isempty(inx)
    eq(inx) = q_s{iw}(inx)-val(inx);
    prob    = mix_s(inx,iw)==0;
    setToZeroRate(inx(prob)) = true;
end

%grat (producer)
inx = getIndicesPerControl(type, status, 'grat');
if ~isempty(inx)
    eq(inx) = q_s{ig}(inx)-val(inx);
    prob    = mix_s(inx,ig)==0;
    setToZeroRate(inx(prob)) = true;
end

%lrat (producer)
inx = getIndicesPerControl(type, status, 'lrat');
if ~isempty(inx)
    eq(inx) = q_s{iw}(inx)+q_s{io}(inx)-val(inx);
    prob    = (mix_s(inx,iw)+mix_s(inx,io))==0;
    setToZeroRate(inx(prob)) = true;
end

%vrat (producer) - special volume rate
inx = getIndicesPerControl(type, status, 'vrat');
if ~isempty(inx)
    eq(inx) = qt_s(inx)-val(inx);
end

%rate is zero or control on phaserate impossible
if ~isempty(setToZeroRate)
    eq(setToZeroRate) = qt_s(setToZeroRate);
end

% In the case when wells are shut, set a trivial condition
if any(disabledWells)
    rateScale = max(mean(abs(double(qt_s))), sqrt(eps));
    bhpScale  = max(mean(abs(bhp)),          1);
    eq(disabledWells) = rateScale*(pBH(disabledWells) - bhp(disabledWells))./bhpScale;
end
end

function inx = getIndicesPerControl(type, active, name)
    inx = find(cellfun(@(x)strcmp(name, x), type) & active);
end
