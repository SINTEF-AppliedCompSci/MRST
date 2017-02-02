function eq = setupWellControlEquationsSingleWell(sol, pBH, q_s, status, mix_s, model)
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
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

actPh = model.getActivePhases();
nPh = numel(q_s);
if model.water
    iw = 1;
end
if model.oil
    io = sum(actPh(1:2));
end
if model.gas
    ig = sum(actPh(1:3));
end


type = sol.type;
val  = sol.val;
bhp  = sol.bhp;

qt_s = 0;
for ph = 1:nPh
    qt_s = qt_s + q_s{ph};
end

setToZeroRate = val == 0 & ~strcmpi(type, 'bhp');
disabledWells  = ~status;

if ~status
    % The well is shut, set a trivial condition
    eq = (pBH(disabledWells) - bhp(disabledWells))/barsa;
elseif setToZeroRate
    % Well equation is simply zero rate - no reasonable control
    eq = qt_s;
else
    switch type
        case 'bhp'
            eq = pBH - val;
        case {'rate', 'vrat'}
            % Controlled by total rate
            eq = qt_s - val;
        case 'orat'
            assert(model.oil, 'Oil phase must be present to control wells on oil rates');
            eq = q_s{io} - val;
            % if mix_s(io) == 0;
            %     % No oil present, singular equation. Set to zero control.
            %     eq = qt_s;
            % end
        case 'wrat'
            assert(model.water, 'Water phase must be present to control wells on oil rates');
            eq = q_s{iw} - val;
            % if mix_s(iw) == 0;
            %     % Set to zero control
            %     eq = qt_s;
            % end
        case 'grat'
            assert(model.gas, 'Gas phase must be present to control wells on oil rates');
            eq = q_s{ig} - val;
            % if mix_s(ig) == 0;
            %     % Set to zero control
            %     eq = qt_s;
            % end
        case 'lrat'
            eq = q_s{iw} + q_s{io} - val;
            % if mix_s(iw) + mix_s(io) == 0;
            %     % Set to zero control
            %     eq = qt_s;
            % end
        otherwise
            error(['Unknown well control ', type, ' for well ', wellSol.name]);
    end
end
end
