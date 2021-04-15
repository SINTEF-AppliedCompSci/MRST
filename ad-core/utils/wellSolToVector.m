function [qWs, qOs, qGs, bhp] = wellSolToVector(wellsols)
%Extract selected summary vectors from cell array of well solutions
%
% SYNOPSIS:
%   [WaterRate, OilRate, GasRate, BHP] = wellSolToVector(wellSols)
%
% PARAMETERS:
%   wellSols - Cell array of well solution structures as produced by
%              runScheduleADI or simulateScheduleADI.  Each solution
%              structure must define the fields 'qWs', 'qOs', 'qGs', and
%              'bhp'.
%
% NOTE:
%   Function wellSolToVector does not support variable number of wells.
%
% RETURNS:
%   In the following 'nw' refers to the number of wells
%   (NUMEL(wellSols{1})) while 'nt' refers to the total number of time
%   steps (NUMEL(wellSols)).
%
%   WaterRate - Numeric array of size nt-by-nw of water rate at surface
%               conditions (unit m^3/s).
%
%   OilRate   - Numeric array of size nt-by-nw of oil rate at surface
%               conditions (unit m^3/s).
%
%   GasRate   - Numeric array of size nt-by-nw of gas rate at surface
%               conditions (unit m^3/s).
%
%   BHP       - Numeric array of size nt-by-nw of well bottom-hole pressure
%               values (unit Pascal).
%
% SEE ALSO:
%   simulateScheduleADI, runScheduleADI.

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

   extract1 = @(wsol, fld) [ wsol.(fld) ];
   extract2 = @( c ) vertcat(c{:});
   extract  = @(fld) extract2(cellfun(@(wsol) extract1(wsol, fld), ...
                                      wellsols, 'UniformOutput', false));

   sgn =        extract('sign');
   bhp =        extract('bhp');
   
   qWs = extract('qWs');
   qOs = extract('qOs');
   qGs = extract('qGs');
   if ~isempty(qWs), qWs = sgn .* qWs; end
   if ~isempty(qOs), qOs = sgn .* qOs; end
   if ~isempty(qGs), qGs = sgn .* qGs; end
end
