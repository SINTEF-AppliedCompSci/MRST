function report = convertUnitsOfReport(report, metric)
%Convert report object to standard metric units
%
% SYNOPSIS:
%   report = convertUnitsOfReport(report, metric)
%
% PARAMETERS:
%   report - Report data structure, as defined by function
%            'addToTimeStruct' for all time steps of a simulation run.
%            Assumed to contain the fields 'TIME', 'WBHP', 'WVPT', 'WOPR',
%            and 'WWPR'.
%
%   metric - Logical value indicating whether or not the output report
%            values should be given in metric or field units.  Typically
%            corresponds to 'grdecl.METRIC'.
%
%            The metric units are
%                - Pressure -> bar (barsa)
%                - Time     -> days
%                - Rate     -> meter^3 / day
%
%            The field units are
%                - Pressure -> psi (psia)
%                - Time     -> days
%                - Rate     -> stb / day
%
% RETURNS:
%   report - Report data structure whose values have been converted as per
%            the caller's request.
%
% SEE ALSO:
%   `wellCalculateProduction`, `addToTimeStruct`.

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


   assert (~isempty(report) && isstruct(report));

   if metric,
      u = struct('time', day, 'press', barsa, 'rate', meter^3/day);
   else
      u = struct('time', day, 'press', psia , 'rate', stb/day);
   end

   report.TIME = convertTo(report.TIME, u.time);
   report.WBHP = convertTo(report.WBHP, u.press);

   for f = {'WWPR', 'WOPR', 'WVPT'},
      report.(f{1}) = convertTo(report.(f{1}), u.rate);
   end
end
