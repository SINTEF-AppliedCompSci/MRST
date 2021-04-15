% RESULTINPUT
%
% Files
%   addFluxesToRestartStates         - Expand Reservoir State to Include Phase Fluxes.
%   convertRestartToStates           - states = convertRestartToStates(fn, G, varargin)
%   convertSummaryToWellSols         - [wellSols, time] = convertSummaryToWellSols(fn, unit)
%   eclOut2mrst                      - Undocumented Utility Function
%   initGridFromEclipseOutput        - Reads eclipse output files and converts to MRST-compatible datastructures
%   processEclipseRestartSpec        - Read unformatted (binary) ECLIPSE restart specification file (*RSSPEC) and 
%   readEclipseOutputFileFmt         - Read formatted (ASCII) ECLIPSE output/result file
%   readEclipseOutputFileUnFmt       - Read unformatted (binary) ECLIPSE output/result file
%   readEclipseRestartFmt            - Read formatted (text/ASCII) ECLIPSE restart data
%   readEclipseRestartUnFmt          - Read unformatted (binary) ECLIPSE unified or multiple restart file(s)
%   readEclipseRestartUnFmt_fallback - Read unformatted (binary) ECLIPSE restart data
%   readEclipseSummaryFmt            - Read formatted (text/ASCII) ECLIPSE summary data
%   readEclipseSummaryUnFmt          - Undocumented utility function
%   readEclipseSummaryUnFmt_fallback - Read unformatted (binary) ECLIPSE summary data

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
