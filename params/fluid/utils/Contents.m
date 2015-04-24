% Files
%   constantProperties.m  - Construct fluid property evaluator from constant data
%   dinterpTable.m        - Compute derivative of one-dimensional interpolant, possibly using splines.
%   evalMultipleRegions.m - Evaluate fluid functions in multiple regions
%   interpRelPermTable.m  - Fast linear interpolation of tabulated relperm function (and derivative)
%   interpTable.m         - Interpolate a one-dimensional table, possibly using splines.
%   sgfn.m                - Construct gas relperm evaluation functions from SGFN table.
%   sgof.m                - Construct gas/oil rel-perm evaluation functions from SGOF table.
%   sof2.m                - Construct oil-water or oil-gas relperm eval. functions from SOF2 table.
%   sof3.m                - Construct oil-water and oil-gas relperm eval. functions from SOF3 table.
%   swfn.m                - Construct water relperm evaluation functions from SWFN table.
%   swof.m                - Construct water/oil rel-perm evaluation functions from SWOF table.
%   tabulatedSatFunc.m    - Construct rel-perm and capillary pressure evaluators from tabulated data

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
