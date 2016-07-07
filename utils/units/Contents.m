% UNITS
%   MRST's implementation of support for different units of measurement
%
% Files
%   atm           - Compute numerical value, in units of Pascal, of one atmosphere.
%   barsa         - Compute numerical value, in units of Pascal, of one bar.
%   centi         - One houndreth prefix.
%   convertFrom   - Convert physical quantity from given unit to equivalent SI.
%   convertTo     - Convert physical quantity from SI to equivalent given unit.
%   darcy         - Compute numerical value, in units of m^2, of the Darcy constant.
%   day           - Give numerical value, in units of seconds, of one day.
%   deci          - One tenth prefix.
%   dyne          - Compute numerical value, in units of Newton of one dyne.
%   ft            - Distance of one foot (in units of meters).
%   gallon        - Compute numerical value, in units of m^3, of one U.S. liquid gallon.
%   getUnitSystem - Define unit conversion factors for input data.
%   giga          - One billion (milliard) prefix.
%   gram          - Mass of one gram, in units of kilogram.
%   hour          - Time span of one hour (in units of seconds).
%   inch          - Distance of one inch (in units of meters).
%   Kelvin        - Temperature of one Kelvin (in units Kelvin)
%   kilo          - One thousand prefix.
%   kilogram      - Mass of one kilogram, in units of kilogram.
%   lbf           - Force excerted by a mass of one avoirdupois pound at Tellus equator.
%   mega          - One million prefix.
%   meter         - Distance of one meter (in units of meters).
%   micro         - One millionth prefix.
%   milli         - One thousandth prefix.
%   minute        - Time span of one minute (in units of seconds).
%   Newton        - Force of one Newton, in units of Newton.
%   Pascal        - Compute numerical value, in units of Pascal, of one Pascal.
%   poise         - Compute numerical value, in units of Pa*s, of one poise (P).
%   pound         - Mass of one avoirdupois pound, in units of kilogram.
%   psia          - Compute numerical value, in units of Pascal, of one Psi.
%   Rankine       - Temperature of one Rankine (in units of Kelvin).
%   second        - Time span of one second (in units of seconds).
%   stb           - Compute numerical value, in units of m^3, of one standard barrel.
%   year          - Give numerical value, in units of seconds, of one year.

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
